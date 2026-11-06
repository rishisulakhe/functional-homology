#!/usr/bin/env python3
"""Generate mean-pooled ESM-2 protein embeddings for a FASTA proteome.

Each protein is embedded by:
  1. truncating the sequence to ``--max-residues`` (default 1022),
  2. tokenising with the model tokenizer (special tokens included),
  3. running the model under ``torch.no_grad()``,
  4. mean-pooling residue-level hidden states, excluding pad/cls/eos tokens.

Outputs ``<outdir>/<name>.npy`` (float32, one row per protein) and
``<outdir>/<name>.ids.txt`` (accessions in identical order), plus a small
``<name>.meta.json`` describing the run.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import numpy as np
import torch
from tqdm import tqdm
from transformers import AutoModel, AutoTokenizer

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fh_common import read_fasta  # noqa: E402

LOG = logging.getLogger("embed")


def mean_pool_masked(hidden: torch.Tensor, token_ids: torch.Tensor,
                     attention_mask: torch.Tensor, special_ids: set[int]) -> torch.Tensor:
    """Mean-pool hidden states over real residues only.

    ``hidden``: (batch, seq_len, dim); ``token_ids``: (batch, seq_len);
    ``attention_mask``: (batch, seq_len).  Special token ids (cls/eos/pad/unk
    as provided) are excluded.  Sequences with no residues fall back to the
    masked mean over attention_mask, then to the first token.
    """
    mask = attention_mask.to(hidden.dtype).clone()
    for tok_id in special_ids:
        if tok_id is not None:
            mask = mask * (token_ids != tok_id).to(hidden.dtype)
    denom = mask.sum(dim=1, keepdim=True)
    pooled = (hidden * mask.unsqueeze(-1)).sum(dim=1) / denom.clamp(min=1.0)
    empty = denom.squeeze(1) == 0
    if empty.any():
        fallback_mask = attention_mask.to(hidden.dtype)
        fallback = (hidden * fallback_mask.unsqueeze(-1)).sum(dim=1) / \
            fallback_mask.sum(dim=1, keepdim=True).clamp(min=1.0)
        pooled[empty] = fallback[empty]
    return pooled


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--fasta", type=Path, required=True, help="input proteome FASTA")
    p.add_argument("--outdir", type=Path, required=True, help="output directory")
    p.add_argument("--name", default=None,
                   help="output prefix (default: input FASTA stem)")
    p.add_argument("--model", default="facebook/esm2_t6_8M_UR50D",
                   help="Hugging Face ESM model id (default: %(default)s)")
    p.add_argument("--batch-size", type=int, default=8)
    p.add_argument("--max-residues", type=int, default=1022,
                   help="truncate sequences longer than this many residues")
    p.add_argument("--device", default="auto", help="auto, cpu, cuda, or mps")
    p.add_argument("--threads", type=int, default=None,
                   help="torch CPU threads (default: leave torch default)")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def pick_device(requested: str) -> torch.device:
    if requested != "auto":
        return torch.device(requested)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if getattr(torch.backends, "mps", None) and torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S")
    if args.threads:
        torch.set_num_threads(args.threads)

    args.outdir.mkdir(parents=True, exist_ok=True)
    name = args.name or args.fasta.stem
    npy_path = args.outdir / f"{name}.npy"
    ids_path = args.outdir / f"{name}.ids.txt"
    meta_path = args.outdir / f"{name}.meta.json"

    records = read_fasta(args.fasta)
    ids = [acc for acc, _ in records]
    seqs = [seq for _, seq in records]
    lengths = np.array([len(s) for s in seqs], dtype=np.int64)
    truncated = [acc for (acc, s) in records if len(s) > args.max_residues]
    seqs = [s[: args.max_residues] for s in seqs]
    LOG.info("Loaded %d proteins from %s (median length %d aa)",
             len(records), args.fasta, int(np.median(lengths)))
    if truncated:
        LOG.info("Truncated %d/%d sequences to %d residues (longest original: %d aa)",
                 len(truncated), len(records), args.max_residues, int(lengths.max()))
    else:
        LOG.info("No sequences exceeded %d residues", args.max_residues)

    device = pick_device(args.device)
    LOG.info("Loading model %s on %s ...", args.model, device)
    t0 = time.time()
    tokenizer = AutoTokenizer.from_pretrained(args.model)
    model = AutoModel.from_pretrained(args.model)
    model.to(device).eval()
    LOG.info("Model loaded in %.1fs (hidden size %d)",
             time.time() - t0, int(model.config.hidden_size))

    special_ids = {tokenizer.cls_token_id, tokenizer.eos_token_id,
                   tokenizer.pad_token_id}
    if getattr(tokenizer, "mask_token_id", None) is not None:
        special_ids.add(tokenizer.mask_token_id)

    out = np.zeros((len(seqs), int(model.config.hidden_size)), dtype=np.float32)
    done = 0
    with torch.no_grad():
        iterator = range(0, len(seqs), args.batch_size)
        for start in tqdm(iterator, desc=f"embedding {name}", unit="batch",
                          disable=not sys.stderr.isatty()):
            batch = seqs[start:start + args.batch_size]
            enc = tokenizer(batch, return_tensors="pt", padding=True,
                            truncation=True, max_length=args.max_residues + 2)
            enc = {k: v.to(device) for k, v in enc.items()}
            hidden = model(**enc).last_hidden_state
            pooled = mean_pool_masked(hidden, enc["input_ids"],
                                      enc["attention_mask"], special_ids)
            out[start:start + len(batch)] = pooled.float().cpu().numpy()
            done += len(batch)
            if not sys.stderr.isatty() and (done == len(seqs) or done % (args.batch_size * 20) == 0):
                LOG.info("  embedded %d/%d proteins", done, len(seqs))

    if not np.isfinite(out).all():
        raise RuntimeError("Non-finite values in embeddings; aborting")
    np.save(npy_path, out)
    ids_path.write_text("\n".join(ids) + "\n", encoding="utf-8")
    meta = {
        "model": args.model,
        "n_proteins": len(ids),
        "embedding_dim": int(out.shape[1]),
        "max_residues": args.max_residues,
        "n_truncated": len(truncated),
        "truncated_ids": truncated,
        "batch_size": args.batch_size,
        "device": str(device),
        "fasta": str(args.fasta),
    }
    meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    LOG.info("Saved %s shape=%s and %s", npy_path, out.shape, ids_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
