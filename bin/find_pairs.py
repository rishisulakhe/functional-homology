#!/usr/bin/env python3
"""Combine ESM-2 embeddings and DIAMOND hits into candidate protein pairs.

For every protein in proteome A, the top-k most cosine-similar proteins in
proteome B are retained.  Each pair is annotated with whether DIAMOND reported
a significant hit (``sequence_hit``).  A random control set, matched in size to
the embedding-only candidate set, is sampled uniformly from the cross-product
of the two proteomes while excluding every top-k pair and every DIAMOND hit.

Outputs (in ``--outdir``):
  pairs.tsv         all top-k embedding pairs, joined with DIAMOND hits
  diamond_pairs.tsv all DIAMOND hits (best HSP per pair)
  random_pairs.tsv  random control pairs
  pairs_summary.json
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fh_common import accession_from_header  # noqa: E402

LOG = logging.getLogger("find_pairs")

PAIR_COLUMNS = ["id_a", "id_b", "cosine", "sequence_hit", "pident", "evalue", "bitscore"]
DIAMOND_COLUMNS = ["id_a", "id_b", "pident", "evalue", "bitscore"]


def load_ids(path: str | Path) -> list[str]:
    ids = [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines()
           if line.strip()]
    if not ids:
        raise ValueError(f"No ids found in {path}")
    if len(set(ids)) != len(ids):
        raise ValueError(f"Duplicate ids in {path}")
    return ids


def normalize_id(seq_id: str) -> str:
    """Normalise a DIAMOND sequence id to a bare UniProt accession."""
    return accession_from_header(seq_id)


def l2_normalize(mat: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    """Row-wise L2 normalisation; zero rows stay zero."""
    mat = np.asarray(mat, dtype=np.float32)
    norms = np.linalg.norm(mat, axis=1, keepdims=True)
    return mat / np.maximum(norms, eps)


def top_k_cosine(emb_a: np.ndarray, emb_b: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (indices, scores) of the top-k B proteins for each A protein."""
    if emb_a.ndim != 2 or emb_b.ndim != 2:
        raise ValueError("Embedding matrices must be 2-D")
    if emb_a.shape[1] != emb_b.shape[1]:
        raise ValueError(f"Embedding dimensions differ: {emb_a.shape[1]} vs {emb_b.shape[1]}")
    if k < 1:
        raise ValueError("--top-k must be >= 1")
    a = l2_normalize(emb_a)
    b = l2_normalize(emb_b)
    sim = a @ b.T
    k_eff = min(k, sim.shape[1])
    part = np.argpartition(-sim, kth=k_eff - 1, axis=1)[:, :k_eff]
    rows = np.arange(sim.shape[0])[:, None]
    vals = sim[rows, part]
    order = np.argsort(-vals, axis=1)
    idx = np.take_along_axis(part, order, axis=1)
    scores = np.take_along_axis(vals, order, axis=1)
    return idx, scores.astype(np.float32)


def parse_diamond_hits(path: str | Path, evalue_cutoff: float,
                       valid_a: set[str] | None = None,
                       valid_b: set[str] | None = None) -> tuple[dict, dict]:
    """Parse DIAMOND outfmt-6 hits, keeping the best HSP per (a, b) pair."""
    best: dict[tuple[str, str], dict] = {}
    stats = {"lines": 0, "malformed": 0, "non_significant": 0,
             "unknown_id": 0, "unique_pairs": 0}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            stats["lines"] += 1
            fields = line.split("\t")
            if len(fields) < 5:
                stats["malformed"] += 1
                continue
            try:
                a = normalize_id(fields[0])
                b = normalize_id(fields[1])
                pident = float(fields[2])
                evalue = float(fields[3])
                bitscore = float(fields[4])
            except ValueError:
                stats["malformed"] += 1
                continue
            if evalue > evalue_cutoff:
                stats["non_significant"] += 1
                continue
            if (valid_a is not None and a not in valid_a) or \
               (valid_b is not None and b not in valid_b):
                stats["unknown_id"] += 1
                continue
            key = (a, b)
            old = best.get(key)
            if old is None or (evalue, -bitscore, -pident) < \
                    (old["evalue"], -old["bitscore"], -old["pident"]):
                best[key] = {"id_a": a, "id_b": b, "pident": pident,
                             "evalue": evalue, "bitscore": bitscore}
    stats["unique_pairs"] = len(best)
    return best, stats


def build_random_control(ids_a: list[str], ids_b: list[str], n: int,
                         excluded: set[tuple[str, str]], seed: int,
                         max_attempts_factor: int = 100) -> list[tuple[str, str]]:
    """Sample ``n`` unique random cross-organism pairs not in ``excluded``.

    Deterministic for a given seed.  If the requested size exceeds the
    available space, returns as many as possible (and logs a warning).
    """
    if n <= 0:
        return []
    total = len(ids_a) * len(ids_b)
    available = total - len(excluded)
    if available < n:
        LOG.warning("Random control requested %d pairs but only %d available; "
                    "returning %d", n, max(available, 0), max(available, 0))
        n = max(available, 0)
    rng = random.Random(seed)
    chosen: set[tuple[str, str]] = set()
    max_attempts = max(1000, n * max_attempts_factor)
    attempts = 0
    while len(chosen) < n and attempts < max_attempts:
        attempts += 1
        pair = (ids_a[rng.randrange(len(ids_a))], ids_b[rng.randrange(len(ids_b))])
        if pair in excluded or pair in chosen:
            continue
        chosen.add(pair)
    if len(chosen) < n:
        LOG.warning("Random sampler reached attempt cap: %d/%d pairs", len(chosen), n)
    return sorted(chosen)


def write_tsv(path: Path, columns: list[str], rows: list[list]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--embeddings-a", type=Path, required=True)
    p.add_argument("--ids-a", type=Path, required=True)
    p.add_argument("--embeddings-b", type=Path, required=True)
    p.add_argument("--ids-b", type=Path, required=True)
    p.add_argument("--diamond-hits", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--top-k", type=int, default=5)
    p.add_argument("--cosine-threshold", type=float, default=0.9)
    p.add_argument("--evalue-cutoff", type=float, default=1e-3)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S")
    args.outdir.mkdir(parents=True, exist_ok=True)

    ids_a = load_ids(args.ids_a)
    ids_b = load_ids(args.ids_b)
    emb_a = np.load(args.embeddings_a)
    emb_b = np.load(args.embeddings_b)
    if emb_a.shape[0] != len(ids_a) or emb_b.shape[0] != len(ids_b):
        raise ValueError(
            f"Embedding/id mismatch: A {emb_a.shape[0]} vs {len(ids_a)} ids; "
            f"B {emb_b.shape[0]} vs {len(ids_b)} ids")
    if not (np.isfinite(emb_a).all() and np.isfinite(emb_b).all()):
        raise ValueError("Embeddings contain non-finite values")
    LOG.info("Embeddings A %s, B %s", emb_a.shape, emb_b.shape)

    idx, scores = top_k_cosine(emb_a, emb_b, args.top_k)
    LOG.info("Computed top-%d matches for %d A proteins", idx.shape[1], idx.shape[0])

    hits, hit_stats = parse_diamond_hits(
        args.diamond_hits, args.evalue_cutoff, set(ids_a), set(ids_b))
    LOG.info("DIAMOND: %d unique significant pairs (from %d lines; skipped %s)",
             len(hits), hit_stats["lines"],
             {k: v for k, v in hit_stats.items() if v and k not in ("lines", "unique_pairs")})

    pair_rows: list[list] = []
    embedding_only: list[list] = []
    topk_keys: set[tuple[str, str]] = set()
    n_seq_hits_topk = 0
    for i, a_id in enumerate(ids_a):
        for j, cos in zip(idx[i], scores[i]):
            b_id = ids_b[int(j)]
            topk_keys.add((a_id, b_id))
            hit = hits.get((a_id, b_id))
            seq_hit = hit is not None
            if seq_hit:
                n_seq_hits_topk += 1
            row = [a_id, b_id, f"{float(cos):.6f}", str(seq_hit),
                   f"{hit['pident']:.2f}" if hit else "",
                   f"{hit['evalue']:.3g}" if hit else "",
                   f"{hit['bitscore']:.3g}" if hit else ""]
            pair_rows.append(row)
            if (not seq_hit) and float(cos) >= args.cosine_threshold:
                embedding_only.append(row)

    diamond_rows = [[h["id_a"], h["id_b"], f"{h['pident']:.2f}",
                     f"{h['evalue']:.3g}", f"{h['bitscore']:.3g}"]
                    for h in hits.values()]
    diamond_rows.sort(key=lambda r: (r[0], r[1]))
    write_tsv(args.outdir / "diamond_pairs.tsv", DIAMOND_COLUMNS, diamond_rows)

    excluded = topk_keys | set(hits.keys())
    random_pairs = build_random_control(ids_a, ids_b, len(embedding_only),
                                        excluded, args.seed)
    random_rows = [[a, b] for a, b in random_pairs]

    write_tsv(args.outdir / "pairs.tsv", PAIR_COLUMNS, pair_rows)
    write_tsv(args.outdir / "random_pairs.tsv", ["id_a", "id_b"], random_rows)
    LOG.info("Wrote %d top-k pairs (%d embedding-only above cos %.3f), "
             "%d DIAMOND pairs, %d random pairs",
             len(pair_rows), len(embedding_only), args.cosine_threshold,
             len(diamond_rows), len(random_rows))

    summary = {
        "params": {
            "top_k": args.top_k,
            "cosine_threshold": args.cosine_threshold,
            "evalue_cutoff": args.evalue_cutoff,
            "seed": args.seed,
        },
        "n_proteins_a": len(ids_a),
        "n_proteins_b": len(ids_b),
        "n_topk_pairs": len(pair_rows),
        "n_embedding_only": len(embedding_only),
        "n_sequence_hits_topk": n_seq_hits_topk,
        "n_diamond_pairs": len(diamond_rows),
        "n_random_pairs": len(random_rows),
        "diamond_parse_stats": hit_stats,
    }
    (args.outdir / "pairs_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
