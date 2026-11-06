#!/usr/bin/env python3
"""Download and prepare the functionalhomology input data.

Sources
-------
* UniProtKB reviewed (Swiss-Prot) FASTA, via the UniProt REST stream endpoint.
* Experimental GO annotations per organism, via the QuickGO annotation search
  API.  We query the ECO term ``ECO:0000269`` ("experimental evidence") so the
  API returns all descendants, then keep only the GO experimental evidence
  codes listed in ``EXPERIMENTAL_CODES``.  ``NOT`` qualifiers are dropped.
* ``go-basic.obo`` from the Gene Ontology Consortium.

The script can also build a small, seeded test dataset (FASTA subset,
annotation subset and an OBO slice closed under is_a/part_of ancestors).

Examples
--------
    python scripts/download_data.py --mode all
    python scripts/download_data.py --mode full --outdir data
    python scripts/download_data.py --mode test --test-size 200 --seed 42
"""

from __future__ import annotations

import argparse
import json
import logging
import random
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import requests
from Bio import SeqIO
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

LOG = logging.getLogger("download_data")

UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"
QUICKGO_SEARCH = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
GO_OBO_URL = "https://current.geneontology.org/ontology/go-basic.obo"

ORGANISMS = {
    "ecoli": {
        "taxid": 83333,
        "label": "Escherichia coli K-12",
        "fasta": "proteome_a.fasta",
        "go": "go_annotations_ecoli.tsv",
        "uniprot_query": "(organism_id:83333) AND (reviewed:true)",
    },
    "bsubtilis": {
        "taxid": 224308,
        "label": "Bacillus subtilis 168",
        "fasta": "proteome_b.fasta",
        "go": "go_annotations_bsubtilis.tsv",
        "uniprot_query": "(organism_id:224308) AND (reviewed:true)",
    },
}

# GO consortium experimental / high-throughput experimental evidence codes.
EXPERIMENTAL_CODES = {
    "EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
    "HTP", "HDA", "HMP", "HGI", "HEP",
}

GO_COLUMNS = [
    "accession", "go_id", "evidence_code", "qualifier", "aspect",
    "reference", "assigned_by", "organism", "taxon_id",
]


def setup_logging(verbose: bool = False) -> None:
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )


def make_session() -> requests.Session:
    sess = requests.Session()
    retry = Retry(
        total=6,
        backoff_factor=1.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset({"GET"}),
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=8, pool_maxsize=8)
    sess.mount("https://", adapter)
    sess.headers.update({"User-Agent": "functionalhomology/1.0 (research pipeline)"})
    return sess


def download_uniprot_fasta(session: requests.Session, spec: dict, dest: Path) -> dict:
    """Download reviewed proteins for one organism. Returns metadata."""
    LOG.info("Downloading %s (taxid %s) from UniProt ...", spec["label"], spec["taxid"])
    params = {"format": "fasta", "query": spec["uniprot_query"]}
    with session.get(UNIPROT_STREAM, params=params, stream=True, timeout=300) as r:
        r.raise_for_status()
        release = r.headers.get("X-UniProt-Release", "unknown")
        tmp = dest.with_suffix(dest.suffix + ".part")
        with tmp.open("wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    fh.write(chunk)
        tmp.replace(dest)
    n = sum(1 for _ in SeqIO.parse(str(dest), "fasta"))
    size = dest.stat().st_size
    LOG.info("  -> %s (%d sequences, %.1f MiB, UniProt release %s)",
             dest, n, size / 1e6, release)
    if n == 0:
        raise RuntimeError(f"No sequences downloaded to {dest}")
    return {"sequences": n, "bytes": size, "uniprot_release": release}


def _quickgo_page(session: requests.Session, taxid: int, page: int, limit: int) -> dict:
    params = {
        "taxonId": taxid,
        "geneProductType": "protein",
        "evidenceCode": "ECO:0000269",
        "limit": limit,
        "page": page,
    }
    r = session.get(QUICKGO_SEARCH, params=params, timeout=180)
    r.raise_for_status()
    return r.json()


def download_quickgo_annotations(
    session: requests.Session,
    spec: dict,
    dest: Path,
    limit: int = 200,
    workers: int = 4,
) -> dict:
    """Download experimental GO annotations for one organism from QuickGO."""
    taxid = spec["taxid"]
    LOG.info("Downloading experimental GO annotations for %s ...", spec["label"])
    first = _quickgo_page(session, taxid, 1, limit)
    total_pages = int(first.get("pageInfo", {}).get("total", 1))
    total_hits = int(first.get("numberOfHits", 0))
    LOG.info("  %d matching annotations over %d pages (limit=%d)",
             total_hits, total_pages, limit)

    pages: dict[int, list] = {1: first.get("results", [])}

    def fetch(page: int) -> tuple[int, list]:
        return page, _quickgo_page(session, taxid, page, limit).get("results", [])

    if total_pages > 1:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(fetch, p): p for p in range(2, total_pages + 1)}
            done = 0
            for fut in as_completed(futures):
                page, results = fut.result()
                pages[page] = results
                done += 1
                if done % 20 == 0 or done == len(futures):
                    LOG.info("  fetched %d/%d pages", done + 1, total_pages)

    rows: set[tuple] = set()
    dropped = {"non_uniprot": 0, "evidence": 0, "not_qualifier": 0, "bad_go": 0}
    for page in range(1, total_pages + 1):
        for rec in pages.get(page, []):
            gp = rec.get("geneProductId") or ""
            if not gp.startswith("UniProtKB:"):
                dropped["non_uniprot"] += 1
                continue
            acc = gp.split(":", 1)[1]
            ev = (rec.get("goEvidence") or "").strip()
            if ev not in EXPERIMENTAL_CODES:
                dropped["evidence"] += 1
                continue
            qual = (rec.get("qualifier") or "").strip()
            if "NOT" in qual.upper().split("|"):
                dropped["not_qualifier"] += 1
                continue
            go_id = (rec.get("goId") or "").strip()
            if not go_id.startswith("GO:"):
                dropped["bad_go"] += 1
                continue
            rows.add((
                acc, go_id, ev, qual, rec.get("goAspect") or "",
                rec.get("reference") or "", rec.get("assignedBy") or "",
                spec["label"], str(taxid),
            ))

    ordered = sorted(rows)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with tmp.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(GO_COLUMNS) + "\n")
        for row in ordered:
            fh.write("\t".join(row) + "\n")
    tmp.replace(dest)
    n_prot = len({r[0] for r in ordered})
    n_terms = len({r[1] for r in ordered})
    LOG.info("  -> %s (%d annotations, %d proteins, %d GO terms; dropped %s)",
             dest, len(ordered), n_prot, n_terms, dropped)
    return {
        "annotations": len(ordered), "proteins": n_prot, "terms": n_terms,
        "dropped": dropped, "quickgo_numberOfHits": total_hits,
        "evidence_codes": sorted({r[2] for r in ordered}),
    }


def download_obo(session: requests.Session, dest: Path) -> dict:
    LOG.info("Downloading go-basic.obo ...")
    with session.get(GO_OBO_URL, stream=True, timeout=300) as r:
        r.raise_for_status()
        tmp = dest.with_suffix(dest.suffix + ".part")
        with tmp.open("wb") as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                if chunk:
                    fh.write(chunk)
        tmp.replace(dest)
    size = dest.stat().st_size
    LOG.info("  -> %s (%.1f MiB)", dest, size / 1e6)
    return {"bytes": size}


def merge_annotations(paths: Iterable[Path], dest: Path) -> int:
    """Merge organism annotation TSVs into one file, de-duplicated."""
    rows: set[tuple] = set()
    for path in paths:
        with path.open(encoding="utf-8") as fh:
            header = fh.readline()
            if not header:
                continue
            for line in fh:
                rows.add(tuple(line.rstrip("\n").split("\t")))
    ordered = sorted(rows)
    with dest.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(GO_COLUMNS) + "\n")
        for row in ordered:
            fh.write("\t".join(row) + "\n")
    LOG.info("Merged %d annotation rows -> %s", len(ordered), dest)
    return len(ordered)


def parse_obo(path: Path):
    """Return (header, term_blocks, parents) from an OBO file.

    parents maps GO id -> set of GO ids connected by is_a or part_of.
    """
    header: list[str] = []
    blocks: dict[str, list[str]] = {}
    parents: dict[str, set[str]] = {}
    current: list[str] | None = None
    current_id: str | None = None
    in_term = False

    def flush() -> None:
        nonlocal current, current_id
        if in_term and current_id:
            blocks[current_id] = current or []
        current, current_id = None, None

    with path.open(encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if line.startswith("[Term]"):
                flush()
                in_term = True
                current = [line]
                continue
            if line.startswith("[") and line != "[Term]":
                flush()
                in_term = False
                current = None
                continue
            if not in_term:
                if not blocks:
                    header.append(line)
                continue
            current = current or []
            current.append(line)
            if line.startswith("id: "):
                current_id = line[4:].strip()
                parents.setdefault(current_id, set())
            elif line.startswith("is_a: "):
                parent = line[6:].split("!", 1)[0].strip()
                if current_id and parent.startswith("GO:"):
                    parents.setdefault(current_id, set()).add(parent)
            elif line.startswith("relationship: part_of "):
                parent = line.split("part_of", 1)[1].split("!", 1)[0].strip()
                if current_id and parent.startswith("GO:"):
                    parents.setdefault(current_id, set()).add(parent)
    flush()
    return header, blocks, parents


def ancestor_closure(terms: Iterable[str], parents: dict[str, set[str]]) -> set[str]:
    seen: set[str] = set()
    stack = list(terms)
    while stack:
        term = stack.pop()
        if term in seen:
            continue
        seen.add(term)
        stack.extend(parents.get(term, ()))
    return seen


def write_mini_obo(full_obo: Path, used_terms: Iterable[str], dest: Path) -> int:
    """Write an OBO slice containing used terms and all their ancestors."""
    header, blocks, parents = parse_obo(full_obo)
    keep = ancestor_closure(used_terms, parents)
    kept = sorted(t for t in keep if t in blocks)
    missing = sorted(t for t in keep if t not in blocks)

    with dest.open("w", encoding="utf-8") as fh:
        fh.write("format-version: 1.2\n")
        fh.write("subset: functionalhomology-test-slice\n")
        fh.write("remark: slice of go-basic.obo closed under is_a and part_of\n\n")
        for term in kept:
            block = blocks[term]
            for line in block:
                fh.write(line + "\n")
            fh.write("\n")
    LOG.info("Wrote OBO slice %s (%d terms; %d referenced terms missing from OBO)",
             dest, len(kept), len(missing))
    return len(kept)


def read_annotations(path: Path) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


def write_annotations(rows: Iterable[dict], dest: Path) -> None:
    rows = list(rows)
    with dest.open("w", encoding="utf-8") as fh:
        fh.write("\t".join(GO_COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")) for c in GO_COLUMNS) + "\n")


def build_test_dataset(
    full_fastas: dict[str, Path],
    full_annotations: Path,
    full_obo: Path,
    testdir: Path,
    test_size: int,
    seed: int,
) -> dict:
    """Create the committed, tiny, seeded test dataset."""
    testdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    meta: dict = {"seed": seed, "requested_test_size": test_size, "organisms": {}}
    test_accessions: set[str] = set()

    for key, spec in ORGANISMS.items():
        records = list(SeqIO.parse(str(full_fastas[key]), "fasta"))
        if not records:
            raise RuntimeError(f"No records in {full_fastas[key]}")
        n = min(test_size, len(records))
        idx = sorted(rng.sample(range(len(records)), n))
        chosen = [records[i] for i in idx]
        dest = testdir / spec["fasta"]
        SeqIO.write(chosen, str(dest), "fasta")
        accs = {rec.id.split("|")[1] if "|" in rec.id else rec.id.split()[0]
                for rec in chosen}
        test_accessions |= accs
        meta["organisms"][key] = {
            "label": spec["label"], "taxid": spec["taxid"],
            "n_proteins": len(chosen), "fasta": dest.name,
        }
        LOG.info("Wrote test FASTA %s (%d proteins)", dest, len(chosen))

    rows = [r for r in read_annotations(full_annotations)
            if r["accession"] in test_accessions]
    dest_ann = testdir / "go_annotations.tsv"
    write_annotations(rows, dest_ann)
    used_terms = {r["go_id"] for r in rows}
    LOG.info("Wrote test annotations %s (%d rows, %d proteins, %d terms)",
             dest_ann, len(rows), len({r["accession"] for r in rows}),
             len(used_terms))

    dest_obo = testdir / "go-basic.mini.obo"
    n_terms = write_mini_obo(full_obo, used_terms, dest_obo)
    meta.update({
        "n_annotations": len(rows),
        "n_annotated_proteins": len({r["accession"] for r in rows}),
        "n_go_terms_annotated": len(used_terms),
        "n_go_terms_in_obo_slice": n_terms,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    })
    (testdir / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    LOG.info("Test dataset ready in %s", testdir)
    return meta


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", choices=["all", "full", "test"], default="all",
                   help="download full data, build test data, or both (default: all)")
    p.add_argument("--outdir", type=Path, default=Path("data"),
                   help="directory for full data (default: data)")
    p.add_argument("--testdir", type=Path, default=Path("testdata"),
                   help="directory for test data (default: testdata)")
    p.add_argument("--test-size", type=int, default=200,
                   help="proteins sampled per organism for the test set (default: 200)")
    p.add_argument("--seed", type=int, default=42, help="test sampling seed (default: 42)")
    p.add_argument("--force", action="store_true",
                   help="re-download even if output files already exist")
    p.add_argument("--workers", type=int, default=4,
                   help="parallel QuickGO page requests (default: 4)")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    setup_logging(args.verbose)
    outdir: Path = args.outdir
    testdir: Path = args.testdir
    outdir.mkdir(parents=True, exist_ok=True)
    session = make_session()
    meta: dict = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sources": {
            "uniprot_stream": UNIPROT_STREAM,
            "quickgo_search": QUICKGO_SEARCH,
            "go_obo": GO_OBO_URL,
        },
        "experimental_evidence_codes": sorted(EXPERIMENTAL_CODES),
        "organisms": {},
    }

    full_fastas: dict[str, Path] = {}
    full_go: dict[str, Path] = {}
    for key, spec in ORGANISMS.items():
        fasta = outdir / spec["fasta"]
        go = outdir / spec["go"]
        full_fastas[key] = fasta
        full_go[key] = go
        if args.force or not fasta.exists():
            meta["organisms"].setdefault(key, {})["fasta"] = download_uniprot_fasta(
                session, spec, fasta)
        else:
            n = sum(1 for _ in SeqIO.parse(str(fasta), "fasta"))
            LOG.info("Using existing %s (%d sequences)", fasta, n)
            meta["organisms"].setdefault(key, {})["fasta"] = {"sequences": n,
                                                              "cached": True}
        if args.mode in ("all", "full") and (args.force or not go.exists()):
            meta["organisms"].setdefault(key, {})["go"] = download_quickgo_annotations(
                session, spec, go, workers=args.workers)
        elif go.exists():
            rows = read_annotations(go)
            LOG.info("Using existing %s (%d annotations)", go, len(rows))
            meta["organisms"].setdefault(key, {})["go"] = {
                "annotations": len(rows),
                "proteins": len({r["accession"] for r in rows}),
                "cached": True,
            }

    obo = outdir / "go-basic.obo"
    if args.mode in ("all", "full"):
        if args.force or not obo.exists():
            meta["go_obo"] = download_obo(session, obo)
        else:
            LOG.info("Using existing %s", obo)

    combined = outdir / "go_annotations.tsv"
    if args.mode in ("all", "full"):
        meta["combined_annotations"] = merge_annotations(full_go.values(), combined)
    elif not combined.exists():
        merge_annotations(full_go.values(), combined)

    (outdir / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")

    if args.mode in ("all", "test"):
        if not all(p.exists() for p in [*full_fastas.values(), combined, obo]):
            LOG.error("Full data missing; run --mode full first")
            return 1
        build_test_dataset(full_fastas, combined, obo, testdir,
                           args.test_size, args.seed)

    LOG.info("Done. Full data in %s, test data in %s", outdir, testdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
