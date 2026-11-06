#!/usr/bin/env python3
"""Evaluate candidate protein pairs with experimental GO annotation similarity.

Groups
------
* ``embedding_only`` : top-k embedding pairs with no DIAMOND hit (cosine
  >= threshold).
* ``diamond_hit``    : all DIAMOND hits passing the e-value cutoff.
* ``random``         : random cross-organism control pairs.

Each protein's experimentally annotated GO terms are propagated to all
``is_a`` and ``part_of`` ancestors.  A pair is scored with Jaccard similarity
of the propagated term sets (and optionally Resnik similarity using
information content estimated from the annotated proteomes).  A pair is only
scored when both proteins have at least one experimental annotation.

Statistics: Mann-Whitney U tests (two-sided) for embedding-only vs random and
DIAMOND vs random, Cliff's delta effect size, and Holm-adjusted p-values.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import math
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fh_common import ancestor_closure, obo_names, parse_obo  # noqa: E402

LOG = logging.getLogger("go_eval")

EXPERIMENTAL_CODES = {
    "EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
    "HTP", "HDA", "HMP", "HGI", "HEP",
}

SCORED_COLUMNS = [
    "group", "id_a", "id_b", "cosine", "pident", "evalue",
    "jaccard", "resnik", "n_terms_a", "n_terms_b",
    "shared_go_ids", "shared_go_names",
]


def jaccard(set_a: set[str], set_b: set[str]) -> float:
    """Jaccard similarity; 0.0 when both sets are empty."""
    union = set_a | set_b
    if not union:
        return 0.0
    return len(set_a & set_b) / len(union)


def read_annotations(path: str | Path,
                     allowed_codes: set[str] | None = None) -> tuple[dict[str, set[str]], dict]:
    """Read GO annotations into accession -> set(GO ids)."""
    allowed = allowed_codes or EXPERIMENTAL_CODES
    ann: dict[str, set[str]] = {}
    stats = {"rows": 0, "kept": 0, "bad_evidence": 0, "not_qualifier": 0}
    with open(path, encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames or "accession" not in reader.fieldnames \
                or "go_id" not in reader.fieldnames:
            raise ValueError(f"{path} must contain 'accession' and 'go_id' columns; "
                             f"found {reader.fieldnames}")
        for row in reader:
            stats["rows"] += 1
            code = (row.get("evidence_code") or "").strip()
            if code and code not in allowed:
                stats["bad_evidence"] += 1
                continue
            qual = (row.get("qualifier") or "").upper()
            if "NOT" in qual.split("|"):
                stats["not_qualifier"] += 1
                continue
            acc = (row.get("accession") or "").strip()
            go_id = (row.get("go_id") or "").strip()
            if not acc or not go_id.startswith("GO:"):
                continue
            ann.setdefault(acc, set()).add(go_id)
            stats["kept"] += 1
    return ann, stats


def propagate_annotations(ann: dict[str, set[str]],
                          parents: dict[str, set[str]]) -> tuple[dict[str, set[str]], int]:
    """Propagate every protein's GO terms to all is_a/part_of ancestors."""
    propagated: dict[str, set[str]] = {}
    missing_terms: set[str] = set()
    for acc, terms in ann.items():
        missing_terms |= {t for t in terms if t not in parents}
        propagated[acc] = ancestor_closure(terms, parents)
    return propagated, len(missing_terms)


def build_ic(propagated: dict[str, set[str]]) -> dict[str, float]:
    """Information content per GO term from the annotated protein corpus."""
    n = len(propagated)
    if n == 0:
        return {}
    counts: Counter[str] = Counter()
    for terms in propagated.values():
        counts.update(terms)
    return {term: -math.log(count / n) for term, count in counts.items()}


def resnik_similarity(terms_a: set[str], terms_b: set[str],
                      ic: dict[str, float]) -> float:
    """Resnik similarity = IC of the most informative common ancestor."""
    common = terms_a & terms_b
    if not common:
        return 0.0
    return max(ic.get(term, 0.0) for term in common)


def read_pairs(path: str | Path, group: str, cosine_threshold: float | None = None,
               only_sequence_hits: bool | None = None,
               has_sequence_flag: bool = True) -> list[dict]:
    """Read a pairs TSV into a list of dicts with a fixed ``group`` label."""
    pairs: list[dict] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"id_a", "id_b"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain columns {sorted(required)}")
        for row in reader:
            row = dict(row)
            row["group"] = group
            if has_sequence_flag and only_sequence_hits is not None:
                flag = str(row.get("sequence_hit", "")).strip().lower() in {"true", "1", "yes"}
                if flag != only_sequence_hits:
                    continue
            if cosine_threshold is not None:
                try:
                    cos = float(row.get("cosine") or "nan")
                except ValueError:
                    cos = float("nan")
                if not (math.isfinite(cos) and cos >= cosine_threshold):
                    continue
            pairs.append(row)
    return pairs


def describe(values: list[float]) -> dict:
    """Summary statistics for a list of scores."""
    if not values:
        return {"n": 0, "mean": None, "median": None, "sd": None,
                "q1": None, "q3": None, "min": None, "max": None}
    import numpy as np
    arr = np.asarray(values, dtype=float)
    return {
        "n": int(arr.size),
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "sd": float(np.std(arr, ddof=1)) if arr.size > 1 else 0.0,
        "q1": float(np.percentile(arr, 25)),
        "q3": float(np.percentile(arr, 75)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


def mann_whitney(x: list[float], y: list[float], label: str) -> dict:
    """Two-sided Mann-Whitney U plus Cliff's delta (x relative to y)."""
    if len(x) == 0 or len(y) == 0:
        return {"comparison": label, "n_x": len(x), "n_y": len(y),
                "u": None, "p_value": None, "p_value_holm": None,
                "cliffs_delta": None, "direction": None,
                "note": "insufficient data (empty group)"}
    res = mannwhitneyu(x, y, alternative="two-sided", method="auto")
    n1, n2 = len(x), len(y)
    delta = 2.0 * float(res.statistic) / (n1 * n2) - 1.0
    direction = "higher" if delta > 0 else ("lower" if delta < 0 else "equal")
    return {
        "comparison": label,
        "n_x": n1,
        "n_y": n2,
        "u": float(res.statistic),
        "p_value": float(res.pvalue),
        "cliffs_delta": delta,
        "direction": direction,
        "note": "Cliff's delta > 0 means x tends to be larger than y",
    }


def holm_adjust(pvalues: list[float | None]) -> list[float | None]:
    """Holm-Bonferroni adjustment, preserving None entries."""
    indexed = [(i, p) for i, p in enumerate(pvalues) if p is not None]
    m = len(indexed)
    adjusted: list[float | None] = [None] * len(pvalues)
    running = 0.0
    for rank, (idx, p) in enumerate(sorted(indexed, key=lambda t: t[1])):
        value = min(1.0, (m - rank) * p)
        running = max(running, value)
        adjusted[idx] = running
    return adjusted


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pairs", type=Path, required=True,
                   help="pairs.tsv from find_pairs.py")
    p.add_argument("--random", type=Path, required=True,
                   help="random_pairs.tsv from find_pairs.py")
    p.add_argument("--diamond", type=Path, default=None,
                   help="diamond_pairs.tsv; if omitted use sequence_hit rows in --pairs")
    p.add_argument("--go-ann", type=Path, required=True,
                   help="merged experimental GO annotations TSV")
    p.add_argument("--go-obo", type=Path, required=True, help="go-basic.obo")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--cosine-threshold", type=float, default=0.9)
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(message)s",
                        datefmt="%H:%M:%S")
    args.outdir.mkdir(parents=True, exist_ok=True)

    blocks, parents = parse_obo(args.go_obo)
    names = obo_names(blocks)
    LOG.info("Loaded GO graph: %d terms, %d with at least one parent",
             len(parents), sum(1 for v in parents.values() if v))

    ann, ann_stats = read_annotations(args.go_ann)
    propagated, n_missing = propagate_annotations(ann, parents)
    ic = build_ic(propagated)
    LOG.info("Annotations: %d rows kept, %d proteins, %d with propagated terms; "
             "%d annotated terms missing from OBO",
             ann_stats["kept"], len(ann), len(propagated), n_missing)
    if ann_stats["bad_evidence"]:
        LOG.warning("Ignored %d rows with non-experimental evidence codes",
                    ann_stats["bad_evidence"])

    groups: dict[str, list[dict]] = {
        "embedding_only": read_pairs(args.pairs, "embedding_only",
                                     cosine_threshold=args.cosine_threshold,
                                     only_sequence_hits=False),
        "random": read_pairs(args.random, "random", has_sequence_flag=False),
    }
    if args.diamond and args.diamond.exists():
        groups["diamond_hit"] = read_pairs(args.diamond, "diamond_hit",
                                           has_sequence_flag=False)
        LOG.info("Using complete DIAMOND hit set for group (b)")
    else:
        groups["diamond_hit"] = read_pairs(args.pairs, "diamond_hit",
                                           only_sequence_hits=True)
        LOG.info("Using sequence_hit rows from pairs.tsv for group (b)")

    scored_rows: list[list] = []
    scored_by_group: dict[str, list[float]] = {g: [] for g in groups}
    counts: dict[str, dict] = {}
    for group, pairs in groups.items():
        n_scored = 0
        n_missing_ann = 0
        for row in pairs:
            a, b = row.get("id_a", ""), row.get("id_b", "")
            set_a = propagated.get(a)
            set_b = propagated.get(b)
            cos = row.get("cosine") or ""
            pident = row.get("pident") or ""
            evalue = row.get("evalue") or ""
            if not set_a or not set_b:
                n_missing_ann += 1
                continue
            jac = jaccard(set_a, set_b)
            res = resnik_similarity(set_a, set_b, ic)
            shared = sorted(set_a & set_b)
            scored_rows.append([
                group, a, b, cos, pident, evalue,
                f"{jac:.6f}", f"{res:.6f}", len(set_a), len(set_b),
                ";".join(shared),
                ";".join(names.get(t, t) for t in shared),
            ])
            scored_by_group[group].append(jac)
            n_scored += 1
        counts[group] = {
            "n_input": len(pairs),
            "n_scored": n_scored,
            "n_unscorable_missing_annotation": n_missing_ann,
            "n_scored_with_shared_terms": int(sum(1 for row in scored_rows
                                                  if row[0] == group and row[10])),
        }
        LOG.info("Group %-15s input=%-6d scored=%-6d unscorable=%-6d",
                 group, len(pairs), n_scored, n_missing_ann)

    with (args.outdir / "pairs_scored.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(SCORED_COLUMNS)
        writer.writerows(scored_rows)

    a = scored_by_group["embedding_only"]
    b = scored_by_group["diamond_hit"]
    c = scored_by_group["random"]
    stats_a = describe(a)
    stats_b = describe(b)
    stats_c = describe(c)
    comparisons = [
        mann_whitney(a, c, "embedding_only_vs_random"),
        mann_whitney(b, c, "diamond_hit_vs_random"),
    ]
    adjusted = holm_adjust([cmp["p_value"] for cmp in comparisons])
    for cmp, p_adj in zip(comparisons, adjusted):
        cmp["p_value_holm"] = p_adj

    def supported(cmp: dict) -> bool | str:
        if cmp.get("p_value") is None:
            return "insufficient_data"
        return bool(cmp["p_value"] < 0.05 and (cmp["cliffs_delta"] or 0) > 0)

    cmp_a = comparisons[0]
    if cmp_a.get("p_value") is None:
        hypothesis = "insufficient_data"
    elif supported(cmp_a):
        hypothesis = "supported"
    else:
        hypothesis = "not_supported"

    summary = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "params": {
            "cosine_threshold": args.cosine_threshold,
            "evidence_codes": sorted(EXPERIMENTAL_CODES),
            "distance": "jaccard on ancestor-propagated term sets (is_a + part_of)",
            "test": "two-sided Mann-Whitney U; Cliff's delta effect size",
        },
        "go_graph": {"n_terms": len(parents), "n_annotated_terms_missing_from_obo": n_missing},
        "annotation_stats": ann_stats,
        "n_annotated_proteins": len(propagated),
        "groups": counts,
        "summary_stats": {
            "embedding_only": stats_a,
            "diamond_hit": stats_b,
            "random": stats_c,
        },
        "comparisons": comparisons,
        "hypothesis": hypothesis,
        "hypothesis_statement": (
            "Group (a) embedding-only pairs have higher GO Jaccard similarity "
            "than group (c) random pairs (two-sided Mann-Whitney U p < 0.05 "
            "and Cliff's delta > 0)."),
        "hypothesis_supported": supported(cmp_a),
    }
    (args.outdir / "eval_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    LOG.info("Wrote %s and %s", args.outdir / "pairs_scored.tsv",
             args.outdir / "eval_summary.json")
    for cmp in comparisons:
        LOG.info("%s: n=%d vs %d, p=%.4g (Holm %.4g), Cliff's delta=%.3f (%s)",
                 cmp["comparison"], cmp["n_x"], cmp["n_y"],
                 cmp["p_value"] if cmp["p_value"] is not None else float("nan"),
                 cmp["p_value_holm"] if cmp["p_value_holm"] is not None else float("nan"),
                 cmp["cliffs_delta"] if cmp["cliffs_delta"] is not None else float("nan"),
                 cmp["direction"])
    LOG.info("Hypothesis: %s", hypothesis)
    return 0


if __name__ == "__main__":
    sys.exit(main())
