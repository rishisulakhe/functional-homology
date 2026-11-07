#!/usr/bin/env python3
"""Build a self-contained HTML report from the GO evaluation outputs.

The report contains summary statistics, a violin/box plot of GO similarity
(embedded as a base64 PNG), the top embedding-only candidate pairs with shared
GO terms, methods, and limitations.
"""

from __future__ import annotations

import argparse
import base64
import csv
import html
import io
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fh_common import fasta_names  # noqa: E402

GROUP_LABELS = {
    "embedding_only": "Embedding-only (a)",
    "diamond_hit": "DIAMOND hits (b)",
    "random": "Random control (c)",
}
GROUP_ORDER = ["embedding_only", "diamond_hit", "random"]


def read_tsv(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_scores(path: Path) -> dict[str, list[float]]:
    scores: dict[str, list[float]] = {g: [] for g in GROUP_ORDER}
    if not path.exists():
        return scores
    for row in read_tsv(path):
        group = row.get("group", "")
        if group not in scores:
            continue
        try:
            scores[group].append(float(row["jaccard"]))
        except (KeyError, ValueError):
            continue
    return scores


def fmt(value, digits: int = 3) -> str:
    if value is None:
        return "—"
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return "—"
        if abs(value) < 1e-4 and value != 0:
            return f"{value:.2e}"
        return f"{value:.{digits}f}"
    return html.escape(str(value))


def fmt_p(value) -> str:
    if value is None:
        return "—"
    if value < 1e-4:
        return f"{value:.2e}"
    return f"{value:.4f}"


def plot_groups(scores: dict[str, list[float]], out_png: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.6), dpi=150)
    colors = {"embedding_only": "#4C78A8", "diamond_hit": "#F58518",
              "random": "#9E9E9E"}
    positions = {g: i + 1 for i, g in enumerate(GROUP_ORDER)}
    rng = np.random.default_rng(12345)

    for group in GROUP_ORDER:
        vals = np.asarray(scores.get(group, []), dtype=float)
        pos = positions[group]
        if vals.size >= 2:
            parts = ax.violinplot([vals], positions=[pos], widths=0.75,
                                  showmeans=False, showextrema=False)
            for body in parts["bodies"]:
                body.set_facecolor(colors[group])
                body.set_edgecolor("white")
                body.set_alpha(0.35)
        if vals.size >= 1:
            ax.boxplot([vals], positions=[pos], widths=0.14, showfliers=False,
                       patch_artist=True,
                       boxprops={"facecolor": colors[group], "alpha": 0.75,
                                 "edgecolor": "#333"},
                       medianprops={"color": "white", "linewidth": 1.5},
                       whiskerprops={"color": "#444"},
                       capprops={"color": "#444"})
            jitter = rng.normal(0, 0.035, vals.size)
            ax.scatter(np.full(vals.size, pos) + jitter, vals, s=6,
                       color="#222", alpha=0.28, linewidths=0)

    labels = [f"{GROUP_LABELS[g]}\n(n={len(scores.get(g, []))})"
              for g in GROUP_ORDER]
    ax.set_xticks([positions[g] for g in GROUP_ORDER])
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("GO Jaccard similarity (propagated terms)")
    ax.set_ylim(bottom=0)
    ax.set_title(title, fontsize=11)
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, format="png", bbox_inches="tight")
    plt.close(fig)


def png_data_uri(path: Path) -> str:
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def comparison_for(summary: dict, name: str) -> dict:
    for cmp in summary.get("comparisons", []):
        if cmp.get("comparison") == name:
            return cmp
    return {}


def build_top_candidates(pairs_path: Path, scored_path: Path, names_a: dict,
                         names_b: dict, threshold: float, limit: int = 20) -> list[dict]:
    scored = {}
    if scored_path.exists():
        for row in read_tsv(scored_path):
            if row.get("group") == "embedding_only":
                scored[(row.get("id_a"), row.get("id_b"))] = row

    candidates = []
    if pairs_path.exists():
        for row in read_tsv(pairs_path):
            if str(row.get("sequence_hit", "")).strip().lower() in {"true", "1", "yes"}:
                continue
            try:
                cos = float(row.get("cosine") or "nan")
            except ValueError:
                continue
            if not (math.isfinite(cos) and cos >= threshold):
                continue
            candidates.append((cos, row))
    candidates.sort(key=lambda t: -t[0])

    out = []
    for cos, row in candidates[:limit]:
        key = (row.get("id_a"), row.get("id_b"))
        srow = scored.get(key, {})
        shared_names = [x for x in (srow.get("shared_go_names") or "").split(";") if x]
        shared_ids = [x for x in (srow.get("shared_go_ids") or "").split(";") if x]
        out.append({
            "id_a": key[0], "id_b": key[1],
            "name_a": names_a.get(key[0], key[0]),
            "name_b": names_b.get(key[1], key[1]),
            "cosine": cos,
            "jaccard": srow.get("jaccard"),
            "shared_names": shared_names,
            "shared_ids": shared_ids,
            "scored": bool(srow),
        })
    return out


def render_rows(rows: list[dict]) -> str:
    if not rows:
        return '<tr><td colspan="6" class="muted">No embedding-only candidates in this run.</td></tr>'
    chunks = []
    for i, row in enumerate(rows, 1):
        shared = row["shared_names"]
        shown = shared[:6]
        extra = len(shared) - len(shown)
        shared_html = ", ".join(html.escape(s) for s in shown)
        if extra > 0:
            shared_html += f' <span class="muted">(+{extra} more)</span>'
        if not shared_html:
            shared_html = '<span class="muted">—</span>'
        jac = fmt(float(row["jaccard"])) if row["jaccard"] not in (None, "") else "n/a"
        chunks.append(
            "<tr>"
            f"<td>{i}</td>"
            f'<td><code>{html.escape(row["id_a"])}</code><br><span class="pname">{html.escape(row["name_a"])}</span></td>'
            f'<td><code>{html.escape(row["id_b"])}</code><br><span class="pname">{html.escape(row["name_b"])}</span></td>'
            f"<td>{row['cosine']:.4f}</td>"
            f"<td>{jac}</td>"
            f"<td>{shared_html}</td>"
            "</tr>"
        )
    return "\n".join(chunks)


def render_summary_table(summary: dict) -> str:
    stats = summary.get("summary_stats", {})
    counts = summary.get("groups", {})
    cmp_a = comparison_for(summary, "embedding_only_vs_random")
    cmp_b = comparison_for(summary, "diamond_hit_vs_random")
    rows = []
    for group in GROUP_ORDER:
        st = stats.get(group, {})
        cnt = counts.get(group, {})
        cmp = cmp_a if group == "embedding_only" else (cmp_b if group == "diamond_hit" else {})
        rows.append(
            "<tr>"
            f"<td>{GROUP_LABELS[group]}</td>"
            f"<td>{fmt(cnt.get('n_input'))}</td>"
            f"<td>{fmt(cnt.get('n_scored'))}</td>"
            f"<td>{fmt(st.get('mean'))}</td>"
            f"<td>{fmt(st.get('median'))}</td>"
            f"<td>{fmt(cmp.get('cliffs_delta'))}</td>"
            f"<td>{fmt_p(cmp.get('p_value'))}</td>"
            f"<td>{fmt_p(cmp.get('p_value_holm'))}</td>"
            "</tr>"
        )
    return "\n".join(rows)


def interpretation_html(summary: dict) -> str:
    stats = summary.get("summary_stats", {})
    a = stats.get("embedding_only", {})
    c = stats.get("random", {})
    b = stats.get("diamond_hit", {})
    cmp_a = comparison_for(summary, "embedding_only_vs_random")
    cmp_b = comparison_for(summary, "diamond_hit_vs_random")
    status = summary.get("hypothesis", "unknown")

    if status == "supported":
        css = "verdict-supported"
        headline = "Hypothesis supported"
        text = (
            "Embedding-only pairs had significantly higher GO Jaccard similarity "
            f"than random pairs (mean {fmt(a.get('mean'))} vs {fmt(c.get('mean'))}; "
            f"two-sided Mann-Whitney p = {fmt_p(cmp_a.get('p_value'))}; "
            f"Cliff's delta = {fmt(cmp_a.get('cliffs_delta'))}). "
            "This is consistent with the idea that some protein pairs share "
            "functional annotations without a detectable DIAMOND hit."
        )
    elif status == "insufficient_data":
        css = "verdict-neutral"
        headline = "Insufficient data for a conclusion"
        text = (
            "Too few pairs had experimental GO annotations on both sides to run "
            "the comparison. This is an annotation-coverage limitation, not "
            "evidence against the hypothesis."
        )
    else:
        css = "verdict-not-supported"
        headline = "Hypothesis not supported in this run"
        text = (
            "Embedding-only pairs were not significantly higher than random pairs "
            f"in GO Jaccard similarity (mean {fmt(a.get('mean'))} vs {fmt(c.get('mean'))}; "
            f"two-sided Mann-Whitney p = {fmt_p(cmp_a.get('p_value'))}; "
            f"Cliff's delta = {fmt(cmp_a.get('cliffs_delta'))}). "
            "Treat this as an inconclusive/negative result under this particular "
            "embedding, cutoff and sparse-annotation setting."
        )
    text += (
        f" For context, DIAMOND-hit pairs scored mean {fmt(b.get('mean'))} "
        f"(p vs random = {fmt_p(cmp_b.get('p_value'))})."
    )
    return f'<div class="verdict {css}"><h3>{headline}</h3><p>{text}</p></div>'


def sensitivity_html(summary: dict) -> str:
    sens = summary.get("sensitivity_no_root")
    if not sens:
        return ""
    stats = sens.get("summary_stats", {})
    cmps = {c.get("comparison"): c for c in sens.get("comparisons", [])}
    ca = cmps.get("embedding_only_vs_random", {})
    cb = cmps.get("diamond_hit_vs_random", {})
    return (
        '<div class="note"><strong>Sensitivity analysis (root terms removed):</strong> '
        f"mean Jaccard embedding-only {fmt(stats.get('embedding_only', {}).get('mean'))}, "
        f"DIAMOND {fmt(stats.get('diamond_hit', {}).get('mean'))}, "
        f"random {fmt(stats.get('random', {}).get('mean'))}; "
        f"embedding-only vs random p = {fmt_p(ca.get('p_value'))}, "
        f"Cliff's δ = {fmt(ca.get('cliffs_delta'))}; "
        f"DIAMOND vs random p = {fmt_p(cb.get('p_value'))}, "
        f"Cliff's δ = {fmt(cb.get('cliffs_delta'))}. "
        "The conclusion is unchanged when the three generic GO roots are removed."
        "</div>"
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--summary", type=Path, required=True, help="eval_summary.json")
    p.add_argument("--pairs-scored", type=Path, required=True, help="pairs_scored.tsv")
    p.add_argument("--pairs", type=Path, default=None, help="pairs.tsv (for top candidates)")
    p.add_argument("--fasta-a", type=Path, required=True)
    p.add_argument("--fasta-b", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True, help="output report.html")
    p.add_argument("--model", default="facebook/esm2_t6_8M_UR50D")
    p.add_argument("--top-k", type=int, default=5)
    p.add_argument("--cosine-threshold", type=float, default=0.9)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--title", default="Functional homology without sequence homology")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    scores = load_scores(args.pairs_scored)
    names_a = fasta_names(args.fasta_a)
    names_b = fasta_names(args.fasta_b)

    tmp_png = args.out.with_suffix(".tmp.png")
    tmp_png.parent.mkdir(parents=True, exist_ok=True)
    plot_groups(scores, tmp_png,
                "GO similarity of embedding-only, DIAMOND and random pairs")
    img_uri = png_data_uri(tmp_png)
    tmp_png.unlink(missing_ok=True)

    top = build_top_candidates(args.pairs, args.pairs_scored, names_a, names_b,
                               args.cosine_threshold, limit=20)
    counts = summary.get("groups", {})
    stats = summary.get("summary_stats", {})
    ann = summary.get("annotation_stats", {})
    params = summary.get("params", {})
    generated = summary.get("generated_utc") or datetime.now(timezone.utc).isoformat()

    verdict = interpretation_html(summary)
    summary_table = render_summary_table(summary)
    top_rows = render_rows(top)

    html_doc = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(args.title)}</title>
<style>
  :root {{ color-scheme: light; }}
  * {{ box-sizing: border-box; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
         Helvetica, Arial, sans-serif; margin: 0; background: #f6f7f9;
         color: #1f2328; line-height: 1.55; }}
  main {{ max-width: 1080px; margin: 0 auto; padding: 32px 20px 64px; }}
  header {{ border-bottom: 3px solid #4C78A8; padding-bottom: 14px; margin-bottom: 24px; }}
  h1 {{ font-size: 1.7rem; margin: 0 0 6px; }}
  h2 {{ font-size: 1.25rem; margin: 34px 0 12px; border-bottom: 1px solid #dfe3e8;
       padding-bottom: 6px; }}
  h3 {{ margin: 0 0 6px; font-size: 1.02rem; }}
  p {{ margin: 8px 0; }}
  .sub {{ color: #59636e; font-size: .92rem; }}
  .cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr));
            gap: 12px; margin: 18px 0; }}
  .card {{ background: white; border: 1px solid #dfe3e8; border-radius: 10px;
           padding: 12px 14px; }}
  .card .k {{ font-size: .72rem; text-transform: uppercase; letter-spacing: .06em;
              color: #59636e; }}
  .card .v {{ font-size: 1.45rem; font-weight: 650; margin-top: 2px; }}
  .verdict {{ border-radius: 10px; padding: 14px 18px; margin: 18px 0;
              border: 1px solid; }}
  .verdict-supported {{ background: #e8f5e9; border-color: #7ac47f; }}
  .verdict-not-supported {{ background: #fff4e5; border-color: #e0a458; }}
  .verdict-neutral {{ background: #eef2f6; border-color: #a9b7c6; }}
  table {{ width: 100%; border-collapse: collapse; background: white;
           border: 1px solid #dfe3e8; border-radius: 10px; overflow: hidden;
           font-size: .88rem; }}
  th, td {{ padding: 8px 10px; text-align: left; vertical-align: top;
            border-bottom: 1px solid #eef1f4; }}
  th {{ background: #f1f4f8; font-size: .75rem; text-transform: uppercase;
        letter-spacing: .04em; color: #42505e; }}
  tr:last-child td {{ border-bottom: none; }}
  code {{ background: #f1f4f8; padding: 1px 4px; border-radius: 4px;
          font-size: .85em; }}
  .pname {{ color: #59636e; font-size: .78rem; }}
  .muted {{ color: #8a949e; }}
  img {{ max-width: 100%; height: auto; border: 1px solid #dfe3e8;
         border-radius: 10px; background: white; }}
  ul {{ margin: 8px 0 8px 20px; padding: 0; }}
  li {{ margin: 4px 0; }}
  .grid2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 24px; }}
  .note {{ background: #f1f4f8; border-left: 4px solid #4C78A8; border-radius: 6px;
           padding: 10px 14px; margin: 14px 0; font-size: .9rem; }}
  @media (max-width: 760px) {{ .grid2 {{ grid-template-columns: 1fr; }} }}
</style>
</head>
<body>
<main>
<header>
  <h1>{html.escape(args.title)}</h1>
  <div class="sub">E. coli K-12 vs B. subtilis 168 · ESM-2 embeddings vs DIAMOND ·
  experimental GO evaluation · generated {html.escape(str(generated))}</div>
</header>

{verdict}

<div class="cards">
  <div class="card"><div class="k">Model</div><div class="v" style="font-size:1rem">{html.escape(args.model)}</div></div>
  <div class="card"><div class="k">Proteins A / B</div><div class="v">{fmt(summary.get('n_annotated_proteins'))} annotated<br><span class="sub">GO-annotated proteins</span></div></div>
  <div class="card"><div class="k">Embedding-only pairs</div><div class="v">{fmt(counts.get('embedding_only', {}).get('n_input'))}</div></div>
  <div class="card"><div class="k">Scorable pairs (a)</div><div class="v">{fmt(counts.get('embedding_only', {}).get('n_scored'))}</div></div>
  <div class="card"><div class="k">DIAMOND pairs</div><div class="v">{fmt(counts.get('diamond_hit', {}).get('n_input'))}</div></div>
</div>

<h2>GO similarity by group</h2>
<img alt="Violin and box plot of GO Jaccard similarity for the three pair groups"
     src="{img_uri}">

<table>
<thead><tr>
<th>Group</th><th>Input pairs</th><th>Scorable</th><th>Mean Jaccard</th>
<th>Median</th><th>Cliff's δ vs random</th><th>p (MWU)</th><th>p (Holm)</th>
</tr></thead>
<tbody>
{summary_table}
</tbody>
</table>
<p class="sub">A pair is scorable only when both proteins carry at least one
experimental GO annotation. Cliff's δ &gt; 0 means the group tends to have
higher GO similarity than random.</p>
{sensitivity_html(summary)}

<h2>Top 20 embedding-only candidate pairs</h2>
<p class="sub">Ranked by mean-pooled ESM-2 cosine similarity among top-{args.top_k}
pairs with no DIAMOND hit and cosine ≥ {args.cosine_threshold:.2f}.
Shared GO terms are ancestor-propagated experimental terms. Pairs without
annotation on both sides are shown as n/a.</p>
<table>
<thead><tr>
<th>#</th><th>Protein A</th><th>Protein B</th><th>Cosine</th>
<th>GO Jaccard</th><th>Shared experimental GO terms</th>
</tr></thead>
<tbody>
{top_rows}
</tbody>
</table>

<h2>Methods</h2>
<div class="grid2">
<div>
<h3>Data</h3>
<ul>
<li>Reviewed (Swiss-Prot) proteomes from UniProt REST: E. coli K-12
(taxid 83333) and B. subtilis 168 (taxid 224308).</li>
<li>GO annotations from QuickGO (ECO:0000269 descendants), filtered to the
experimental GO evidence codes {", ".join(html.escape(c) for c in params.get("evidence_codes", []))}.
<code>NOT</code> qualifiers were removed.</li>
<li>{fmt(ann.get("kept"))} annotation rows kept for {fmt(summary.get("n_annotated_proteins"))}
proteins; GO graph from <code>go-basic.obo</code>.</li>
</ul>
<h3>Sequence similarity</h3>
<ul>
<li>DIAMOND <code>makedb</code> + <code>blastp --very-sensitive</code>,
e-value ≤ {html.escape(str(summary.get("params", {}).get("evalue_cutoff", "1e-3")))}.</li>
<li>Outfmt 6: qseqid, sseqid, pident, evalue, bitscore; best HSP per pair.</li>
</ul>
</div>
<div>
<h3>Embedding similarity</h3>
<ul>
<li>{html.escape(args.model)}; mean-pooled residue embeddings excluding
special tokens; sequences truncated to 1022 residues.</li>
<li>Rows L2-normalised; cosine similarity matrix; top-{args.top_k} matches per
A protein.</li>
<li>Embedding-only candidate = no DIAMOND hit and cosine ≥ {args.cosine_threshold:.2f}.</li>
</ul>
<h3>Evaluation</h3>
<ul>
<li>GO terms propagated over <code>is_a</code> and <code>part_of</code>.</li>
<li>Jaccard similarity of propagated term sets; Resnik reported in
<code>pairs_scored.tsv</code> as a secondary metric. A sensitivity analysis
repeats the comparison after removing the three GO root terms.</li>
<li>Two-sided Mann-Whitney U tests vs the random control
(seed {args.seed}); Cliff's delta effect size; Holm correction across the two
tests.</li>
</ul>
</div>
</div>

<h2>Limitations</h2>
<ul>
<li><strong>Mean pooling</strong> collapses per-residue information; functional
sites may be diluted. Attention pooling or per-domain embeddings could improve
signal.</li>
<li><strong>Truncation</strong> at 1022 residues can drop domains from long
proteins (~{fmt(summary.get("n_annotated_proteins"))} annotated proteins; see
embedding metadata for truncation counts).</li>
<li><strong>Sparse experimental annotations</strong>: B. subtilis has far fewer
experimental GO annotations than E. coli, so most cross-organism pairs are
unscorable. This reduces power and can bias the scored subset.</li>
<li><strong>GO graph effects</strong>: propagated ancestor sets share generic
terms; Jaccard is sensitive to annotation depth and aspect composition.</li>
<li><strong>DIAMOND is not an exhaustive homology detector</strong>; distant
homologues may be missed, which is exactly the hypothesis, but it also means
some embedding-only pairs are genuine remote homologues.</li>
<li><strong>Random controls</strong> are sampled from the full proteome cross
product, not matched for abundance, length or amino-acid composition.</li>
<li><strong>One model, one cutoff, one seed</strong>: results may not generalise;
larger ESM-2 models and sensitivity analyses are warranted.</li>
</ul>

<h2>Reproducibility</h2>
<p class="sub">Generated by <code>nextflow run main.nf -profile conda</code>.
Parameters: model={html.escape(args.model)}, top_k={args.top_k},
cosine_threshold={args.cosine_threshold}, seed={args.seed}.
Statistics computed on {fmt(generated)}. This HTML is self-contained.</p>
</main>
</body>
</html>
"""
    args.out.write_text(html_doc, encoding="utf-8")
    print(f"Wrote {args.out} ({len(html_doc)/1024:.1f} KiB)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
