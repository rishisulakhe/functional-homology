# functionalhomology

Do two proteins from different organisms ever share a function without
detectable sequence homology? This pipeline compares the proteomes of
**Escherichia coli K-12** and **Bacillus subtilis 168** using two independent
views of protein similarity:

1. **Sequence similarity** — DIAMOND `blastp --very-sensitive` (the baseline).
2. **Embedding similarity** — mean-pooled ESM-2 protein language-model
   embeddings and cosine similarity.

Candidate pairs are embedding-similar pairs for which DIAMOND reports **no**
hit. Those candidates are then evaluated with **experimental Gene Ontology
(GO) annotations** and compared against DIAMOND-hit pairs and a size-matched
random-pair negative control.

> **Result in one line:** on the full proteomes, embedding-only pairs have
> significantly higher experimental GO similarity than random pairs
> (mean Jaccard 0.125 vs 0.081, Mann-Whitney p = 6.9×10⁻³⁶, Cliff's δ = 0.235),
> so the hypothesis is **supported at the group level**, though the effect is
> modest and DIAMOND-hit pairs remain the strongest group. See
> [Results](#results-full-run) for the full and honest interpretation.

---

## Pipeline at a glance

```
UniProt FASTA A ──► EMBED (ESM-2) ──► embeddings A ─┐
UniProt FASTA B ──► EMBED (ESM-2) ──► embeddings B ─┤
                                                     ├─► FIND_PAIRS ─► GO_EVAL ─► REPORT
UniProt FASTA A+B ► DIAMOND makedb + blastp ─────────┘
                                │
                                └─► diamond_hits.tsv
```

`DIAMOND` and both `EMBED` tasks run in parallel; `FIND_PAIRS` joins top-k
cosine matches with DIAMOND hits and samples the random control; `GO_EVAL`
propagates GO terms and runs the statistics; `REPORT` writes a self-contained
HTML report.

---

## Requirements

- Linux or macOS (tested on Arch Linux)
- **Java 17–26** for Nextflow (Nextflow 26.10.5 was used; system Java 27 is too
  new for this Nextflow release)
- **Nextflow ≥ 23** (tested with 26.10.5)
- **Conda / Miniforge** for the `conda` profile. Nextflow's conda activation
  requires a real `conda` executable; **micromamba alone is not sufficient**.
- ~4 GB free disk for the full data and embedding outputs
- Internet access on first run (UniProt, QuickGO, GO, and Hugging Face model
  download)

Python 3.11, DIAMOND, and all Python packages are pinned in
[`environment.yml`](environment.yml). Nextflow creates that environment
automatically with `-profile conda`. You do not need to install Python packages
on the host.

### Install Nextflow (example)

```bash
mkdir -p "$HOME/.local/bin"
cd /tmp
NXF_VER=26.10.5 bash -c 'curl -s https://get.nextflow.io | bash'
mv nextflow "$HOME/.local/bin/"
export PATH="$HOME/.local/bin:$PATH"
nextflow -version
```

### Java and conda

With Miniforge/Miniconda already installed:

```bash
conda create -n fh-java -c conda-forge openjdk=17 -y
conda activate fh-java
```

If you prefer to use an already-created project environment instead of letting
Nextflow create one:

```bash
conda env create -f environment.yml
conda activate functionalhomology
nextflow run main.nf -profile test,local
```

### Docker (optional)

```bash
docker build -t functionalhomology:latest .
nextflow run main.nf -profile test,docker
```

The Dockerfile builds the same pinned environment with micromamba. Both the
`conda` and `docker` test profiles have been run end to end; the conda profile is
the primary path.

---

## Quick start: tiny test dataset

The repository ships a committed, seeded test dataset in `testdata/`: 200
proteins per organism (seed 42), the matching experimental GO annotation
subset, and an OBO slice closed under `is_a`/`part_of` ancestors. No download
is needed for the test profile.

```bash
git clone <this-repo> functionalhomology
cd functionalhomology
nextflow run main.nf -profile test,conda
```

Expected outputs:

```
results/test/report.html
results/test/go_eval/eval_summary.json
results/test/go_eval/pairs_scored.tsv
...
```

The test profile is intentionally small (only ~200 scorable pairs), so its
statistics are not meaningful; it exists to verify every step end to end.

---

## Full analysis

### 1. Download data

The download script needs the Python packages from `environment.yml`. Create
the named environment once (it is also handy for `pytest`):

```bash
conda env create -f environment.yml
conda activate functionalhomology
python scripts/download_data.py --mode full
```

If you do not want a named environment, you can instead let Nextflow create its
own conda environment for the pipeline and run the download script with any
Python 3.11 that has `requests` and `biopython` installed. The download itself
only uses those two packages.

This downloads, into `data/` (git-ignored):

| File | Contents |
|---|---|
| `proteome_a.fasta` | Reviewed (Swiss-Prot) *E. coli* K-12, taxid 83333 |
| `proteome_b.fasta` | Reviewed (Swiss-Prot) *B. subtilis* 168, taxid 224308 |
| `go_annotations_ecoli.tsv` | Experimental GO annotations, *E. coli* |
| `go_annotations_bsubtilis.tsv` | Experimental GO annotations, *B. subtilis* |
| `go_annotations.tsv` | Merged annotation table used by the pipeline |
| `go-basic.obo` | GO ontology |
| `metadata.json` | URLs, UniProt release, counts, timestamp |

`--mode all` also regenerates `testdata/` from the full data. The script is
idempotent; use `--force` to re-download.

**Annotation source and evidence filter.** Annotations come from the
[QuickGO annotation search API](https://www.ebi.ac.uk/QuickGO/services/annotation/search)
with `taxonId` and `evidenceCode=ECO:0000269` (the ECO ancestor of experimental
evidence). The JSON response is filtered to the GO experimental evidence codes:

```
EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP
```

`NOT` qualifiers are removed. The API's bulk GAF download is capped at 10,000
annotations, so the script pages through the JSON search endpoint instead.

We use **experimental annotations only** to avoid circular evaluation:
IEA/ISS/ISO and many other electronic annotations are themselves inferred from
sequence similarity, so scoring sequence-based predictions against them would
partly ask whether DIAMOND agrees with DIAMOND. The high-throughput
experimental codes (HTP/HDA/HMP/HGI/HEP) are included because they are
experimental; this is documented as a possible source of noise.

### 2. Run

With the named environment active (recommended, avoids creating a second env):

```bash
nextflow run main.nf -profile full,local
```

Or let Nextflow create its own conda environment from `environment.yml`:

```bash
nextflow run main.nf -profile full,conda
```

Equivalent explicit command (works with either environment):

```bash
nextflow run main.nf -profile conda \
  --proteome_a data/proteome_a.fasta \
  --proteome_b data/proteome_b.fasta \
  --go_ann data/go_annotations.tsv \
  --go_obo data/go-basic.obo \
  --outdir results/full
```

Use `-resume` to reuse completed tasks after a failure.

### 3. Larger ESM-2 model

```bash
nextflow run main.nf -profile full,conda \
  --esm_model facebook/esm2_t30_150M_UR50D \
  --outdir results/full_150m
```

Larger models need more RAM and CPU time. The default
`facebook/esm2_t6_8M_UR50D` is CPU-friendly; the embeddings are 320-dimensional.
The model is downloaded once into the Hugging Face cache.

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--proteome_a` / `--proteome_b` | — | Input FASTA files (required unless a profile supplies them) |
| `--go_ann` | — | Merged experimental annotation TSV |
| `--go_obo` | — | `go-basic.obo` or a closed slice |
| `--outdir` | `results` | Output directory |
| `--esm_model` | `facebook/esm2_t6_8M_UR50D` | Hugging Face ESM-2 model id |
| `--embed_batch_size` | `8` | Proteins per forward pass |
| `--max_residues` | `1022` | Truncation length per sequence |
| `--top_k` | `5` | Embedding neighbours retained per protein in A |
| `--cosine_threshold` | `0.9` | Minimum cosine for embedding-only candidates |
| `--seed` | `42` | Seed for the random control and test subsampling |
| `--evalue_cutoff` | `1e-3` | DIAMOND significance cutoff |

Profiles: `test`, `full`, `conda`, `docker`, `local`.

---

## Outputs

| Path | Description |
|---|---|
| `results/<out>/diamond/diamond_hits.tsv` | DIAMOND outfmt 6: qseqid, sseqid, pident, evalue, bitscore |
| `results/<out>/embeddings/proteome_{a,b}.npy` | Mean-pooled embeddings, one row per protein |
| `.../proteome_{a,b}.ids.txt` | Accessions in embedding row order |
| `.../proteome_{a,b}.meta.json` | Model, dimension, truncation count/ids |
| `results/<out>/pairs/pairs.tsv` | Top-k embedding pairs joined with DIAMOND hits |
| `.../diamond_pairs.tsv` | Best HSP per DIAMOND pair |
| `.../random_pairs.tsv` | Size-matched random control |
| `.../pairs_summary.json` | Pair counts and parameters |
| `results/<out>/go_eval/pairs_scored.tsv` | Scored pairs with Jaccard, no-root Jaccard, Resnik, shared terms |
| `.../eval_summary.json` | Group statistics, tests, effect sizes, sensitivity analysis |
| `results/<out>/report.html` | Self-contained HTML report (plot embedded as base64) |
| `docs/full_run_report.html` | Archived report from the full run documented here |
| `results/<out>/pipeline_info/` | Nextflow trace and report |

---

## Evaluation design

### Pair groups

- **(a) Embedding-only:** top-k ESM-2 cosine neighbours with **no** DIAMOND hit
  passing e-value 1e-3, and cosine ≥ 0.9.
- **(b) DIAMOND hits:** all DIAMOND pairs passing e-value 1e-3 (not restricted
  to the embedding top-k).
- **(c) Random control:** uniformly sampled cross-organism pairs, matched in
  size to group (a), excluding every top-k pair and every DIAMOND hit; fixed
  seed.

A pair is **scorable** only when both proteins have at least one experimental
GO annotation. Unscorable pairs are counted and reported, not silently dropped.

### GO similarity

- Each protein's experimental GO terms are propagated to all `is_a` and
  `part_of` ancestors.
- Primary metric: **Jaccard similarity** of the propagated term sets.
- Secondary metric: **Resnik similarity**, with information content estimated
  from the annotated proteins in this run.
- Sensitivity metric: Jaccard after removing the three GO roots
  (`GO:0003674`, `GO:0005575`, `GO:0008150`), so the result cannot be driven only
  by shared generic terms.

### Statistics

- Two-sided **Mann-Whitney U** for (a) vs (c) and (b) vs (c).
- **Cliff's delta** effect size (`2U/(n₁n₂) − 1`); positive means the first group
  tends to score higher.
- **Holm** adjustment across the two comparisons.
- The hypothesis is called *supported* when (a) vs (c) has p < 0.05 and
  Cliff's δ > 0. Effect size and scorable counts are always reported alongside
  p-values, because large *n* can make tiny differences significant.

---

## Results (full run)

Run date: 2026-11-06. UniProt release 2026_04. Nextflow 26.10.5, Python 3.11.18,
DIAMOND 2.2.8, `facebook/esm2_t6_8M_UR50D`, top-k = 5, cosine ≥ 0.9, seed 42.

### Data and candidate counts

| Quantity | Value |
|---|---|
| *E. coli* K-12 reviewed proteins | 4,531 |
| *B. subtilis* 168 reviewed proteins | 4,191 |
| Experimental annotation rows | 24,438 |
| Proteins with ≥1 experimental annotation | 4,123 |
| DIAMOND pairs (e ≤ 1e-3) | 15,219 |
| Top-k embedding pairs | 22,655 |
| Embedding-only candidates (a) | 16,751 |
| Random controls (c) | 16,751 |
| Truncated sequences (> 1022 aa) | 57 (*E. coli*) + 45 (*B. subtilis*) |

### GO similarity by group

| Group | Input pairs | Scorable | Mean Jaccard | Median | Cliff's δ vs random | p (MWU) | p (Holm) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Embedding-only (a) | 16,751 | 1,828 | **0.1246** | 0.0833 | **0.235** | **6.9×10⁻³⁶** | 6.9×10⁻³⁶ |
| DIAMOND hits (b) | 15,219 | 1,843 | **0.1644** | 0.1053 | **0.356** | **1.7×10⁻⁸⁰** | 3.5×10⁻⁸⁰ |
| Random (c) | 16,751 | 1,947 | 0.0807 | 0.0526 | — | — | — |

**Sensitivity (root terms removed):** embedding-only mean 0.0976 vs random
0.0551 (p = 1.7×10⁻⁴¹, δ = 0.250); DIAMOND mean 0.1445 vs random 0.0551
(p = 2.4×10⁻¹⁰⁴, δ = 0.402). The conclusion is unchanged.

### Interpretation

- **The hypothesis is supported at the group level.** Embedding-only pairs have
  higher experimental GO similarity than random pairs, with a very small
  p-value and a **small-to-moderate** effect size (Cliff's δ = 0.235). The
  no-root sensitivity analysis is slightly stronger, so the result is not an
  artefact of shared root terms.
- **Sequence homology remains the stronger signal.** DIAMOND-hit pairs have the
  highest mean GO similarity (0.164) and the largest effect size (δ = 0.356).
  Embeddings add a weaker, complementary signal, not a replacement for sequence
  search.
- **The absolute overlap is low.** Mean Jaccard is only 0.125 for group (a),
  and many high-cosine pairs share only broad terms such as
  "catalytic activity" or "molecular_function". Finding a pair in the top-k is
  a hypothesis, not evidence that the pair is a functional analogue.
- **Annotation coverage is the binding constraint.** 14,923 of 16,751
  embedding-only candidates (89%) could not be scored because at least one
  protein lacked experimental annotation (mostly the *B. subtilis* side:
  1,736 annotation rows / 631 proteins versus 22,702 rows / 3,492 proteins for
  *E. coli*). The statistical result applies to the 1,828 scorable candidates,
  which are not a random sample of all candidates.
- **Illustrative pair:** *E. coli* P69826 (PTS mannitol-specific cryptic EIICB)
  and *B. subtilis* P05306 (PTS sucrose-specific EIIBC) have cosine 0.993 and no
  DIAMOND hit — a biologically plausible functional analogy — yet Jaccard is 0
  because their sparse experimental annotations do not overlap. This captures
  both the promise and the measurement problem.

### Honest caveats

The result does **not** show that any particular embedding-only pair is a true
functional analogue, nor that function can be transferred without homology. It
shows that the embedding neighbourhood is enriched for shared experimental GO
terms relative to random pairs, under one model, one cutoff and one annotation
snapshot. The effect is modest; the annotation sample is sparse and
taxonomically imbalanced.

---

## Known limitations

- **Mean pooling** collapses residue-level information; active sites and
  domain architecture are diluted. Attention pooling, per-domain embeddings, or
  alignments of embeddings could improve signal.
- **Truncation** at 1022 residues drops domains from 102 long proteins in this
  run; 57 *E. coli* and 45 *B. subtilis* sequences were affected.
- **Sparse, imbalanced experimental annotations.** *B. subtilis* has far fewer
  experimental GO terms than *E. coli*; most cross-organism pairs are
  unscorable, reducing power and potentially biasing the scorable subset.
- **GO graph effects.** Propagated ancestor sets share generic terms; Jaccard is
  sensitive to annotation depth and aspect composition. The no-root sensitivity
  analysis mitigates but does not eliminate this.
- **DIAMOND is not exhaustive.** Very divergent homologues can be missed, which
  is the phenomenon under study, but it also means some "embedding-only" pairs
  are genuine remote homologues.
- **Random controls are not matched** for protein length, amino-acid
  composition, abundance, or annotation depth. Embedding similarity can capture
  composition/length effects that correlate with GO similarity.
- **One model, one cutoff, one seed.** Results may not generalise; larger
  ESM-2 models and sensitivity sweeps over `top_k`, `cosine_threshold` and
  seeds are the natural next step.

---

## Repository layout

```
functionalhomology/
├── main.nf                  # Nextflow DSL2 workflow
├── nextflow.config          # params, profiles (test/full/conda/docker/local), resources
├── environment.yml          # pinned Python 3.11 / DIAMOND / PyTorch CPU stack
├── Dockerfile
├── bin/
│   ├── fh_common.py         # shared FASTA/OBO helpers
│   ├── embed.py             # ESM-2 mean-pooled embeddings
│   ├── find_pairs.py        # cosine top-k, DIAMOND join, random control
│   ├── go_eval.py           # GO propagation, Jaccard/Resnik, statistics
│   └── make_report.py       # self-contained HTML report
├── scripts/
│   └── download_data.py     # UniProt + QuickGO + GO download, test subset
├── testdata/                # committed tiny seeded dataset
├── tests/                   # pytest unit tests
├── docs/
│   ├── methods.md           # detailed methods
│   ├── dev_notes.md         # problems encountered and fixes
│   ├── full_run_report.html # archived HTML report from the full run
│   ├── full_run_eval_summary.json
│   └── full_run_pairs_summary.json
├── PLAN.md
└── LICENSE                  # MIT
```

## Tests

```bash
conda run -n functionalhomology pytest -q tests
# or, with the environment active:
pytest -q tests
```

Tests cover L2 normalisation, cosine top-k selection, DIAMOND parsing,
random-control sampling, OBO parsing, ancestor propagation, Jaccard, Resnik,
and evidence-code filtering. They use tiny synthetic data and do not load the
language model.

## Reproducibility notes

- All versions are pinned in `environment.yml`; Nextflow, Java and Nextflow
  version are documented above.
- Random seeds: test subsampling seed 42; random control seed 42.
- Full data is downloaded live from UniProt/QuickGO/GO into `data/` and is
  git-ignored. Only the tiny test dataset is committed.
- The HTML report embeds its plot and needs no external assets.
- `results/`, `work/`, `data/` and `tools/` are git-ignored.

## License

MIT — see [LICENSE](LICENSE).
