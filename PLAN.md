# Plan and assumptions

## Goal
Test whether protein pairs across two organisms can share function without
detectable sequence homology, by comparing ESM-2 embedding similarity with
DIAMOND sequence similarity and scoring both against experimental GO
annotation similarity.

## Plan (milestones)
1. **M1** Environment check, repo skeleton, reproducible data download
   (UniProt FASTA, QuickGO experimental annotations, go-basic.obo), tiny
   seeded test dataset.
2. **M2** `embed.py` (ESM-2 mean-pooled embeddings) and a manual DIAMOND run.
3. **M3** `find_pairs.py` and `go_eval.py` working manually on test data; pytest
   unit tests.
4. **M4** Wrap all steps in a Nextflow DSL2 workflow; `-profile test,conda`
   runs end to end.
5. **M5** Self-contained HTML report.
6. **M6** Full E. coli K-12 vs B. subtilis 168 run, README, results and
   limitations, cleanup.

## Assumptions
- "Experimental GO evidence" = EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP,
  HGI, HEP (the GO consortium experimental/high-throughput experimental set).
  QuickGO `evidenceCode=ECO:0000269` returns descendants of *experimental
  evidence*; we additionally filter the returned evidence-code string against
  this exact allow-list, and drop `NOT` qualifiers.
- Orthologues/identical proteins are not filtered out. DIAMOND hits therefore
  include genuine homologues and the analysis asks only whether *sequence
  hits* have higher GO similarity than random.
- Mean-pooled, L2-normalised ESM-2 embeddings are used as the functional
  representation. Residues beyond 1022 are truncated.
- DIAMOND e-value cutoff is 1e-3 (`--very-sensitive`).
- Cross-organism random pairs are sampled uniformly from the Cartesian product
  of the two proteomes, excluding all top-k embedding pairs and all DIAMOND
  hits, with a fixed seed.
- Group (b) in GO evaluation is the complete set of DIAMOND hits passing the
  e-value cutoff (not only those that also appear in the embedding top-k).
- GO term similarity uses Jaccard on ancestor-propagated term sets (is_a and
  part_of). Resnik is reported as a secondary metric using information content
  estimated from the annotated proteins in this dataset.
- Only pairs where BOTH proteins have at least one experimental annotation are
  scored; scorable counts per group are reported.

## Environment notes (this machine)
- Arch Linux, 12 CPUs, 15 GiB RAM, Docker available, system Java 27 (too new
  for Nextflow 26.x), Python 3.14 (too new for pinned stack).
- Installed locally under `tools/` (git-ignored): micromamba, a Python 3.11
  environment with OpenJDK 17, DIAMOND 2.2.8, PyTorch CPU, transformers,
  goatools; standalone Nextflow 26.10.5.
- Everything else is documented in README.md for a fresh clone.
