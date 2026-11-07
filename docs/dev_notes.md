# Development notes: problems and fixes

This is a running log of the issues encountered while building and running the
pipeline, and how each was resolved.

## Environment

**System state at start:** Arch Linux, 12 CPUs, 15 GiB RAM, 207 GiB free,
Docker available, network available. System Java 27, Python 3.14. No Nextflow,
no conda/mamba, no DIAMOND.

1. **Java 27 is too new for Nextflow 26.10.5.**
   - Symptom: `nextflow -version` → "Cannot find Java or it's a wrong version
     -- please make sure that Java 17 or later (up to 26) is installed".
   - Fix: installed OpenJDK 17 into a local micromamba environment and set
     `JAVA_HOME` / prepended its `bin` to `PATH` for Nextflow. The same OpenJDK
     17 is pinned in `environment.yml`.

2. **System Python 3.14 does not match the pinned scientific stack.**
   - Fix: created a local Python 3.11.18 environment with micromamba and
     installed CPU PyTorch/transformers via pip. The same stack is pinned in
     `environment.yml`, which Nextflow installs with the `conda` profile.

3. **No Nextflow / DIAMOND.**
   - Fix: installed Nextflow standalone (`get.nextflow.io`, version 26.10.5) and
     DIAMOND 2.2.8 from bioconda into the local environment. Users only need
     Java + Nextflow + conda; DIAMOND and Python come from `environment.yml`.

## Data acquisition

4. **QuickGO `evidenceCode` does not accept GO codes such as `EXP` or `IDA`.**
   - Symptom: HTTP 400 "The 'Evidence code identifier' parameter contains
     invalid values: IDA".
   - Cause: the current QuickGO annotation API expects ECO identifiers.
   - Fix: query with `evidenceCode=ECO:0000269` (experimental evidence
     ancestor), then filter the returned `goEvidence` values against the GO
     experimental code allow-list. Also filter `NOT` qualifiers.

5. **QuickGO bulk GAF download is silently capped at 10,000 annotations.**
   - Symptom: `downloadSearch` with `Accept: text/gaf` returned exactly 10,000
     annotation rows for *E. coli*, while the search API reported 23,668
     matching annotations. `limit`/`page` were ignored by the download
     endpoint.
   - Fix: use the paged JSON search endpoint with `limit=200` and parallel page
     fetch (4–6 workers). `pageInfo.total` is the number of pages, not the
     number of rows; the script handles that.

6. **`go-basic.obo` PURL redirects.**
   - `http://purl.obolibrary.org/obo/go/go-basic.obo` returns 302 to
     `https://current.geneontology.org/ontology/go-basic.obo`. `requests`
     follows redirects, and the final URL is recorded in the code/metadata.

7. **Test data would have required the full 32 MiB OBO.**
   - Fix: generate `testdata/go-basic.mini.obo`, containing only annotated
     terms plus all `is_a`/`part_of` ancestors. The slice is closed under the
     same edges used by the evaluator, so test and full propagation agree.

## DIAMOND

8. **Default DIAMOND cap would truncate group (b).**
   - `diamond blastp` defaults to `--max-target-seqs 25` in DIAMOND 2.2.8.
   - Verified with a synthetic one-query/40-target test: default produced 25
     hits, `--max-target-seqs 0` produced 40.
   - Fix: set `--max-target-seqs 0` in the pipeline. `--max-hsps 1` is kept and
     the Python parser additionally keeps the best HSP per pair.

9. **DIAMOND ids are full UniProt headers, not accessions.**
   - `qseqid`/`sseqid` are `sp|P12345|NAME_ORGANISM`.
   - Fix: normalise with the same `accession_from_header` helper used for the
     FASTA and embeddings before joining.

## Nextflow

10. **Top-level `log.info` / `if` statements are illegal in a DSL2 script with
    an explicit workflow.**
    - Symptom: "Statements cannot be mixed with script declarations".
    - Fix: moved the banner and parameter validation inside `workflow { ... }`.

11. **Calling the same process twice (`EMBED`) was rejected.**
    - Symptom: "Process 'EMBED' was called twice -- include it with a different
      alias".
    - Fix: one `EMBED` invocation receives a two-element channel
      (`proteome_a`, `proteome_b`); the two tasks still run in parallel. The
      outputs are filtered by the `name` field and re-tupled for `FIND_PAIRS`.

12. **Scripts were not executable.**
    - Symptom: REPORT exit 126, `/…/bin/make_report.py: Permission denied`.
    - Fix: `chmod +x bin/*.py scripts/*.py`; the executable bit is committed.

13. **Nextflow conda profile failed with micromamba.**
    - The env was created, but `.command.run` sourced
      `$(conda info --json | awk ...)/bin/activate`, which resolved to
      `/bin/activate` because micromamba's JSON lacks a `conda_prefix` key and
      micromamba envs lack `bin/activate`.
    - Fix: installed Miniforge locally and pointed the `conda` wrapper at real
      `conda` (26.9.4). The README warns that micromamba alone is insufficient.

14. **Conda env creation prompted for confirmation under Nextflow.**
    - Nextflow invokes `conda env create --prefix … --file …` without `-y`.
    - Fix (local development only): a wrapper appended `--yes` to
      create/install/update/remove commands. End users with a normal conda
      should not need this; Nextflow's conda integration normally passes
      `-y`. If it prompts, run `conda config --set always_yes true`.

15. **Reruns warned that `trace.tsv`/report already existed.**
    - Fix: `trace.overwrite = true`, `report.overwrite = true`.

## Analysis / report

16. **Cosine threshold 0.8 was too permissive.**
    - On the test set, 937/1000 top-k pairs exceeded 0.8 (~94%).
    - Fix: default `cosine_threshold = 0.9`; it retains ~34% of test top-k
      pairs and produces a manageable candidate set. It is a parameter.

17. **Shared GO terms were sorted by GO id, so generic ancestors dominated.**
    - Fix: sort shared terms by decreasing information content (most specific
      first) before writing `pairs_scored.tsv`; the report displays the first
      six.

18. **Jaccard including GO roots could be dismissed as a root-term artefact.**
    - Fix: added `jaccard_no_root` and a full sensitivity analysis in
      `eval_summary.json` and the HTML report. The conclusion is unchanged
      (actually slightly stronger without roots).

19. **Report needed protein names, not just accessions.**
    - Fix: parse `sp|ACC|ENTRY Description` headers from both FASTAs and show
      entry name + description in the top-candidate table.

## Verification performed

- `pytest -q tests`: 15 passed.
- `nextflow run main.nf -profile test,conda`: succeeds end to end (6 processes).
- Full run: `nextflow run main.nf -profile full,conda` (or explicit `--proteome_a` etc.); 6/6 processes succeeded.
- Sanity checks: embedding shapes (4531×320, 4191×320), 22,655 top-k pairs,
  15,219 DIAMOND pairs, 16,751 embedding-only and 16,751 random pairs,
  1,828/1,843/1,947 scorable pairs. Supporting row counts were checked against
  the generated files.
- Robustness: reran the comparison after removing GO roots; conclusion
  unchanged.

## Docker

20. **Container ran as `mambauser` (UID 57439) and could not write to the
    host work directory.**
    - Symptom: `touch: cannot touch '.command.trace': Permission denied`.
    - Fix: added `USER root` to the Dockerfile and kept
      `docker.fixOwnership = true` in the docker profile.

21. **Nextflow task metrics failed because `ps` was missing.**
    - Symptom: "Command 'ps' required by nextflow to collect task metrics
      cannot be found".
    - Fix: install `procps-ng` into the image with micromamba.

22. **Docker profile verified.**
    - `docker build -t functionalhomology:latest .` succeeds.
    - `nextflow run main.nf -profile test,docker` completes all 6 processes.
