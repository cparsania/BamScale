# BamScale benchmark harness

This folder holds two complementary, budget-aware benchmark drivers and the
assets used to summarise their results. Nothing here is loaded at package run
time — it exists to reproduce the numbers reported in the benchmark article.

- **`run_server_benchmark.R`** — read-pattern micro-benchmarks (`step1` /
  `galignments` / `seqqual`) against `Rsamtools` and `GenomicAlignments`.
- **`run_workflow_benchmark.R`** — end-to-end workflow benchmarks
  (coverage → bigWig, ATAC fragment-size QC) reported as a read / compute /
  write phase decomposition so the Amdahl bound on the end-to-end speedup is
  explicit.
- **`bench_common.R`** — shared helpers + the phase-aware timing engine sourced
  by the workflow driver (defines functions only; no side effects on source).
- **`download_atac_data.R`** — fetch + index a public paired-end ATAC-seq
  dataset; run once before the ATAC workflow (no ATAC data ships with the host).

Results are written as machine-readable CSVs to `--outdir` (e.g. a git-ignored
`benchmark_results/`); none are committed. The narrative summary and figures live
in `vignettes/benchmark-results.Rmd` (the pkgdown "Benchmarks" article).

## Headline results

Intel Xeon Gold 6252 (96 cores), warm page cache, median of 5 iterations. All
BamScale output is verified byte-identical to the standard tool.

| Workload | Comparator | Speedup |
| --- | --- | :--: |
| Core alignment fields (`step1`) | `Rsamtools::scanBam` | **2.3×** |
| `GAlignments` object | `GenomicAlignments::readGAlignments` | **3.2×** |
| Sequence + base quality (compatible) | `Rsamtools::scanBam` | **2.8×** |
| ATAC fragment-size QC (end-to-end) | `ATACseqQC::fragSizeDist` seam | **3.7×** |
| Coverage → `RleList` (end-to-end) | `readGAlignments` + `coverage()` | **2.5×** |
| Coverage → bigWig (end-to-end) | + `export.bw()` (write-bound) | **1.2×** |

The end-to-end workflow speedup tracks how read-bound the workflow is (Amdahl's
law): near-pure-read ATAC QC keeps almost the full read speedup, while the
shared, single-threaded bigWig writer bounds the coverage→bigWig gain.

## Read-pattern micro-benchmarks — `run_server_benchmark.R`

```bash
Rscript inst/benchmarks/run_server_benchmark.R \
  --profile=balanced --budget-threads=48 --max-threads=48 \
  --bam-dir=/path/to/bam_directory --n-files=12 \
  --bp-backend=snow --iterations=5 --outdir=benchmark_results
```

If no `--bam-dir` is given, BAMs are resolved from `chipseqDBData::H3K9acData()`.

**Tracks.** `fair` = comparator-safe runs used for BamScale vs
Rsamtools/GenomicAlignments comparisons; `optimized` = BamScale-only runs
(compact `seq/qual`). Report sections that claim cross-package comparison use
only the `fair` track.

**Profiles.** `bamscale_showcase` (default, `step1`-focused), `balanced`
(broader workloads for fair comparison), `full` (widest grid).

**Budget model.** `--budget-threads` is the total compute budget
(`min(48, cores)` for showcase/balanced, `cores` for `full`); `--max-threads`
is the per-run thread ceiling. Multi-file BamScale uses a balanced split:
`threads_each = floor(max_threads / workers)`, so total ≈ `workers * threads_each`.

**Useful options:** `--profile`, `--budget-threads`, `--max-threads`,
`--threads` (single-file grid), `--workers` (multi-file grid),
`--single-workloads` / `--multi-workloads=step1,galignments,seqqual`,
`--seqqual-compact=true|false`, `--include-rsamtools`, `--include-galignments`,
`--ensure-index`, `--allow-repeat-files`, `--bp-backend`, `--iterations`.

## End-to-end workflow benchmarks — `run_workflow_benchmark.R`

BamScale replaces exactly the BAM-read step of each workflow; the compute and
write steps are byte-identical work on both arms. Single-file cases sweep the
OpenMP thread axis with full per-phase timing; multi-file cases sweep the
`BiocParallel` worker axis at a fixed core budget.

```bash
# Coverage workflow (uses the same chipseqDBData BAMs by default)
Rscript inst/benchmarks/run_workflow_benchmark.R \
  --include-coverage=true --include-atacqc=false \
  --threads=1,4,12,24,48 --workers=1,2,4,8,12 \
  --n-files=12 --budget-threads=48 --max-threads=48 \
  --bp-backend=snow --iterations=5 --write-bigwig=true \
  --outdir=benchmark_results

# ATAC fragment-size QC (fetch data first, once)
Rscript inst/benchmarks/download_atac_data.R \
  --experiment=ENCSRxxxxxxx --outdir=/path/to/atac/bams
Rscript inst/benchmarks/run_workflow_benchmark.R \
  --include-coverage=false --include-atacqc=true \
  --atac-bam-dir=/path/to/atac/bams \
  --threads=1,4,12,24,48 --iterations=5 --outdir=benchmark_results
```

`download_atac_data.R` accepts `--urls=`, `--experiment=ENCSR…` (resolved via the
ENCODE REST API), or `--files=ENCFF…`; it indexes each BAM and verifies paired-end
reads.

## Outputs

Read-pattern runs: `summary.csv` (per-case stats), `iterations.csv`,
`files.csv`, `config.txt`, `sessionInfo.txt`, and (with `ggplot2`)
`plot_single_scaling.png` / `plot_multi_scaling.png`.

Workflow runs add per-phase columns (`read_s`, `compute_s`, `write_s`,
`total_s`, plus per-phase CPU) and a correctness table recording the
`identical()` / byte-identical checks against the standard tools.

## Render the report

```bash
Rscript -e "rmarkdown::render('vignettes/benchmark-results.Rmd')"
```
