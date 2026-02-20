# BamScale Server Benchmark Workflow (budget-aware)

This folder contains a server-first benchmark workflow for BamScale with configurable compute budget and profiles.

## Benchmark tracks

- `fair`: comparator-safe runs used for BamScale vs Rsamtools/GenomicAlignments comparisons.
- `optimized`: BamScale-only optimization runs (currently compact `seq/qual` mode).

Report sections that claim cross-package comparison use only the `fair` track.

## Profiles

- `bamscale_showcase` (default): step1-focused, designed to highlight BamScale strengths.
- `balanced`: includes broader workloads for mixed/fair comparison.
- `full`: widest workload/grid coverage.

## Budget model

- `--budget-threads`: target total compute budget (default: `min(48, detected_cores)`).
- `--max-threads`: effective per-run thread ceiling (defaults to `--budget-threads`).
- `--seqqual-compact=true|false`: include BamScale compact `seq/qual` runs (default: `TRUE` for `balanced/full`, `FALSE` for `bamscale_showcase`).
- Multi-file BamScale uses balanced split:
  - `threads_each = floor(max_threads / workers)`
  - approx total threads: `workers * threads_each`.

## 1) Run benchmark on server (recommended)

Example (showcase profile, favorable to BamScale):

```bash
Rscript inst/benchmarks/run_server_benchmark.R \
  --profile=bamscale_showcase \
  --budget-threads=48 \
  --max-threads=48 \
  --bam-dir=/path/to/bam_directory \
  --n-files=12 \
  --bp-backend=snow \
  --iterations=3 \
  --outdir=/path/to/benchmark_results
```

Example (balanced profile):

```bash
Rscript inst/benchmarks/run_server_benchmark.R \
  --profile=balanced \
  --budget-threads=48 \
  --max-threads=48 \
  --bam-dir=/path/to/bam_directory \
  --n-files=12 \
  --bp-backend=snow \
  --iterations=3 \
  --outdir=/path/to/benchmark_results
```

If no BAM input is provided, the script resolves BAMs from `chipseqDBData::H3K9acData()`.

## 2) Render report

```bash
quarto render inst/benchmarks/benchmark_report.qmd \
  -P results_dir:/path/to/benchmark_results/run_YYYYMMDD_HHMMSS
```

## Main outputs

- `summary.csv`: per-case summary statistics
- `iterations.csv`: per-iteration elapsed times
- `files.csv`: selected BAM files and sizes
- `config.txt`: parsed and effective config
- `sessionInfo.txt`: reproducibility snapshot
- `plot_single_scaling.png`, `plot_multi_scaling.png` (if `ggplot2` installed)

## Useful options

- `--profile=bamscale_showcase|balanced|full`
- `--budget-threads=<int>`
- `--max-threads=<int>`
- `--threads=...` (override single-file grid)
- `--workers=...` (override multi-file worker grid)
- `--multi-workloads=step1,seqqual,galignments`
- `--seqqual-compact=true|false`
- `--include-rsamtools=true|false`
- `--include-galignments=true|false`
- `--ensure-index=true|false`
