# BamScale Server Benchmark Workflow (budget-aware)

This folder contains a server-first benchmark workflow for BamScale with configurable compute budget and profiles.

## Current benchmark assets

- Final `step1` + `galignments` run:
  - `inst/benchmarks/run_20260320_133141/`
- Final `seqqual` run:
  - `inst/benchmarks/run_20260320_162359/`
- Pkgdown article source:
  - `vignettes/benchmark-results.Rmd`

## Benchmark tracks

- `fair`: comparator-safe runs used for BamScale vs Rsamtools/GenomicAlignments comparisons.
- `optimized`: BamScale-only optimization runs (currently compact `seq/qual` mode).

Report sections that claim cross-package comparison use only the `fair` track.

## Profiles

- `bamscale_showcase` (default): step1-focused, designed to highlight BamScale strengths.
- `balanced`: includes broader workloads for mixed/fair comparison.
- `full`: widest workload/grid coverage.

## Budget model

- `--budget-threads`: target total compute budget.
  Default: `min(48, detected_cores)` for `bamscale_showcase`/`balanced`, and `detected_cores` for `full`.
- `--max-threads`: effective per-run thread ceiling (defaults to `--budget-threads`).
- `--seqqual-compact=true|false`: include BamScale compact `seq/qual` runs (default: `TRUE` for `balanced/full`, `FALSE` for `bamscale_showcase`).
- Multi-file BamScale uses balanced split:
  - `threads_each = floor(max_threads / workers)`
  - approx total threads: `workers * threads_each`.

## Scaling interpretation note

- In multi-file runs, increasing `workers` while reducing per-worker threads can be slower even at similar `total_threads`.
- Main causes: process/socket serialization overhead (`SnowParam`), object allocation pressure, and shared disk/memory bandwidth contention.
- This pattern can appear for BamScale, Rsamtools, and GenomicAlignments in multi-process mode.
- The report includes a dedicated **Worker/Thread Tradeoff Note** table to make this explicit.

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

## 2) Render article/report

Pkgdown article source:

```bash
/usr/local/bin/Rscript -e "rmarkdown::render('vignettes/benchmark-results.Rmd')"
```

## Main outputs

- `summary.csv`: per-case summary statistics
- `iterations.csv`: per-iteration elapsed times
- `files.csv`: selected BAM files and sizes
- `config.txt`: parsed and effective config
- `sessionInfo.txt`: reproducibility snapshot
- `plot_single_scaling.png`, `plot_multi_scaling.png` (if `ggplot2` installed)

## Current benchmark interpretation

- `step1`: strong BamScale win in single-file mode and best-case multi-file mode; best multi-file plan favored low `W`, high `T`
- `galignments`: strong BamScale win in both single-file and multi-file settings
- `seqqual`:
  - fair compatibility path remains slower than the best Rsamtools multi-file result
  - compact mode substantially improves BamScale and produces the strongest single-file `seqqual` result

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
