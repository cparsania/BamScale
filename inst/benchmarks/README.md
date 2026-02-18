# BamScale Server Benchmark Workflow

This folder contains a server-first benchmark workflow for BamScale.

## Prerequisites

- `BamScale` must be installed and loadable on the server (`library(BamScale)`).
- Optional comparator packages: `Rsamtools`, `GenomicAlignments`, `BiocParallel`.
- Optional report/plots: `quarto`, `knitr`, `ggplot2`.

## 1) Run Benchmarks (recommended)

Use the plain R script for execution on server/HPC nodes.

```bash
Rscript inst/benchmarks/run_server_benchmark.R \
  --bam-dir=/path/to/bam_directory \
  --n-files=24 \
  --threads=1,2,4,8,12,18,24,36,48,72 \
  --workers=1,2,4,6,8,12,18,24,36,72 \
  --max-threads=72 \
  --iterations=3 \
  --multi-workloads=step1,seqqual,galignments \
  --outdir=/path/to/benchmark_results
```

Key options:

- `--bam-dir`: directory of BAM files (recursive by default)
- `--bam-files`: comma-separated explicit BAM paths
- `--bam-file-list`: text file with one BAM path per line
- `--threads`: single-file thread grid
- `--workers`: multi-file worker grid
- `--max-threads`: requested ceiling (for `auto_threads` and balanced runs)
- `--n-files`: number of files used for multi-file tests
- `--multi-workloads`: any of `step1`, `seqqual`, `galignments`
- `--ensure-index=true|false`: create `.bai` when missing for comparator methods
- `--allow-repeat-files=true|false`: allow repeated files when fewer than `n-files`

If no BAM input is provided, the script tries `chipseqDBData::H3K9acData()`.

Outputs are written to a timestamped run directory:

- `summary.csv`
- `iterations.csv`
- `files.csv`
- `config.txt`
- `sessionInfo.txt`
- plot PNGs (if `ggplot2` is installed)

## 2) Render Optional Report (Quarto)

Use Quarto only for reporting from already-generated CSV outputs.

```bash
quarto render inst/benchmarks/benchmark_report.qmd \
  -P results_dir:/path/to/benchmark_results/run_YYYYMMDD_HHMMSS
```

## Recommendation

- **Execution**: plain `Rscript` (robust, scheduler-friendly, easier retries/logging).
- **Presentation**: Quarto report generated after execution from saved CSV results.
