# BamScale

<img src="man/figures/BamScale-logo.png" align="right" height="170" alt="BamScale logo" />

[![R-CMD-check](https://github.com/cparsania/BamScale/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cparsania/BamScale/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/cparsania/BamScale/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/cparsania/BamScale/actions/workflows/pkgdown.yaml)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-submission_in_progress-orange)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

BamScale is a multithreaded BAM processing package for R built on top of the `ompBAM` C++ engine. It is designed for Bioconductor users who need high-throughput BAM parsing while preserving familiar `scanBam`/`readGAlignments`-style workflows.

## Why BamScale in Bioconductor Workflows?

BamScale focuses on three goals:

- speed on modern multi-core systems,
- compatibility with common Bioconductor input/output contracts,
- transparent benchmarking and reproducibility.

Key capabilities:

- OpenMP-enabled per-file parallelism via `threads`,
- optional multi-file parallelism via `BiocParallel` (`BPPARAM`),
- `ScanBamParam`-like filtering (`mapqFilter`, `flag`, `which`, `what`, `tag`),
- multiple output modes:
  - `data.frame` and `S4Vectors::DataFrame`,
  - `GenomicAlignments::GAlignments`,
  - `GenomicAlignments::GAlignmentPairs`,
  - `scanBam`-shaped list output (`as = "scanBam"`).

## Installation

### Prerequisites

- R with a C++17 toolchain
- OpenMP-capable compiler/runtime
- `ompBAM` available in your R library

### Current install route (pre-Bioconductor release)

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ompBAM")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("cparsania/BamScale")
```

### After Bioconductor acceptance

```r
BiocManager::install("BamScale")
```

## Quick Start

```r
library(BamScale)

bam <- ompBAM::example_BAM("Unsorted")

# 1) Step1-style extraction
x <- bam_read(
  file = bam,
  what = c("qname", "flag", "rname", "pos", "mapq", "cigar"),
  threads = 4
)

# 2) Seq/qual in comparator-compatible mode
sq <- bam_read(
  file = bam,
  what = c("qname", "seq", "qual"),
  as = "data.frame",
  seqqual_mode = "compatible",
  threads = 4
)

# 3) Fast chromosome-level counts
cnt <- bam_count(file = bam, threads = 4)
```

## Parallelism Model

BamScale can parallelize on two axes:

- across files via `BPPARAM` workers,
- within each file via OpenMP `threads`.

Approximate effective concurrency:

`min(length(file), bpnworkers(BPPARAM)) * threads`

When `auto_threads = TRUE`, BamScale can cap per-file threads under multi-worker execution to reduce oversubscription.

## Reproducibility and Output Equivalence

A dedicated vignette validates output comparability with Rsamtools:

- `vignettes/output-comparison.Rmd`
- compares `step1` fields and `seq+qual` fields in compatible mode

Render locally:

```r
rmarkdown::render("vignettes/output-comparison.Rmd")
```

## Benchmarking

Benchmark assets are in `inst/benchmarks/`:

- runner: `inst/benchmarks/run_server_benchmark.R`
- report: `inst/benchmarks/benchmark_report.qmd`
- protocol: `inst/benchmarks/README.md`

The benchmark design separates two tracks:

- `fair`: cross-package comparisons (BamScale vs Rsamtools/GenomicAlignments)
- `optimized`: BamScale-only optimizations (for example compact `seq/qual`)

This separation is intentional and avoids mixing optimization-only paths into cross-package claims.

## Interpreting Worker vs Thread Scaling

More workers do not always mean faster runtime at fixed total CPU budget.

For multi-file BAM workloads, performance can drop when workers increase and per-worker threads shrink because of:

- worker/process overhead,
- memory and serialization pressure,
- shared storage bandwidth contention.

This behavior is not unique to BamScale and can appear with Rsamtools/GenomicAlignments multi-process patterns as well. Benchmark reporting should always state worker/thread allocation explicitly.

## Current Limitations

- `param$which` is currently implemented as sequential filtering rather than indexed random-access jumps.
- `seqqual_mode = "compact"` is optimization-oriented and not intended for strict cross-package output-equivalence comparisons.
- `GAlignments` and `GAlignmentPairs` outputs exclude unmapped records by design.

## Community and Support

- Bioconductor support site (recommended for user-facing questions):
  - https://support.bioconductor.org
- Development issues and feature requests:
  - https://github.com/cparsania/BamScale/issues

When posting performance reports, include:

- package versions,
- hardware/storage context,
- exact benchmark command and profile,
- `threads`/`BPPARAM` settings.

## Citation

```r
citation("BamScale")
```

If BamScale contributes to performance claims, please also cite `ompBAM`.

## Contributing

Pull requests are welcome. Please include:

- a short motivation,
- tests for behavior changes,
- benchmark evidence for performance claims.

Before opening a PR, run:

```bash
R CMD build .
R CMD check --as-cran BamScale_*.tar.gz
Rscript -e "BiocCheck::BiocCheck('.')"
```

## License

MIT (`LICENSE`).
