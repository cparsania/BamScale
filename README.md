
# BamScale <img src="man/figures/BamScale-logo.png" align="right" height="240"/>

[![R-CMD-check](https://github.com/cparsania/BamScale/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/cparsania/BamScale/actions/workflows/R-CMD-check.yaml)
[![Benchmark Report](https://img.shields.io/badge/Benchmarks-Quarto-blue)](inst/benchmarks/)
[![codecov](https://img.shields.io/badge/Coverage-Codecov-lightgrey)](https://app.codecov.io/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-targeted-green)](https://bioconductor.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)


BamScale is a Bioconductor-friendly, multithreaded BAM processing package for R built on top of the `ompBAM` C++ engine.

It is designed for users who want `scanBam`/`readGAlignments`-style ergonomics with explicit control of modern CPU parallelism.

## Why BamScale?

BAM analysis pipelines in R often face two practical constraints:

-   high per-file latency for large BAMs when using single-thread readers,
-   workflow friction when switching between fast low-level tools and Bioconductor-native data structures.

BamScale addresses both by combining:

-   an OpenMP-enabled C++ backend (`ompBAM`),
-   a Bioconductor-compatible R interface (`file`, `param`, `what`, `tag`, `BPPARAM`),
-   direct output modes for downstream ecosystem interoperability.

## ompBAM vs BamScale

`ompBAM` and BamScale are complementary layers.

-   `ompBAM` provides the high-performance C++ BAM decode/scan core.
-   BamScale provides the user-facing R contract for Bioconductor workflows.

### What ompBAM already solves

-   fast sequential BAM parsing and decompression in C++,
-   thread-level parallelism via OpenMP,
-   low-level access to alignment fields.

### What BamScale adds

-   Bioconductor-style API:
    -   `bam_read()` for read extraction,
    -   `bam_count()` for chromosome-level counting,
-   compatibility with `ScanBamParam` semantics (`mapqFilter`, `flag`, `which`, `what`, `tag`),
-   output modes tailored to R users:
    -   `data.frame`, `S4Vectors::DataFrame`,
    -   `GenomicAlignments::GAlignments`, `GenomicAlignments::GAlignmentPairs`,
    -   scan-like `list` output (`as = "scanBam"`),
-   nested parallelism control (`threads`, `BPPARAM`, `auto_threads`) with cpuset-aware thread capping,
-   reproducible benchmark framework with strict `fair` vs `optimized` tracks.

## Main Functions

-   `bam_read()`
    -   Inputs: `character`, `Rsamtools::BamFile`, `Rsamtools::BamFileList`
    -   Supported fields: `qname`, `flag`, `rname`, `strand`, `pos`, `qwidth`, `mapq`, `cigar`, `mrnm`, `mpos`, `isize`, `seq`, `qual`
    -   Output modes: `"DataFrame"`, `"data.frame"`, `"GAlignments"`, `"GAlignmentPairs"`, `"scanBam"`
    -   `seqqual_mode`:
        -   `"compatible"` (default): character seq/qual for Bioconductor compatibility
        -   `"compact"`: raw list-columns for high-throughput seq/qual extraction
-   `bam_count()`
    -   fast chromosome-level counts with the same filtering semantics (`mapqFilter`, `flag`, `which`)

## Installation

BamScale requires:

-   R with C++17 toolchain,
-   OpenMP-capable compiler/runtime,
-   `ompBAM` installed and discoverable by R.

Install from source (inside this repository):

``` bash
R CMD INSTALL .
```

## Quick Start

``` r
library(BamScale)

bam <- ompBAM::example_BAM("Unsorted")

# 1) Bioconductor-familiar extraction
x <- bam_read(
  file = bam,
  what = c("qname", "flag", "rname", "pos", "mapq", "cigar"),
  threads = 4
)

# 2) Seq/qual extraction (default compatible mode)
sq <- bam_read(
  file = bam,
  what = c("qname", "seq", "qual"),
  as = "data.frame",
  seqqual_mode = "compatible",
  threads = 4
)

# 3) High-throughput compact seq/qual path
sq_compact <- bam_read(
  file = bam,
  what = c("qname", "seq", "qual"),
  as = "data.frame",
  seqqual_mode = "compact",
  threads = 4
)

# 4) scanBam-shaped output
scan_like <- bam_read(
  file = bam,
  what = c("qname", "flag"),
  tag = c("NM"),
  as = "scanBam",
  threads = 4
)

# 5) Fast count summary
counts <- bam_count(file = bam, threads = 4)
```

## Parallelism Model

BamScale supports two parallel axes:

-   across files: `BPPARAM` workers,
-   within file: OpenMP `threads`.

Approximate total concurrency:

`min(length(file), BiocParallel::bpnworkers(BPPARAM)) * threads`

If `auto_threads = TRUE`, per-file OpenMP threads are capped to reduce oversubscription under multi-worker runs.

## Benchmarking and Reproducibility

Benchmark assets live under `inst/benchmarks/`.

-   Runner: `inst/benchmarks/run_server_benchmark.R`
-   Report: `inst/benchmarks/benchmark_report.qmd`
-   Protocol: `inst/benchmarks/README.md`

The benchmark framework separates:

-   `fair` track: strict BamScale vs comparator comparison (`Rsamtools`, `GenomicAlignments`),
-   `optimized` track: BamScale-only optimized configurations (for example compact `seq/qual`).

This prevents mixing optimization-only paths into cross-package claims.

## Current Limitations (Transparent)

-   `param$which` is implemented as sequential filtering, not indexed random-access jumps.
-   `seqqual_mode = "compact"` is currently supported for `as = "data.frame"` and `as = "DataFrame"`.
-   `GAlignments` / `GAlignmentPairs` outputs drop unmapped records by design.
-   Multi-file performance does not necessarily improve monotonically with more workers.
    -   At fixed total thread budget, higher worker counts can lose to lower worker counts because of process overhead, serialization, memory pressure, and storage contention.
-   Observed speedups are workload-, storage-, and hardware-dependent.

## What Has Been Improved Recently

-   Added compact seq/qual extraction path for lower-overhead high-throughput runs.
-   Optimized core C++ hot paths (decode/merge/allocation behavior).
-   Added cpuset-aware core detection for more reliable thread budgeting in container/HPC environments.
-   Redesigned benchmark workflow/report to make `fair` vs `optimized` interpretation explicit.

## License

MIT (`LICENSE`).
