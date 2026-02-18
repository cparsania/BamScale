# BamScale

BamScale is a multithreaded BAM reader/scanner for R built on top of `ompBAM`, with an argument surface designed to feel familiar to Bioconductor users.

## Goals

- Faster sequential BAM read/scan throughput via OpenMP.
- Argument compatibility with common Bioc patterns (`file`, `param`, `what`, `tag`, `BPPARAM`).
- Explicit nested parallelism control via `threads`, `BPPARAM`, and `auto_threads`.
- Output compatibility for downstream pipelines (`data.frame`, `S4Vectors::DataFrame`, `GenomicAlignments::GAlignments`, `GenomicAlignments::GAlignmentPairs`, scan-like list output).

## Main Functions

- `bam_read()`
  - Supports `character`, `BamFile`, `BamFileList` inputs.
  - Honors key `ScanBamParam` fields: `mapqFilter`, `flag`, `which`, `what`, `tag`.
  - Supports `as = "DataFrame" | "data.frame" | "GAlignments" | "GAlignmentPairs" | "scanBam"`.

- `bam_count()`
  - Fast chromosome-level counts with same `param` filtering semantics.

## Usage

```r
library(BamScale)

bam <- ompBAM::example_BAM("Unsorted")

# scanBam/readGAlignments-like usage
x <- bam_read(
  file = bam,
  what = c("qname", "flag", "rname", "pos", "mapq", "cigar"),
  threads = 4
)

# Tag extraction
x_tag <- bam_read(
  file = bam,
  what = c("qname", "rname", "pos"),
  tag = c("NH", "NM"),
  threads = 4
)

# scanBam()-shaped list output
x_scan <- bam_read(
  file = bam,
  what = c("qname", "flag", "rname", "pos"),
  tag = c("NM"),
  as = "scanBam",
  threads = 4
)

# Strict scanBam() range batching (one element per range label, including empty ranges)
x_scan_ranges <- bam_read(
  file = bam,
  param = list(
    which = data.frame(
      seqname = c("1", "1"),
      start = c(1L, 250000000L),
      end = c(100000L, 250010000L),
      label = c("region_hit", "region_empty")
    )
  ),
  what = c("qname", "flag"),
  as = "scanBam",
  threads = 4
)

# Paired-end output (when mates are available)
x_pairs <- bam_read(
  file = bam,
  what = c("qname", "flag", "rname", "pos", "cigar", "strand"),
  as = "GAlignmentPairs",
  include_unmapped = FALSE,
  threads = 4
)

# Count summary
counts <- bam_count(file = bam, threads = 4)
```

## Parallelism (`threads` vs `BPPARAM`)

- `BPPARAM` parallelizes across BAM files.
- `threads` parallelizes within each BAM file through OpenMP.
- Approximate total concurrency is
  `min(length(file), BiocParallel::bpnworkers(BPPARAM)) * threads`.
- With `auto_threads = TRUE`, per-file OpenMP threads are automatically capped to
  `max(1, min(threads, floor(available_cores / workers_eff)))`, where
  `workers_eff = min(length(file), BiocParallel::bpnworkers(BPPARAM))`.

```r
# Multi-file run with automatic nested parallelism balancing
if (requireNamespace("BiocParallel", quietly = TRUE)) {
  files <- c(sample1 = bam, sample2 = bam)
  bp <- BiocParallel::SerialParam()

  out <- bam_read(
    file = files,
    what = c("qname", "flag", "rname", "pos"),
    threads = 8,
    BPPARAM = bp,
    auto_threads = TRUE
  )
}
```

## Notes

- Region filtering via `param$which` is supported as a sequential filter (not random index jumps).
- For `as = "GAlignments"` and `as = "GAlignmentPairs"`, unmapped records are dropped.
- `as = "scanBam"` mirrors `scanBam()` list shape: no `which` gives one unnamed batch; `which` gives one named batch per range label, including empty ranges.
- In `as = "scanBam"` mode, `seq` and `qual` are returned as `DNAStringSet` and `PhredQuality` when Biostrings is available.
- Multi-file inputs return a named list; optionally parallelized with `BPPARAM`.
