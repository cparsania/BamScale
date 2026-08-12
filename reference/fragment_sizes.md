# Fragment-size distribution, aggregated in C++

Computes the paired-end fragment-size (insert-size) distribution – a
table of `abs(TLEN)` counts – by accumulating the histogram inside the
multithreaded C++ reader, without ever materialising the per-read insert
sizes as an R vector. This is the fast path for ATAC-seq / paired-end
fragment-size QC on large BAMs: the histogram is folded together inside
the OpenMP threads and only the compact `(fragment_size, count)` table
crosses back into R (cf. `scanBam(what = "isize")` followed by
`table(abs(isize))`, which materialises every insert size in R first).

## Usage

``` r
fragment_sizes(
  file,
  param = NULL,
  threads = 1L,
  BPPARAM = BiocParallel::bpparam(),
  auto_threads = FALSE,
  max_fragment = 100000L,
  drop_mate_unmapped = TRUE
)
```

## Arguments

- file:

  A BAM file path,
  [Rsamtools::BamFile](https://rdrr.io/pkg/Rsamtools/man/BamFile-class.html),
  or a vector / `BamFileList` of them.

- param:

  Optional
  [Rsamtools::ScanBamParam](https://rdrr.io/pkg/Rsamtools/man/ScanBamParam-class.html)
  (or a compatible list). Only its `flag` and `mapqFilter` are honoured,
  applied per read in C++ – e.g.
  `Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag( isSecondaryAlignment = FALSE, isUnmappedQuery = FALSE, isNotPassingQualityControls = FALSE))`
  to reproduce the standard ATAC fragment-size filter.

- threads:

  Number of OpenMP threads for within-file decoding.

- BPPARAM:

  A
  [BiocParallel::BiocParallelParam](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  for across-file parallelism when `file` has more than one element.

- auto_threads:

  Logical; if `TRUE`, balance the thread / worker split.

- max_fragment:

  Integer; fragment sizes up to this value are tallied in a dense
  histogram and larger ones in an overflow map. Every observed size is
  returned exactly regardless of this value (default `1e5`).

- drop_mate_unmapped:

  Logical; if `TRUE` (default), reads whose mate is unmapped (SAM flag
  `0x8`) are excluded – matching
  [`Rsamtools::scanBam`](https://rdrr.io/pkg/Rsamtools/man/scanBam.html),
  which reports their `isize` as `NA` (no defined template length), so
  that `table(abs(isize))` drops them. Set `FALSE` to tally `abs(TLEN)`
  over all flag-passing reads.

## Value

A `data.frame` with columns `fragment_size` (integer `abs(TLEN)`) and
`count` (numeric), one row per observed fragment size in ascending
order. For multiple input files, a named list of such data.frames.

## See also

[`bam_read()`](https://cparsania.github.io/BamScale/reference/bam_read.md)
for reading records,
[`bam_count()`](https://cparsania.github.io/BamScale/reference/bam_count.md)
for per-chromosome counts.
