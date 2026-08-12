# Mapping-quality (MAPQ) distribution, aggregated in C++

Computes the MAPQ distribution – a table of per-value alignment-count –
by accumulating the histogram inside the multithreaded C++ reader,
without ever materialising the per-read MAPQ values as an R vector. This
is the fast path for mapping-quality QC on large BAMs: each thread
tallies MAPQ into a private 256-slot histogram inside the OpenMP region
and only the compact `(mapq, count)` table crosses back into R (cf.
`scanBam(what = "mapq")` followed by `table(mapq)`, which materialises
every MAPQ in R first).

## Usage

``` r
mapq_dist(
  file,
  param = NULL,
  threads = 1L,
  BPPARAM = BiocParallel::bpparam(),
  auto_threads = FALSE
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
  applied per read in C++. Note that a non-zero `mapqFilter` removes
  reads below that MAPQ from the distribution.

- threads:

  Number of OpenMP threads for within-file decoding.

- BPPARAM:

  A
  [BiocParallel::BiocParallelParam](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  for across-file parallelism when `file` has more than one element.

- auto_threads:

  Logical; if `TRUE`, balance the thread / worker split.

## Value

A `data.frame` with columns `mapq` (integer, 0-255) and `count`
(numeric), one row per observed MAPQ value in ascending order. For
multiple input files, a named list of such data.frames.

## See also

[`fragment_sizes()`](https://cparsania.github.io/BamScale/reference/fragment_sizes.md)
for the fragment-size distribution,
[`bam_count()`](https://cparsania.github.io/BamScale/reference/bam_count.md)
for per-chromosome counts.
