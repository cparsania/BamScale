# Write per-base coverage to a bigWig, computed in C++

Computes per-base coverage (identical to
[`bam_coverage()`](https://cparsania.github.io/BamScale/reference/bam_coverage.md))
and writes it directly to a bigWig file entirely in C++ via the bundled
libBigWig, without ever materialising the coverage as an R object. This
is the fast path for generating coverage tracks from large BAMs: the
usual route (`rtracklayer::export.bw(bam_coverage(file), out)`)
round-trips a large `RleList` through R, and the serialisation
dominates. Only non-zero coverage runs are written (uncovered bases are
implicitly zero, the bigWig convention); importing the file reconstructs
full-length coverage from the header lengths.

## Usage

``` r
bam_coverage_bigwig(
  file,
  outfile,
  param = NULL,
  threads = 1L,
  BPPARAM = BiocParallel::bpparam(),
  auto_threads = FALSE,
  n_zooms = 10L,
  compress_level = -1L,
  parallel = TRUE,
  verbose = FALSE
)
```

## Arguments

- file:

  A BAM file path,
  [Rsamtools::BamFile](https://rdrr.io/pkg/Rsamtools/man/BamFile-class.html),
  or a vector / `BamFileList` of them.

- outfile:

  Output bigWig path(s); must have the same length as `file`.

- param:

  Optional
  [Rsamtools::ScanBamParam](https://rdrr.io/pkg/Rsamtools/man/ScanBamParam-class.html)
  (or a compatible list). Only its `flag` and `mapqFilter` are honoured,
  applied per read in C++. With the default `NULL` the filter matches
  `readGAlignments()` (drop only unmapped reads).

- threads:

  Number of OpenMP threads for within-file decoding.

- BPPARAM:

  A
  [BiocParallel::BiocParallelParam](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  for across-file parallelism when `file` has more than one element.

- auto_threads:

  Logical; if `TRUE`, balance the thread / worker split.

- n_zooms:

  Maximum number of bigWig zoom levels to generate (default `10`).

- compress_level:

  Integer zlib deflate level for the bigWig blocks: `-1` (default) uses
  zlib's default (level 6, matching other bigWig writers); `1` to `9`
  trade file size for speed (`1` is the fastest write and gives the
  largest file). Coverage values are unaffected – only file size and
  write time change.

- parallel:

  Logical; if `TRUE` (default), compress the bigWig data and zoom blocks
  across the OpenMP `threads` (deferred parallel compression). The
  output is byte-identical to serial compression. Set `FALSE` for the
  classic single-threaded write.

- verbose:

  Logical; if `TRUE`, print a per-phase timing breakdown (coverage
  compute / interval write / zoom + index finalize) to stderr.

## Value

The output path(s), invisibly.

## Details

The written values are identical to
`GenomicAlignments::coverage(GenomicAlignments::readGAlignments(file))`
(`M`/`=`/`X`/`D` contribute, `N` is a gap, `I`/`S`/`H`/`P` are ignored;
only unmapped reads are dropped by default). Coverage values are
integers stored as the bigWig `float` type, which is exact for depths up
to `2^24`.

## See also

[`bam_coverage()`](https://cparsania.github.io/BamScale/reference/bam_coverage.md)
for the in-memory `RleList`.
