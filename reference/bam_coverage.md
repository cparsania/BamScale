# Per-base coverage, aggregated in C++

Computes per-base read coverage by folding a per-contig difference array
inside the multithreaded C++ reader and run-length encoding the result,
without ever materialising the alignments as an R `GAlignments` object.
This is the fast path for coverage / bigWig generation on large BAMs:
the standard route
(`GenomicAlignments::coverage(readGAlignments(file))`) builds a full
`GAlignments` in R first, which is the dominant cost on 100M+-read
files.

## Usage

``` r
bam_coverage(
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
  applied per read in C++. With the default `NULL` the filter matches
  `readGAlignments()` (drop only unmapped reads). Note that `which`
  (region) restriction is not applied.

- threads:

  Number of OpenMP threads for within-file decoding.

- BPPARAM:

  A
  [BiocParallel::BiocParallelParam](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
  for across-file parallelism when `file` has more than one element.

- auto_threads:

  Logical; if `TRUE`, balance the thread / worker split.

## Value

An
[IRanges::RleList](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
(`SimpleRleList`) of integer coverage, one element per BAM-header
sequence in header order, each of length equal to the header sequence
length (contigs with no reads are all-zero runs). For multiple input
files, a named list of such `RleList`s.

## Details

The result is **byte-identical**
([`identical()`](https://rdrr.io/r/base/identical.html)) to
`GenomicAlignments::coverage(GenomicAlignments::readGAlignments(file))`:
CIGAR operations `M`, `=`, `X` and `D` contribute to coverage, `N`
introduces a gap, and `I`/`S`/`H`/`P` are ignored; only unmapped reads
(SAM flag `0x4`) are dropped by default – secondary, supplementary,
duplicate and QC-fail reads are counted, exactly as `readGAlignments()`
does.

## See also

[`fragment_sizes()`](https://cparsania.github.io/BamScale/reference/fragment_sizes.md),
[`mapq_dist()`](https://cparsania.github.io/BamScale/reference/mapq_dist.md),
[`bam_read()`](https://cparsania.github.io/BamScale/reference/bam_read.md).
