# Count BAM records with Bioconductor-compatible filtering

`bam_count()` provides a fast chromosome-level count summary, honoring
key filtering fields from `ScanBamParam` (`mapqFilter`, `flag`, and
`which`).

## Usage

``` r
bam_count(
  file,
  param = NULL,
  threads = 1L,
  BPPARAM = NULL,
  auto_threads = FALSE,
  include_unmapped = TRUE
)
```

## Arguments

- file:

  BAM input (`character`, `BamFile`, or `BamFileList`).

- param:

  Optional
  [`Rsamtools::ScanBamParam`](https://rdrr.io/pkg/Rsamtools/man/ScanBamParam-class.html)
  (or compatible list).

- threads:

  Requested number of OpenMP threads. May be capped when
  `auto_threads = TRUE`.

- BPPARAM:

  Optional `BiocParallel` parameter for multi-file operation.

- auto_threads:

  Logical; when `TRUE` and `BPPARAM` has multiple workers, BamScale
  automatically caps per-file OpenMP threads to avoid oversubscription.

- include_unmapped:

  Whether to include an extra `*` row for unmapped records.

## Value

For one file: a `data.frame` with columns `seqname`, `seqlength`,
`count`. For multiple files: named list of such `data.frame`s.

## Details

Parallelism behavior matches
[`bam_read()`](https://cparsania.github.io/BamScale/reference/bam_read.md):
`BPPARAM` distributes work across BAM files, while `threads` controls
OpenMP work within each file. If `auto_threads = TRUE` and `BPPARAM` has
multiple workers, per-file OpenMP threads are capped using
`max(1, min(threads, floor(available_cores / workers_eff)))`, where
`workers_eff = min(length(file), BiocParallel::bpnworkers(BPPARAM))`.

## Examples

``` r
if (requireNamespace("ompBAM", quietly = TRUE)) {
  bam <- ompBAM::example_BAM("Unsorted")
  bam_count(bam, threads = 2)
}
#>    seqname seqlength count
#> 1        1 248956422  1308
#> 2       10 133797422   308
#> 3       11 135086622   600
#> 4       12 133275309   648
#> 5       13 114364328   162
#> 6       14 107043718   230
#> 7       15 101991189   294
#> 8       16  90338345   334
#> 9       17  83257441   570
#> 10      18  80373285    78
#> 11      19  58617616   600
#> 12       2 242193529   586
#> 13      20  64444167   168
#> 14      21  46709983    76
#> 15      22  50818468   200
#> 16       3 198295559   454
#> 17       4 190214555   252
#> 18       5 181538259   516
#> 19       6 170805979   824
#> 20       7 159345973   412
#> 21       8 145138636   452
#> 22       9 138394717   372
#> 23      MT     16569   288
#> 24       X 156040895   254
#> 25       Y  57227415    14
#> 26       *        NA     0
```
