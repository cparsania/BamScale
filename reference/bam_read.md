# Fast BAM reading with Bioconductor-compatible arguments

`bam_read()` is a multithreaded sequential BAM reader built on top of
`ompBAM`. The interface is designed to be familiar to users of
[`Rsamtools::scanBam()`](https://rdrr.io/pkg/Rsamtools/man/scanBam.html),
[`GenomicAlignments::readGAlignments()`](https://rdrr.io/pkg/GenomicAlignments/man/readGAlignments.html),
and
[`GenomicAlignments::readGAlignmentPairs()`](https://rdrr.io/pkg/GenomicAlignments/man/readGAlignments.html).

## Usage

``` r
bam_read(
  file,
  param = NULL,
  what = NULL,
  tag = NULL,
  as = c("DataFrame", "data.frame", "GAlignments", "GAlignmentPairs", "scanBam"),
  seqqual_mode = c("compatible", "compact"),
  threads = 1L,
  BPPARAM = BiocParallel::bpparam(),
  auto_threads = FALSE,
  use.names = FALSE,
  with.which_label = FALSE,
  include_unmapped = TRUE
)
```

## Arguments

- file:

  A BAM input. Supported values are:

  - a single BAM path (`character(1)`) or multiple BAM paths,

  - a
    [`Rsamtools::BamFile`](https://rdrr.io/pkg/Rsamtools/man/BamFile-class.html),

  - a
    [`Rsamtools::BamFileList`](https://rdrr.io/pkg/Rsamtools/man/BamFile-class.html).

- param:

  Optional
  [`Rsamtools::ScanBamParam`](https://rdrr.io/pkg/Rsamtools/man/ScanBamParam-class.html)
  (or a compatible list for lightweight use). The following fields are
  honored: `mapqFilter`, `flag`, `which`, `what`, and `tag`.

- what:

  Character vector of fields to return, similar to `scanBam(what=...)`.
  Supported fields are `qname`, `flag`, `rname`, `strand`, `pos`,
  `qwidth`, `mapq`, `cigar`, `mrnm`, `mpos`, `isize`, `seq`, `qual`.

- tag:

  Character vector of 2-letter tag names to extract.

- as:

  Output format:

  - `"DataFrame"`: returns
    [`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
    (default),

  - `"data.frame"`: returns base `data.frame`,

  - `"GAlignments"`: returns
    [`GenomicAlignments::GAlignments`](https://rdrr.io/pkg/GenomicAlignments/man/GAlignments-class.html),

  - `"GAlignmentPairs"`: returns
    [`GenomicAlignments::GAlignmentPairs`](https://rdrr.io/pkg/GenomicAlignments/man/GAlignmentPairs-class.html),

  - `"scanBam"`: returns a `scanBam()`-shaped list-of-lists.

- seqqual_mode:

  Controls representation of `seq`/`qual` when those fields are
  requested:

  - `"compatible"` (default): return character vectors matching
    `scanBam`-style expectations,

  - `"compact"`: return lower-level raw list-columns for
    faster/lower-overhead extraction. In compact mode, `seq` is returned
    as one `raw` vector per read containing BAM-native packed sequence
    bytes (two bases per byte), and `qual` is returned as one `raw`
    vector per read containing per-base Phred bytes. These are not plain
    character strings. `qwidth` is needed to decode compact `seq` back
    to base letters, and `qual` values of `255` correspond to missing
    quality values. This mode is currently supported for
    `as = "data.frame"` or `as = "DataFrame"`.

- threads:

  Requested number of OpenMP threads used for reading/decompression. May
  be capped when `auto_threads = TRUE`.

- BPPARAM:

  `BiocParallel` parameter used when `file` contains more than one BAM.
  Defaults to
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html).
  Set to `NULL` to force serial file processing.

- auto_threads:

  Logical; when `TRUE` and `BPPARAM` has multiple workers, BamScale
  adaptively avoids oversubscription by preserving higher per-file
  OpenMP thread counts when possible and reducing the number of
  concurrently active file workers before shrinking per-file threads.

- use.names:

  Passed to alignment object conversion. When `TRUE`, read names
  (`qname`) are used as object names.

- with.which_label:

  Logical; if `TRUE` and `param` includes `which`, an extra
  `which_label` column is returned.

- include_unmapped:

  Logical; whether unmapped records are retained (subject to
  `param$flag` constraints).

## Value

If `file` is length 1: one object in the format specified by `as`. If
`file` has length \> 1 (or is a `BamFileList`): a named list of outputs,
one per BAM file.

## Details

`bam_read()` is intentionally column-compatible with common BAM fields
used by Bioconductor workflows and can be used as a fast drop-in reader
before conversion to downstream classes.

Parallelism model:

- `BPPARAM` parallelizes across files (one file per BiocParallel
  worker).

- `threads` parallelizes within each file via OpenMP.

- Effective total concurrency is approximately
  `min(length(file), BiocParallel::bpnworkers(BPPARAM)) * threads`.

- If `auto_threads = TRUE` and `BPPARAM` has multiple workers, BamScale
  first limits the number of concurrently active workers to preserve the
  requested per-file thread count within the detected core budget, then
  caps per-file OpenMP threads only if a single file would still
  oversubscribe the machine.

Compatibility notes:

- Region filtering via `param$which` is supported as a sequential filter
  (not index-jump random access).

- Flag filtering uses `ScanBamFlag` semantics by converting logical flag
  requirements into required-set and required-unset bit masks.

- Tag values are returned as character columns. Scalar tags are scalar
  strings; `B` tags are comma-separated vectors.

- `seqqual_mode = "compact"` is optimized for throughput-oriented
  benchmarking and returns raw list-columns for `seq`/`qual`, not
  ordinary sequence or quality strings. In this representation, `seq`
  contains BAM-packed nucleotide bytes and `qual` contains raw Phred
  bytes. Compact output is intended for users who want to defer or avoid
  full string-materialization costs; use
  [`decode_compact_seq()`](https://cparsania.github.io/BamScale/reference/decode_compact_seq.md),
  [`decode_compact_qual()`](https://cparsania.github.io/BamScale/reference/decode_compact_qual.md),
  or
  [`decode_seqqual_compact()`](https://cparsania.github.io/BamScale/reference/decode_seqqual_compact.md)
  to decode compact output back to standard string form when needed.

- `"GAlignments"` and `"GAlignmentPairs"` output exclude unmapped
  records.

- `as = "scanBam"` returns a strict scan-like list-of-lists: without
  `param$which`, it returns one unnamed batch; with `param$which`, it
  returns one batch per range label (including empty ranges), with
  requested `what` fields and `tag` values under `$tag`. In this output
  mode, `seq` and `qual` are returned as
  [`Biostrings::DNAStringSet`](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  and
  [`Biostrings::PhredQuality`](https://rdrr.io/pkg/Biostrings/man/XStringQuality-class.html)
  for closer `scanBam()` compatibility.

## Examples

``` r
bam <- ompBAM::example_BAM("Unsorted")

# Familiar scanBam-like field selection
x <- bam_read(bam, what = c("qname", "flag", "rname", "pos", "cigar"))

# Include sequence + quality
y <- bam_read(bam, what = c("qname", "seq", "qual"), threads = 2)

# scanBam-shaped output
z <- bam_read(bam, what = c("qname", "flag"), tag = "NM", as = "scanBam")
```
