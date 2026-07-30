# Changelog

## BamScale 0.99.13

### New features

- On the compatible read path, `seq` and `qual` are now built directly
  in C as `DNAStringSet` / `PhredQuality` objects from the packed BAM
  byte buffer (via the IRanges/XVector C-callable constructors),
  avoiding a character round-trip. Output remains byte-identical to
  [`Rsamtools::scanBam`](https://rdrr.io/pkg/Rsamtools/man/scanBam.html).
- Added an end-to-end workflow benchmark harness under
  `inst/benchmarks/` (coverage → bigWig and ATAC fragment-size QC) and a
  self-contained `benchmark-results` vignette.

### Bug fixes

- `GAlignments` output now carries the BAM-header sequence lengths in
  its `Seqinfo`. Previously the seqlengths were `NA`, so `coverage()`
  extended each seqlevel only to its maximum end and `export.bw()` wrote
  truncated chromosome sizes, diverging from
  [`GenomicAlignments::readGAlignments`](https://rdrr.io/pkg/GenomicAlignments/man/readGAlignments.html).
- Force-load the C-callable providers (`S4Vectors`, `IRanges`,
  `XVector`, `Biostrings`) and `GenomicAlignments` in `.onLoad`, so the
  compatible seq/qual and `GAlignments` paths resolve their constructors
  and classes on fresh `SnowParam` SOCK workers.

### Documentation

- Refreshed the README with the latest benchmark results (read-pattern
  and end-to-end workflow speedups) and Bioconductor installation
  instructions.
