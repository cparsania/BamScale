# Changelog

## BamScale 0.99.14

### New features

- Four new exported functions compute common BAM summaries entirely
  inside the multithreaded C++ reader (per-thread accumulators folded in
  the OpenMP region; only the compact result crosses into R – no
  per-read R objects are ever materialised), each verified
  byte-identical to its standard Bioconductor equivalent:
  - [`fragment_sizes()`](https://cparsania.github.io/BamScale/reference/fragment_sizes.md):
    paired-end fragment-size (insert-size) distribution; identical to
    `table(abs(scanBam(isize)))`. `drop_mate_unmapped = TRUE` (default)
    replicates `scanBam`’s `isize = NA` rule for mate-unmapped reads.
  - [`mapq_dist()`](https://cparsania.github.io/BamScale/reference/mapq_dist.md):
    mapping-quality distribution; identical to `table(scanBam(mapq))`
    (MAPQ 255 reported as 255, where scanBam uses `NA`).
  - [`bam_coverage()`](https://cparsania.github.io/BamScale/reference/bam_coverage.md):
    per-base coverage as a `SimpleRleList`;
    [`identical()`](https://rdrr.io/r/base/identical.html) to
    `GenomicAlignments::coverage(readGAlignments(...))`, including the
    CIGAR conventions (`M/=/X/D` cover, `N` gaps,
    `drop.D.ranges = FALSE`) and the readGAlignments record filter (only
    unmapped reads dropped).
  - [`bam_coverage_bigwig()`](https://cparsania.github.io/BamScale/reference/bam_coverage_bigwig.md):
    single-pass BAM -\> coverage bigWig written natively via the bundled
    libBigWig – no R-side coverage object, no `rtracklayer::export.bw`
    round-trip. Data and zoom blocks are compressed across the OpenMP
    threads (`parallel = TRUE`, default); the parallel output is
    bit-for-bit identical to the serial writer. A `compress_level`
    argument exposes the zlib deflate level, and `verbose = TRUE` prints
    a per-phase timing breakdown.
- `bam_read(as = "GAlignments")` now uses a dedicated fast path that
  assembles the `GAlignments` slots in C++: seqnames/strand accumulated
  as run-length pairs inside the parallel region, the cigar column
  materialised through a CHARSXP cache that exploits CIGAR-string
  redundancy, and qname decoded only when actually requested. Output is
  [`identical()`](https://rdrr.io/r/base/identical.html) to
  [`GenomicAlignments::readGAlignments()`](https://rdrr.io/pkg/GenomicAlignments/man/readGAlignments.html)
  (verified at 226M-read scale) and to the previous construction path.
  `options(BamScale.ga_fastpath = FALSE)` restores the old path.
- Vendored libBigWig 0.4.8 (MIT, Devon Ryan) under `src/libBigWig/`
  (built with `-DNOCURL`, local files only). Local modifications, each
  commented in `bwWrite.c`: a tunable zlib deflate level, deferred
  parallel block compression, and zoom-level construction from in-memory
  coverage runs instead of re-reading the just-written data section.

### Performance

- Fixed a quadratic reallocation pattern in the reader’s batch merge:
  per-batch exact-capacity `reserve()` calls copied all accumulated data
  on every batch. Buffers now grow geometrically. This affects every
  read path – e.g. on a 229M-read BAM at 48 threads,
  `as = "GAlignments"` dropped from ~300 s to ~56 s, and plain
  data-frame reads improved 2.8-4.4x.
- CIGAR strings are formatted with an allocation-free digit writer
  instead of per-op `std::to_string`.
- Representative large-file results (229M-read ATAC BAM, 48 threads,
  warm cache):
  [`bam_coverage()`](https://cparsania.github.io/BamScale/reference/bam_coverage.md)
  7.2x vs `coverage(readGAlignments())`;
  [`bam_coverage_bigwig()`](https://cparsania.github.io/BamScale/reference/bam_coverage_bigwig.md)
  2.6x vs megadepth at equal compression;
  [`fragment_sizes()`](https://cparsania.github.io/BamScale/reference/fragment_sizes.md)
  5.4x vs `scanBam` + [`table()`](https://rdrr.io/r/base/table.html);
  `as = "GAlignments"` 4.8x vs `readGAlignments()` (and faster even
  single-threaded).

### Bug fixes

- Removed a long-standing `^src/Makevars$` entry from `.Rbuildignore`
  (dating to the pkgdown setup): built tarballs – including the
  Bioconductor builds – shipped without `src/Makevars` and therefore
  compiled WITHOUT OpenMP (and without the explicit `-lz`), silently
  running single-threaded when installed from a tarball. In-place
  `R CMD INSTALL .` was unaffected, which is why local installs always
  threaded. Tarball builds now carry the full build configuration
  (OpenMP, zlib, `-DNOCURL`, and the bundled libBigWig objects).

### Benchmarks

- `run_workflow_benchmark.R` gains first-class arms for the four new
  functions (`--include-fastcov/-bigwig/-fragsize/-mapq`), each behind
  an [`identical()`](https://rdrr.io/r/base/identical.html) correctness
  gate, with an optional megadepth context arm (`--megadepth-bin`).
- `run_server_benchmark.R`: BiocParallel PSOCK clusters are now started
  and their namespaces pre-loaded outside the timed region in all
  multi-file arms (previously cluster spin-up was charged to every
  iteration, diluting comparator ratios at high worker counts); added
  `--include-single` so single and multi scenarios can run in separate
  stages.
- New reproducible manuscript-benchmark tooling: `select_encode_atac.R`
  (ENCODE portal selection -\> accession/md5 manifest), manifest-driven
  `download_atac_data.R` with md5 verification, `samtools_reference.sh`,
  `merge_runs.R`, and the staged `run_manuscript_final.sh` orchestrator.

### Tests

- New testthat gates: byte-identity of
  [`fragment_sizes()`](https://cparsania.github.io/BamScale/reference/fragment_sizes.md)
  /
  [`mapq_dist()`](https://cparsania.github.io/BamScale/reference/mapq_dist.md)
  /
  [`bam_coverage()`](https://cparsania.github.io/BamScale/reference/bam_coverage.md)
  against their Rsamtools/GenomicAlignments equivalents;
  [`bam_coverage_bigwig()`](https://cparsania.github.io/BamScale/reference/bam_coverage_bigwig.md)
  value-identity after re-import plus bit-for-bit parallel-vs-serial
  output; six `GAlignments` fast-path identity gates (threads 1/4, mcols
  order, `use.names`, flag/mapq filters, fast == slow, `gctorture`).

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
