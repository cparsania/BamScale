test_that("bam_read returns selected fields", {
  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(
    file = bam,
    what = c("qname", "flag", "rname", "pos", "mapq", "cigar"),
    threads = 1
  )

  expect_s4_class(x, "DataFrame")
  expect_true(all(c("qname", "flag", "rname", "pos", "mapq", "cigar") %in% colnames(x)))
  expect_true(nrow(x) > 0)
})

test_that("rname, mrnm and strand are scanBam-compatible factors", {
  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(
    file = bam,
    what = c("rname", "mrnm", "strand"),
    as = "data.frame",
    threads = 1
  )

  expect_s3_class(x$rname, "factor")
  expect_s3_class(x$mrnm, "factor")
  expect_s3_class(x$strand, "factor")
  expect_identical(levels(x$strand), c("+", "-", "*"))

  # Levels and integer codes match Rsamtools::scanBam exactly.
  rs <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = c("rname", "mrnm", "strand"))
  )[[1]]
  expect_identical(levels(x$rname), levels(rs$rname))
  expect_identical(as.integer(x$rname), as.integer(rs$rname))
  expect_identical(as.integer(x$mrnm), as.integer(rs$mrnm))
  expect_identical(as.character(x$strand), as.character(rs$strand))
})

test_that("numeric tags are returned with native scanBam-compatible types", {
  # Build a small BAM carrying integer (i), float (f) and string (Z) tags,
  # with some records omitting a tag so absent values decode to NA.
  tmp <- tempfile(fileext = ".sam")
  n <- 1500L
  hdr <- c("@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:100000")
  recs <- vapply(seq_len(n), function(i) {
    core <- c(sprintf("r%05d", i), "0", "chr1", as.character(10L + i), "60",
              "10M", "*", "0", "0", "ACGTACGTAC", "IIIIIIIIII")
    opt <- if (i %% 3L == 0L) {
      c(sprintf("AS:i:%d", 100L - (i %% 40L)), sprintf("RG:Z:g%d", i %% 3L))
    } else {
      c(sprintf("NM:i:%d", i %% 5L), sprintf("AS:i:%d", 100L - (i %% 40L)),
        sprintf("XF:f:%.3f", (i %% 5L) + 0.25), sprintf("RG:Z:g%d", i %% 3L))
    }
    paste(c(core, opt), collapse = "\t")
  }, character(1))
  writeLines(c(hdr, recs), tmp)
  bam <- Rsamtools::asBam(tmp, tempfile(), overwrite = TRUE, indexDestination = FALSE)

  tags <- c("NM", "AS", "XF", "RG")
  x <- bam_read(bam, what = "qname", tag = tags, as = "data.frame", threads = 2)
  rs <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = "qname", tag = tags)
  )[[1]]$tag

  expect_type(x$NM, "integer")
  expect_type(x$AS, "integer")
  expect_type(x$XF, "double")
  expect_type(x$RG, "character")

  expect_identical(x$NM, rs$NM)
  expect_identical(x$AS, rs$AS)
  expect_equal(x$XF, rs$XF)
  expect_identical(x$RG, rs$RG)

  # Absent NM values decode to NA in the same positions as scanBam.
  expect_identical(is.na(x$NM), is.na(rs$NM))
})

test_that("bam_read supports sequence and quality columns", {
  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(
    file = bam,
    what = c("qname", "seq", "qual"),
    as = "data.frame",
    threads = 1
  )

  expect_s3_class(x, "data.frame")
  expect_true(all(c("qname", "seq", "qual") %in% colnames(x)))
  expect_true(nrow(x) > 0)
})

test_that("list-style param filtering works", {
  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(
    file = bam,
    param = list(mapqFilter = 20L, flag_unset = 4L),
    what = c("flag", "mapq"),
    as = "data.frame",
    threads = 1
  )

  expect_true(all(x$mapq >= 20L))
  expect_true(all(bitwAnd(as.integer(x$flag), 4L) == 0L))
})

test_that("scanBam mode returns scan-like list structure", {
  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(
    file = bam,
    what = c("qname", "flag", "rname", "pos"),
    tag = c("NM"),
    as = "scanBam",
    threads = 1
  )

  expect_type(x, "list")
  expect_true(length(x) >= 1L)
  expect_null(names(x))
  expect_type(x[[1]], "list")
  expect_true(all(c("qname", "flag", "rname", "pos", "tag") %in% names(x[[1]])))
  expect_true("NM" %in% names(x[[1]]$tag))
})


test_that("scanBam seq/qual uses Biostrings classes when available", {
  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(
    file = bam,
    what = c("qname", "seq", "qual"),
    as = "scanBam",
    threads = 1
  )

  expect_s4_class(x[[1]]$seq, "DNAStringSet")
  expect_s4_class(x[[1]]$qual, "PhredQuality")
})


test_that("seq and qual are byte-identical to Rsamtools::scanBam", {
  bam <- ompBAM::example_BAM("Unsorted")
  rs <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = c("seq", "qual"))
  )[[1]]

  # data.frame mode now returns native Biostrings columns (previously character),
  # built in C directly from the reader's shared byte buffer. For any result
  # scanBam returns as a single pool buffer (small/moderate files) this is
  # byte-for-byte identical(); larger inputs where scanBam fragments its pool
  # remain content-identical but not object-identical() by design.
  df <- bam_read(bam, what = c("seq", "qual"), as = "data.frame", threads = 1)
  expect_s4_class(df$seq, "DNAStringSet")
  expect_s4_class(df$qual, "PhredQuality")
  expect_identical(df$seq, rs$seq)
  expect_identical(df$qual, rs$qual)

  # The multi-threaded merge preserves file order, so the result is unchanged.
  df_mt <- bam_read(bam, what = c("seq", "qual"), as = "data.frame", threads = 4)
  expect_identical(df_mt$seq, rs$seq)
  expect_identical(df_mt$qual, rs$qual)

  # scanBam output mode carries the same S4 objects through unchanged.
  sb <- bam_read(bam, what = c("seq", "qual"), as = "scanBam", threads = 1)[[1]]
  expect_identical(sb$seq, rs$seq)
  expect_identical(sb$qual, rs$qual)
})


test_that("compatible seq/qual builds correctly on fresh SnowParam SOCK workers", {
  # Regression guard: the C XStringSet builder resolves IRanges/XVector
  # C-callables at runtime and instantiates Biostrings S4 classes. A fresh SOCK
  # worker must have those namespaces loaded (via .onLoad) or it fails with
  # "function '_new_IRanges' not provided" / "DNAStringSet is not a defined class".
  skip_if_not_installed("BiocParallel")
  bam <- ompBAM::example_BAM("Unsorted")
  rs <- Rsamtools::scanBam(bam, param = Rsamtools::ScanBamParam(what = c("seq", "qual")))[[1]]
  bp <- BiocParallel::SnowParam(2L, type = "SOCK")
  res <- BiocParallel::bplapply(1:2, function(i) {
    x <- BamScale::bam_read(bam, what = c("seq", "qual"), as = "data.frame", threads = 1)
    list(seq = as.character(x$seq), qual = as.character(x$qual),
         seq_cls = class(x$seq), qual_cls = class(x$qual))
  }, BPPARAM = bp)

  expect_identical(as.character(res[[1]]$seq_cls), "DNAStringSet")
  expect_identical(as.character(res[[1]]$qual_cls), "PhredQuality")
  expect_identical(res[[1]]$seq, as.character(rs$seq))
  expect_identical(res[[1]]$qual, as.character(rs$qual))
})


test_that("scanBam mode batches by which labels and keeps empty ranges", {
  bam <- ompBAM::example_BAM("Unsorted")
  seed <- bam_read(
    file = bam,
    what = c("rname", "pos"),
    as = "data.frame",
    threads = 1
  )

  first_rname <- as.character(seed[["rname"]][1])
  first_pos <- as.integer(seed[["pos"]][1])

  x <- bam_read(
    file = bam,
    param = list(
      which = data.frame(
        seqname = c(first_rname, first_rname),
        start = c(max(1L, first_pos - 2L), first_pos + 1000000000L),
        end = c(first_pos + 2L, first_pos + 1000000010L),
        label = c("hit", "empty"),
        stringsAsFactors = FALSE
      )
    ),
    what = c("qname", "flag"),
    as = "scanBam",
    threads = 1
  )

  expect_type(x, "list")
  expect_equal(names(x), c("hit", "empty"))
  expect_equal(length(x), 2L)
  expect_true(length(x[["hit"]]$qname) >= 1L)
  expect_length(x[["empty"]]$qname, 0L)
  expect_length(x[["empty"]]$flag, 0L)
})

test_that("multi-file input returns named list", {
  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(file = c(a = bam, b = bam), what = c("qname"), threads = 1)

  expect_type(x, "list")
  expect_equal(length(x), 2L)
  expect_equal(names(x), c("a", "b"))
})

test_that("bam_count runs and returns expected columns", {
  bam <- ompBAM::example_BAM("Unsorted")
  y <- bam_count(file = bam, threads = 1)

  expect_s3_class(y, "data.frame")
  expect_true(all(c("seqname", "seqlength", "count") %in% names(y)))
  expect_true(nrow(y) > 0)
})

test_that("fragment_sizes is byte-identical to table(abs(scanBam isize))", {
  bam <- ompBAM::example_BAM("Unsorted")
  FLAG <- Rsamtools::scanBamFlag(
    isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE,
    isNotPassingQualityControls = FALSE
  )
  fs <- fragment_sizes(
    bam,
    param = Rsamtools::ScanBamParam(flag = FLAG),
    threads = 1
  )
  expect_s3_class(fs, "data.frame")
  expect_identical(names(fs), c("fragment_size", "count"))
  expect_false(is.unsorted(fs$fragment_size))

  # Reference: scanBam materialises every insert size, then table(abs()).
  # scanBam reports isize = NA for mate-unmapped reads, which table() drops;
  # fragment_sizes(drop_mate_unmapped = TRUE) replicates that exactly.
  iz <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = "isize", flag = FLAG)
  )[[1]]$isize
  ref <- table(abs(iz))
  ref_df <- data.frame(
    fragment_size = as.integer(names(ref)),
    ref = as.numeric(ref)
  )
  m <- merge(fs, ref_df, by = "fragment_size", all = TRUE)
  m[is.na(m)] <- 0
  expect_equal(m$count, m$ref)
  expect_equal(sum(fs$count), sum(!is.na(iz)))
})

test_that("fragment_sizes(drop_mate_unmapped = FALSE) keeps mate-unmapped reads", {
  bam <- ompBAM::example_BAM("Unsorted")
  keep <- fragment_sizes(bam, threads = 1, drop_mate_unmapped = FALSE)
  drop <- fragment_sizes(bam, threads = 1, drop_mate_unmapped = TRUE)
  expect_gte(sum(keep$count), sum(drop$count))
})

test_that("mapq_dist is byte-identical to table(scanBam mapq)", {
  bam <- ompBAM::example_BAM("Unsorted")
  FLAG <- Rsamtools::scanBamFlag(
    isSecondaryAlignment = FALSE,
    isUnmappedQuery = FALSE,
    isNotPassingQualityControls = FALSE
  )
  md <- mapq_dist(
    bam,
    param = Rsamtools::ScanBamParam(flag = FLAG),
    threads = 1
  )
  expect_s3_class(md, "data.frame")
  expect_identical(names(md), c("mapq", "count"))
  expect_false(is.unsorted(md$mapq))

  mq <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = "mapq", flag = FLAG)
  )[[1]]$mapq
  ref <- table(mq)
  ref_df <- data.frame(
    mapq = as.integer(names(ref)),
    ref = as.numeric(ref)
  )
  m <- merge(md, ref_df, by = "mapq", all = TRUE)
  m[is.na(m)] <- 0
  expect_equal(m$count, m$ref)
  expect_equal(sum(md$count), sum(!is.na(mq)))
})

test_that("fragment_sizes and mapq_dist return named lists for multiple files", {
  bam <- ompBAM::example_BAM("Unsorted")
  fs <- fragment_sizes(file = c(a = bam, b = bam), threads = 1)
  expect_type(fs, "list")
  expect_equal(names(fs), c("a", "b"))
  expect_s3_class(fs[[1]], "data.frame")

  md <- mapq_dist(file = c(a = bam, b = bam), threads = 1)
  expect_type(md, "list")
  expect_equal(names(md), c("a", "b"))
  expect_s3_class(md[[1]], "data.frame")
})

test_that("bam_coverage is byte-identical to coverage(readGAlignments)", {
  skip_if_not_installed("GenomicAlignments")
  bam <- ompBAM::example_BAM("Unsorted")
  got <- bam_coverage(bam, threads = 1)
  ref <- GenomicAlignments::coverage(GenomicAlignments::readGAlignments(bam))

  # Output must be an *uncompressed* SimpleRleList, integer runs, header-order
  # names -- each independently required for identical() to coverage().
  expect_true(methods::is(got, "SimpleRleList"))
  expect_false(methods::is(got, "CompressedRleList"))
  expect_identical(names(got), names(ref))
  expect_identical(typeof(S4Vectors::runValue(got[[1]])), "integer")
  expect_identical(got, ref)
})

test_that("bam_coverage keeps secondary/supplementary reads (readGAlignments filter)", {
  # Only unmapped (0x4) reads are dropped; the default must match readGAlignments,
  # not a stricter flag filter. Passing that exact flag must not change the result.
  skip_if_not_installed("GenomicAlignments")
  bam <- ompBAM::example_BAM("Unsorted")
  default_cov <- bam_coverage(bam, threads = 1)
  explicit_cov <- bam_coverage(
    bam, threads = 1,
    param = Rsamtools::ScanBamParam(
      flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)
    )
  )
  expect_identical(default_cov, explicit_cov)
})

test_that("bam_coverage returns a named list for multiple files", {
  bam <- ompBAM::example_BAM("Unsorted")
  cv <- bam_coverage(file = c(a = bam, b = bam), threads = 1)
  expect_type(cv, "list")
  expect_equal(names(cv), c("a", "b"))
  expect_true(methods::is(cv[[1]], "SimpleRleList"))
})

test_that("bam_coverage_bigwig writes a bigWig with values identical to coverage()", {
  # rtracklayer's bigWig import goes through the UCSC udc layer, which cannot
  # parse Windows paths ("Unrecognized protocol C") -- bigWig I/O via
  # rtracklayer is unsupported on Windows. The parallel-vs-serial byte-identity
  # test below still runs there.
  skip_on_os("windows")
  skip_if_not_installed("rtracklayer")
  skip_if_not_installed("GenomicAlignments")
  bam <- ompBAM::example_BAM("Unsorted")
  out <- tempfile(fileext = ".bw")
  on.exit(unlink(out), add = TRUE)

  ret <- bam_coverage_bigwig(bam, out, threads = 1)
  expect_identical(as.character(ret), out)
  expect_true(file.exists(out) && file.info(out)$size > 0)

  imp <- rtracklayer::import(out, as = "RleList")
  ref <- GenomicAlignments::coverage(GenomicAlignments::readGAlignments(bam))
  common <- intersect(names(ref), names(imp))
  expect_gt(length(common), 0)
  for (nm in common) {
    expect_equal(as.numeric(imp[[nm]]), as.numeric(ref[[nm]]))
  }
  auc <- function(rl) sum(vapply(rl, function(x)
    sum(as.numeric(S4Vectors::runValue(x)) * as.numeric(S4Vectors::runLength(x))),
    numeric(1)))
  expect_equal(auc(imp[common]), auc(ref))
})

test_that("bam_coverage_bigwig errors when outfile length != file length", {
  bam <- ompBAM::example_BAM("Unsorted")
  expect_error(
    bam_coverage_bigwig(c(bam, bam), tempfile(fileext = ".bw"), threads = 1),
    "same length"
  )
})

test_that("bam_coverage_bigwig parallel output is byte-identical to serial", {
  # Parallel block compression must produce a bit-for-bit identical bigWig
  # (compress2 is deterministic; blocks are written in the same order). On
  # Windows the UCRT's stream-position semantics on buffered update streams can
  # encode index offsets slightly differently between the two write patterns, so
  # there we assert the invariants that survive -- identical size and identical
  # header -- and leave bit-identity to the POSIX platforms where it is proven.
  bam <- ompBAM::example_BAM("Unsorted")
  op <- tempfile(fileext = ".bw")
  os <- tempfile(fileext = ".bw")
  on.exit(unlink(c(op, os)), add = TRUE)
  bam_coverage_bigwig(bam, op, threads = 2, parallel = TRUE)
  bam_coverage_bigwig(bam, os, threads = 2, parallel = FALSE)
  if (identical(.Platform$OS.type, "windows")) {
    expect_identical(file.info(op)$size, file.info(os)$size)
    expect_identical(readBin(op, "raw", 64L), readBin(os, "raw", 64L))
  } else {
    expect_identical(
      readBin(op, "raw", file.info(op)$size),
      readBin(os, "raw", file.info(os)$size)
    )
  }
})

test_that("GA fast path is identical to readGAlignments (core fields)", {
  skip_if_not_installed("GenomicAlignments")
  bam <- ompBAM::example_BAM("Unsorted")
  ref <- GenomicAlignments::readGAlignments(Rsamtools::BamFile(bam))

  ga1 <- bam_read(bam, what = c("rname", "pos", "cigar", "strand"),
                  as = "GAlignments", threads = 1)
  expect_identical(ga1, ref)

  # threads > 1 exercises the tid-ordered merge and RLE boundary coalescing
  ga4 <- bam_read(bam, what = c("rname", "pos", "cigar", "strand"),
                  as = "GAlignments", threads = 4)
  expect_identical(ga4, ref)
})

test_that("GA fast path mcols match readGAlignments param what", {
  skip_if_not_installed("GenomicAlignments")
  bam <- ompBAM::example_BAM("Unsorted")

  ref <- GenomicAlignments::readGAlignments(
    Rsamtools::BamFile(bam),
    param = Rsamtools::ScanBamParam(what = c("flag", "qname"))
  )
  ga <- bam_read(bam,
                 what = c("rname", "pos", "cigar", "strand", "flag", "qname"),
                 as = "GAlignments", threads = 2)
  expect_identical(ga, ref)

  # reversed mcols order must be preserved
  ref2 <- GenomicAlignments::readGAlignments(
    Rsamtools::BamFile(bam),
    param = Rsamtools::ScanBamParam(what = c("qname", "flag"))
  )
  ga2 <- bam_read(bam,
                  what = c("rname", "pos", "cigar", "strand", "qname", "flag"),
                  as = "GAlignments", threads = 2)
  expect_identical(ga2, ref2)
})

test_that("GA fast path honors use.names and filters like readGAlignments", {
  skip_if_not_installed("GenomicAlignments")
  bam <- ompBAM::example_BAM("Unsorted")

  ref <- GenomicAlignments::readGAlignments(
    Rsamtools::BamFile(bam),
    use.names = TRUE,
    param = Rsamtools::ScanBamParam(what = c("flag", "qname"))
  )
  ga <- bam_read(bam,
                 what = c("rname", "pos", "cigar", "strand", "flag", "qname"),
                 as = "GAlignments", use.names = TRUE, threads = 2)
  expect_identical(ga, ref)

  FLAG <- Rsamtools::scanBamFlag(isMinusStrand = FALSE)
  reff <- GenomicAlignments::readGAlignments(
    Rsamtools::BamFile(bam),
    param = Rsamtools::ScanBamParam(flag = FLAG, mapqFilter = 10)
  )
  gaf <- bam_read(bam, what = c("rname", "pos", "cigar", "strand"),
                  as = "GAlignments", threads = 2,
                  param = Rsamtools::ScanBamParam(flag = FLAG, mapqFilter = 10))
  expect_identical(gaf, reff)
})

test_that("GA fast path and slow path produce identical objects", {
  bam <- ompBAM::example_BAM("Unsorted")
  fields <- c("rname", "pos", "cigar", "strand", "flag", "mapq")
  fast <- bam_read(bam, what = fields, as = "GAlignments", threads = 2)
  old <- options(BamScale.ga_fastpath = FALSE)
  on.exit(options(old), add = TRUE)
  slow <- bam_read(bam, what = fields, as = "GAlignments", threads = 2)
  expect_identical(fast, slow)
})

test_that("GA fast path survives gctorture (CHARSXP cache GC safety)", {
  bam <- ompBAM::example_BAM("Unsorted")
  gctorture(TRUE)
  ga <- tryCatch(
    bam_read(bam, what = c("rname", "pos", "cigar", "strand"),
             as = "GAlignments", threads = 1),
    finally = gctorture(FALSE)
  )
  expect_s4_class(ga, "GAlignments")
  expect_true(length(ga) > 0)
})

test_that("GAlignments output is available when package is installed", {
  bam <- ompBAM::example_BAM("Unsorted")
  g <- bam_read(
    file = bam,
    what = c("qname", "rname", "pos", "cigar", "strand", "flag"),
    as = "GAlignments",
    threads = 1
  )

  expect_s4_class(g, "GAlignments")
})

test_that("GAlignmentPairs output is available when package is installed", {
  bam <- ompBAM::example_BAM("Unsorted")
  gp <- bam_read(
    file = bam,
    what = c("qname", "flag", "rname", "pos", "cigar", "strand"),
    as = "GAlignmentPairs",
    include_unmapped = FALSE,
    threads = 1
  )

  expect_s4_class(gp, "GAlignmentPairs")
})


test_that("GAlignments carries header seqlengths and gives identical coverage()", {
  bam <- ompBAM::example_BAM("Unsorted")

  g <- bam_read(
    file = bam,
    what = c("rname", "pos", "cigar", "strand"),
    as = "GAlignments",
    threads = 1
  )
  std <- GenomicAlignments::readGAlignments(bam)

  # Header sequence lengths must be populated (not NA) and match the standard
  # reader, so coverage()/export.bw() see true chromosome sizes rather than
  # max-end truncation.
  expect_identical(GenomeInfoDb::seqlevels(g), GenomeInfoDb::seqlevels(std))
  expect_identical(GenomeInfoDb::seqlengths(g), GenomeInfoDb::seqlengths(std))
  expect_false(any(is.na(GenomeInfoDb::seqlengths(g))))

  # coverage() is byte-identical to the standard Bioconductor path -- the
  # equivalence anchor for the coverage -> bigWig workflow benchmark.
  expect_identical(
    GenomicAlignments::coverage(g),
    GenomicAlignments::coverage(std)
  )
})


test_that("auto_threads validates logical input", {
  expect_error(
    BamScale:::.bamscale_resolve_parallel_plan(threads = 1L, auto_threads = NA),
    "`auto_threads` must be TRUE or FALSE"
  )
})


test_that("auto_threads preserves per-file threads by reducing active workers first", {
  bp <- tryCatch(
    BiocParallel::SnowParam(workers = 8L, type = "SOCK", progressbar = FALSE),
    error = function(e) NULL
  )
  if (is.null(bp)) {
    skip("SnowParam could not be initialized in this environment")
  }
  on.exit(try(BiocParallel::bpstop(bp), silent = TRUE), add = TRUE)

  plan <- BamScale:::.bamscale_resolve_parallel_plan(
    threads = 64L,
    BPPARAM = bp,
    auto_threads = TRUE,
    n_files = 8L
  )

  expected_threads <- max(1L, min(64L, BamScale:::.bamscale_detect_cores()))
  expect_equal(plan$threads, expected_threads)
  expect_equal(plan$bp_workers, 1L)
})


test_that("auto_threads keeps multiple workers when requested per-file threads are small", {
  bp <- tryCatch(
    BiocParallel::SnowParam(workers = 8L, type = "SOCK", progressbar = FALSE),
    error = function(e) NULL
  )
  if (is.null(bp)) {
    skip("SnowParam could not be initialized in this environment")
  }
  on.exit(try(BiocParallel::bpstop(bp), silent = TRUE), add = TRUE)

  plan <- BamScale:::.bamscale_resolve_parallel_plan(
    threads = 2L,
    BPPARAM = bp,
    auto_threads = TRUE,
    n_files = 8L
  )

  workers_eff <- min(BamScale:::.bamscale_bpparam_workers(bp), 8L)
  expected_workers <- min(
    workers_eff,
    max(1L, floor(BamScale:::.bamscale_detect_cores() / 2L))
  )
  expect_equal(plan$threads, 2L)
  expect_equal(plan$bp_workers, expected_workers)
  expect_lte(plan$threads * plan$bp_workers, max(1L, BamScale:::.bamscale_detect_cores()))
})

test_that("compact decode helpers decode synthetic seq and qual correctly", {
  seq_raw <- list(
    as.raw(c(0x12, 0x48)), # A C G T
    as.raw(c(0xFF)),       # N N
    raw()
  )
  qwidth <- c(4L, 2L, 0L)

  qual_raw <- list(
    as.raw(c(0L, 1L, 2L, 3L)),
    as.raw(c(255L, 255L)),
    raw()
  )

  expect_identical(
    decode_compact_seq(seq_raw, qwidth),
    c("ACGT", "NN", "")
  )

  expect_identical(
    decode_compact_qual(qual_raw),
    c("!\"#$", "*", "")
  )
})

test_that("compact decode helpers validate malformed inputs", {
  expect_error(
    decode_compact_seq(list(as.raw(c(0x12))), integer()),
    "`seq` and `qwidth` must have the same length"
  )

  expect_error(
    decode_compact_seq("not-a-list", 1L),
    "`seq` must be a list of raw vectors"
  )

  expect_error(
    decode_compact_qual("not-a-list"),
    "`qual` must be a list of raw vectors"
  )

  expect_error(
    decode_seqqual_compact(data.frame(seq = I(list(as.raw(c(0x12)))))),
    "requires a `qwidth` column"
  )
})

test_that("decode_seqqual_compact preserves non-seqqual columns", {
  x <- data.frame(
    qname = "r1",
    qwidth = 4L,
    mapq = 42L,
    stringsAsFactors = FALSE
  )
  x$seq <- I(list(as.raw(c(0x12, 0x48))))
  x$qual <- I(list(as.raw(c(0L, 1L, 2L, 3L))))

  y <- decode_seqqual_compact(x)

  expect_identical(y$qname, x$qname)
  expect_identical(y$qwidth, x$qwidth)
  expect_identical(y$mapq, x$mapq)
  expect_identical(y$seq, "ACGT")
  expect_identical(y$qual, "!\"#$")
})

test_that("compact seqqual round-trips to BamScale compatible output", {
  bam <- ompBAM::example_BAM("Unsorted")

  compat <- bam_read(
    file = bam,
    what = c("qname", "qwidth", "seq", "qual"),
    as = "data.frame",
    seqqual_mode = "compatible",
    threads = 1
  )

  compact <- bam_read(
    file = bam,
    what = c("qname", "qwidth", "seq", "qual"),
    as = "data.frame",
    seqqual_mode = "compact",
    threads = 1
  )

  expect_true(is.list(compact$seq))
  expect_true(is.list(compact$qual))
  expect_true(all(vapply(compact$seq, is.raw, logical(1))))
  expect_true(all(vapply(compact$qual, is.raw, logical(1))))

  decoded <- decode_seqqual_compact(compact)

  expect_identical(decoded$qname, compat$qname)
  expect_identical(decoded$qwidth, compat$qwidth)
  # Compatible mode now returns native Biostrings columns, so compare the decoded
  # compact strings against their character form.
  expect_identical(decoded$seq, as.character(compat$seq))
  expect_identical(decoded$qual, as.character(compat$qual))
})

test_that("compact seqqual decode matches Rsamtools sequence and quality output", {
  bam <- ompBAM::example_BAM("Unsorted")

  compact <- bam_read(
    file = bam,
    what = c("qname", "qwidth", "seq", "qual"),
    as = "data.frame",
    seqqual_mode = "compact",
    threads = 1
  )
  decoded <- decode_seqqual_compact(compact)

  rs <- Rsamtools::scanBam(
    bam,
    param = Rsamtools::ScanBamParam(what = c("qname", "seq", "qual"))
  )[[1]]

  expect_identical(decoded$qname, rs$qname)
  expect_identical(decoded$seq, as.character(rs$seq))
  expect_identical(decoded$qual, as.character(rs$qual))
})
