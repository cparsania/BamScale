test_that("bam_read returns selected fields", {
  skip_if_not_installed("ompBAM")

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

test_that("bam_read supports sequence and quality columns", {
  skip_if_not_installed("ompBAM")

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
  skip_if_not_installed("ompBAM")

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
  skip_if_not_installed("ompBAM")

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
  skip_if_not_installed("ompBAM")

  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(
    file = bam,
    what = c("qname", "seq", "qual"),
    as = "scanBam",
    threads = 1
  )

  if (requireNamespace("Biostrings", quietly = TRUE)) {
    expect_s4_class(x[[1]]$seq, "DNAStringSet")
    expect_s4_class(x[[1]]$qual, "PhredQuality")
  } else {
    expect_type(x[[1]]$seq, "character")
    expect_type(x[[1]]$qual, "character")
  }
})


test_that("scanBam mode batches by which labels and keeps empty ranges", {
  skip_if_not_installed("ompBAM")

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
  skip_if_not_installed("ompBAM")

  bam <- ompBAM::example_BAM("Unsorted")
  x <- bam_read(file = c(a = bam, b = bam), what = c("qname"), threads = 1)

  expect_type(x, "list")
  expect_equal(length(x), 2L)
  expect_equal(names(x), c("a", "b"))
})

test_that("bam_count runs and returns expected columns", {
  skip_if_not_installed("ompBAM")

  bam <- ompBAM::example_BAM("Unsorted")
  y <- bam_count(file = bam, threads = 1)

  expect_s3_class(y, "data.frame")
  expect_true(all(c("seqname", "seqlength", "count") %in% names(y)))
  expect_true(nrow(y) > 0)
})

test_that("GAlignments output is available when package is installed", {
  skip_if_not_installed("ompBAM")
  skip_if_not_installed("GenomicAlignments")

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
  skip_if_not_installed("ompBAM")
  skip_if_not_installed("GenomicAlignments")

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


test_that("auto_threads validates logical input", {
  expect_error(
    BamScale:::.bamscale_resolve_threads(threads = 1L, auto_threads = NA),
    "`auto_threads` must be TRUE or FALSE"
  )
})


test_that("auto_threads preserves per-file threads by reducing active workers first", {
  skip_if_not_installed("BiocParallel")

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
  skip_if_not_installed("BiocParallel")

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

  expected_workers <- min(8L, max(1L, floor(BamScale:::.bamscale_detect_cores() / 2L)))
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
  skip_if_not_installed("ompBAM")

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
  expect_identical(decoded$seq, compat$seq)
  expect_identical(decoded$qual, compat$qual)
})

test_that("compact seqqual decode matches Rsamtools sequence and quality output", {
  skip_if_not_installed("ompBAM")
  skip_if_not_installed("Rsamtools")
  skip_if_not_installed("Biostrings")

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
