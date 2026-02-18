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


test_that("auto_threads caps per-file OpenMP threads using BPPARAM workers", {
  skip_if_not_installed("BiocParallel")

  bp <- tryCatch(
    BiocParallel::SnowParam(workers = 2L, type = "SOCK", progressbar = FALSE),
    error = function(e) NULL
  )
  if (is.null(bp)) {
    skip("SnowParam could not be initialized in this environment")
  }
  on.exit(try(BiocParallel::bpstop(bp), silent = TRUE), add = TRUE)

  resolved <- BamScale:::.bamscale_resolve_threads(
    threads = 64L,
    BPPARAM = bp,
    auto_threads = TRUE,
    n_files = 4L
  )

  expected <- max(1L, min(64L, floor(BamScale:::.bamscale_detect_cores() / 2L)))
  expect_equal(resolved, expected)
})
