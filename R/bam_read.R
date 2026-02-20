.bamscale_default_what <- c("qname", "flag", "rname", "pos", "mapq", "cigar")
.bamscale_supported_what <- c(
    "qname", "flag", "rname", "strand", "pos", "qwidth", "mapq", "cigar",
    "mrnm", "mpos", "isize", "seq", "qual"
)
.bamscale_tag_na_sentinel <- "__BAMSCALE_TAG_NA__"

.bamscale_field_bits <- c(
    qname = bitwShiftL(1L, 0L),
    flag = bitwShiftL(1L, 1L),
    rname = bitwShiftL(1L, 2L),
    strand = bitwShiftL(1L, 3L),
    pos = bitwShiftL(1L, 4L),
    qwidth = bitwShiftL(1L, 5L),
    mapq = bitwShiftL(1L, 6L),
    cigar = bitwShiftL(1L, 7L),
    mrnm = bitwShiftL(1L, 8L),
    mpos = bitwShiftL(1L, 9L),
    isize = bitwShiftL(1L, 10L),
    seq = bitwShiftL(1L, 11L),
    qual = bitwShiftL(1L, 12L)
)


#' Fast BAM reading with Bioconductor-compatible arguments
#'
#' `bam_read()` is a multithreaded sequential BAM reader built on top of
#' `ompBAM`. The interface is designed to be familiar to users of
#' `Rsamtools::scanBam()`, `GenomicAlignments::readGAlignments()`, and
#' `GenomicAlignments::readGAlignmentPairs()`.
#'
#' @param file A BAM input. Supported values are:
#' - a single BAM path (`character(1)`) or multiple BAM paths,
#' - a `Rsamtools::BamFile`,
#' - a `Rsamtools::BamFileList`.
#' @param param Optional `Rsamtools::ScanBamParam` (or a compatible list for
#'   lightweight use). The following fields are honored:
#'   `mapqFilter`, `flag`, `which`, `what`, and `tag`.
#' @param what Character vector of fields to return, similar to
#'   `scanBam(what=...)`. Supported fields are
#'   `qname`, `flag`, `rname`, `strand`, `pos`, `qwidth`, `mapq`, `cigar`,
#'   `mrnm`, `mpos`, `isize`, `seq`, `qual`.
#' @param tag Character vector of 2-letter tag names to extract.
#' @param as Output format:
#' - `"DataFrame"`: returns `S4Vectors::DataFrame` (default),
#' - `"data.frame"`: returns base `data.frame`,
#' - `"GAlignments"`: returns `GenomicAlignments::GAlignments`,
#' - `"GAlignmentPairs"`: returns `GenomicAlignments::GAlignmentPairs`,
#' - `"scanBam"`: returns a `scanBam()`-shaped list-of-lists.
#' @param seqqual_mode Controls representation of `seq`/`qual` when those
#'   fields are requested:
#' - `"compatible"` (default): return character vectors matching
#'   `scanBam`-style expectations,
#' - `"compact"`: return raw list-columns for faster/lower-overhead
#'   extraction. This mode is currently supported for
#'   `as = "data.frame"` or `as = "DataFrame"`.
#' @param threads Requested number of OpenMP threads used for
#'   reading/decompression. May be capped when `auto_threads = TRUE`.
#' @param BPPARAM Optional `BiocParallel` parameter used when `file` contains
#'   more than one BAM. If `NULL`, files are processed serially.
#' @param auto_threads Logical; when `TRUE` and `BPPARAM` has multiple workers,
#'   BamScale automatically caps per-file OpenMP threads to avoid
#'   oversubscription.
#' @param use.names Passed to alignment object conversion. When `TRUE`, read names
#'   (`qname`) are used as object names.
#' @param with.which_label Logical; if `TRUE` and `param` includes `which`,
#'   an extra `which_label` column is returned.
#' @param include_unmapped Logical; whether unmapped records are retained
#'   (subject to `param$flag` constraints).
#'
#' @return
#' If `file` is length 1: one object in the format specified by `as`.
#' If `file` has length > 1 (or is a `BamFileList`): a named list of outputs,
#' one per BAM file.
#'
#' @details
#' `bam_read()` is intentionally column-compatible with common BAM fields used by
#' Bioconductor workflows and can be used as a fast drop-in reader before
#' conversion to downstream classes.
#'
#' Parallelism model:
#' - `BPPARAM` parallelizes across files (one file per BiocParallel worker).
#' - `threads` parallelizes within each file via OpenMP.
#' - Effective total concurrency is approximately
#'   `min(length(file), BiocParallel::bpnworkers(BPPARAM)) * threads`.
#' - If `auto_threads = TRUE` and `BPPARAM` has multiple workers, per-file
#'   OpenMP threads are set to
#'   `max(1, min(threads, floor(available_cores / workers_eff)))`, where
#'   `workers_eff = min(length(file), BiocParallel::bpnworkers(BPPARAM))`.
#'
#' Compatibility notes:
#' - Region filtering via `param$which` is supported as a sequential filter
#'   (not index-jump random access).
#' - Flag filtering uses `ScanBamFlag` semantics by converting logical flag
#'   requirements into required-set and required-unset bit masks.
#' - Tag values are returned as character columns. Scalar tags are scalar
#'   strings; `B` tags are comma-separated vectors.
#' - `seqqual_mode = "compact"` is optimized for throughput-oriented
#'   benchmarking and returns raw list-columns for `seq`/`qual`.
#' - `"GAlignments"` and `"GAlignmentPairs"` output exclude unmapped records.
#' - `as = "scanBam"` returns a strict scan-like list-of-lists:
#'   without `param$which`, it returns one unnamed batch; with `param$which`,
#'   it returns one batch per range label (including empty ranges), with
#'   requested `what` fields and `tag` values under `$tag`.
#'   If Biostrings is installed, `seq` and `qual` are returned as
#'   `DNAStringSet` and `PhredQuality` for closer `scanBam()` compatibility.
#'
#' @examples
#' if (requireNamespace("ompBAM", quietly = TRUE)) {
#'   bam <- ompBAM::example_BAM("Unsorted")
#'
#'   # Familiar scanBam-like field selection
#'   x <- bam_read(bam, what = c("qname", "flag", "rname", "pos", "cigar"))
#'
#'   # Include sequence + quality
#'   y <- bam_read(bam, what = c("qname", "seq", "qual"), threads = 2)
#'
#'   # scanBam-shaped output
#'   z <- bam_read(bam, what = c("qname", "flag"), tag = "NM", as = "scanBam")
#' }
#' @export
bam_read <- function(
    file,
    param = NULL,
    what = NULL,
    tag = NULL,
    as = c("DataFrame", "data.frame", "GAlignments", "GAlignmentPairs", "scanBam"),
    seqqual_mode = c("compatible", "compact"),
    threads = 1L,
    BPPARAM = NULL,
    auto_threads = FALSE,
    use.names = FALSE,
    with.which_label = FALSE,
    include_unmapped = TRUE
) {
    as <- match.arg(as)
    seqqual_mode <- match.arg(seqqual_mode)

    files <- .bamscale_normalize_files(file)
    threads <- .bamscale_resolve_threads(
        threads = threads,
        BPPARAM = BPPARAM,
        auto_threads = auto_threads,
        n_files = length(files)
    )
    parsed <- .bamscale_parse_param(param)

    what_final <- what
    if (is.null(what_final)) {
        what_final <- parsed$what
    }
    if (length(what_final) == 0L) {
        what_final <- .bamscale_default_what
    }
    what_final <- unique(as.character(what_final))

    unsupported <- setdiff(what_final, .bamscale_supported_what)
    if (length(unsupported) > 0L) {
        stop(
            "Unsupported `what` fields: ",
            paste(unsupported, collapse = ", "),
            call. = FALSE
        )
    }

    tag_final <- tag
    if (is.null(tag_final)) {
        tag_final <- parsed$tag
    }
    tag_final <- unique(as.character(tag_final))
    tag_final <- tag_final[nzchar(tag_final)]

    split_scan_by_which <- identical(as, "scanBam") && nrow(parsed$which) > 0L
    internal_with_which_label <- isTRUE(with.which_label) || split_scan_by_which

    needed_fields <- unique(what_final)
    if (as %in% c("GAlignments", "GAlignmentPairs")) {
        needed_fields <- unique(c(needed_fields, "rname", "pos", "cigar", "strand"))
    }
    if (as == "GAlignmentPairs") {
        needed_fields <- unique(c(needed_fields, "qname", "flag"))
    }

    include_seq <- "seq" %in% needed_fields
    include_qual <- "qual" %in% needed_fields
    compact_seqqual <- identical(seqqual_mode, "compact") && (include_seq || include_qual)

    if (isTRUE(compact_seqqual) && !as %in% c("data.frame", "DataFrame")) {
        warning(
            "`seqqual_mode='compact'` is currently supported only for `as='data.frame'` ",
            "or `as='DataFrame'`; falling back to `seqqual_mode='compatible'`."
        )
        compact_seqqual <- FALSE
    }

    field_mask <- .bamscale_make_field_mask(needed_fields)

    worker <- function(path_one) {
        raw_df <- .Call(
            `_BamScale_read_bam_cpp`,
            path_one,
            threads,
            as.integer(parsed$min_mapq),
            as.logical(include_unmapped),
            as.logical(include_seq),
            as.logical(include_qual),
            as.logical(compact_seqqual),
            as.integer(parsed$flag_set),
            as.integer(parsed$flag_unset),
            as.character(tag_final),
            as.character(parsed$which$seqname),
            as.integer(parsed$which$start),
            as.integer(parsed$which$end),
            as.character(parsed$which$label),
            as.logical(internal_with_which_label),
            as.integer(field_mask)
        )

        raw_df <- .bamscale_decode_tag_sentinel(raw_df, tag_final)

        keep <- unique(c(what_final, tag_final))
        if (as %in% c("GAlignments", "GAlignmentPairs")) {
            keep <- unique(c(keep, "rname", "pos", "cigar", "strand"))
        }
        if (as == "GAlignmentPairs") {
            keep <- unique(c(keep, "qname", "flag"))
        }
        if (isTRUE(internal_with_which_label)) {
            keep <- unique(c(keep, "which_label"))
        }
        keep <- intersect(keep, colnames(raw_df))

        out <- raw_df[, keep, drop = FALSE]
        .bamscale_format_output(
            out,
            raw_df,
            as = as,
            use.names = use.names,
            what = what_final,
            tag = tag_final,
            with.which_label = with.which_label,
            split_scan_by_which = split_scan_by_which,
            which_df = parsed$which
        )
    }

    .bamscale_apply_files(files, worker, BPPARAM = BPPARAM)
}

#' Count BAM records with Bioconductor-compatible filtering
#'
#' `bam_count()` provides a fast chromosome-level count summary, honoring key
#' filtering fields from `ScanBamParam` (`mapqFilter`, `flag`, and `which`).
#'
#' @param file BAM input (`character`, `BamFile`, or `BamFileList`).
#' @param param Optional `Rsamtools::ScanBamParam` (or compatible list).
#' @param threads Requested number of OpenMP threads. May be capped when
#'   `auto_threads = TRUE`.
#' @param BPPARAM Optional `BiocParallel` parameter for multi-file operation.
#' @param auto_threads Logical; when `TRUE` and `BPPARAM` has multiple workers,
#'   BamScale automatically caps per-file OpenMP threads to avoid
#'   oversubscription.
#' @param include_unmapped Whether to include an extra `*` row for unmapped
#'   records.
#'
#' @return
#' For one file: a `data.frame` with columns `seqname`, `seqlength`, `count`.
#' For multiple files: named list of such `data.frame`s.
#'
#' @details
#' Parallelism behavior matches `bam_read()`: `BPPARAM` distributes work across
#' BAM files, while `threads` controls OpenMP work within each file. If
#' `auto_threads = TRUE` and `BPPARAM` has multiple workers, per-file OpenMP
#' threads are capped using
#' `max(1, min(threads, floor(available_cores / workers_eff)))`, where
#' `workers_eff = min(length(file), BiocParallel::bpnworkers(BPPARAM))`.
#'
#' @examples
#' if (requireNamespace("ompBAM", quietly = TRUE)) {
#'   bam <- ompBAM::example_BAM("Unsorted")
#'   bam_count(bam, threads = 2)
#' }
#' @export
bam_count <- function(
    file,
    param = NULL,
    threads = 1L,
    BPPARAM = NULL,
    auto_threads = FALSE,
    include_unmapped = TRUE
) {
    files <- .bamscale_normalize_files(file)
    threads <- .bamscale_resolve_threads(
        threads = threads,
        BPPARAM = BPPARAM,
        auto_threads = auto_threads,
        n_files = length(files)
    )
    parsed <- .bamscale_parse_param(param)

    worker <- function(path_one) {
        .Call(
            `_BamScale_count_bam_cpp`,
            path_one,
            threads,
            as.integer(parsed$min_mapq),
            as.logical(include_unmapped),
            as.integer(parsed$flag_set),
            as.integer(parsed$flag_unset),
            as.character(parsed$which$seqname),
            as.integer(parsed$which$start),
            as.integer(parsed$which$end)
        )
    }

    .bamscale_apply_files(files, worker, BPPARAM = BPPARAM)
}


.bamscale_parse_cpu_list <- function(x) {
    if (is.null(x) || !nzchar(x)) return(NA_integer_)

    parts <- strsplit(gsub("\\s+", "", x), ",", fixed = TRUE)[[1L]]
    total <- 0L

    for (p in parts) {
        if (!nzchar(p)) next

        if (grepl("-", p, fixed = TRUE)) {
            bounds <- strsplit(p, "-", fixed = TRUE)[[1L]]
            if (length(bounds) == 2L) {
                lo <- suppressWarnings(as.integer(bounds[[1L]]))
                hi <- suppressWarnings(as.integer(bounds[[2L]]))
                if (!is.na(lo) && !is.na(hi) && hi >= lo) {
                    total <- total + (hi - lo + 1L)
                }
            }
        } else {
            v <- suppressWarnings(as.integer(p))
            if (!is.na(v)) total <- total + 1L
        }
    }

    if (total > 0L) total else NA_integer_
}

.bamscale_detect_cores <- function() {
    # Honor Linux cpuset limits when available (containers, schedulers).
    if (.Platform$OS.type == "unix" && file.exists("/proc/self/status")) {
        status <- tryCatch(readLines("/proc/self/status", warn = FALSE), error = function(e) character())
        hit <- grep("^Cpus_allowed_list\\s*:", status, value = TRUE)
        if (length(hit) > 0L) {
            allowed <- sub("^Cpus_allowed_list\\s*:\\s*", "", hit[[1L]])
            n_allowed <- .bamscale_parse_cpu_list(allowed)
            if (!is.na(n_allowed) && n_allowed >= 1L) return(as.integer(n_allowed))
        }
    }

    physical <- suppressWarnings(as.integer(parallel::detectCores(logical = FALSE)))
    if (length(physical) == 1L && !is.na(physical) && physical >= 1L) {
        return(physical)
    }

    logical_cores <- suppressWarnings(as.integer(parallel::detectCores(logical = TRUE)))
    if (length(logical_cores) == 1L && !is.na(logical_cores) && logical_cores >= 1L) {
        return(logical_cores)
    }

    slurm <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")))
    if (length(slurm) == 1L && !is.na(slurm) && slurm >= 1L) {
        return(slurm)
    }

    1L
}

.bamscale_bpparam_workers <- function(BPPARAM) {
    workers <- tryCatch(
        as.integer(BiocParallel::bpnworkers(BPPARAM)),
        error = function(e) NA_integer_
    )

    if (length(workers) != 1L || is.na(workers) || workers < 1L) {
        return(1L)
    }

    workers
}

.bamscale_resolve_threads <- function(threads, BPPARAM = NULL, auto_threads = FALSE, n_files = 1L) {
    threads <- as.integer(threads)
    if (length(threads) != 1L || is.na(threads) || threads < 1L) {
        stop("`threads` must be a positive integer", call. = FALSE)
    }

    if (length(auto_threads) != 1L || is.na(auto_threads)) {
        stop("`auto_threads` must be TRUE or FALSE", call. = FALSE)
    }

    n_files <- as.integer(n_files)
    if (length(n_files) != 1L || is.na(n_files) || n_files < 1L) {
        n_files <- 1L
    }

    if (!isTRUE(auto_threads) || is.null(BPPARAM)) {
        return(threads)
    }

    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
        stop("`auto_threads = TRUE` with `BPPARAM` requires BiocParallel", call. = FALSE)
    }

    workers <- min(.bamscale_bpparam_workers(BPPARAM), n_files)
    if (workers <= 1L) {
        return(threads)
    }

    available_cores <- .bamscale_detect_cores()
    budget <- max(1L, floor(available_cores / workers))

    as.integer(max(1L, min(threads, budget)))
}


.bamscale_make_field_mask <- function(fields) {
    fields <- unique(as.character(fields))
    bad <- setdiff(fields, names(.bamscale_field_bits))
    if (length(bad) > 0L) {
        stop("Unsupported fields in internal mask: ", paste(bad, collapse = ", "), call. = FALSE)
    }

    mask <- 0L
    for (f in fields) {
        mask <- bitwOr(mask, unname(.bamscale_field_bits[[f]]))
    }
    as.integer(mask)
}

.bamscale_apply_files <- function(files, FUN, BPPARAM = NULL) {
    if (!is.null(BPPARAM)) {
        if (!requireNamespace("BiocParallel", quietly = TRUE)) {
            stop("`BPPARAM` was provided but BiocParallel is not installed", call. = FALSE)
        }
        out <- BiocParallel::bplapply(files, FUN, BPPARAM = BPPARAM)
    } else {
        out <- lapply(files, FUN)
    }

    if (length(files) == 1L) {
        return(out[[1L]])
    }

    names(out) <- names(files)
    out
}


.bamscale_bamfile_path <- function(x) {
    p <- tryCatch(
        {
            if (requireNamespace("BiocGenerics", quietly = TRUE)) {
                BiocGenerics::path(x)
            } else if ("path" %in% names(x)) {
                x$path
            } else {
                stop("No path accessor available")
            }
        },
        error = function(e) {
            if ("path" %in% names(x)) {
                x$path
            } else {
                stop(
                    "Failed to extract BAM path from BamFile: ",
                    conditionMessage(e),
                    call. = FALSE
                )
            }
        }
    )

    p <- as.character(p)
    if (length(p) != 1L || !nzchar(p)) {
        stop("Invalid BAM path extracted from BamFile", call. = FALSE)
    }
    p
}

.bamscale_normalize_files <- function(file) {
    if (is.character(file)) {
        if (length(file) < 1L || any(!nzchar(file))) {
            stop("`file` must contain non-empty BAM paths", call. = FALSE)
        }
        if (any(!file.exists(file))) {
            missing <- file[!file.exists(file)]
            stop("BAM file(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
        }
        nm <- names(file)
        if (is.null(nm)) nm <- basename(file)
        names(file) <- nm
        return(file)
    }

    if (!requireNamespace("Rsamtools", quietly = TRUE)) {
        stop(
            "Non-character BAM inputs require Rsamtools (BamFile/BamFileList)",
            call. = FALSE
        )
    }

    if (methods::is(file, "BamFile")) {
        p <- .bamscale_bamfile_path(file)
        names(p) <- if (!is.null(names(file)) && nzchar(names(file)[1L])) names(file)[1L] else basename(p)
        if (!file.exists(p)) stop("BAM file not found: ", p, call. = FALSE)
        return(p)
    }

    if (methods::is(file, "BamFileList")) {
        p <- vapply(file, .bamscale_bamfile_path, FUN.VALUE = character(1))
        nm <- names(file)
        if (is.null(nm) || any(!nzchar(nm))) nm <- basename(p)
        names(p) <- nm
        if (any(!file.exists(p))) {
            missing <- p[!file.exists(p)]
            stop("BAM file(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
        }
        return(p)
    }

    stop("`file` must be character, BamFile, or BamFileList", call. = FALSE)
}
.bamscale_empty_which <- function() {
    data.frame(
        seqname = character(),
        start = integer(),
        end = integer(),
        label = character(),
        stringsAsFactors = FALSE
    )
}

.bamscale_parse_param <- function(param) {
    out <- list(
        min_mapq = 0L,
        flag_set = 0L,
        flag_unset = 0L,
        which = .bamscale_empty_which(),
        what = character(),
        tag = character()
    )

    if (is.null(param)) {
        return(out)
    }

    if (is.list(param) && !methods::is(param, "ScanBamParam")) {
        if (!is.null(param$mapqFilter)) out$min_mapq <- as.integer(param$mapqFilter)
        if (!is.null(param$min_mapq)) out$min_mapq <- as.integer(param$min_mapq)

        if (!is.null(param$flag_set)) out$flag_set <- as.integer(param$flag_set)
        if (!is.null(param$flag_unset)) out$flag_unset <- as.integer(param$flag_unset)
        if (!is.null(param$flag) && methods::is(param$flag, "ScanBamFlag")) {
            masks <- .bamscale_scan_flag_masks(param$flag)
            out$flag_set <- masks$set
            out$flag_unset <- masks$unset
        }

        if (!is.null(param$what)) out$what <- as.character(param$what)
        if (!is.null(param$tag)) out$tag <- as.character(param$tag)
        if (!is.null(param$which)) out$which <- .bamscale_which_to_df(param$which)
        return(out)
    }

    if (!requireNamespace("Rsamtools", quietly = TRUE)) {
        stop("`param` as ScanBamParam requires Rsamtools", call. = FALSE)
    }

    if (!methods::is(param, "ScanBamParam")) {
        stop("`param` must be NULL, ScanBamParam, or a compatible list", call. = FALSE)
    }

    out$min_mapq <- as.integer(methods::slot(param, "mapqFilter"))
    out$what <- as.character(methods::slot(param, "what"))
    out$tag <- as.character(methods::slot(param, "tag"))

    flag <- methods::slot(param, "flag")
    masks <- .bamscale_scan_flag_masks(flag)
    out$flag_set <- masks$set
    out$flag_unset <- masks$unset

    out$which <- .bamscale_which_to_df(methods::slot(param, "which"))
    out
}

.bamscale_scan_flag_masks <- function(flag_obj) {
    bit_map <- c(
        isPaired = 0x1L,
        isProperPair = 0x2L,
        isUnmappedQuery = 0x4L,
        hasUnmappedMate = 0x8L,
        isMinusStrand = 0x10L,
        isMateMinusStrand = 0x20L,
        isFirstMateRead = 0x40L,
        isSecondMateRead = 0x80L,
        isNotPrimaryRead = 0x100L,
        isSecondaryAlignment = 0x100L,
        isNotPassingQualityControls = 0x200L,
        isDuplicateRead = 0x400L,
        isDuplicate = 0x400L,
        isSupplementaryAlignment = 0x800L
    )

    set_mask <- 0L
    unset_mask <- 0L

    for (nm in intersect(methods::slotNames(flag_obj), names(bit_map))) {
        val <- methods::slot(flag_obj, nm)
        if (length(val) != 1L || is.na(val)) next

        if (isTRUE(val)) {
            set_mask <- bitwOr(set_mask, bit_map[[nm]])
        } else if (identical(val, FALSE)) {
            unset_mask <- bitwOr(unset_mask, bit_map[[nm]])
        }
    }

    list(set = as.integer(set_mask), unset = as.integer(unset_mask))
}

.bamscale_which_to_df <- function(which_obj) {
    if (is.null(which_obj) || length(which_obj) == 0L) {
        return(.bamscale_empty_which())
    }

    if (is.data.frame(which_obj)) {
        req <- c("seqname", "start", "end")
        if (!all(req %in% names(which_obj))) {
            stop("`which` data.frame must contain seqname/start/end columns", call. = FALSE)
        }
        out <- data.frame(
            seqname = as.character(which_obj$seqname),
            start = as.integer(which_obj$start),
            end = as.integer(which_obj$end),
            label = if ("label" %in% names(which_obj)) as.character(which_obj$label) else NA_character_,
            stringsAsFactors = FALSE
        )
        missing_label <- is.na(out$label) | !nzchar(out$label)
        out$label[missing_label] <- paste0(out$seqname[missing_label], ":", out$start[missing_label], "-", out$end[missing_label])
        return(out)
    }

    if (methods::is(which_obj, "GenomicRanges")) {
        seqname <- as.character(GenomicRanges::seqnames(which_obj))
        start <- as.integer(IRanges::start(which_obj))
        end <- as.integer(IRanges::end(which_obj))
        label <- names(which_obj)
        if (is.null(label)) label <- rep(NA_character_, length(seqname))

        out <- data.frame(
            seqname = seqname,
            start = start,
            end = end,
            label = as.character(label),
            stringsAsFactors = FALSE
        )
        missing_label <- is.na(out$label) | !nzchar(out$label)
        out$label[missing_label] <- paste0(out$seqname[missing_label], ":", out$start[missing_label], "-", out$end[missing_label])
        return(out)
    }

    if (methods::is(which_obj, "IntegerRangesList")) {
        unlisted <- BiocGenerics::unlist(which_obj, use.names = FALSE)
        seqname <- rep(names(which_obj), times = as.integer(IRanges::elementNROWS(which_obj)))
        start <- as.integer(IRanges::start(unlisted))
        end <- as.integer(IRanges::end(unlisted))
        label <- paste0(seqname, ":", start, "-", end)

        return(data.frame(
            seqname = seqname,
            start = start,
            end = end,
            label = label,
            stringsAsFactors = FALSE
        ))
    }

    stop("Unsupported `which` type. Use GRanges, IntegerRangesList, or data.frame", call. = FALSE)
}

.bamscale_decode_tag_sentinel <- function(df, tag_names) {
    if (length(tag_names) == 0L) return(df)
    for (tg in intersect(tag_names, names(df))) {
        values <- as.character(df[[tg]])
        values[values == .bamscale_tag_na_sentinel] <- NA_character_
        df[[tg]] <- values
    }
    df
}

.bamscale_format_output <- function(
    selected_df,
    raw_df,
    as = "DataFrame",
    use.names = FALSE,
    what = character(),
    tag = character(),
    with.which_label = FALSE,
    split_scan_by_which = FALSE,
    which_df = .bamscale_empty_which()
) {
    if (identical(as, "data.frame")) {
        return(selected_df)
    }

    if (identical(as, "DataFrame")) {
        if (!requireNamespace("S4Vectors", quietly = TRUE)) {
            warning("S4Vectors not installed; returning data.frame")
            return(selected_df)
        }
        return(S4Vectors::DataFrame(selected_df, check.names = FALSE))
    }

    if (identical(as, "GAlignments")) {
        return(.bamscale_as_galignments(selected_df, raw_df, use.names = use.names))
    }

    if (identical(as, "GAlignmentPairs")) {
        return(.bamscale_as_galignment_pairs(selected_df, raw_df, use.names = use.names))
    }

    if (identical(as, "scanBam")) {
        return(.bamscale_as_scanbam(
            selected_df,
            what = what,
            tag = tag,
            split_by_which = split_scan_by_which,
            which_df = which_df
        ))
    }

    stop("Unsupported `as`: ", as, call. = FALSE)
}

.bamscale_as_scanbam <- function(selected_df, what, tag, split_by_which = FALSE, which_df = .bamscale_empty_which()) {
    df <- as.data.frame(selected_df, stringsAsFactors = FALSE)

    if (isTRUE(split_by_which)) {
        labels <- as.character(which_df$label)
        if (length(labels) == 0L) {
            labels <- character()
        }

        if ("which_label" %in% names(df) && nrow(df) > 0L) {
            by_label <- split(seq_len(nrow(df)), df$which_label, drop = TRUE)
        } else {
            by_label <- list()
        }

        batches <- vector("list", length(labels))
        for (i in seq_along(labels)) {
            idx <- by_label[[labels[[i]]]]
            if (is.null(idx)) idx <- integer()
            batches[[i]] <- idx
        }
        names(batches) <- labels
    } else if (nrow(df) == 0L) {
        batches <- list(integer())
    } else {
        batches <- list(seq_len(nrow(df)))
    }

    fields <- unique(as.character(what))
    tags <- unique(as.character(tag))

    out <- vector("list", length(batches))
    nms <- names(batches)
    for (i in seq_along(batches)) {
        idx <- batches[[i]]
        rec <- list()

        for (field in fields) {
            if (field %in% names(df)) {
                rec[[field]] <- df[[field]][idx]
            } else {
                rec[[field]] <- .bamscale_empty_field(field, length(idx))
            }
        }

        if (length(tags) > 0L) {
            tag_list <- stats::setNames(vector("list", length(tags)), tags)
            for (tg in tags) {
                if (tg %in% names(df)) {
                    tag_list[[tg]] <- df[[tg]][idx]
                } else {
                    tag_list[[tg]] <- rep(NA_character_, length(idx))
                }
            }
            rec$tag <- tag_list
        }

        rec <- .bamscale_scanbam_biostrings(rec)
        out[[i]] <- rec
    }

    if (is.null(nms)) {
        return(out)
    }
    names(out) <- nms
    out
}


.bamscale_scanbam_biostrings <- function(rec) {
    if (!requireNamespace("Biostrings", quietly = TRUE)) {
        return(rec)
    }

    if ("seq" %in% names(rec)) {
        seq_vals <- as.character(rec$seq)
        seq_vals[is.na(seq_vals)] <- "N"
        rec$seq <- Biostrings::DNAStringSet(seq_vals)
    }

    if ("qual" %in% names(rec)) {
        qual_vals <- as.character(rec$qual)
        qual_vals[is.na(qual_vals)] <- "*"
        rec$qual <- Biostrings::PhredQuality(qual_vals)
    }

    rec
}

.bamscale_empty_field <- function(field, n) {
    char_fields <- c("qname", "rname", "strand", "cigar", "mrnm", "seq", "qual")
    int_fields <- c("flag", "pos", "qwidth", "mapq", "mpos", "isize")

    if (field %in% char_fields) {
        return(rep(NA_character_, n))
    }
    if (field %in% int_fields) {
        return(rep(NA_integer_, n))
    }
    rep(NA, n)
}

.bamscale_as_galignments <- function(selected_df, raw_df, use.names = FALSE) {
    if (!requireNamespace("GenomicAlignments", quietly = TRUE)) {
        stop("`as = 'GAlignments'` requires GenomicAlignments", call. = FALSE)
    }

    required <- c("rname", "pos", "cigar", "strand")
    if (!all(required %in% names(selected_df))) {
        stop(
            "`as='GAlignments'` requires fields: ",
            paste(required, collapse = ", "),
            ". Add them to `what`.",
            call. = FALSE
        )
    }

    mapped <- !is.na(selected_df$pos) & selected_df$rname != "*"
    if (any(!mapped)) {
        warning("Dropping unmapped records for GAlignments output")
    }
    x <- selected_df[mapped, , drop = FALSE]

    if (nrow(x) == 0L) {
        return(GenomicAlignments::GAlignments())
    }

    extra_cols <- x[, setdiff(names(x), c("rname", "pos", "cigar", "strand")), drop = FALSE]

    header_seqnames <- attr(raw_df, "seqnames_header")
    if (!is.null(header_seqnames)) {
        seq_factor <- factor(x$rname, levels = as.character(header_seqnames))
    } else {
        seq_factor <- factor(x$rname)
    }

    ga_args <- c(
        list(
            seqnames = S4Vectors::Rle(seq_factor),
            pos = as.integer(x$pos),
            cigar = as.character(x$cigar),
            strand = as.character(x$strand),
            names = NULL
        ),
        as.list(extra_cols)
    )
    out <- do.call(GenomicAlignments::GAlignments, ga_args)

    if (isTRUE(use.names) && "qname" %in% names(x)) {
        names(out) <- as.character(x$qname)
    }

    if (!is.null(attr(raw_df, "seqnames_header")) && !is.null(attr(raw_df, "seqlengths_header")) &&
        requireNamespace("GenomeInfoDb", quietly = TRUE)) {
        seqinfo <- GenomeInfoDb::Seqinfo(
            seqnames = as.character(attr(raw_df, "seqnames_header")),
            seqlengths = as.integer(attr(raw_df, "seqlengths_header"))
        )
        suppressWarnings(GenomeInfoDb::seqinfo(out) <- seqinfo[GenomeInfoDb::seqlevels(out)])
    }

    out
}

.bamscale_as_galignment_pairs <- function(selected_df, raw_df, use.names = FALSE) {
    if (!requireNamespace("GenomicAlignments", quietly = TRUE)) {
        stop("`as = 'GAlignmentPairs'` requires GenomicAlignments", call. = FALSE)
    }

    required <- c("qname", "flag", "rname", "pos", "cigar", "strand")
    if (!all(required %in% names(selected_df))) {
        stop(
            "`as='GAlignmentPairs'` requires fields: ",
            paste(required, collapse = ", "),
            ". Add them to `what`.",
            call. = FALSE
        )
    }

    flag <- as.integer(selected_df$flag)
    mapped <- !is.na(selected_df$pos) & selected_df$rname != "*"
    paired <- bitwAnd(flag, 0x1L) != 0L
    keep <- mapped & paired

    if (any(!keep)) {
        warning("Dropping non-paired or unmapped records for GAlignmentPairs output")
    }

    x <- selected_df[keep, , drop = FALSE]
    if (nrow(x) == 0L) {
        return(GenomicAlignments::GAlignmentPairs())
    }

    ga <- .bamscale_as_galignments(x, raw_df, use.names = FALSE)

    flag <- as.integer(x$flag)
    qname <- as.character(x$qname)
    idx_first <- which(bitwAnd(flag, 0x40L) != 0L)
    idx_second <- which(bitwAnd(flag, 0x80L) != 0L)

    if (length(idx_first) == 0L || length(idx_second) == 0L) {
        warning("No complete first/second mate pairs detected")
        return(GenomicAlignments::GAlignmentPairs())
    }

    split_first <- split(idx_first, qname[idx_first])
    split_second <- split(idx_second, qname[idx_second])
    keys <- intersect(names(split_first), names(split_second))

    n_keys <- length(keys)
    p1_parts <- vector("list", n_keys)
    p2_parts <- vector("list", n_keys)
    name_parts <- vector("list", n_keys)
    dropped <- 0L

    if (n_keys > 0L) {
        for (i in seq_along(keys)) {
            k <- keys[[i]]
            f <- split_first[[k]]
            s <- split_second[[k]]
            n <- min(length(f), length(s))

            if (n > 0L) {
                p1_parts[[i]] <- f[seq_len(n)]
                p2_parts[[i]] <- s[seq_len(n)]
                name_parts[[i]] <- rep(k, n)
            }

            dropped <- dropped + abs(length(f) - length(s))
        }

        keep <- vapply(p1_parts, length, integer(1)) > 0L
        p1 <- unlist(p1_parts[keep], use.names = FALSE)
        p2 <- unlist(p2_parts[keep], use.names = FALSE)
        pair_names <- unlist(name_parts[keep], use.names = FALSE)
    } else {
        p1 <- integer()
        p2 <- integer()
        pair_names <- character()
    }

    if (length(p1) == 0L) {
        warning("No pairable first/second mate records remained after filtering")
        return(GenomicAlignments::GAlignmentPairs())
    }

    first <- ga[p1]
    second <- ga[p2]

    pairs <- tryCatch(
        GenomicAlignments::GAlignmentPairs(first, second),
        error = function(e) GenomicAlignments::makeGAlignmentPairs(first, second)
    )

    if (isTRUE(use.names)) {
        names(pairs) <- pair_names
    }

    if (dropped > 0L) {
        warning("Dropped ", dropped, " unpaired mate records while forming GAlignmentPairs")
    }

    pairs
}
