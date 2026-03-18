#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(BamScale)
})


.bamscale_parse_cpu_list <- function(x) {
  if (is.null(x) || !nzchar(x)) return(NA_integer_)
  parts <- strsplit(gsub("\\s+", "", x), ",", fixed = TRUE)[[1L]]
  total <- 0L
  for (p in parts) {
    if (!nzchar(p)) next
    if (grepl("-", p, fixed = TRUE)) {
      ends <- strsplit(p, "-", fixed = TRUE)[[1L]]
      if (length(ends) == 2L) {
        lo <- suppressWarnings(as.integer(ends[[1L]]))
        hi <- suppressWarnings(as.integer(ends[[2L]]))
        if (!is.na(lo) && !is.na(hi) && hi >= lo) total <- total + (hi - lo + 1L)
      }
    } else {
      v <- suppressWarnings(as.integer(p))
      if (!is.na(v)) total <- total + 1L
    }
  }
  if (total > 0L) total else NA_integer_
}

.bamscale_detect_cores <- function() {
  if (.Platform$OS.type == "unix" && file.exists("/proc/self/status")) {
    s <- tryCatch(readLines("/proc/self/status", warn = FALSE), error = function(e) character())
    hit <- grep("^Cpus_allowed_list\\s*:", s, value = TRUE)
    if (length(hit) > 0L) {
      allowed <- sub("^Cpus_allowed_list\\s*:\\s*", "", hit[[1L]])
      n_allowed <- .bamscale_parse_cpu_list(allowed)
      if (!is.na(n_allowed) && n_allowed >= 1L) return(as.integer(n_allowed))
    }
  }

  n <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) NA_integer_)
  n <- suppressWarnings(as.integer(n))
  if (!is.na(n) && n >= 1L) return(n)

  slurm <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "")))
  if (!is.na(slurm) && slurm >= 1L) return(slurm)

  1L
}



assert_bamscale_symbols <- function() {
  required <- c("_BamScale_read_bam_cpp", "_BamScale_count_bam_cpp")
  missing <- required[!vapply(
    required,
    function(sym) {
      isTRUE(tryCatch({
        getNativeSymbolInfo(sym, PACKAGE = "BamScale")
        TRUE
      }, error = function(e) FALSE))
    },
    logical(1)
  )]

  if (length(missing) > 0L) {
    stop(
      paste0(
        "BamScale native symbols are not loaded: ", paste(missing, collapse = ", "), "\n",
        "Reinstall BamScale from source and re-run."
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

assert_bamscale_symbols()

parse_args <- function(args) {
  out <- list()
  out$.positional <- character()

  if (length(args) == 0L) return(out)

  for (arg in args) {
    if (grepl("^--[^=]+=", arg)) {
      key <- sub("^--([^=]+)=.*$", "\\1", arg)
      val <- sub("^--[^=]+=", "", arg)
      out[[key]] <- val
    } else if (grepl("^--", arg)) {
      key <- sub("^--", "", arg)
      out[[key]] <- "TRUE"
    } else {
      out$.positional <- c(out$.positional, arg)
    }
  }

  out
}

split_csv <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character())
  vals <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
  vals <- trimws(vals)
  vals[nzchar(vals)]
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0L || !nzchar(x)) return(default)
  y <- tolower(as.character(x[[1L]]))
  if (y %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (y %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop("Invalid logical value: ", x)
}

as_int_scalar <- function(x, default) {
  if (is.null(x) || length(x) == 0L || !nzchar(x)) return(as.integer(default))
  out <- suppressWarnings(as.integer(x[[1L]]))
  if (length(out) != 1L || is.na(out)) stop("Invalid integer value: ", x)
  out
}

as_num_scalar <- function(x, default) {
  if (is.null(x) || length(x) == 0L || !nzchar(x)) return(as.numeric(default))
  out <- suppressWarnings(as.numeric(x[[1L]]))
  if (length(out) != 1L || is.na(out) || !is.finite(out)) stop("Invalid numeric value: ", x)
  out
}

as_int_vec <- function(x, default) {
  if (is.null(x) || length(x) == 0L || !nzchar(x)) return(as.integer(default))
  vals <- split_csv(x)
  out <- suppressWarnings(as.integer(vals))
  out <- out[!is.na(out) & out >= 1L]
  if (length(out) == 0L) stop("Invalid integer vector value: ", x)
  unique(sort(out))
}

as_chr_vec <- function(x, default = character()) {
  vals <- split_csv(x)
  if (length(vals) == 0L) return(default)
  unique(vals)
}

require_or_stop <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required package: ", pkg, call. = FALSE)
  }
}

bamfile_path <- function(x) {
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
        stop("Unable to extract BAM path: ", conditionMessage(e), call. = FALSE)
      }
    }
  )

  p <- as.character(p)
  p[[1L]]
}

fetch_chipseqdbdata_bams <- function(n_target = 8L) {
  require_or_stop("chipseqDBData")

  tbl <- chipseqDBData::H3K9acData()
  if (nrow(tbl) == 0L || !("Path" %in% colnames(tbl))) {
    stop("chipseqDBData::H3K9acData() returned no BAM paths", call. = FALSE)
  }

  reads_num <- if ("Reads" %in% colnames(tbl)) {
    suppressWarnings(as.numeric(gsub("[^0-9.]+", "", as.character(tbl$Reads))))
  } else {
    rep(NA_real_, nrow(tbl))
  }

  score <- ifelse(is.na(reads_num), -Inf, reads_num)
  ord <- order(score, decreasing = TRUE)

  paths <- vapply(tbl$Path, bamfile_path, FUN.VALUE = character(1))
  paths <- paths[ord]
  paths <- unique(paths[file.exists(paths)])

  if (length(paths) == 0L) {
    stop("No existing BAM files were resolved from chipseqDBData", call. = FALSE)
  }

  head(paths, n = as.integer(n_target))
}

bam_index_exists <- function(bam) {
  file.exists(paste0(bam, ".bai")) || file.exists(sub("\\.bam$", ".bai", bam, ignore.case = TRUE))
}

ensure_bam_index <- function(bam) {
  if (bam_index_exists(bam)) return(invisible(FALSE))
  require_or_stop("Rsamtools")
  Rsamtools::indexBam(bam)
  invisible(TRUE)
}

select_n_files <- function(paths, n_files, allow_repeat = FALSE) {
  if (length(paths) >= n_files) {
    return(paths[seq_len(n_files)])
  }

  if (!allow_repeat) {
    stop(
      "Requested ", n_files, " BAM files but only ", length(paths),
      " are available. Provide more files or use --allow-repeat-files=true.",
      call. = FALSE
    )
  }

  rep(paths, length.out = n_files)
}

make_bpparam <- function(workers, backend = c("snow", "multicore")) {
  workers <- as.integer(workers)
  if (workers <= 1L) return(NULL)

  require_or_stop("BiocParallel")
  backend <- match.arg(backend)

  if (backend == "multicore" && .Platform$OS.type == "unix") {
    BiocParallel::MulticoreParam(
      workers = workers,
      progressbar = FALSE,
      stop.on.error = TRUE
    )
  } else {
    BiocParallel::SnowParam(
      workers = workers,
      type = "SOCK",
      progressbar = FALSE,
      stop.on.error = TRUE
    )
  }
}

stop_bpparam <- function(bp) {
  if (is.null(bp)) return(invisible(NULL))
  try(BiocParallel::bpstop(bp), silent = TRUE)
  invisible(NULL)
}

count_scan_batch <- function(batch) {
  if (!is.list(batch) || length(batch) == 0L) return(0)

  field_names <- setdiff(names(batch), "tag")
  if (length(field_names) > 0L && !is.null(batch[[field_names[[1L]]]])) {
    return(length(batch[[field_names[[1L]]]]))
  }

  if ("tag" %in% names(batch) && length(batch$tag) > 0L) {
    return(length(batch$tag[[1L]]))
  }

  0
}

count_records <- function(x) {
  if (is.null(x)) return(0)

  if (inherits(x, "GAlignments") || inherits(x, "GAlignmentPairs")) {
    return(as.numeric(length(x)))
  }

  if (is.data.frame(x) || inherits(x, "DataFrame")) {
    return(as.numeric(nrow(x)))
  }

  if (is.list(x)) {
    if (length(x) == 0L) return(0)

    if (!is.null(names(x)) && any(names(x) %in% c("qname", "flag", "rname", "pos", "seq", "qual", "tag"))) {
      return(as.numeric(count_scan_batch(x)))
    }

    return(sum(vapply(x, count_records, numeric(1)), na.rm = TRUE))
  }

  NA_real_
}

bench_times <- function(fun, iterations = 3L, warmup = TRUE) {
  elapsed <- rep(NA_real_, iterations)
  errors <- rep(NA_character_, iterations)

  if (isTRUE(warmup)) {
    message("  warmup: start")
    warm_status <- tryCatch(
      {
        invisible(fun())
        "ok"
      },
      error = function(e) {
        errors[[1L]] <<- paste("warmup:", conditionMessage(e))
        "error"
      }
    )
    message("  warmup: ", warm_status)
    if (warm_status == "error") return(list(elapsed = elapsed, errors = errors))
  }

  for (i in seq_len(iterations)) {
    gc(verbose = FALSE)
    message(sprintf("  iteration %d/%d: start", i, iterations))

    t0 <- proc.time()[["elapsed"]]
    iter_ok <- tryCatch(
      {
        invisible(fun())
        TRUE
      },
      error = function(e) {
        errors[[i]] <<- conditionMessage(e)
        FALSE
      }
    )
    t1 <- proc.time()[["elapsed"]]

    if (isTRUE(iter_ok)) {
      elapsed[[i]] <- t1 - t0
      message(sprintf("  iteration %d/%d: %.3f s", i, iterations, elapsed[[i]]))
    } else {
      message(sprintf("  iteration %d/%d: ERROR (%s)", i, iterations, errors[[i]]))
    }
  }

  list(elapsed = elapsed, errors = errors)
}

summarize_elapsed <- function(valid) {
  min_s <- if (length(valid) > 0L) min(valid) else NA_real_
  median_s <- if (length(valid) > 0L) median(valid) else NA_real_
  max_s <- if (length(valid) > 0L) max(valid) else NA_real_
  mean_s <- if (length(valid) > 0L) mean(valid) else NA_real_
  sd_s <- if (length(valid) > 1L) stats::sd(valid) else NA_real_
  sem_s <- if (length(valid) > 1L && is.finite(sd_s)) sd_s / sqrt(length(valid)) else NA_real_
  mad_s <- if (length(valid) > 0L) stats::mad(valid) else NA_real_
  p25_s <- if (length(valid) > 0L) unname(stats::quantile(valid, probs = 0.25, names = FALSE, type = 7)) else NA_real_
  p75_s <- if (length(valid) > 0L) unname(stats::quantile(valid, probs = 0.75, names = FALSE, type = 7)) else NA_real_
  iqr_s <- if (is.finite(p25_s) && is.finite(p75_s)) p75_s - p25_s else NA_real_
  cv_s <- if (is.finite(mean_s) && mean_s > 0 && is.finite(sd_s)) sd_s / mean_s else NA_real_

  list(
    min_s = min_s,
    median_s = median_s,
    max_s = max_s,
    mean_s = mean_s,
    sd_s = sd_s,
    sem_s = sem_s,
    mad_s = mad_s,
    p25_s = p25_s,
    p75_s = p75_s,
    iqr_s = iqr_s,
    cv_s = cv_s
  )
}

run_case <- function(
    scenario,
    workload,
    method,
    method_family,
    threads_requested,
    threads_effective,
    bp_workers_requested,
    bp_workers,
    n_files,
    n_records,
    total_mb,
    iterations,
    warmup,
    fun,
    case_id,
    execution_rank = case_id,
    outer_run = getOption("bamscale.benchmark.outer_run", 1L),
    comparison_track = "fair",
    seqqual_mode = NA_character_,
    auto_threads = FALSE,
    bp_backend = NA_character_
) {
  message(
    sprintf(
      "[case %03d] %s | %s | %s | req_t=%s eff_t=%s workers=%s files=%s",
      case_id,
      scenario,
      workload,
      method,
      threads_requested,
      threads_effective,
      bp_workers,
      n_files
    )
  )

  out <- bench_times(fun, iterations = iterations, warmup = warmup)
  elapsed <- out$elapsed
  errors <- out$errors

  iter_df <- data.frame(
    case_id = rep(case_id, iterations),
    case_execution_id = rep(sprintf("%03d_%02d", as.integer(case_id), as.integer(outer_run)), iterations),
    outer_run = rep(as.integer(outer_run), iterations),
    execution_rank = rep(as.integer(execution_rank), iterations),
    scenario = rep(scenario, iterations),
    workload = rep(workload, iterations),
    method = rep(method, iterations),
    method_family = rep(method_family, iterations),
    comparison_track = rep(comparison_track, iterations),
    seqqual_mode = rep(if (is.na(seqqual_mode)) NA_character_ else as.character(seqqual_mode), iterations),
    auto_threads = rep(as.logical(auto_threads), iterations),
    threads_requested = rep(as.integer(threads_requested), iterations),
    threads_effective = rep(as.integer(threads_effective), iterations),
    bp_workers_requested = rep(as.integer(bp_workers_requested), iterations),
    bp_workers_effective = rep(as.integer(bp_workers), iterations),
    total_threads_requested = rep(as.integer(max(1L, bp_workers_requested) * max(1L, threads_requested)), iterations),
    total_threads_effective = rep(as.integer(max(1L, bp_workers) * max(1L, threads_effective)), iterations),
    bp_backend = rep(if (is.na(bp_backend)) NA_character_ else as.character(bp_backend), iterations),
    iteration = seq_len(iterations),
    elapsed_s = elapsed,
    status = ifelse(is.finite(elapsed), "ok", "error"),
    error_message = ifelse(is.finite(elapsed), NA_character_, errors),
    stringsAsFactors = FALSE
  )

  valid <- elapsed[is.finite(elapsed)]
  status <- if (length(valid) > 0L) "ok" else "error"
  error_message <- if (status == "error") paste(unique(na.omit(errors)), collapse = " | ") else NA_character_
  elapsed_stats <- summarize_elapsed(valid)

  records_per_s <- if (!is.na(n_records) && !is.na(elapsed_stats$median_s) && elapsed_stats$median_s > 0) n_records / elapsed_stats$median_s else NA_real_
  mb_per_s <- if (!is.na(total_mb) && !is.na(elapsed_stats$median_s) && elapsed_stats$median_s > 0) total_mb / elapsed_stats$median_s else NA_real_
  records_per_s_mean <- if (!is.na(n_records) && !is.na(elapsed_stats$mean_s) && elapsed_stats$mean_s > 0) n_records / elapsed_stats$mean_s else NA_real_
  mb_per_s_mean <- if (!is.na(total_mb) && !is.na(elapsed_stats$mean_s) && elapsed_stats$mean_s > 0) total_mb / elapsed_stats$mean_s else NA_real_

  comparison_track_scalar <- if (
    is.null(comparison_track) || length(comparison_track) == 0L ||
      is.na(comparison_track[[1L]]) || !nzchar(as.character(comparison_track[[1L]]))
  ) {
    "fair"
  } else {
    as.character(comparison_track[[1L]])
  }

  seqqual_mode_scalar <- if (
    is.null(seqqual_mode) || length(seqqual_mode) == 0L ||
      is.na(seqqual_mode[[1L]]) || !nzchar(as.character(seqqual_mode[[1L]]))
  ) {
    NA_character_
  } else {
    as.character(seqqual_mode[[1L]])
  }

  summary_df <- data.frame(
    case_id = case_id,
    case_execution_id = sprintf("%03d_%02d", as.integer(case_id), as.integer(outer_run)),
    outer_run = as.integer(outer_run),
    execution_rank = as.integer(execution_rank),
    scenario = scenario,
    workload = workload,
    method = method,
    method_family = method_family,
    comparison_track = comparison_track_scalar,
    seqqual_mode = seqqual_mode_scalar,
    auto_threads = as.logical(auto_threads),
    threads_requested = as.integer(threads_requested),
    threads_effective = as.integer(threads_effective),
    bp_workers_requested = as.integer(bp_workers_requested),
    bp_workers = as.integer(bp_workers),
    bp_workers_effective = as.integer(bp_workers),
    total_threads_requested = as.integer(max(1L, bp_workers_requested) * max(1L, threads_requested)),
    total_threads = as.integer(max(1L, bp_workers) * max(1L, threads_effective)),
    total_threads_effective = as.integer(max(1L, bp_workers) * max(1L, threads_effective)),
    bp_backend = if (is.na(bp_backend)) NA_character_ else as.character(bp_backend),
    n_files = as.integer(n_files),
    n_records = as.numeric(n_records),
    total_mb = as.numeric(total_mb),
    iterations_planned = as.integer(iterations),
    iterations_completed = as.integer(length(valid)),
    status = status,
    error_message = error_message,
    min_s = elapsed_stats$min_s,
    median_s = elapsed_stats$median_s,
    mean_s = elapsed_stats$mean_s,
    max_s = elapsed_stats$max_s,
    sd_s = elapsed_stats$sd_s,
    sem_s = elapsed_stats$sem_s,
    mad_s = elapsed_stats$mad_s,
    p25_s = elapsed_stats$p25_s,
    p75_s = elapsed_stats$p75_s,
    iqr_s = elapsed_stats$iqr_s,
    cv_s = elapsed_stats$cv_s,
    records_per_s = records_per_s,
    records_per_s_mean = records_per_s_mean,
    mb_per_s = mb_per_s,
    mb_per_s_mean = mb_per_s_mean,
    stringsAsFactors = FALSE
  )

  list(summary = summary_df, iterations = iter_df)
}

run_rsamtools_multi <- function(files, what, bp = NULL) {
  require_or_stop("Rsamtools")
  param <- Rsamtools::ScanBamParam(what = as.character(what))
  fun_one <- function(path_one) {
    Rsamtools::scanBam(path_one, param = param)
  }

  if (is.null(bp)) {
    lapply(files, fun_one)
  } else {
    BiocParallel::bplapply(files, fun_one, BPPARAM = bp)
  }
}

run_galignments_multi <- function(files, what, bp = NULL) {
  require_or_stop("Rsamtools")
  require_or_stop("GenomicAlignments")
  param <- Rsamtools::ScanBamParam(what = as.character(what))
  fun_one <- function(path_one) {
    GenomicAlignments::readGAlignments(path_one, param = param)
  }

  if (is.null(bp)) {
    lapply(files, fun_one)
  } else {
    BiocParallel::bplapply(files, fun_one, BPPARAM = bp)
  }
}

plain_vector <- function(x) {
  if (is.null(x)) return(x)
  if (is.factor(x)) return(as.character(x))
  if (inherits(x, "Rle")) return(as.vector(x))
  if (inherits(x, "XStringSet") || inherits(x, "BStringSet") || inherits(x, "PhredQuality")) {
    return(as.character(x))
  }
  if (isS4(x)) {
    as_vec <- tryCatch(as.vector(x), error = function(e) NULL)
    if (!is.null(as_vec)) return(as_vec)
    as_chr <- tryCatch(as.character(x), error = function(e) NULL)
    if (!is.null(as_chr)) return(as_chr)
  }
  x
}

plain_df <- function(x) {
  if (inherits(x, "DataFrame")) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }
  if (!is.data.frame(x)) stop("Expected a data.frame-like object", call. = FALSE)
  x[] <- lapply(x, plain_vector)
  rownames(x) <- NULL
  x
}

scanbam_batch_to_df <- function(batch, fields) {
  out <- setNames(vector("list", length(fields)), fields)
  for (nm in fields) {
    out[[nm]] <- plain_vector(batch[[nm]])
  }
  plain_df(as.data.frame(out, stringsAsFactors = FALSE))
}

galignments_to_df <- function(x, fields = c("qname", "rname", "pos", "cigar", "strand", "flag")) {
  out <- list(
    qname = if ("qname" %in% names(S4Vectors::mcols(x))) plain_vector(S4Vectors::mcols(x)$qname) else rep(NA_character_, length(x)),
    rname = as.character(GenomicRanges::seqnames(x)),
    pos = as.integer(GenomicRanges::start(x)),
    cigar = as.character(GenomicAlignments::cigar(x)),
    strand = as.character(BiocGenerics::strand(x)),
    flag = if ("flag" %in% names(S4Vectors::mcols(x))) as.integer(S4Vectors::mcols(x)$flag) else rep(NA_integer_, length(x))
  )
  plain_df(as.data.frame(out[fields], stringsAsFactors = FALSE))
}

sort_plain_df <- function(df) {
  df <- plain_df(df)
  if (nrow(df) <= 1L || ncol(df) == 0L) return(df)
  ord_cols <- lapply(df, function(col) ifelse(is.na(col), "__NA__", as.character(plain_vector(col))))
  ord <- do.call(order, c(ord_cols, list(na.last = TRUE)))
  df[ord, , drop = FALSE]
}

compare_plain_df <- function(lhs, rhs) {
  lhs <- sort_plain_df(lhs)
  rhs <- sort_plain_df(rhs)

  cols <- intersect(colnames(lhs), colnames(rhs))
  lhs <- lhs[, cols, drop = FALSE]
  rhs <- rhs[, cols, drop = FALSE]

  same_nrow <- nrow(lhs) == nrow(rhs)
  identical_data <- same_nrow && identical(lhs, rhs)

  list(
    same_nrow = same_nrow,
    identical_data = identical_data,
    columns_compared = paste(cols, collapse = ","),
    bam_n = nrow(lhs),
    ref_n = nrow(rhs)
  )
}

make_correctness_region <- function(bam, window_bp = 1000L) {
  seed <- BamScale::bam_read(
    file = bam,
    what = c("rname", "pos", "flag"),
    as = "data.frame",
    threads = 1L,
    include_unmapped = FALSE,
    auto_threads = FALSE
  )
  seed <- plain_df(seed)
  keep <- which(!is.na(seed$rname) & nzchar(as.character(seed$rname)) & !is.na(seed$pos) & seed$pos > 0L)
  if (length(keep) == 0L) {
    stop("Unable to derive a correctness region from the representative BAM", call. = FALSE)
  }
  idx <- keep[[1L]]
  start <- max(1L, as.integer(seed$pos[[idx]]))
  end <- start + as.integer(window_bp) - 1L
  data.frame(
    seqname = as.character(seed$rname[[idx]]),
    start = start,
    end = end,
    stringsAsFactors = FALSE
  )
}

run_correctness_preflight <- function(single_file, workloads, active_workloads, cfg, out_dir) {
  if (!isTRUE(cfg$correctness_preflight)) {
    df <- data.frame(
      workload = character(),
      comparator = character(),
      status = character(),
      passed = logical(),
      note = character(),
      stringsAsFactors = FALSE
    )
    write.csv(df, file = file.path(out_dir, "correctness_preflight.csv"), row.names = FALSE)
    return(df)
  }

  comparator_for <- function(w) {
    if (w == "galignments") "GenomicAlignments::readGAlignments" else "Rsamtools::scanBam"
  }

  if (!bam_index_exists(single_file)) {
    df <- data.frame(
      workload = active_workloads,
      comparator = vapply(active_workloads, comparator_for, character(1)),
      status = rep("skipped", length(active_workloads)),
      passed = rep(NA, length(active_workloads)),
      seqname = rep(NA_character_, length(active_workloads)),
      start = rep(NA_integer_, length(active_workloads)),
      end = rep(NA_integer_, length(active_workloads)),
      columns_compared = rep(NA_character_, length(active_workloads)),
      bam_n = rep(NA_integer_, length(active_workloads)),
      ref_n = rep(NA_integer_, length(active_workloads)),
      same_nrow = rep(NA, length(active_workloads)),
      identical_data = rep(NA, length(active_workloads)),
      note = rep("Skipping correctness preflight because the representative BAM has no index.", length(active_workloads)),
      stringsAsFactors = FALSE
    )
    write.csv(df, file = file.path(out_dir, "correctness_preflight.csv"), row.names = FALSE)
    return(df)
  }

  region <- tryCatch(make_correctness_region(single_file, window_bp = cfg$correctness_window_bp), error = function(e) e)
  if (inherits(region, "error")) {
    df <- data.frame(
      workload = active_workloads,
      comparator = vapply(active_workloads, comparator_for, character(1)),
      status = rep("error", length(active_workloads)),
      passed = rep(FALSE, length(active_workloads)),
      seqname = rep(NA_character_, length(active_workloads)),
      start = rep(NA_integer_, length(active_workloads)),
      end = rep(NA_integer_, length(active_workloads)),
      columns_compared = rep(NA_character_, length(active_workloads)),
      bam_n = rep(NA_integer_, length(active_workloads)),
      ref_n = rep(NA_integer_, length(active_workloads)),
      same_nrow = rep(NA, length(active_workloads)),
      identical_data = rep(NA, length(active_workloads)),
      note = rep(conditionMessage(region), length(active_workloads)),
      stringsAsFactors = FALSE
    )
    write.csv(df, file = file.path(out_dir, "correctness_preflight.csv"), row.names = FALSE)
    return(df)
  }

  region_gr <- GenomicRanges::GRanges(
    seqnames = region$seqname,
    ranges = IRanges::IRanges(start = region$start, end = region$end)
  )
  param_df <- data.frame(
    seqname = region$seqname,
    start = region$start,
    end = region$end,
    label = "correctness",
    stringsAsFactors = FALSE
  )

  out <- lapply(active_workloads, function(w) {
    comparator <- comparator_for(w)

    if (w %in% c("step1", "seqqual") && !requireNamespace("Rsamtools", quietly = TRUE)) {
      return(data.frame(
        workload = w,
        comparator = comparator,
        status = "skipped",
        passed = NA,
        seqname = region$seqname,
        start = region$start,
        end = region$end,
        columns_compared = NA_character_,
        bam_n = NA_integer_,
        ref_n = NA_integer_,
        same_nrow = NA,
        identical_data = NA,
        note = "Rsamtools is not installed.",
        stringsAsFactors = FALSE
      ))
    }
    if (w == "galignments" && (!requireNamespace("Rsamtools", quietly = TRUE) || !requireNamespace("GenomicAlignments", quietly = TRUE))) {
      return(data.frame(
        workload = w,
        comparator = comparator,
        status = "skipped",
        passed = NA,
        seqname = region$seqname,
        start = region$start,
        end = region$end,
        columns_compared = NA_character_,
        bam_n = NA_integer_,
        ref_n = NA_integer_,
        same_nrow = NA,
        identical_data = NA,
        note = "GenomicAlignments and/or Rsamtools is not installed.",
        stringsAsFactors = FALSE
      ))
    }

    check <- tryCatch(
      {
        if (w %in% c("step1", "seqqual")) {
          fields <- workloads[[w]]$what
          bam_df <- plain_df(BamScale::bam_read(
            file = single_file,
            what = fields,
            as = "data.frame",
            param = list(which = param_df),
            threads = 1L,
            auto_threads = FALSE
          ))
          ref_batch <- Rsamtools::scanBam(
            single_file,
            param = Rsamtools::ScanBamParam(what = as.character(fields), which = region_gr)
          )[[1L]]
          ref_df <- scanbam_batch_to_df(ref_batch, fields)
          compare_plain_df(bam_df, ref_df)
        } else {
          fields <- c("qname", "rname", "pos", "cigar", "strand", "flag")
          bam_ga <- BamScale::bam_read(
            file = single_file,
            what = fields,
            as = "GAlignments",
            param = list(which = param_df),
            threads = 1L,
            auto_threads = FALSE
          )
          ref_ga <- GenomicAlignments::readGAlignments(
            single_file,
            param = Rsamtools::ScanBamParam(what = as.character(fields), which = region_gr)
          )
          compare_plain_df(galignments_to_df(bam_ga, fields = fields), galignments_to_df(ref_ga, fields = fields))
        }
      },
      error = function(e) e
    )

    if (inherits(check, "error")) {
      data.frame(
        workload = w,
        comparator = comparator,
        status = "error",
        passed = FALSE,
        seqname = region$seqname,
        start = region$start,
        end = region$end,
        columns_compared = NA_character_,
        bam_n = NA_integer_,
        ref_n = NA_integer_,
        same_nrow = NA,
        identical_data = NA,
        note = conditionMessage(check),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        workload = w,
        comparator = comparator,
        status = if (isTRUE(check$identical_data)) "ok" else "mismatch",
        passed = isTRUE(check$identical_data),
        seqname = region$seqname,
        start = region$start,
        end = region$end,
        columns_compared = check$columns_compared,
        bam_n = as.integer(check$bam_n),
        ref_n = as.integer(check$ref_n),
        same_nrow = isTRUE(check$same_nrow),
        identical_data = isTRUE(check$identical_data),
        note = if (isTRUE(check$identical_data)) "" else "Normalized subset comparison did not match comparator output.",
        stringsAsFactors = FALSE
      )
    }
  })

  out <- do.call(rbind, out)
  write.csv(out, file = file.path(out_dir, "correctness_preflight.csv"), row.names = FALSE)

  if (isTRUE(cfg$correctness_stop_on_fail) && any(out$status %in% c("error", "mismatch"))) {
    stop("Correctness preflight failed; review correctness_preflight.csv", call. = FALSE)
  }

  out
}

write_session_info <- function(path_out) {
  con <- file(path_out, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(capture.output(sessionInfo()), con = con)
}

read_linux_mem_bytes <- function() {
  if (!file.exists("/proc/meminfo")) return(NA_real_)
  s <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character())
  hit <- grep("^MemTotal:\\s+", s, value = TRUE)
  if (length(hit) == 0L) return(NA_real_)
  kb <- suppressWarnings(as.numeric(sub("^MemTotal:\\s+([0-9]+)\\s+kB$", "\\1", hit[[1L]])))
  if (is.na(kb)) return(NA_real_)
  kb * 1024
}

read_darwin_sysctl_num <- function(name) {
  if (!identical(Sys.info()[["sysname"]], "Darwin")) return(NA_real_)
  out <- tryCatch(
    suppressWarnings(system2("sysctl", c("-n", name), stdout = TRUE, stderr = FALSE)),
    error = function(e) character()
  )
  if (length(out) == 0L) return(NA_real_)
  val <- suppressWarnings(as.numeric(out[[1L]]))
  if (is.na(val)) NA_real_ else val
}

detect_total_mem_bytes <- function() {
  linux <- read_linux_mem_bytes()
  if (!is.na(linux) && linux > 0) return(linux)

  darwin <- read_darwin_sysctl_num("hw.memsize")
  if (!is.na(darwin) && darwin > 0) return(darwin)

  NA_real_
}

detect_cpu_model <- function() {
  sysname <- Sys.info()[["sysname"]]

  if (identical(sysname, "Linux") && file.exists("/proc/cpuinfo")) {
    s <- tryCatch(readLines("/proc/cpuinfo", warn = FALSE), error = function(e) character())
    hit <- grep("^model name\\s*:", s, value = TRUE)
    if (length(hit) > 0L) {
      return(trimws(sub("^[^:]+:", "", hit[[1L]])))
    }
  }

  if (identical(sysname, "Darwin")) {
    out <- tryCatch(
      suppressWarnings(system2("sysctl", c("-n", "machdep.cpu.brand_string"), stdout = TRUE, stderr = FALSE)),
      error = function(e) character()
    )
    if (length(out) > 0L && nzchar(out[[1L]])) {
      return(trimws(out[[1L]]))
    }
  }

  NA_character_
}

write_host_info <- function(path_out, cfg) {
  si <- Sys.info()
  logical_cores <- suppressWarnings(as.integer(parallel::detectCores(logical = TRUE)))
  physical_cores <- suppressWarnings(as.integer(parallel::detectCores(logical = FALSE)))
  mem_bytes <- detect_total_mem_bytes()

  df <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    nodename = unname(si[["nodename"]]),
    sysname = unname(si[["sysname"]]),
    release = unname(si[["release"]]),
    version = unname(si[["version"]]),
    machine = unname(si[["machine"]]),
    cpu_model = detect_cpu_model(),
    logical_cores = if (is.na(logical_cores)) NA_integer_ else as.integer(logical_cores),
    physical_cores = if (is.na(physical_cores)) NA_integer_ else as.integer(physical_cores),
    detected_cores = as.integer(cfg$detected_cores),
    total_mem_bytes = if (is.na(mem_bytes)) NA_real_ else as.numeric(mem_bytes),
    total_mem_gb = if (is.na(mem_bytes)) NA_real_ else round(as.numeric(mem_bytes) / 1024^3, 3),
    profile = as.character(cfg$profile),
    bp_backend = as.character(cfg$bp_backend),
    stringsAsFactors = FALSE
  )

  write.csv(df, file = path_out, row.names = FALSE)
  invisible(df)
}

write_env_vars <- function(path_out, vars = c(
  "OMP_NUM_THREADS", "OMP_PROC_BIND", "OMP_PLACES", "OMP_SCHEDULE",
  "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "BLIS_NUM_THREADS",
  "VECLIB_MAXIMUM_THREADS", "SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE",
  "SLURM_JOB_ID", "SLURM_NODELIST"
)) {
  vals <- Sys.getenv(vars, unset = NA_character_)
  df <- data.frame(
    variable = vars,
    value = unname(vals),
    stringsAsFactors = FALSE
  )
  write.csv(df, file = path_out, row.names = FALSE)
  invisible(df)
}

write_config <- function(path_out, cfg, opts) {
  con <- file(path_out, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("# Parsed options", con)
  for (nm in names(opts)) {
    writeLines(sprintf("%s: %s", nm, paste(opts[[nm]], collapse = ",")), con)
  }

  writeLines("\n# Effective configuration", con)
  for (nm in names(cfg)) {
    writeLines(sprintf("%s: %s", nm, paste(cfg[[nm]], collapse = ",")), con)
  }
}

write_checkpoint <- function(summary_rows, iter_rows, out_dir) {
  if (length(summary_rows) == 0L || length(iter_rows) == 0L) {
    return(invisible(NULL))
  }

  summary_df <- do.call(rbind, summary_rows)
  iter_df <- do.call(rbind, iter_rows)

  summary_df <- summary_df[order(
    summary_df$outer_run,
    summary_df$scenario,
    summary_df$workload,
    summary_df$comparison_track,
    summary_df$method_family,
    summary_df$seqqual_mode,
    summary_df$threads_effective,
    summary_df$bp_workers
  ), ]
  iter_df <- iter_df[order(iter_df$outer_run, iter_df$case_id, iter_df$iteration), ]

  write.csv(summary_df, file = file.path(out_dir, "summary.csv"), row.names = FALSE)
  write.csv(iter_df, file = file.path(out_dir, "iterations.csv"), row.names = FALSE)

  invisible(list(summary = summary_df, iterations = iter_df))
}

build_outer_run_summary <- function(summary_df, outer_runs_planned = 1L) {
  if (is.null(summary_df) || nrow(summary_df) == 0L) return(summary_df)

  keep_cols <- c(
    "case_id", "scenario", "workload", "method", "method_family", "comparison_track",
    "seqqual_mode", "auto_threads", "threads_requested", "threads_effective",
    "bp_workers_requested", "bp_workers", "bp_workers_effective",
    "total_threads_requested", "total_threads", "total_threads_effective",
    "bp_backend", "n_files", "n_records", "total_mb", "iterations_planned"
  )

  out <- lapply(sort(unique(summary_df$case_id)), function(case_id_one) {
    rows <- summary_df[summary_df$case_id == case_id_one, , drop = FALSE]
    valid <- rows$median_s[is.finite(rows$median_s)]
    elapsed_stats <- summarize_elapsed(valid)
    completed <- length(valid)
    status <- if (completed == outer_runs_planned) "ok" else if (completed > 0L) "partial" else "error"
    error_message <- paste(unique(na.omit(rows$error_message)), collapse = " | ")
    if (!nzchar(error_message)) error_message <- NA_character_

    base <- rows[1L, keep_cols, drop = FALSE]
    base$outer_runs_planned <- as.integer(outer_runs_planned)
    base$outer_runs_completed <- as.integer(completed)
    base$status <- status
    base$error_message <- error_message
    base$time_stat_basis <- if (outer_runs_planned > 1L) "outer_run_median_s" else "single_run_median_s"
    base$min_s <- elapsed_stats$min_s
    base$median_s <- elapsed_stats$median_s
    base$mean_s <- elapsed_stats$mean_s
    base$max_s <- elapsed_stats$max_s
    base$sd_s <- elapsed_stats$sd_s
    base$sem_s <- elapsed_stats$sem_s
    base$mad_s <- elapsed_stats$mad_s
    base$p25_s <- elapsed_stats$p25_s
    base$p75_s <- elapsed_stats$p75_s
    base$iqr_s <- elapsed_stats$iqr_s
    base$cv_s <- elapsed_stats$cv_s
    base$records_per_s <- if (!is.na(base$n_records[[1L]]) && is.finite(base$median_s[[1L]]) && base$median_s[[1L]] > 0) base$n_records[[1L]] / base$median_s[[1L]] else NA_real_
    base$records_per_s_mean <- if (!is.na(base$n_records[[1L]]) && is.finite(base$mean_s[[1L]]) && base$mean_s[[1L]] > 0) base$n_records[[1L]] / base$mean_s[[1L]] else NA_real_
    base$mb_per_s <- if (!is.na(base$total_mb[[1L]]) && is.finite(base$median_s[[1L]]) && base$median_s[[1L]] > 0) base$total_mb[[1L]] / base$median_s[[1L]] else NA_real_
    base$mb_per_s_mean <- if (!is.na(base$total_mb[[1L]]) && is.finite(base$mean_s[[1L]]) && base$mean_s[[1L]] > 0) base$total_mb[[1L]] / base$mean_s[[1L]] else NA_real_
    base
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

write_resolved_plans <- function(case_plan_df, out_dir) {
  resolved_df <- case_plan_df[, c(
    "outer_run", "case_id", "global_execution_rank", "scenario", "workload",
    "method", "comparison_track", "auto_threads", "threads_requested",
    "threads_effective", "bp_workers_requested", "bp_workers_effective",
    "total_threads_requested", "total_threads_effective", "bp_backend"
  )]
  resolved_df$effective_thread_fraction <- ifelse(
    is.finite(resolved_df$total_threads_requested) & resolved_df$total_threads_requested > 0,
    resolved_df$total_threads_effective / resolved_df$total_threads_requested,
    NA_real_
  )
  write.csv(resolved_df, file = file.path(out_dir, "resolved_plans.csv"), row.names = FALSE)
  invisible(resolved_df)
}

value_key <- function(x, missing = "none") {
  out <- as.character(x)
  out[is.na(out) | !nzchar(out)] <- missing
  out
}

build_publication_summary <- function(summary_df) {
  if (is.null(summary_df) || nrow(summary_df) == 0L) return(summary_df)

  df <- summary_df
  rank_proxy <- if ("execution_rank" %in% colnames(df)) df$execution_rank else df$case_id
  df$seqqual_mode_label <- value_key(df$seqqual_mode, missing = "none")
  df$bp_backend_label <- value_key(df$bp_backend, missing = "serial")
  df$resource_plan <- paste0(df$bp_workers_effective, "x", df$threads_effective)
  df$throughput_gb_per_s <- ifelse(is.finite(df$mb_per_s), df$mb_per_s / 1024, NA_real_)
  df$throughput_gb_per_s_mean <- ifelse(is.finite(df$mb_per_s_mean), df$mb_per_s_mean / 1024, NA_real_)
  df$effective_thread_fraction <- ifelse(
    is.finite(df$total_threads_requested) & df$total_threads_requested > 0,
    df$total_threads_effective / df$total_threads_requested,
    NA_real_
  )

  valid <- df$status == "ok" & is.finite(df$median_s)
  df$baseline_case_id_same_method <- NA_integer_
  df$baseline_total_threads_same_method <- NA_integer_
  df$baseline_median_s_same_method <- NA_real_
  df$speedup_vs_same_method_baseline <- NA_real_
  df$parallel_efficiency_vs_same_method_baseline <- NA_real_

  method_group_key <- interaction(
    df$scenario,
    df$workload,
    df$method,
    df$method_family,
    df$comparison_track,
    value_key(df$seqqual_mode, missing = "none"),
    df$auto_threads,
    df$n_files,
    value_key(df$bp_backend, missing = "serial"),
    drop = TRUE,
    lex.order = TRUE
  )

  for (key in unique(method_group_key[valid])) {
    idx <- which(method_group_key == key & valid)
    idx <- idx[order(
      df$total_threads_effective[idx],
      df$bp_workers_effective[idx],
      df$threads_effective[idx],
      rank_proxy[idx]
    )]
    baseline_idx <- idx[[1L]]
    baseline_threads <- df$total_threads_effective[[baseline_idx]]
    baseline_median <- df$median_s[[baseline_idx]]

    df$baseline_case_id_same_method[idx] <- df$case_id[[baseline_idx]]
    df$baseline_total_threads_same_method[idx] <- baseline_threads
    df$baseline_median_s_same_method[idx] <- baseline_median
    df$speedup_vs_same_method_baseline[idx] <- baseline_median / df$median_s[idx]
    df$parallel_efficiency_vs_same_method_baseline[idx] <- df$speedup_vs_same_method_baseline[idx] /
      (df$total_threads_effective[idx] / baseline_threads)
  }

  df$best_case_id_matched_resources <- NA_integer_
  df$best_median_s_matched_resources <- NA_real_
  df$slowdown_vs_best_matched_resources <- NA_real_
  df$rank_matched_resources <- NA_integer_

  matched_group_key <- interaction(
    df$scenario,
    df$workload,
    df$n_files,
    df$total_threads_effective,
    df$bp_workers_effective,
    value_key(df$seqqual_mode, missing = "none"),
    value_key(df$bp_backend, missing = "serial"),
    drop = TRUE,
    lex.order = TRUE
  )

  for (key in unique(matched_group_key[valid])) {
    idx <- which(matched_group_key == key & valid)
    idx <- idx[order(df$median_s[idx], rank_proxy[idx])]
    best_idx <- idx[[1L]]
    df$best_case_id_matched_resources[idx] <- df$case_id[[best_idx]]
    df$best_median_s_matched_resources[idx] <- df$median_s[[best_idx]]
    df$slowdown_vs_best_matched_resources[idx] <- df$median_s[idx] / df$median_s[[best_idx]]
    df$rank_matched_resources[idx] <- seq_along(idx)
  }

  df
}

normalize_benchmark_summary <- function(df) {
  if (is.null(df) || nrow(df) == 0L) return(df)

  if (!"case_id" %in% colnames(df)) df$case_id <- seq_len(nrow(df))
  if (!"method_family" %in% colnames(df)) df$method_family <- df$method
  if (!"comparison_track" %in% colnames(df)) df$comparison_track <- "fair"
  if (!"seqqual_mode" %in% colnames(df)) df$seqqual_mode <- NA_character_
  if (!"auto_threads" %in% colnames(df)) df$auto_threads <- FALSE
  if (!"threads_requested" %in% colnames(df) && "threads_effective" %in% colnames(df)) df$threads_requested <- df$threads_effective
  if (!"threads_effective" %in% colnames(df) && "threads_requested" %in% colnames(df)) df$threads_effective <- df$threads_requested
  if (!"bp_workers_requested" %in% colnames(df) && "bp_workers" %in% colnames(df)) df$bp_workers_requested <- df$bp_workers
  if (!"bp_workers_effective" %in% colnames(df) && "bp_workers" %in% colnames(df)) df$bp_workers_effective <- df$bp_workers
  if (!"bp_workers" %in% colnames(df) && "bp_workers_effective" %in% colnames(df)) df$bp_workers <- df$bp_workers_effective
  if (!"total_threads_requested" %in% colnames(df) && "total_threads" %in% colnames(df)) df$total_threads_requested <- df$total_threads
  if (!"total_threads_effective" %in% colnames(df) && "total_threads" %in% colnames(df)) df$total_threads_effective <- df$total_threads
  if (!"total_threads" %in% colnames(df) && "total_threads_effective" %in% colnames(df)) df$total_threads <- df$total_threads_effective
  if (!"status" %in% colnames(df)) df$status <- "ok"
  if (!"median_s" %in% colnames(df)) stop("Benchmark summary is missing median_s", call. = FALSE)
  if (!"n_files" %in% colnames(df)) df$n_files <- 1L
  if (!"n_records" %in% colnames(df)) df$n_records <- NA_real_
  if (!"total_mb" %in% colnames(df)) df$total_mb <- NA_real_

  df$comparison_track <- ifelse(is.na(df$comparison_track) | !nzchar(df$comparison_track), "fair", as.character(df$comparison_track))
  df$seqqual_mode <- ifelse(is.na(df$seqqual_mode) | !nzchar(df$seqqual_mode), NA_character_, as.character(df$seqqual_mode))
  df$threads_effective <- as.integer(df$threads_effective)
  df$bp_workers_effective <- as.integer(df$bp_workers_effective)
  df$n_files <- as.integer(df$n_files)
  df$median_s <- as.numeric(df$median_s)
  df$n_records <- as.numeric(df$n_records)
  df$total_mb <- as.numeric(df$total_mb)
  df
}

read_benchmark_summary <- function(results_dir) {
  candidates <- c(
    file.path(results_dir, "summary_qmd.csv"),
    file.path(results_dir, "outer_run_summary.csv"),
    file.path(results_dir, "summary.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop("No summary_qmd.csv, outer_run_summary.csv, or summary.csv found in baseline results dir", call. = FALSE)
  }

  df <- utils::read.csv(hit[[1L]], stringsAsFactors = FALSE, check.names = FALSE)
  df <- normalize_benchmark_summary(df)
  attr(df, "source_path") <- hit[[1L]]
  df
}

make_benchmark_match_key <- function(df) {
  paste(
    df$scenario,
    df$workload,
    df$method,
    df$comparison_track,
    ifelse(is.na(df$seqqual_mode), "", df$seqqual_mode),
    df$threads_effective,
    df$bp_workers_effective,
    df$n_files,
    sep = "||"
  )
}

compare_against_baseline <- function(current_df, baseline_df, slowdown_threshold_pct = 5, scope = c("bamscale", "all")) {
  scope <- match.arg(scope)
  current_df <- normalize_benchmark_summary(current_df)
  baseline_df <- normalize_benchmark_summary(baseline_df)

  current_df$key_match <- make_benchmark_match_key(current_df)
  baseline_df$key_match <- make_benchmark_match_key(baseline_df)

  baseline_counts <- table(baseline_df$key_match)
  baseline_idx <- match(current_df$key_match, baseline_df$key_match)
  baseline_found <- !is.na(baseline_idx)

  out <- current_df[, c(
    "case_id", "scenario", "workload", "method", "method_family", "comparison_track",
    "seqqual_mode", "auto_threads", "threads_effective", "bp_workers_effective",
    "n_files", "n_records", "total_mb", "status", "median_s"
  ), drop = FALSE]
  names(out)[names(out) == "median_s"] <- "current_median_s"
  names(out)[names(out) == "status"] <- "current_status"

  out$baseline_case_id <- ifelse(baseline_found, as.integer(baseline_df$case_id[baseline_idx]), NA_integer_)
  out$baseline_status <- ifelse(baseline_found, as.character(baseline_df$status[baseline_idx]), NA_character_)
  out$baseline_median_s <- ifelse(baseline_found, as.numeric(baseline_df$median_s[baseline_idx]), NA_real_)
  out$baseline_n_records <- ifelse(baseline_found, as.numeric(baseline_df$n_records[baseline_idx]), NA_real_)
  out$baseline_total_mb <- ifelse(baseline_found, as.numeric(baseline_df$total_mb[baseline_idx]), NA_real_)
  out$baseline_key_duplicates <- as.integer(ifelse(baseline_found, baseline_counts[current_df$key_match], 0L))
  out$baseline_match_found <- baseline_found
  out$delta_s <- out$current_median_s - out$baseline_median_s
  out$delta_pct <- ifelse(
    is.finite(out$baseline_median_s) & out$baseline_median_s > 0 & is.finite(out$current_median_s),
    100 * (out$current_median_s - out$baseline_median_s) / out$baseline_median_s,
    NA_real_
  )
  out$speedup_vs_baseline <- ifelse(
    is.finite(out$current_median_s) & out$current_median_s > 0 & is.finite(out$baseline_median_s),
    out$baseline_median_s / out$current_median_s,
    NA_real_
  )
  out$same_n_records <- ifelse(
    baseline_found & is.finite(out$n_records) & is.finite(out$baseline_n_records),
    out$n_records == out$baseline_n_records,
    NA
  )
  out$same_total_mb <- ifelse(
    baseline_found & is.finite(out$total_mb) & is.finite(out$baseline_total_mb),
    abs(out$total_mb - out$baseline_total_mb) < 1e-9,
    NA
  )
  out$input_match <- ifelse(
    baseline_found & !is.na(out$same_n_records) & !is.na(out$same_total_mb),
    out$same_n_records & out$same_total_mb,
    NA
  )
  out$regression_threshold_pct <- as.numeric(slowdown_threshold_pct)
  out$comparison_status <- ifelse(
    !out$baseline_match_found,
    "baseline_missing",
    ifelse(
      out$baseline_key_duplicates > 1L,
      "baseline_duplicate",
      ifelse(
        !is.finite(out$current_median_s) | out$current_status != "ok",
        "current_failed",
        ifelse(
          !is.finite(out$baseline_median_s) | out$baseline_status != "ok",
          "baseline_failed",
          ifelse(
            out$delta_pct > slowdown_threshold_pct,
            "regressed",
            ifelse(out$delta_pct < -slowdown_threshold_pct, "improved", "flat")
          )
        )
      )
    )
  )
  out$gate_scope <- if (scope == "all") TRUE else out$method_family == "BamScale"
  out$gate_regression <- out$gate_scope & out$comparison_status == "regressed"
  out[order(out$scenario, out$workload, out$method_family, out$method, out$bp_workers_effective, out$threads_effective), , drop = FALSE]
}

write_baseline_comparison <- function(current_df, cfg, out_dir) {
  empty_df <- data.frame(
    case_id = integer(),
    scenario = character(),
    workload = character(),
    method = character(),
    method_family = character(),
    comparison_track = character(),
    seqqual_mode = character(),
    auto_threads = logical(),
    threads_effective = integer(),
    bp_workers_effective = integer(),
    n_files = integer(),
    n_records = numeric(),
    total_mb = numeric(),
    current_status = character(),
    current_median_s = numeric(),
    baseline_case_id = integer(),
    baseline_status = character(),
    baseline_median_s = numeric(),
    baseline_n_records = numeric(),
    baseline_total_mb = numeric(),
    baseline_key_duplicates = integer(),
    baseline_match_found = logical(),
    delta_s = numeric(),
    delta_pct = numeric(),
    speedup_vs_baseline = numeric(),
    same_n_records = logical(),
    same_total_mb = logical(),
    input_match = logical(),
    regression_threshold_pct = numeric(),
    comparison_status = character(),
    gate_scope = logical(),
    gate_regression = logical(),
    stringsAsFactors = FALSE
  )

  if (is.null(cfg$baseline_results_dir) || !nzchar(cfg$baseline_results_dir)) {
    write.csv(empty_df, file = file.path(out_dir, "baseline_comparison.csv"), row.names = FALSE)
    write.csv(empty_df, file = file.path(out_dir, "baseline_regressions.csv"), row.names = FALSE)
    return(list(comparison = empty_df, regressions = empty_df, baseline_source = NA_character_))
  }

  baseline_df <- read_benchmark_summary(cfg$baseline_results_dir)
  comparison_df <- compare_against_baseline(
    current_df = current_df,
    baseline_df = baseline_df,
    slowdown_threshold_pct = cfg$regression_threshold_pct,
    scope = cfg$regression_scope
  )
  regressions_df <- comparison_df[comparison_df$gate_regression, , drop = FALSE]

  write.csv(comparison_df, file = file.path(out_dir, "baseline_comparison.csv"), row.names = FALSE)
  write.csv(regressions_df, file = file.path(out_dir, "baseline_regressions.csv"), row.names = FALSE)

  list(
    comparison = comparison_df,
    regressions = regressions_df,
    baseline_source = attr(baseline_df, "source_path")
  )
}

write_publication_tables <- function(summary_df, iter_df, out_dir) {
  outer_runs_planned <- if ("outer_run" %in% colnames(summary_df)) max(1L, suppressWarnings(as.integer(max(summary_df$outer_run, na.rm = TRUE)))) else 1L
  outer_summary_df <- build_outer_run_summary(summary_df, outer_runs_planned = outer_runs_planned)
  write.csv(outer_summary_df, file = file.path(out_dir, "outer_run_summary.csv"), row.names = FALSE)

  summary_qmd <- build_publication_summary(outer_summary_df)
  write.csv(summary_qmd, file = file.path(out_dir, "summary_qmd.csv"), row.names = FALSE)

  counts_df <- unique(summary_qmd[, c("scenario", "workload", "n_files", "n_records", "total_mb")])
  counts_df <- counts_df[order(counts_df$scenario, counts_df$workload, counts_df$n_files), , drop = FALSE]
  rownames(counts_df) <- NULL
  write.csv(counts_df, file = file.path(out_dir, "reference_counts.csv"), row.names = FALSE)

  iter_ok <- iter_df[iter_df$status == "ok" & is.finite(iter_df$elapsed_s), , drop = FALSE]
  if (nrow(iter_ok) > 0L) {
    case_keys <- unique(iter_ok[, c("case_id", "outer_run"), drop = FALSE])
    iteration_stats <- do.call(rbind, lapply(seq_len(nrow(case_keys)), function(i) {
      case_id_one <- case_keys$case_id[[i]]
      outer_run_one <- case_keys$outer_run[[i]]
      vals <- iter_ok$elapsed_s[iter_ok$case_id == case_id_one & iter_ok$outer_run == outer_run_one]
      data.frame(
        case_id = as.integer(case_id_one),
        outer_run = as.integer(outer_run_one),
        iteration_mean_s = mean(vals),
        iteration_sd_s = if (length(vals) > 1L) stats::sd(vals) else NA_real_,
        iteration_sem_s = if (length(vals) > 1L) stats::sd(vals) / sqrt(length(vals)) else NA_real_,
        iteration_min_s = min(vals),
        iteration_median_s = median(vals),
        iteration_max_s = max(vals),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    iteration_stats <- data.frame(
      case_id = integer(),
      outer_run = integer(),
      iteration_mean_s = numeric(),
      iteration_sd_s = numeric(),
      iteration_sem_s = numeric(),
      iteration_min_s = numeric(),
      iteration_median_s = numeric(),
      iteration_max_s = numeric(),
      stringsAsFactors = FALSE
    )
  }
  write.csv(iteration_stats, file = file.path(out_dir, "iteration_stats.csv"), row.names = FALSE)

  invisible(list(
    outer_run_summary = outer_summary_df,
    summary_qmd = summary_qmd,
    reference_counts = counts_df,
    iteration_stats = iteration_stats
  ))
}

write_artifact_manifest <- function(out_dir) {
  entries <- list(
    list(file = "summary.csv", category = "results", description = "Per-case benchmark summary with timing and throughput statistics."),
    list(file = "outer_run_summary.csv", category = "results", description = "Per-case aggregate summary across independent outer runs."),
    list(file = "summary_qmd.csv", category = "results", description = "QMD-ready summary with derived scaling, speedup, and matched-resource ranking fields."),
    list(file = "iterations.csv", category = "results", description = "Per-iteration elapsed times and status for every case."),
    list(file = "iteration_stats.csv", category = "results", description = "Iteration-level aggregate statistics per case."),
    list(file = "baseline_comparison.csv", category = "validation", description = "Matched-case comparison against a reference benchmark run, including slowdown percentages and input-drift fields."),
    list(file = "baseline_regressions.csv", category = "validation", description = "Subset of matched baseline comparisons that exceeded the configured slowdown threshold within the selected gate scope."),
    list(file = "case_plan.csv", category = "design", description = "Benchmark execution plan after optional shuffling, including requested and effective resources."),
    list(file = "resolved_plans.csv", category = "design", description = "Explicit requested-versus-effective resource plans, including adaptive auto_threads resolutions."),
    list(file = "reference_counts.csv", category = "design", description = "Reference record counts and input sizes used for throughput normalization."),
    list(file = "correctness_preflight.csv", category = "validation", description = "Comparator-based correctness preflight on a normalized subset of the representative file."),
    list(file = "files.csv", category = "inputs", description = "Benchmark file manifest with source, size, and selection flags."),
    list(file = "host_info.csv", category = "metadata", description = "Host machine and runtime metadata for the run."),
    list(file = "env_vars.csv", category = "metadata", description = "Threading-related environment variables captured for the run."),
    list(file = "config.txt", category = "metadata", description = "Parsed CLI options and effective benchmark configuration."),
    list(file = "sessionInfo.txt", category = "metadata", description = "R session information for package and platform provenance."),
    list(file = "plot_single_scaling.png", category = "plots", description = "Auto-generated single-file BamScale scaling plot, if ggplot2 is available."),
    list(file = "plot_multi_scaling.png", category = "plots", description = "Auto-generated multi-file scaling plot, if ggplot2 is available.")
  )

  manifest_df <- do.call(rbind, lapply(entries, function(x) {
    path_one <- file.path(out_dir, x$file)
    data.frame(
      file = x$file,
      category = x$category,
      exists = file.exists(path_one),
      description = x$description,
      stringsAsFactors = FALSE
    )
  }))

  write.csv(manifest_df, file = file.path(out_dir, "artifact_manifest.csv"), row.names = FALSE)
  invisible(manifest_df)
}

plot_if_possible <- function(summary_df, out_dir) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) return(invisible(FALSE))

  ok <- summary_df[summary_df$status == "ok", , drop = FALSE]
  if (nrow(ok) == 0L) return(invisible(FALSE))

  p1_df <- ok[ok$scenario == "single" & ok$method_family == "BamScale", , drop = FALSE]
  if (nrow(p1_df) > 0L) {
    p <- ggplot2::ggplot(
      p1_df,
      ggplot2::aes(x = threads_effective, y = median_s, color = workload)
    ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2) +
      ggplot2::scale_x_continuous(breaks = sort(unique(p1_df$threads_effective))) +
      ggplot2::labs(
        title = "Single-file scaling (BamScale)",
        x = "Threads",
        y = "Median elapsed (s)",
        color = "Workload"
      ) +
      ggplot2::theme_minimal(base_size = 12)

    ggplot2::ggsave(
      filename = file.path(out_dir, "plot_single_scaling.png"),
      plot = p,
      width = 10,
      height = 6,
      dpi = 160
    )
  }

  p2_df <- ok[ok$scenario == "multi", , drop = FALSE]
  if (nrow(p2_df) > 0L) {
    p <- ggplot2::ggplot(
      p2_df,
      ggplot2::aes(x = bp_workers, y = median_s, color = method)
    ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 1.8) +
      ggplot2::facet_wrap(~workload, scales = "free_y") +
      ggplot2::scale_x_continuous(breaks = sort(unique(p2_df$bp_workers))) +
      ggplot2::labs(
        title = "Multi-file scaling",
        x = "BiocParallel workers",
        y = "Median elapsed (s)",
        color = "Method"
      ) +
      ggplot2::theme_minimal(base_size = 12)

    ggplot2::ggsave(
      filename = file.path(out_dir, "plot_multi_scaling.png"),
      plot = p,
      width = 11,
      height = 7,
      dpi = 160
    )
  }

  invisible(TRUE)
}

args <- commandArgs(trailingOnly = TRUE)
opts <- parse_args(args)

detected_cores <- as.integer(.bamscale_detect_cores())
if (length(detected_cores) != 1L || is.na(detected_cores) || detected_cores < 1L) {
  detected_cores <- 1L
}

profile <- if (!is.null(opts[["profile"]])) tolower(opts[["profile"]]) else "bamscale_showcase"
if (!profile %in% c("bamscale_showcase", "balanced", "full")) {
  stop("--profile must be one of: bamscale_showcase, balanced, full", call. = FALSE)
}

if (profile == "bamscale_showcase") {
  default_n_files <- 12L
  default_threads <- c(1L, 2L, 4L, 8L, 12L, 16L, 24L, 32L, 48L)
  default_workers <- c(1L, 2L, 4L, 8L, 12L, 16L)
  default_include_seqqual <- FALSE
  default_seqqual_compact <- FALSE
  default_include_galignments <- FALSE
  default_include_auto_threads <- FALSE
  default_include_rsamtools <- TRUE
  default_multi_workloads <- c("step1")
} else if (profile == "balanced") {
  default_n_files <- 12L
  default_threads <- c(1L, 2L, 4L, 8L, 12L, 16L, 24L, 32L, 48L)
  default_workers <- c(1L, 2L, 4L, 8L, 12L, 16L)
  default_include_seqqual <- TRUE
  default_seqqual_compact <- TRUE
  default_include_galignments <- TRUE
  default_include_auto_threads <- TRUE
  default_include_rsamtools <- TRUE
  default_multi_workloads <- c("step1", "seqqual")
} else {
  default_n_files <- 16L
  default_threads <- c(1L, 2L, 4L, 8L, 12L, 16L, 24L, 32L, 48L, 64L, 96L)
  default_workers <- c(1L, 2L, 4L, 8L, 12L, 16L, 24L, 32L, 48L, 64L, 96L)
  default_include_seqqual <- TRUE
  default_seqqual_compact <- TRUE
  default_include_galignments <- TRUE
  default_include_auto_threads <- TRUE
  default_include_rsamtools <- TRUE
  default_multi_workloads <- c("step1", "seqqual", "galignments")
}

default_budget_threads <- if (profile == "full") {
  as.integer(max(1L, detected_cores))
} else {
  as.integer(min(48L, detected_cores))
}
requested_budget_threads <- as_int_scalar(opts[["budget-threads"]], default_budget_threads)
requested_max_threads <- as_int_scalar(opts[["max-threads"]], requested_budget_threads)
max_threads <- as.integer(max(1L, min(requested_max_threads, detected_cores)))
if (requested_max_threads > detected_cores) {
  message("Capping --max-threads to detected cores: ", detected_cores)
}

cfg <- list(
  outdir = if (!is.null(opts$outdir)) opts$outdir else "benchmark_results",
  profile = profile,
  detected_cores = as.integer(detected_cores),
  budget_threads = as.integer(requested_budget_threads),
  n_files = as_int_scalar(opts[["n-files"]], default_n_files),
  outer_runs = as_int_scalar(opts[["outer-runs"]], 1L),
  iterations = as_int_scalar(opts$iterations, 3L),
  warmup = as_bool(opts$warmup, TRUE),
  max_threads = as.integer(max_threads),
  threads = as_int_vec(opts$threads, default_threads),
  workers = as_int_vec(opts$workers, default_workers),
  bp_backend = if (!is.null(opts[["bp-backend"]])) tolower(opts[["bp-backend"]]) else "snow",
  include_seqqual = as_bool(opts[["include-seqqual"]], default_include_seqqual),
  seqqual_compact = as_bool(opts[["seqqual-compact"]], default_seqqual_compact),
  include_galignments = as_bool(opts[["include-galignments"]], default_include_galignments),
  include_auto_threads = as_bool(opts[["include-auto-threads"]], default_include_auto_threads),
  include_multi = as_bool(opts[["include-multi"]], TRUE),
  include_rsamtools = as_bool(opts[["include-rsamtools"]], default_include_rsamtools),
  ensure_index = as_bool(opts[["ensure-index"]], TRUE),
  allow_repeat_files = as_bool(opts[["allow-repeat-files"]], FALSE),
  recursive = as_bool(opts$recursive, TRUE),
  multi_workloads = as_chr_vec(opts[["multi-workloads"]], default_multi_workloads),
  compute_records = as_bool(opts[["compute-records"]], TRUE),
  shuffle_cases = as_bool(opts[["shuffle-cases"]], TRUE),
  shuffle_seed = as_int_scalar(opts[["shuffle-seed"]], 1L),
  correctness_preflight = as_bool(opts[["correctness-preflight"]], TRUE),
  correctness_window_bp = as_int_scalar(opts[["correctness-window-bp"]], 1000L),
  correctness_stop_on_fail = as_bool(opts[["correctness-stop-on-fail"]], FALSE),
  baseline_results_dir = if (!is.null(opts[["baseline-results-dir"]])) opts[["baseline-results-dir"]] else "",
  regression_threshold_pct = as_num_scalar(opts[["regression-threshold-pct"]], 5),
  regression_stop_on_fail = as_bool(opts[["regression-stop-on-fail"]], FALSE),
  regression_scope = if (!is.null(opts[["regression-scope"]])) tolower(opts[["regression-scope"]]) else "bamscale"
)
if (!cfg$bp_backend %in% c("snow", "multicore")) {
  stop("--bp-backend must be one of: snow,multicore", call. = FALSE)
}
if (!cfg$regression_scope %in% c("bamscale", "all")) {
  stop("--regression-scope must be one of: bamscale,all", call. = FALSE)
}
if (!is.finite(cfg$regression_threshold_pct) || cfg$regression_threshold_pct < 0) {
  stop("--regression-threshold-pct must be a non-negative number", call. = FALSE)
}
if (nzchar(cfg$baseline_results_dir)) {
  cfg$baseline_results_dir <- normalizePath(cfg$baseline_results_dir, winslash = "/", mustWork = TRUE)
}

cfg$threads <- cfg$threads[cfg$threads >= 1L & cfg$threads <= cfg$max_threads]
if (length(cfg$threads) == 0L) {
  stop("No valid --threads values after applying max thread cap", call. = FALSE)
}

cfg$workers <- cfg$workers[cfg$workers >= 1L & cfg$workers <= cfg$max_threads]
if (length(cfg$workers) == 0L) cfg$workers <- 1L
if (!isTRUE(cfg$include_seqqual)) cfg$seqqual_compact <- FALSE

run_stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(cfg$outdir, paste0("run_", run_stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Output directory: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))

manual_files <- character()
if (!is.null(opts[["bam-files"]])) {
  manual_files <- c(manual_files, split_csv(opts[["bam-files"]]))
}
if (!is.null(opts[["bam-file-list"]])) {
  lst <- readLines(opts[["bam-file-list"]], warn = FALSE)
  manual_files <- c(manual_files, trimws(lst))
}
if (!is.null(opts[["bam-dir"]])) {
  from_dir <- list.files(
    path = opts[["bam-dir"]],
    pattern = "\\.bam$",
    full.names = TRUE,
    recursive = cfg$recursive,
    ignore.case = TRUE
  )
  manual_files <- c(manual_files, from_dir)
}

manual_files <- unique(manual_files[nzchar(manual_files)])

if (length(manual_files) > 0L) {
  file_source <- "manual"
  bam_paths <- normalizePath(manual_files, winslash = "/", mustWork = TRUE)
} else {
  message("No BAM files provided, resolving BAMs from chipseqDBData...")
  file_source <- "chipseqDBData"
  bam_paths <- fetch_chipseqdbdata_bams(n_target = max(cfg$n_files, 2L))
  bam_paths <- normalizePath(bam_paths, winslash = "/", mustWork = TRUE)
}

bam_paths <- unique(bam_paths[file.exists(bam_paths)])
if (length(bam_paths) == 0L) {
  stop("No BAM files available for benchmarking", call. = FALSE)
}

sizes <- file.size(bam_paths)
ord <- order(sizes, decreasing = TRUE)
bam_paths <- bam_paths[ord]

single_file <- bam_paths[[1L]]
multi_files <- select_n_files(
  bam_paths,
  n_files = cfg$n_files,
  allow_repeat = cfg$allow_repeat_files
)

cfg$workers <- cfg$workers[cfg$workers <= length(multi_files)]
if (length(cfg$workers) == 0L) cfg$workers <- 1L

selected_files <- unique(c(single_file, multi_files))
files_df <- data.frame(
  file = selected_files,
  source = file_source,
  selected_for_single = selected_files == single_file,
  selected_for_multi = selected_files %in% multi_files,
  selection_rank = match(selected_files, bam_paths),
  size_bytes = as.numeric(file.size(selected_files)),
  size_mb = round(as.numeric(file.size(selected_files)) / 1024^2, 3),
  has_index = vapply(selected_files, bam_index_exists, logical(1)),
  stringsAsFactors = FALSE
)
write.csv(files_df, file = file.path(out_dir, "files.csv"), row.names = FALSE)

if ((cfg$include_rsamtools || cfg$correctness_preflight) && cfg$ensure_index) {
  message("Ensuring BAM index files for comparator methods and correctness preflight...")
  for (path_one in unique(c(single_file, multi_files))) {
    if (!bam_index_exists(path_one)) {
      message("Indexing: ", path_one)
      try(ensure_bam_index(path_one), silent = FALSE)
    }
  }
}

workloads <- list(
  step1 = list(
    what = c("qname", "flag", "rname", "pos", "mapq", "cigar"),
    as = "data.frame",
    seqqual_mode = "compatible"
  ),
  seqqual = list(
    what = c("qname", "seq", "qual"),
    as = "data.frame",
    seqqual_mode = "compatible"
  ),
  galignments = list(
    what = c("qname", "rname", "pos", "cigar", "strand", "flag"),
    as = "GAlignments",
    seqqual_mode = "compatible"
  )
)

active_workloads <- c("step1")
if (isTRUE(cfg$include_seqqual)) active_workloads <- c(active_workloads, "seqqual")
if (isTRUE(cfg$include_galignments)) active_workloads <- c(active_workloads, "galignments")

cfg$multi_workloads <- unique(intersect(cfg$multi_workloads, active_workloads))
if (length(cfg$multi_workloads) == 0L) cfg$multi_workloads <- "step1"

message("Single-file benchmark target: ", single_file)
message("Multi-file count: ", length(multi_files))
message(paste0("Threads grid (<= max_threads=", cfg$max_threads, "): ", paste(cfg$threads, collapse = ", ")))
message("Workers grid: ", paste(cfg$workers, collapse = ", "))
message("Profile: ", cfg$profile)
message("Budget threads: ", cfg$budget_threads, " | Detected cores: ", cfg$detected_cores)
message("BiocParallel backend: ", cfg$bp_backend)
message("Seq/qual compact track: ", cfg$seqqual_compact)
message("Adaptive auto_threads track: ", cfg$include_auto_threads)
message("Outer runs: ", cfg$outer_runs)
message("Shuffle cases: ", cfg$shuffle_cases, " (seed=", cfg$shuffle_seed, ")")
message("Correctness preflight: ", cfg$correctness_preflight, " (window_bp=", cfg$correctness_window_bp, ")")
if (nzchar(cfg$baseline_results_dir)) {
  message("Baseline comparison: ", cfg$baseline_results_dir)
  message("Regression gate: ", cfg$regression_scope, " | threshold=", cfg$regression_threshold_pct, "% | stop_on_fail=", cfg$regression_stop_on_fail)
}

single_mb <- as.numeric(file.size(single_file) / 1024^2)
multi_mb <- as.numeric(sum(file.size(multi_files)) / 1024^2)

single_records <- setNames(rep(NA_real_, length(active_workloads)), active_workloads)
if (isTRUE(cfg$compute_records)) {
  message("Computing single-file record counts...")
  for (w in active_workloads) {
    spec <- workloads[[w]]
    ref <- BamScale::bam_read(
      file = single_file,
      what = spec$what,
      as = spec$as,
      seqqual_mode = spec$seqqual_mode,
      threads = 1L,
      auto_threads = FALSE
    )
    single_records[[w]] <- count_records(ref)
  }
}

multi_records <- setNames(rep(NA_real_, length(cfg$multi_workloads)), cfg$multi_workloads)
if (isTRUE(cfg$compute_records)) {
  message("Computing multi-file record counts...")
  for (w in cfg$multi_workloads) {
    spec <- workloads[[w]]
    ref <- BamScale::bam_read(
      file = multi_files,
      what = spec$what,
      as = spec$as,
      seqqual_mode = spec$seqqual_mode,
      threads = 1L,
      BPPARAM = NULL,
      auto_threads = FALSE
    )
    multi_records[[w]] <- count_records(ref)
  }
}

message("Running correctness preflight...")
correctness_preflight <- run_correctness_preflight(
  single_file = single_file,
  workloads = workloads,
  active_workloads = active_workloads,
  cfg = cfg,
  out_dir = out_dir
)

summary_rows <- list()
iter_rows <- list()
case_specs <- list()
build_order <- 0L
cfg$file_source <- file_source

make_case_meta <- function(
    scenario,
    workload,
    method,
    method_family,
    comparison_track,
    seqqual_mode,
    auto_threads,
    threads_requested,
    threads_effective,
    bp_workers_requested,
    bp_workers_effective,
    n_files,
    n_records,
    total_mb,
    bp_backend
) {
  data.frame(
    scenario = scenario,
    workload = workload,
    method = method,
    method_family = method_family,
    comparison_track = comparison_track,
    seqqual_mode = if (is.na(seqqual_mode)) NA_character_ else as.character(seqqual_mode),
    auto_threads = as.logical(auto_threads),
    threads_requested = as.integer(threads_requested),
    threads_effective = as.integer(threads_effective),
    bp_workers_requested = as.integer(bp_workers_requested),
    bp_workers_effective = as.integer(bp_workers_effective),
    total_threads_requested = as.integer(max(1L, bp_workers_requested) * max(1L, threads_requested)),
    total_threads_effective = as.integer(max(1L, bp_workers_effective) * max(1L, threads_effective)),
    n_files = as.integer(n_files),
    n_records = as.numeric(n_records),
    total_mb = as.numeric(total_mb),
    bp_backend = if (is.na(bp_backend)) NA_character_ else as.character(bp_backend),
    stringsAsFactors = FALSE
  )
}

append_case <- function(meta, runner) {
  build_order <<- build_order + 1L
  meta$build_order <- as.integer(build_order)
  case_specs[[length(case_specs) + 1L]] <<- list(meta = meta, runner = runner)
  invisible(NULL)
}

for (w in active_workloads) {
  spec <- workloads[[w]]

  for (t in cfg$threads) {
    append_case(
      make_case_meta(
        scenario = "single",
        workload = w,
        method = "BamScale",
        method_family = "BamScale",
        comparison_track = "fair",
        seqqual_mode = if (w == "seqqual") "compatible" else NA_character_,
        auto_threads = FALSE,
        threads_requested = t,
        threads_effective = t,
        bp_workers_requested = 1L,
        bp_workers_effective = 1L,
        n_files = 1L,
        n_records = single_records[[w]],
        total_mb = single_mb,
        bp_backend = NA_character_
      ),
      local({
        w_local <- w
        spec_local <- spec
        t_local <- t
        n_records_local <- single_records[[w]]
        total_mb_local <- single_mb
        function(case_id, execution_rank) {
          run_case(
            scenario = "single",
            workload = w_local,
            method = "BamScale",
            method_family = "BamScale",
            threads_requested = t_local,
            threads_effective = t_local,
            bp_workers_requested = 1L,
            bp_workers = 1L,
            n_files = 1L,
            n_records = n_records_local,
            total_mb = total_mb_local,
            iterations = cfg$iterations,
            warmup = cfg$warmup,
            fun = function() {
              BamScale::bam_read(
                file = single_file,
                what = spec_local$what,
                as = spec_local$as,
                seqqual_mode = spec_local$seqqual_mode,
                threads = t_local,
                BPPARAM = NULL,
                auto_threads = FALSE
              )
            },
            case_id = case_id,
            execution_rank = execution_rank,
            comparison_track = "fair",
            seqqual_mode = if (w_local == "seqqual") "compatible" else NA_character_,
            auto_threads = FALSE,
            bp_backend = NA_character_
          )
        }
      })
    )

    if (w == "seqqual" && isTRUE(cfg$seqqual_compact)) {
      append_case(
        make_case_meta(
          scenario = "single",
          workload = w,
          method = "BamScale (compact seqqual)",
          method_family = "BamScale",
          comparison_track = "optimized",
          seqqual_mode = "compact",
          auto_threads = FALSE,
          threads_requested = t,
          threads_effective = t,
          bp_workers_requested = 1L,
          bp_workers_effective = 1L,
          n_files = 1L,
          n_records = single_records[[w]],
          total_mb = single_mb,
          bp_backend = NA_character_
        ),
        local({
          spec_local <- spec
          t_local <- t
          n_records_local <- single_records[[w]]
          total_mb_local <- single_mb
          function(case_id, execution_rank) {
            run_case(
              scenario = "single",
              workload = "seqqual",
              method = "BamScale (compact seqqual)",
              method_family = "BamScale",
              threads_requested = t_local,
              threads_effective = t_local,
              bp_workers_requested = 1L,
              bp_workers = 1L,
              n_files = 1L,
              n_records = n_records_local,
              total_mb = total_mb_local,
              iterations = cfg$iterations,
              warmup = cfg$warmup,
              fun = function() {
                BamScale::bam_read(
                  file = single_file,
                  what = spec_local$what,
                  as = spec_local$as,
                  seqqual_mode = "compact",
                  threads = t_local,
                  BPPARAM = NULL,
                  auto_threads = FALSE
                )
              },
              case_id = case_id,
              execution_rank = execution_rank,
              comparison_track = "optimized",
              seqqual_mode = "compact",
              auto_threads = FALSE,
              bp_backend = NA_character_
            )
          }
        })
      )
    }
  }

  if (cfg$include_rsamtools && w %in% c("step1", "seqqual") && requireNamespace("Rsamtools", quietly = TRUE)) {
    append_case(
      make_case_meta(
        scenario = "single",
        workload = w,
        method = "Rsamtools::scanBam",
        method_family = "Rsamtools",
        comparison_track = "fair",
        seqqual_mode = if (w == "seqqual") "compatible" else NA_character_,
        auto_threads = FALSE,
        threads_requested = 1L,
        threads_effective = 1L,
        bp_workers_requested = 1L,
        bp_workers_effective = 1L,
        n_files = 1L,
        n_records = single_records[[w]],
        total_mb = single_mb,
        bp_backend = NA_character_
      ),
      local({
        w_local <- w
        spec_local <- spec
        n_records_local <- single_records[[w]]
        total_mb_local <- single_mb
        function(case_id, execution_rank) {
          run_case(
            scenario = "single",
            workload = w_local,
            method = "Rsamtools::scanBam",
            method_family = "Rsamtools",
            threads_requested = 1L,
            threads_effective = 1L,
            bp_workers_requested = 1L,
            bp_workers = 1L,
            n_files = 1L,
            n_records = n_records_local,
            total_mb = total_mb_local,
            iterations = cfg$iterations,
            warmup = cfg$warmup,
            fun = function() {
              param <- Rsamtools::ScanBamParam(what = as.character(spec_local$what))
              Rsamtools::scanBam(single_file, param = param)
            },
            case_id = case_id,
            execution_rank = execution_rank,
            comparison_track = "fair",
            seqqual_mode = if (w_local == "seqqual") "compatible" else NA_character_,
            auto_threads = FALSE,
            bp_backend = NA_character_
          )
        }
      })
    )
  }

  if (w == "galignments" && requireNamespace("GenomicAlignments", quietly = TRUE) && requireNamespace("Rsamtools", quietly = TRUE)) {
    append_case(
      make_case_meta(
        scenario = "single",
        workload = w,
        method = "GenomicAlignments::readGAlignments",
        method_family = "GenomicAlignments",
        comparison_track = "fair",
        seqqual_mode = NA_character_,
        auto_threads = FALSE,
        threads_requested = 1L,
        threads_effective = 1L,
        bp_workers_requested = 1L,
        bp_workers_effective = 1L,
        n_files = 1L,
        n_records = single_records[[w]],
        total_mb = single_mb,
        bp_backend = NA_character_
      ),
      local({
        spec_local <- spec
        n_records_local <- single_records[[w]]
        total_mb_local <- single_mb
        function(case_id, execution_rank) {
          run_case(
            scenario = "single",
            workload = "galignments",
            method = "GenomicAlignments::readGAlignments",
            method_family = "GenomicAlignments",
            threads_requested = 1L,
            threads_effective = 1L,
            bp_workers_requested = 1L,
            bp_workers = 1L,
            n_files = 1L,
            n_records = n_records_local,
            total_mb = total_mb_local,
            iterations = cfg$iterations,
            warmup = cfg$warmup,
            fun = function() {
              param <- Rsamtools::ScanBamParam(what = as.character(spec_local$what))
              GenomicAlignments::readGAlignments(single_file, param = param)
            },
            case_id = case_id,
            execution_rank = execution_rank,
            comparison_track = "fair",
            seqqual_mode = NA_character_,
            auto_threads = FALSE,
            bp_backend = NA_character_
          )
        }
      })
    )
  }
}

if (isTRUE(cfg$include_multi)) {
  for (w in cfg$multi_workloads) {
    spec <- workloads[[w]]

    for (workers in cfg$workers) {
      bp_template <- tryCatch(
        make_bpparam(workers, backend = cfg$bp_backend),
        error = function(e) e
      )

      if (inherits(bp_template, "error")) {
        message("Skipping workers=", workers, " because BPPARAM init failed: ", conditionMessage(bp_template))
        next
      }

      threads_each <- max(1L, floor(cfg$max_threads / workers))
      auto_plan <- if (isTRUE(cfg$include_auto_threads) && workers > 1L) {
        BamScale:::.bamscale_resolve_parallel_plan(
          threads = threads_each,
          BPPARAM = bp_template,
          auto_threads = TRUE,
          n_files = length(multi_files)
        )
      } else {
        NULL
      }
      stop_bpparam(bp_template)

      append_case(
        make_case_meta(
          scenario = "multi",
          workload = w,
          method = "BamScale (balanced budget)",
          method_family = "BamScale",
          comparison_track = "fair",
          seqqual_mode = if (w == "seqqual") "compatible" else NA_character_,
          auto_threads = FALSE,
          threads_requested = threads_each,
          threads_effective = threads_each,
          bp_workers_requested = workers,
          bp_workers_effective = workers,
          n_files = length(multi_files),
          n_records = multi_records[[w]],
          total_mb = multi_mb,
          bp_backend = cfg$bp_backend
        ),
        local({
          w_local <- w
          spec_local <- spec
          workers_local <- workers
          threads_local <- threads_each
          backend_local <- cfg$bp_backend
          n_files_local <- length(multi_files)
          n_records_local <- multi_records[[w]]
          total_mb_local <- multi_mb
          function(case_id, execution_rank) {
            bp <- make_bpparam(workers_local, backend = backend_local)
            on.exit(stop_bpparam(bp), add = TRUE)
            run_case(
              scenario = "multi",
              workload = w_local,
              method = "BamScale (balanced budget)",
              method_family = "BamScale",
              threads_requested = threads_local,
              threads_effective = threads_local,
              bp_workers_requested = workers_local,
              bp_workers = workers_local,
              n_files = n_files_local,
              n_records = n_records_local,
              total_mb = total_mb_local,
              iterations = cfg$iterations,
              warmup = cfg$warmup,
              fun = function() {
                BamScale::bam_read(
                  file = multi_files,
                  what = spec_local$what,
                  as = spec_local$as,
                  seqqual_mode = spec_local$seqqual_mode,
                  threads = threads_local,
                  BPPARAM = bp,
                  auto_threads = FALSE
                )
              },
              case_id = case_id,
              execution_rank = execution_rank,
              comparison_track = "fair",
              seqqual_mode = if (w_local == "seqqual") "compatible" else NA_character_,
              auto_threads = FALSE,
              bp_backend = backend_local
            )
          }
        })
      )

      if (!is.null(auto_plan)) {
        append_case(
          make_case_meta(
            scenario = "multi",
            workload = w,
            method = "BamScale (adaptive auto_threads)",
            method_family = "BamScale",
            comparison_track = "adaptive",
            seqqual_mode = if (w == "seqqual") "compatible" else NA_character_,
            auto_threads = TRUE,
            threads_requested = threads_each,
            threads_effective = auto_plan$threads,
            bp_workers_requested = workers,
            bp_workers_effective = auto_plan$bp_workers,
            n_files = length(multi_files),
            n_records = multi_records[[w]],
            total_mb = multi_mb,
            bp_backend = cfg$bp_backend
          ),
          local({
            w_local <- w
            spec_local <- spec
            workers_local <- workers
            threads_local <- threads_each
            plan_local <- auto_plan
            backend_local <- cfg$bp_backend
            n_files_local <- length(multi_files)
            n_records_local <- multi_records[[w]]
            total_mb_local <- multi_mb
            function(case_id, execution_rank) {
              bp <- make_bpparam(workers_local, backend = backend_local)
              on.exit(stop_bpparam(bp), add = TRUE)
              run_case(
                scenario = "multi",
                workload = w_local,
                method = "BamScale (adaptive auto_threads)",
                method_family = "BamScale",
                threads_requested = threads_local,
                threads_effective = plan_local$threads,
                bp_workers_requested = workers_local,
                bp_workers = plan_local$bp_workers,
                n_files = n_files_local,
                n_records = n_records_local,
                total_mb = total_mb_local,
                iterations = cfg$iterations,
                warmup = cfg$warmup,
                fun = function() {
                  BamScale::bam_read(
                    file = multi_files,
                    what = spec_local$what,
                    as = spec_local$as,
                    seqqual_mode = spec_local$seqqual_mode,
                    threads = threads_local,
                    BPPARAM = bp,
                    auto_threads = TRUE
                  )
                },
                case_id = case_id,
                execution_rank = execution_rank,
                comparison_track = "adaptive",
                seqqual_mode = if (w_local == "seqqual") "compatible" else NA_character_,
                auto_threads = TRUE,
                bp_backend = backend_local
              )
            }
          })
        )
      }

      if (w == "seqqual" && isTRUE(cfg$seqqual_compact)) {
        append_case(
          make_case_meta(
            scenario = "multi",
            workload = w,
            method = "BamScale (compact seqqual budget)",
            method_family = "BamScale",
            comparison_track = "optimized",
            seqqual_mode = "compact",
            auto_threads = FALSE,
            threads_requested = threads_each,
            threads_effective = threads_each,
            bp_workers_requested = workers,
            bp_workers_effective = workers,
            n_files = length(multi_files),
            n_records = multi_records[[w]],
            total_mb = multi_mb,
            bp_backend = cfg$bp_backend
          ),
          local({
            workers_local <- workers
            threads_local <- threads_each
            backend_local <- cfg$bp_backend
            n_files_local <- length(multi_files)
            n_records_local <- multi_records[[w]]
            total_mb_local <- multi_mb
            spec_local <- spec
            function(case_id, execution_rank) {
              bp <- make_bpparam(workers_local, backend = backend_local)
              on.exit(stop_bpparam(bp), add = TRUE)
              run_case(
                scenario = "multi",
                workload = "seqqual",
                method = "BamScale (compact seqqual budget)",
                method_family = "BamScale",
                threads_requested = threads_local,
                threads_effective = threads_local,
                bp_workers_requested = workers_local,
                bp_workers = workers_local,
                n_files = n_files_local,
                n_records = n_records_local,
                total_mb = total_mb_local,
                iterations = cfg$iterations,
                warmup = cfg$warmup,
                fun = function() {
                  BamScale::bam_read(
                    file = multi_files,
                    what = spec_local$what,
                    as = spec_local$as,
                    seqqual_mode = "compact",
                    threads = threads_local,
                    BPPARAM = bp,
                    auto_threads = FALSE
                  )
                },
                case_id = case_id,
                execution_rank = execution_rank,
                comparison_track = "optimized",
                seqqual_mode = "compact",
                auto_threads = FALSE,
                bp_backend = backend_local
              )
            }
          })
        )
      }

      if (cfg$include_rsamtools && w %in% c("step1", "seqqual") && requireNamespace("Rsamtools", quietly = TRUE)) {
        append_case(
          make_case_meta(
            scenario = "multi",
            workload = w,
            method = "Rsamtools::scanBam + BiocParallel",
            method_family = "Rsamtools",
            comparison_track = "fair",
            seqqual_mode = if (w == "seqqual") "compatible" else NA_character_,
            auto_threads = FALSE,
            threads_requested = 1L,
            threads_effective = 1L,
            bp_workers_requested = workers,
            bp_workers_effective = workers,
            n_files = length(multi_files),
            n_records = multi_records[[w]],
            total_mb = multi_mb,
            bp_backend = cfg$bp_backend
          ),
          local({
            w_local <- w
            spec_local <- spec
            workers_local <- workers
            backend_local <- cfg$bp_backend
            n_files_local <- length(multi_files)
            n_records_local <- multi_records[[w]]
            total_mb_local <- multi_mb
            function(case_id, execution_rank) {
              bp <- make_bpparam(workers_local, backend = backend_local)
              on.exit(stop_bpparam(bp), add = TRUE)
              run_case(
                scenario = "multi",
                workload = w_local,
                method = "Rsamtools::scanBam + BiocParallel",
                method_family = "Rsamtools",
                threads_requested = 1L,
                threads_effective = 1L,
                bp_workers_requested = workers_local,
                bp_workers = workers_local,
                n_files = n_files_local,
                n_records = n_records_local,
                total_mb = total_mb_local,
                iterations = cfg$iterations,
                warmup = cfg$warmup,
                fun = function() {
                  run_rsamtools_multi(multi_files, what = spec_local$what, bp = bp)
                },
                case_id = case_id,
                execution_rank = execution_rank,
                comparison_track = "fair",
                seqqual_mode = if (w_local == "seqqual") "compatible" else NA_character_,
                auto_threads = FALSE,
                bp_backend = backend_local
              )
            }
          })
        )
      }

      if (w == "galignments" && requireNamespace("GenomicAlignments", quietly = TRUE) && requireNamespace("Rsamtools", quietly = TRUE)) {
        append_case(
          make_case_meta(
            scenario = "multi",
            workload = w,
            method = "GenomicAlignments::readGAlignments + BiocParallel",
            method_family = "GenomicAlignments",
            comparison_track = "fair",
            seqqual_mode = NA_character_,
            auto_threads = FALSE,
            threads_requested = 1L,
            threads_effective = 1L,
            bp_workers_requested = workers,
            bp_workers_effective = workers,
            n_files = length(multi_files),
            n_records = multi_records[[w]],
            total_mb = multi_mb,
            bp_backend = cfg$bp_backend
          ),
          local({
            workers_local <- workers
            backend_local <- cfg$bp_backend
            n_files_local <- length(multi_files)
            n_records_local <- multi_records[[w]]
            total_mb_local <- multi_mb
            spec_local <- spec
            function(case_id, execution_rank) {
              bp <- make_bpparam(workers_local, backend = backend_local)
              on.exit(stop_bpparam(bp), add = TRUE)
              run_case(
                scenario = "multi",
                workload = "galignments",
                method = "GenomicAlignments::readGAlignments + BiocParallel",
                method_family = "GenomicAlignments",
                threads_requested = 1L,
                threads_effective = 1L,
                bp_workers_requested = workers_local,
                bp_workers = workers_local,
                n_files = n_files_local,
                n_records = n_records_local,
                total_mb = total_mb_local,
                iterations = cfg$iterations,
                warmup = cfg$warmup,
                fun = function() {
                  run_galignments_multi(multi_files, what = spec_local$what, bp = bp)
                },
                case_id = case_id,
                execution_rank = execution_rank,
                comparison_track = "fair",
                seqqual_mode = NA_character_,
                auto_threads = FALSE,
                bp_backend = backend_local
              )
            }
          })
        )
      }
    }
  }
}

if (length(case_specs) == 0L) {
  stop("No benchmark cases configured", call. = FALSE)
}

case_plan_rows <- list()
global_execution_rank <- 0L
for (outer_run in seq_len(cfg$outer_runs)) {
  case_order <- seq_along(case_specs)
  if (isTRUE(cfg$shuffle_cases) && length(case_specs) > 1L) {
    set.seed(cfg$shuffle_seed + outer_run - 1L)
    case_order <- sample(case_order, length(case_order), replace = FALSE)
  }

  for (exec_rank in seq_along(case_order)) {
    case_id <- case_order[[exec_rank]]
    meta <- case_specs[[case_id]]$meta
    global_execution_rank <- global_execution_rank + 1L
    meta$case_id <- as.integer(case_id)
    meta$outer_run <- as.integer(outer_run)
    meta$execution_rank <- as.integer(exec_rank)
    meta$global_execution_rank <- as.integer(global_execution_rank)
    case_plan_rows[[length(case_plan_rows) + 1L]] <- meta
  }
}

case_plan_df <- do.call(rbind, case_plan_rows)
case_plan_df <- case_plan_df[, c(
  "case_id", "outer_run", "execution_rank", "global_execution_rank", "build_order", "scenario", "workload",
  "method", "method_family", "comparison_track", "seqqual_mode",
  "auto_threads", "threads_requested", "threads_effective",
  "bp_workers_requested", "bp_workers_effective",
  "total_threads_requested", "total_threads_effective",
  "n_files", "n_records", "total_mb", "bp_backend"
)]
write.csv(case_plan_df, file = file.path(out_dir, "case_plan.csv"), row.names = FALSE)
write_resolved_plans(case_plan_df, out_dir)

old_outer_run_option <- getOption("bamscale.benchmark.outer_run")
on.exit(options(bamscale.benchmark.outer_run = old_outer_run_option), add = TRUE)

for (i in seq_len(nrow(case_plan_df))) {
  options(bamscale.benchmark.outer_run = case_plan_df$outer_run[[i]])
  out <- case_specs[[case_plan_df$case_id[[i]]]]$runner(
    case_id = case_plan_df$case_id[[i]],
    execution_rank = case_plan_df$execution_rank[[i]]
  )
  summary_rows[[length(summary_rows) + 1L]] <- out$summary
  iter_rows[[length(iter_rows) + 1L]] <- out$iterations
  write_checkpoint(summary_rows, iter_rows, out_dir)
}

final <- write_checkpoint(summary_rows, iter_rows, out_dir)
summary_df <- final$summary
iter_df <- final$iterations

write_config(file.path(out_dir, "config.txt"), cfg = cfg, opts = opts)
write_session_info(file.path(out_dir, "sessionInfo.txt"))
write_host_info(file.path(out_dir, "host_info.csv"), cfg = cfg)
write_env_vars(file.path(out_dir, "env_vars.csv"))
publication_tables <- write_publication_tables(summary_df, iter_df, out_dir)
plot_if_possible(publication_tables$summary_qmd, out_dir)
baseline_check <- write_baseline_comparison(publication_tables$outer_run_summary, cfg = cfg, out_dir = out_dir)
write_artifact_manifest(out_dir)

ok_n <- sum(summary_df$status == "ok")
err_n <- sum(summary_df$status != "ok")

message("Benchmark complete.")
message("Successful cases: ", ok_n)
message("Failed cases: ", err_n)
message("Results: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))

if (err_n > 0L) {
  message("Some cases failed. Review summary.csv (status/error_message).")
}

if (nzchar(cfg$baseline_results_dir)) {
  matched_n <- sum(baseline_check$comparison$baseline_match_found, na.rm = TRUE)
  regress_n <- nrow(baseline_check$regressions)
  message("Baseline source: ", baseline_check$baseline_source)
  message("Baseline-matched cases: ", matched_n)
  message("Regression-gated cases above threshold: ", regress_n)
  if (regress_n > 0L) {
    message("Review baseline_regressions.csv for matched cases that regressed beyond the configured threshold.")
  }
  if (isTRUE(cfg$regression_stop_on_fail) && regress_n > 0L) {
    stop("Regression gate failed; review baseline_regressions.csv", call. = FALSE)
  }
}
