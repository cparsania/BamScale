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

run_case <- function(
    scenario,
    workload,
    method,
    method_family,
    threads_requested,
    threads_effective,
    bp_workers,
    n_files,
    n_records,
    total_mb,
    iterations,
    warmup,
    fun,
    case_id,
    comparison_track = "fair",
    seqqual_mode = NA_character_
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
    iteration = seq_len(iterations),
    elapsed_s = elapsed,
    status = ifelse(is.finite(elapsed), "ok", "error"),
    error_message = ifelse(is.finite(elapsed), NA_character_, errors),
    stringsAsFactors = FALSE
  )

  valid <- elapsed[is.finite(elapsed)]
  status <- if (length(valid) > 0L) "ok" else "error"
  error_message <- if (status == "error") paste(unique(na.omit(errors)), collapse = " | ") else NA_character_

  min_s <- if (length(valid) > 0L) min(valid) else NA_real_
  median_s <- if (length(valid) > 0L) median(valid) else NA_real_
  max_s <- if (length(valid) > 0L) max(valid) else NA_real_

  records_per_s <- if (!is.na(n_records) && !is.na(median_s) && median_s > 0) n_records / median_s else NA_real_
  mb_per_s <- if (!is.na(total_mb) && !is.na(median_s) && median_s > 0) total_mb / median_s else NA_real_

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
    scenario = scenario,
    workload = workload,
    method = method,
    method_family = method_family,
    comparison_track = comparison_track_scalar,
    seqqual_mode = seqqual_mode_scalar,
    threads_requested = as.integer(threads_requested),
    threads_effective = as.integer(threads_effective),
    bp_workers = as.integer(bp_workers),
    total_threads = as.integer(max(1L, bp_workers) * max(1L, threads_effective)),
    n_files = as.integer(n_files),
    n_records = as.numeric(n_records),
    total_mb = as.numeric(total_mb),
    iterations_planned = as.integer(iterations),
    iterations_completed = as.integer(length(valid)),
    status = status,
    error_message = error_message,
    min_s = min_s,
    median_s = median_s,
    max_s = max_s,
    records_per_s = records_per_s,
    mb_per_s = mb_per_s,
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

write_session_info <- function(path_out) {
  con <- file(path_out, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(capture.output(sessionInfo()), con = con)
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
    summary_df$scenario,
    summary_df$workload,
    summary_df$comparison_track,
    summary_df$method_family,
    summary_df$seqqual_mode,
    summary_df$threads_effective,
    summary_df$bp_workers
  ), ]
  iter_df <- iter_df[order(iter_df$case_id, iter_df$iteration), ]

  write.csv(summary_df, file = file.path(out_dir, "summary.csv"), row.names = FALSE)
  write.csv(iter_df, file = file.path(out_dir, "iterations.csv"), row.names = FALSE)

  invisible(list(summary = summary_df, iterations = iter_df))
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
  default_include_rsamtools <- TRUE
  default_multi_workloads <- c("step1")
} else if (profile == "balanced") {
  default_n_files <- 12L
  default_threads <- c(1L, 2L, 4L, 8L, 12L, 16L, 24L, 32L, 48L)
  default_workers <- c(1L, 2L, 4L, 8L, 12L, 16L)
  default_include_seqqual <- TRUE
  default_seqqual_compact <- TRUE
  default_include_galignments <- TRUE
  default_include_rsamtools <- TRUE
  default_multi_workloads <- c("step1", "seqqual")
} else {
  default_n_files <- 16L
  default_threads <- c(1L, 2L, 4L, 8L, 12L, 16L, 24L, 32L, 48L, 64L, 96L)
  default_workers <- c(1L, 2L, 4L, 8L, 12L, 16L, 24L, 32L, 48L, 64L, 96L)
  default_include_seqqual <- TRUE
  default_seqqual_compact <- TRUE
  default_include_galignments <- TRUE
  default_include_rsamtools <- TRUE
  default_multi_workloads <- c("step1", "seqqual", "galignments")
}

default_budget_threads <- as.integer(min(48L, detected_cores))
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
  iterations = as_int_scalar(opts$iterations, 3L),
  warmup = as_bool(opts$warmup, TRUE),
  max_threads = as.integer(max_threads),
  threads = as_int_vec(opts$threads, default_threads),
  workers = as_int_vec(opts$workers, default_workers),
  bp_backend = if (!is.null(opts[["bp-backend"]])) tolower(opts[["bp-backend"]]) else "snow",
  include_seqqual = as_bool(opts[["include-seqqual"]], default_include_seqqual),
  seqqual_compact = as_bool(opts[["seqqual-compact"]], default_seqqual_compact),
  include_galignments = as_bool(opts[["include-galignments"]], default_include_galignments),
  include_multi = as_bool(opts[["include-multi"]], TRUE),
  include_rsamtools = as_bool(opts[["include-rsamtools"]], default_include_rsamtools),
  ensure_index = as_bool(opts[["ensure-index"]], TRUE),
  allow_repeat_files = as_bool(opts[["allow-repeat-files"]], FALSE),
  recursive = as_bool(opts$recursive, TRUE),
  multi_workloads = as_chr_vec(opts[["multi-workloads"]], default_multi_workloads),
  compute_records = as_bool(opts[["compute-records"]], TRUE)
)
if (!cfg$bp_backend %in% c("snow", "multicore")) {
  stop("--bp-backend must be one of: snow,multicore", call. = FALSE)
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
  bam_paths <- normalizePath(manual_files, winslash = "/", mustWork = TRUE)
} else {
  message("No BAM files provided, resolving BAMs from chipseqDBData...")
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

files_df <- data.frame(
  file = multi_files,
  size_mb = round(file.size(multi_files) / 1024^2, 3),
  stringsAsFactors = FALSE
)
write.csv(files_df, file = file.path(out_dir, "files.csv"), row.names = FALSE)

if (cfg$include_rsamtools && cfg$ensure_index) {
  message("Ensuring BAM index files for comparator methods...")
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

summary_rows <- list()
iter_rows <- list()
case_id <- 0L

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

for (w in active_workloads) {
  spec <- workloads[[w]]

  for (t in cfg$threads) {
    case_id <- case_id + 1L
    out <- run_case(
      scenario = "single",
      workload = w,
      method = "BamScale",
      method_family = "BamScale",
      threads_requested = t,
      threads_effective = t,
      bp_workers = 1L,
      n_files = 1L,
      n_records = single_records[[w]],
      total_mb = single_mb,
      iterations = cfg$iterations,
      warmup = cfg$warmup,
      fun = function() {
        BamScale::bam_read(
          file = single_file,
          what = spec$what,
          as = spec$as,
          seqqual_mode = spec$seqqual_mode,
          threads = t,
          BPPARAM = NULL,
          auto_threads = FALSE
        )
      },
      case_id = case_id,
      comparison_track = "fair",
      seqqual_mode = if (w == "seqqual") "compatible" else NA_character_
    )

    summary_rows[[length(summary_rows) + 1L]] <- out$summary
    iter_rows[[length(iter_rows) + 1L]] <- out$iterations
    write_checkpoint(summary_rows, iter_rows, out_dir)

    if (w == "seqqual" && isTRUE(cfg$seqqual_compact)) {
      case_id <- case_id + 1L
      out <- run_case(
        scenario = "single",
        workload = w,
        method = "BamScale (compact seqqual)",
        method_family = "BamScale",
        threads_requested = t,
        threads_effective = t,
        bp_workers = 1L,
        n_files = 1L,
        n_records = single_records[[w]],
        total_mb = single_mb,
        iterations = cfg$iterations,
        warmup = cfg$warmup,
        fun = function() {
          BamScale::bam_read(
            file = single_file,
            what = spec$what,
            as = spec$as,
            seqqual_mode = "compact",
            threads = t,
            BPPARAM = NULL,
            auto_threads = FALSE
          )
        },
        case_id = case_id,
        comparison_track = "optimized",
        seqqual_mode = "compact"
      )

      summary_rows[[length(summary_rows) + 1L]] <- out$summary
      iter_rows[[length(iter_rows) + 1L]] <- out$iterations
      write_checkpoint(summary_rows, iter_rows, out_dir)
    }
  }

  if (cfg$include_rsamtools && w %in% c("step1", "seqqual") && requireNamespace("Rsamtools", quietly = TRUE)) {
    case_id <- case_id + 1L
    out <- run_case(
      scenario = "single",
      workload = w,
      method = "Rsamtools::scanBam",
      method_family = "Rsamtools",
      threads_requested = 1L,
      threads_effective = 1L,
      bp_workers = 1L,
      n_files = 1L,
      n_records = single_records[[w]],
      total_mb = single_mb,
      iterations = cfg$iterations,
      warmup = cfg$warmup,
      fun = function() {
        param <- Rsamtools::ScanBamParam(what = as.character(spec$what))
        Rsamtools::scanBam(single_file, param = param)
      },
      case_id = case_id
    )

    summary_rows[[length(summary_rows) + 1L]] <- out$summary
    iter_rows[[length(iter_rows) + 1L]] <- out$iterations
    write_checkpoint(summary_rows, iter_rows, out_dir)
  }

  if (w == "galignments" && requireNamespace("GenomicAlignments", quietly = TRUE) && requireNamespace("Rsamtools", quietly = TRUE)) {
    case_id <- case_id + 1L
    out <- run_case(
      scenario = "single",
      workload = w,
      method = "GenomicAlignments::readGAlignments",
      method_family = "GenomicAlignments",
      threads_requested = 1L,
      threads_effective = 1L,
      bp_workers = 1L,
      n_files = 1L,
      n_records = single_records[[w]],
      total_mb = single_mb,
      iterations = cfg$iterations,
      warmup = cfg$warmup,
      fun = function() {
        param <- Rsamtools::ScanBamParam(what = as.character(spec$what))
        GenomicAlignments::readGAlignments(single_file, param = param)
      },
      case_id = case_id
    )

    summary_rows[[length(summary_rows) + 1L]] <- out$summary
    iter_rows[[length(iter_rows) + 1L]] <- out$iterations
    write_checkpoint(summary_rows, iter_rows, out_dir)
  }
}

if (isTRUE(cfg$include_multi)) {
  for (w in cfg$multi_workloads) {
    spec <- workloads[[w]]

    for (workers in cfg$workers) {
      bp <- tryCatch(
        make_bpparam(workers, backend = cfg$bp_backend),
        error = function(e) e
      )

      if (inherits(bp, "error")) {
        message("Skipping workers=", workers, " because BPPARAM init failed: ", conditionMessage(bp))
        next
      }

      threads_each <- max(1L, floor(cfg$max_threads / workers))

      case_id <- case_id + 1L
      out <- run_case(
        scenario = "multi",
        workload = w,
        method = "BamScale (balanced budget)",
        method_family = "BamScale",
        threads_requested = threads_each,
        threads_effective = threads_each,
        bp_workers = workers,
        n_files = length(multi_files),
        n_records = multi_records[[w]],
        total_mb = multi_mb,
        iterations = cfg$iterations,
        warmup = cfg$warmup,
        fun = function() {
          BamScale::bam_read(
            file = multi_files,
            what = spec$what,
            as = spec$as,
            seqqual_mode = spec$seqqual_mode,
            threads = threads_each,
            BPPARAM = bp,
            auto_threads = FALSE
          )
        },
        case_id = case_id,
        comparison_track = "fair",
        seqqual_mode = if (w == "seqqual") "compatible" else NA_character_
      )

      summary_rows[[length(summary_rows) + 1L]] <- out$summary
      iter_rows[[length(iter_rows) + 1L]] <- out$iterations
      write_checkpoint(summary_rows, iter_rows, out_dir)

      if (w == "seqqual" && isTRUE(cfg$seqqual_compact)) {
        case_id <- case_id + 1L
        out <- run_case(
          scenario = "multi",
          workload = w,
          method = "BamScale (compact seqqual budget)",
          method_family = "BamScale",
          threads_requested = threads_each,
          threads_effective = threads_each,
          bp_workers = workers,
          n_files = length(multi_files),
          n_records = multi_records[[w]],
          total_mb = multi_mb,
          iterations = cfg$iterations,
          warmup = cfg$warmup,
          fun = function() {
            BamScale::bam_read(
              file = multi_files,
              what = spec$what,
              as = spec$as,
              seqqual_mode = "compact",
              threads = threads_each,
              BPPARAM = bp,
              auto_threads = FALSE
            )
          },
          case_id = case_id,
          comparison_track = "optimized",
          seqqual_mode = "compact"
        )

        summary_rows[[length(summary_rows) + 1L]] <- out$summary
        iter_rows[[length(iter_rows) + 1L]] <- out$iterations
        write_checkpoint(summary_rows, iter_rows, out_dir)
      }

      if (cfg$include_rsamtools && w %in% c("step1", "seqqual") && requireNamespace("Rsamtools", quietly = TRUE)) {
        case_id <- case_id + 1L
        out <- run_case(
          scenario = "multi",
          workload = w,
          method = "Rsamtools::scanBam + BiocParallel",
          method_family = "Rsamtools",
          threads_requested = 1L,
          threads_effective = 1L,
          bp_workers = workers,
          n_files = length(multi_files),
          n_records = multi_records[[w]],
          total_mb = multi_mb,
          iterations = cfg$iterations,
          warmup = cfg$warmup,
          fun = function() {
            run_rsamtools_multi(multi_files, what = spec$what, bp = bp)
          },
          case_id = case_id
        )

        summary_rows[[length(summary_rows) + 1L]] <- out$summary
        iter_rows[[length(iter_rows) + 1L]] <- out$iterations
        write_checkpoint(summary_rows, iter_rows, out_dir)
      }

      if (w == "galignments" && requireNamespace("GenomicAlignments", quietly = TRUE) && requireNamespace("Rsamtools", quietly = TRUE)) {
        case_id <- case_id + 1L
        out <- run_case(
          scenario = "multi",
          workload = w,
          method = "GenomicAlignments::readGAlignments + BiocParallel",
          method_family = "GenomicAlignments",
          threads_requested = 1L,
          threads_effective = 1L,
          bp_workers = workers,
          n_files = length(multi_files),
          n_records = multi_records[[w]],
          total_mb = multi_mb,
          iterations = cfg$iterations,
          warmup = cfg$warmup,
          fun = function() {
            run_galignments_multi(multi_files, what = spec$what, bp = bp)
          },
          case_id = case_id
        )

        summary_rows[[length(summary_rows) + 1L]] <- out$summary
        iter_rows[[length(iter_rows) + 1L]] <- out$iterations
        write_checkpoint(summary_rows, iter_rows, out_dir)
      }

      stop_bpparam(bp)
    }
  }
}

final <- write_checkpoint(summary_rows, iter_rows, out_dir)
summary_df <- final$summary

write_config(file.path(out_dir, "config.txt"), cfg = cfg, opts = opts)
write_session_info(file.path(out_dir, "sessionInfo.txt"))
plot_if_possible(summary_df, out_dir)

ok_n <- sum(summary_df$status == "ok")
err_n <- sum(summary_df$status != "ok")

message("Benchmark complete.")
message("Successful cases: ", ok_n)
message("Failed cases: ", err_n)
message("Results: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))

if (err_n > 0L) {
  message("Some cases failed. Review summary.csv (status/error_message).")
}
