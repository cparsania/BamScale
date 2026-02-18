#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(BamScale)
})

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
        "This usually means BamScale was not installed/loaded with compiled code.\n",
        "Fix on server:\n",
        "1) Restart R session\n",
        "2) Reinstall from source: install.packages('/path/to/BamScale', repos = NULL, type = 'source')\n",
        "3) Validate: getNativeSymbolInfo('_BamScale_read_bam_cpp', PACKAGE='BamScale')"
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

fetch_chipseqdbdata_bams <- function(n_target = 12L) {
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
      " available. Provide more files or set --allow-repeat-files=true.",
      call. = FALSE
    )
  }

  rep(paths, length.out = n_files)
}

make_bpparam <- function(workers) {
  workers <- as.integer(workers)
  if (workers <= 1L) return(NULL)

  require_or_stop("BiocParallel")

  if (.Platform$OS.type == "unix") {
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

resolve_effective_threads <- function(requested_threads, BPPARAM, auto_threads, n_files) {
  fn <- getFromNamespace(".bamscale_resolve_threads", "BamScale")
  fn(
    threads = as.integer(requested_threads),
    BPPARAM = BPPARAM,
    auto_threads = as.logical(auto_threads),
    n_files = as.integer(n_files)
  )
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
  if (isTRUE(warmup)) {
    invisible(fun())
  }

  elapsed <- numeric(iterations)
  for (i in seq_len(iterations)) {
    gc(verbose = FALSE)
    elapsed[[i]] <- system.time(fun())[["elapsed"]]
  }

  elapsed
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
    case_id
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

  iter_df <- data.frame(
    case_id = integer(),
    iteration = integer(),
    elapsed_s = numeric(),
    status = character(),
    error_message = character(),
    stringsAsFactors = FALSE
  )

  error_message <- NA_character_
  status <- "ok"

  elapsed <- tryCatch(
    bench_times(fun, iterations = iterations, warmup = warmup),
    error = function(e) {
      status <<- "error"
      error_message <<- conditionMessage(e)
      rep(NA_real_, iterations)
    }
  )

  for (i in seq_len(iterations)) {
    iter_df <- rbind(
      iter_df,
      data.frame(
        case_id = case_id,
        iteration = i,
        elapsed_s = elapsed[[i]],
        status = if (is.finite(elapsed[[i]])) "ok" else "error",
        error_message = if (is.finite(elapsed[[i]])) NA_character_ else error_message,
        stringsAsFactors = FALSE
      )
    )
  }

  valid <- elapsed[is.finite(elapsed)]
  if (length(valid) == 0L) {
    status <- "error"
  }

  min_s <- if (length(valid) > 0L) min(valid) else NA_real_
  median_s <- if (length(valid) > 0L) median(valid) else NA_real_
  max_s <- if (length(valid) > 0L) max(valid) else NA_real_

  records_per_s <- if (!is.na(n_records) && !is.na(median_s) && median_s > 0) n_records / median_s else NA_real_
  mb_per_s <- if (!is.na(total_mb) && !is.na(median_s) && median_s > 0) total_mb / median_s else NA_real_

  summary_df <- data.frame(
    case_id = case_id,
    scenario = scenario,
    workload = workload,
    method = method,
    method_family = method_family,
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

  p2_df <- ok[ok$scenario == "multi" & ok$method_family == "BamScale", , drop = FALSE]
  if (nrow(p2_df) > 0L) {
    p <- ggplot2::ggplot(
      p2_df,
      ggplot2::aes(x = bp_workers, y = median_s, color = method)
    ) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2) +
      ggplot2::facet_wrap(~workload, scales = "free_y") +
      ggplot2::scale_x_continuous(breaks = sort(unique(p2_df$bp_workers))) +
      ggplot2::labs(
        title = "Multi-file scaling (BamScale)",
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

cfg <- list(
  outdir = if (!is.null(opts$outdir)) opts$outdir else "benchmark_results",
  n_files = as_int_scalar(opts[["n-files"]], 12L),
  iterations = as_int_scalar(opts$iterations, 3L),
  warmup = as_bool(opts$warmup, TRUE),
  threads = as_int_vec(opts$threads, c(1L, 2L, 4L, 8L, 12L, 18L, 24L, 36L, 48L, 72L)),
  workers = as_int_vec(opts$workers, c(1L, 2L, 4L, 6L, 8L, 12L, 18L, 24L, 36L, 72L)),
  max_threads = as_int_scalar(opts[["max-threads"]], 72L),
  include_seqqual = as_bool(opts[["include-seqqual"]], TRUE),
  include_galignments = as_bool(opts[["include-galignments"]], TRUE),
  include_multi = as_bool(opts[["include-multi"]], TRUE),
  include_rsamtools = as_bool(opts[["include-rsamtools"]], TRUE),
  ensure_index = as_bool(opts[["ensure-index"]], TRUE),
  allow_repeat_files = as_bool(opts[["allow-repeat-files"]], FALSE),
  recursive = as_bool(opts$recursive, TRUE),
  multi_workloads = as_chr_vec(opts[["multi-workloads"]], c("step1")),
  compute_records = as_bool(opts[["compute-records"]], TRUE)
)

cfg$threads <- cfg$threads[cfg$threads >= 1L & cfg$threads <= cfg$max_threads]
if (length(cfg$threads) == 0L) {
  stop("No valid --threads values after filtering by --max-threads", call. = FALSE)
}

cfg$workers <- cfg$workers[cfg$workers >= 1L & cfg$workers <= cfg$max_threads]
if (length(cfg$workers) == 0L) {
  cfg$workers <- 1L
}

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
sizes <- sizes[ord]

single_file <- bam_paths[[1L]]
multi_files <- select_n_files(
  bam_paths,
  n_files = cfg$n_files,
  allow_repeat = cfg$allow_repeat_files
)

effective_workers_cap <- min(length(multi_files), cfg$max_threads)
cfg$workers <- cfg$workers[cfg$workers <= effective_workers_cap]
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
    as = "data.frame"
  ),
  seqqual = list(
    what = c("qname", "seq", "qual"),
    as = "data.frame"
  ),
  galignments = list(
    what = c("qname", "rname", "pos", "cigar", "strand", "flag"),
    as = "GAlignments"
  )
)

active_workloads <- c("step1")
if (isTRUE(cfg$include_seqqual)) active_workloads <- c(active_workloads, "seqqual")
if (isTRUE(cfg$include_galignments)) active_workloads <- c(active_workloads, "galignments")

cfg$multi_workloads <- unique(intersect(cfg$multi_workloads, names(workloads)))
if (length(cfg$multi_workloads) == 0L) cfg$multi_workloads <- "step1"

message("Single-file benchmark target: ", single_file)
message("Multi-file count: ", length(multi_files))
message("Threads grid: ", paste(cfg$threads, collapse = ", "))
message("Workers grid: ", paste(cfg$workers, collapse = ", "))

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
      threads = 1L
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
          threads = t,
          BPPARAM = NULL,
          auto_threads = FALSE
        )
      },
      case_id = case_id
    )

    summary_rows[[length(summary_rows) + 1L]] <- out$summary
    iter_rows[[length(iter_rows) + 1L]] <- out$iterations
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
  }
}

if (isTRUE(cfg$include_multi)) {
  for (w in cfg$multi_workloads) {
    spec <- workloads[[w]]

    for (workers in cfg$workers) {
      bp <- tryCatch(make_bpparam(workers), error = function(e) e)
      if (inherits(bp, "error")) {
        message("Skipping workers=", workers, " because BPPARAM init failed: ", conditionMessage(bp))
        next
      }

      on.exit(stop_bpparam(bp), add = TRUE)

      eff_auto <- resolve_effective_threads(
        requested_threads = cfg$max_threads,
        BPPARAM = bp,
        auto_threads = TRUE,
        n_files = length(multi_files)
      )

      case_id <- case_id + 1L
      out <- run_case(
        scenario = "multi",
        workload = w,
        method = "BamScale (auto_threads)",
        method_family = "BamScale",
        threads_requested = cfg$max_threads,
        threads_effective = eff_auto,
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
            threads = cfg$max_threads,
            BPPARAM = bp,
            auto_threads = TRUE
          )
        },
        case_id = case_id
      )

      summary_rows[[length(summary_rows) + 1L]] <- out$summary
      iter_rows[[length(iter_rows) + 1L]] <- out$iterations

      balanced_threads <- max(1L, floor(cfg$max_threads / workers))
      case_id <- case_id + 1L
      out <- run_case(
        scenario = "multi",
        workload = w,
        method = "BamScale (manual balanced)",
        method_family = "BamScale",
        threads_requested = balanced_threads,
        threads_effective = balanced_threads,
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
            threads = balanced_threads,
            BPPARAM = bp,
            auto_threads = FALSE
          )
        },
        case_id = case_id
      )

      summary_rows[[length(summary_rows) + 1L]] <- out$summary
      iter_rows[[length(iter_rows) + 1L]] <- out$iterations

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
      }

      stop_bpparam(bp)
    }
  }
}

summary_df <- do.call(rbind, summary_rows)
iter_df <- do.call(rbind, iter_rows)

summary_df <- summary_df[order(summary_df$scenario, summary_df$workload, summary_df$method_family, summary_df$threads_effective, summary_df$bp_workers), ]
iter_df <- iter_df[order(iter_df$case_id, iter_df$iteration), ]

write.csv(summary_df, file = file.path(out_dir, "summary.csv"), row.names = FALSE)
write.csv(iter_df, file = file.path(out_dir, "iterations.csv"), row.names = FALSE)
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
