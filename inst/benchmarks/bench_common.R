# bench_common.R -- shared helpers for the BamScale workflow benchmarks.
#
# This module is sourced by run_workflow_benchmark.R. It intentionally defines
# functions only (no side effects on source) so it can be reused by other
# drivers. It mirrors the timing/output conventions of run_server_benchmark.R
# (same summary column names where they overlap) but adds a phase-aware timing
# engine that separates the read / compute / write stages of an end-to-end
# workflow -- the decomposition the Amdahl argument in the report depends on.

# ---------------------------------------------------------------------------
# Argument parsing (mirrors run_server_benchmark.R's --key=value convention)
# ---------------------------------------------------------------------------

.wf_parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
    out <- list()
    for (a in argv) {
        if (!grepl("^--", a)) next
        a <- sub("^--", "", a)
        if (grepl("=", a, fixed = TRUE)) {
            k <- sub("=.*$", "", a)
            v <- sub("^[^=]*=", "", a)
        } else {
            k <- a
            v <- "true"
        }
        out[[gsub("-", "_", k)]] <- v
    }
    out
}

.wf_as_bool <- function(x, default = FALSE) {
    if (is.null(x)) return(default)
    tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

.wf_as_int <- function(x, default = NA_integer_) {
    if (is.null(x)) return(default)
    v <- suppressWarnings(as.integer(x))
    if (length(v) != 1L || is.na(v)) default else v
}

.wf_int_vec <- function(x, default = integer()) {
    if (is.null(x)) return(default)
    v <- suppressWarnings(as.integer(strsplit(as.character(x), ",", fixed = TRUE)[[1L]]))
    v <- v[!is.na(v)]
    if (!length(v)) default else v
}

.wf_split_csv <- function(x, default = character()) {
    if (is.null(x)) return(default)
    v <- trimws(strsplit(as.character(x), ",", fixed = TRUE)[[1L]])
    v <- v[nzchar(v)]
    if (!length(v)) default else v
}

.wf_log <- function(...) {
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(..., collapse = "")))
    utils::flush.console()
}

# ---------------------------------------------------------------------------
# Cores / parallel backend
# ---------------------------------------------------------------------------

.wf_detect_cores <- function() {
    # Honor Linux cpuset limits (containers/schedulers) before falling back to
    # detectCores(), matching BamScale's own core detection.
    if (.Platform$OS.type == "unix" && file.exists("/proc/self/status")) {
        status <- tryCatch(readLines("/proc/self/status", warn = FALSE),
                           error = function(e) character())
        hit <- grep("^Cpus_allowed_list", status, value = TRUE)
        if (length(hit)) {
            allowed <- sub("^Cpus_allowed_list\\s*:\\s*", "", hit[[1L]])
            parts <- strsplit(gsub("\\s+", "", allowed), ",", fixed = TRUE)[[1L]]
            total <- 0L
            for (p in parts) {
                if (!nzchar(p)) next
                if (grepl("-", p, fixed = TRUE)) {
                    b <- as.integer(strsplit(p, "-", fixed = TRUE)[[1L]])
                    if (length(b) == 2L && !anyNA(b) && b[2] >= b[1]) total <- total + (b[2] - b[1] + 1L)
                } else if (!is.na(suppressWarnings(as.integer(p)))) {
                    total <- total + 1L
                }
            }
            if (total > 0L) return(as.integer(total))
        }
    }
    phys <- suppressWarnings(as.integer(parallel::detectCores(logical = FALSE)))
    if (length(phys) == 1L && !is.na(phys) && phys >= 1L) return(phys)
    log <- suppressWarnings(as.integer(parallel::detectCores(logical = TRUE)))
    if (length(log) == 1L && !is.na(log) && log >= 1L) return(log)
    1L
}

.wf_make_bpparam <- function(backend = "snow", workers = 1L) {
    workers <- max(1L, as.integer(workers))
    if (workers <= 1L) return(NULL)
    if (identical(backend, "multicore")) {
        BiocParallel::MulticoreParam(workers = workers, progressbar = FALSE)
    } else {
        BiocParallel::SnowParam(workers = workers, type = "SOCK", progressbar = FALSE)
    }
}

# ---------------------------------------------------------------------------
# Phase-aware timing engine
# ---------------------------------------------------------------------------

# Time a single expression, returning its value plus wall / user / sys seconds.
# user/sys include waited-for children (fork/multicore); PSOCK worker CPU is not
# captured at the master -- documented in the report's methods.
.wf_time_phase <- function(expr) {
    t0 <- proc.time()
    value <- force(expr)
    d <- proc.time() - t0
    list(
        value = value,
        elapsed = unname(d[["elapsed"]]),
        user = unname(d[["user.self"]] + d[["user.child"]]),
        sys = unname(d[["sys.self"]] + d[["sys.child"]])
    )
}

# Run `run_once` (a 0-arg closure) once as warm-up (optional) then `iterations`
# timed times, gc() before each. `run_once` must return a named list whose
# elements are numeric c(elapsed=, user=, sys=) vectors, one per phase name; it
# may also return a `total` element to override the summed-phase total (used by
# the multi-file case where phases overlap across workers).
.wf_bench_phases <- function(run_once, iterations = 3L, warmup = TRUE,
                             phase_names = c("read", "compute", "write")) {
    errors <- character()

    if (isTRUE(warmup)) {
        ok <- tryCatch({ run_once(); TRUE },
                       error = function(e) { errors <<- c(errors, paste0("warmup: ", conditionMessage(e))); FALSE })
        if (!isTRUE(ok)) return(list(status = "error", error_message = errors[[1]], iterations_completed = 0L))
    }

    rows <- vector("list", iterations)
    for (i in seq_len(iterations)) {
        gc(FALSE)
        rows[[i]] <- tryCatch(run_once(),
                              error = function(e) { errors <<- c(errors, sprintf("iter%d: %s", i, conditionMessage(e))); NULL })
    }
    ok <- Filter(Negate(is.null), rows)
    if (!length(ok)) {
        return(list(status = "error",
                    error_message = if (length(errors)) errors[[1]] else "all iterations failed",
                    iterations_completed = 0L))
    }

    phase_field <- function(r, ph, field) {
        v <- r[[ph]]
        if (is.null(v) || is.na(v[[field]])) NA_real_ else as.numeric(v[[field]])
    }

    # Dispersion (median/min/max/sd) so small ratios can be shown to sit above
    # noise; multithreaded runs have asymmetric run-to-run variance vs the
    # single-threaded baseline, so the spread must be reported, not just a median.
    disp <- function(x) {
        x <- x[is.finite(x)]
        if (!length(x)) return(c(med = NA_real_, min = NA_real_, max = NA_real_, sd = NA_real_))
        c(med = stats::median(x), min = min(x), max = max(x),
          sd = if (length(x) > 1L) stats::sd(x) else 0)
    }

    agg <- list()
    for (ph in phase_names) {
        d <- disp(vapply(ok, phase_field, numeric(1), ph, "elapsed"))
        agg[[paste0(ph, "_s")]]     <- unname(d[["med"]])
        agg[[paste0(ph, "_min_s")]] <- unname(d[["min"]])
        agg[[paste0(ph, "_max_s")]] <- unname(d[["max"]])
        agg[[paste0(ph, "_sd_s")]]  <- unname(d[["sd"]])
        agg[[paste0(ph, "_user_s")]] <- stats::median(vapply(ok, phase_field, numeric(1), ph, "user"), na.rm = TRUE)
        agg[[paste0(ph, "_sys_s")]]  <- stats::median(vapply(ok, phase_field, numeric(1), ph, "sys"),  na.rm = TRUE)
    }

    if (!is.null(ok[[1]][["total"]])) {
        tot_el <- vapply(ok, function(r) as.numeric(r[["total"]][["elapsed"]]), numeric(1))
        tot_us <- vapply(ok, function(r) as.numeric(r[["total"]][["user"]]),    numeric(1))
        tot_sy <- vapply(ok, function(r) as.numeric(r[["total"]][["sys"]]),     numeric(1))
    } else {
        sum_field <- function(r, field) sum(vapply(phase_names, function(ph) {
            v <- phase_field(r, ph, field); if (is.na(v)) 0 else v
        }, numeric(1)))
        tot_el <- vapply(ok, sum_field, numeric(1), "elapsed")
        tot_us <- vapply(ok, sum_field, numeric(1), "user")
        tot_sy <- vapply(ok, sum_field, numeric(1), "sys")
    }

    dt <- disp(tot_el)
    agg$total_s      <- unname(dt[["med"]])
    agg$total_min_s  <- unname(dt[["min"]])
    agg$total_max_s  <- unname(dt[["max"]])
    agg$total_sd_s   <- unname(dt[["sd"]])
    agg$total_user_s <- stats::median(tot_us)
    agg$total_sys_s  <- stats::median(tot_sy)
    agg$iterations_completed <- length(ok)
    agg$status <- "ok"
    agg$error_message <- if (length(errors)) paste(errors, collapse = "; ") else NA_character_
    agg
}

# rbind a list of data.frames with differing columns (union columns, NA-fill).
.wf_rbind_fill <- function(rows) {
    rows <- Filter(function(r) !is.null(r) && nrow(r) > 0L, rows)
    if (!length(rows)) return(NULL)
    cols <- unique(unlist(lapply(rows, names)))
    aligned <- lapply(rows, function(r) {
        for (m in setdiff(cols, names(r))) r[[m]] <- NA
        r[cols]
    })
    do.call(rbind, aligned)
}

# ---------------------------------------------------------------------------
# BAM inputs
# ---------------------------------------------------------------------------

.wf_file_size_mb <- function(paths) as.numeric(file.size(paths)) / 1024^2

# Resolve a chipseqDBData `Path` element (a BamFile) to a plain path string.
# Uses BiocGenerics::path (the generic BamFile honours) and takes the first
# element; a BamFile is an S4 object over an environment, so as.character() on
# the object itself must be avoided.
.wf_bamfile_path <- function(x) {
    p <- tryCatch(
        BiocGenerics::path(x),
        error = function(e) if ("path" %in% names(x)) x$path
                            else stop("Unable to extract BAM path: ", conditionMessage(e), call. = FALSE)
    )
    as.character(p)[[1L]]
}

.wf_fetch_chipseqdbdata <- function(n_target) {
    if (!requireNamespace("chipseqDBData", quietly = TRUE)) {
        stop("chipseqDBData is required to resolve default BAMs; supply --bam-dir instead.", call. = FALSE)
    }
    tbl <- chipseqDBData::H3K9acData()
    if (nrow(tbl) == 0L || !("Path" %in% colnames(tbl))) {
        stop("chipseqDBData::H3K9acData() returned no BAM paths.", call. = FALSE)
    }
    # Order by library size (Reads column) descending; strip any non-numeric text.
    reads_num <- if ("Reads" %in% colnames(tbl)) {
        suppressWarnings(as.numeric(gsub("[^0-9.]+", "", as.character(tbl$Reads))))
    } else rep(NA_real_, nrow(tbl))
    ord <- order(ifelse(is.na(reads_num), -Inf, reads_num), decreasing = TRUE)
    paths <- vapply(tbl$Path, .wf_bamfile_path, character(1))[ord]
    paths <- unique(paths[file.exists(paths)])
    if (!length(paths)) stop("No existing chipseqDBData BAM paths resolved.", call. = FALSE)
    utils::head(paths, max(1L, n_target))
}

.wf_select_n_files <- function(paths, n_files, allow_repeat = FALSE) {
    if (length(paths) >= n_files) return(paths[seq_len(n_files)])
    if (!allow_repeat) {
        stop(sprintf("Requested %d files but only %d available; pass --allow-repeat-files=true.",
                     n_files, length(paths)), call. = FALSE)
    }
    rep(paths, length.out = n_files)
}

.wf_list_bam_dir <- function(dir, recursive = FALSE) {
    if (is.null(dir) || !nzchar(dir)) return(character())
    if (!dir.exists(dir)) stop("--bam-dir does not exist: ", dir, call. = FALSE)
    list.files(dir, pattern = "\\.bam$", recursive = recursive,
               ignore.case = TRUE, full.names = TRUE)
}

.wf_ensure_index <- function(bam) {
    bai <- paste0(bam, ".bai")
    if (file.exists(bai) || file.exists(sub("\\.bam$", ".bai", bam))) return(invisible(TRUE))
    tryCatch({ Rsamtools::indexBam(bam); invisible(TRUE) },
             error = function(e) { .wf_log("index failed for ", basename(bam), ": ", conditionMessage(e)); invisible(FALSE) })
}

# Count records once per file (no index needed; sequential pass) for throughput.
.wf_count_records <- function(bam) {
    tryCatch(as.numeric(Rsamtools::countBam(bam)$records), error = function(e) NA_real_)
}

# ---------------------------------------------------------------------------
# Output writers (schema overlaps run_server_benchmark.R + phase columns)
# ---------------------------------------------------------------------------

.wf_write_csv <- function(df, path) {
    utils::write.csv(df, path, row.names = FALSE)
    invisible(path)
}

.wf_host_info <- function() {
    si <- Sys.info()
    cpu <- tryCatch({
        m <- grep("model name", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)
        if (length(m)) trimws(sub("^.*:", "", m[[1]])) else NA_character_
    }, error = function(e) NA_character_)
    memgb <- tryCatch({
        m <- grep("^MemTotal", readLines("/proc/meminfo", warn = FALSE), value = TRUE)
        if (length(m)) as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", m[[1]])) / 1024^2 else NA_real_
    }, error = function(e) NA_real_)
    data.frame(
        nodename = unname(si[["nodename"]]),
        sysname = unname(si[["sysname"]]),
        release = unname(si[["release"]]),
        cpu_model = cpu,
        logical_cores = parallel::detectCores(logical = TRUE),
        physical_cores = parallel::detectCores(logical = FALSE),
        detected_cores = .wf_detect_cores(),
        total_mem_gb = round(memgb, 1),
        r_version = paste(R.version$major, R.version$minor, sep = "."),
        omp_num_threads = Sys.getenv("OMP_NUM_THREADS", unset = NA),
        omp_proc_bind = Sys.getenv("OMP_PROC_BIND", unset = NA),
        stringsAsFactors = FALSE
    )
}

.wf_pkg_versions <- function(pkgs = c("BamScale", "Rsamtools", "GenomicAlignments",
                                      "rtracklayer", "ompBAM", "BiocParallel", "ATACseqQC")) {
    v <- vapply(pkgs, function(p) tryCatch(as.character(utils::packageVersion(p)),
                                           error = function(e) NA_character_), character(1))
    data.frame(package = pkgs, version = unname(v), stringsAsFactors = FALSE)
}
