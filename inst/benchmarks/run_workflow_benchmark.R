#!/usr/bin/env Rscript
# run_workflow_benchmark.R
#
# End-to-end workflow benchmarks for BamScale, complementing the read-pattern
# micro-benchmarks in run_server_benchmark.R. Two workflows are measured as a
# read / compute / write phase decomposition so the report can state the Amdahl
# bound on the end-to-end speedup explicitly:
#
#   1. coverage  : readGAlignments/bam_read -> coverage() [RleList] -> export.bw()
#   2. atacqc    : scanBam/bam_read isize    -> table(abs(isize)) fragment sizes
#
# In each workflow BamScale replaces exactly the BAM-read step; the compute and
# write steps are byte-identical work on both arms. Single-file cases sweep the
# OpenMP thread axis and report per-phase timing; multi-file cases sweep the
# BiocParallel worker axis (fixed core budget) and report end-to-end wall-clock
# throughput.
#
# Usage (defaults mirror the read benchmark's balanced profile):
#   Rscript inst/benchmarks/run_workflow_benchmark.R \
#     --threads=1,4,12,24,48 --workers=1,2,4,8,12 \
#     --n-files=12 --budget-threads=48 --max-threads=48 --bp-backend=snow \
#     --iterations=3 --warmup=true \
#     --include-coverage=true --include-atacqc=true \
#     --atac-bam-dir=/path/to/atac/bams --write-bigwig=true \
#     --outdir=benchmark_results

suppressWarnings(suppressMessages({
    ok <- requireNamespace("BamScale", quietly = TRUE) &&
        requireNamespace("Rsamtools", quietly = TRUE) &&
        requireNamespace("GenomicAlignments", quietly = TRUE) &&
        requireNamespace("rtracklayer", quietly = TRUE) &&
        requireNamespace("BiocParallel", quietly = TRUE)
}))
if (!ok) stop("Required packages missing (BamScale, Rsamtools, GenomicAlignments, rtracklayer, BiocParallel).")

.wf_script_dir <- function() {
    ca <- commandArgs(FALSE)
    f <- ca[grep("^--file=", ca)]
    if (length(f)) return(normalizePath(dirname(sub("^--file=", "", f[[1]])), mustWork = FALSE))
    getwd()
}
source(file.path(.wf_script_dir(), "bench_common.R"))

# ---------------------------------------------------------------------------
# Workflow pipelines (BamScale replaces only the read phase)
# ---------------------------------------------------------------------------

.wf_atac_flag <- function() {
    Rsamtools::scanBamFlag(isSecondaryAlignment = FALSE,
                           isUnmappedQuery = FALSE,
                           isNotPassingQualityControls = FALSE)
}

.wf_cov_single <- function(path, engine, threads, write_bigwig) {
    # Both arms explicitly exclude unmapped reads at the C level (readGAlignments
    # does this by default via isUnmappedQuery=FALSE), so neither decodes records
    # the other skips. Read exactly the 4 fields readGAlignments uses; omitting
    # `what` would default to 7 incl. the packed qname (an unfair extra decode).
    flt <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
    reader <- if (identical(engine, "BamScale")) {
        function() BamScale::bam_read(path, what = c("rname", "pos", "cigar", "strand"),
                                      as = "GAlignments", threads = threads, BPPARAM = NULL, param = flt)
    } else {
        function() GenomicAlignments::readGAlignments(path, param = flt)
    }
    rd <- .wf_time_phase(reader());              ga  <- rd$value; invisible(length(ga))
    cp <- .wf_time_phase(GenomicAlignments::coverage(ga)); cvg <- cp$value; invisible(length(cvg))
    if (isTRUE(write_bigwig)) {
        bw <- tempfile(fileext = ".bw")
        wr <- .wf_time_phase(rtracklayer::export.bw(cvg, bw))
        unlink(bw)
    } else {
        wr <- list(elapsed = 0, user = 0, sys = 0)
    }
    list(read    = c(elapsed = rd$elapsed, user = rd$user, sys = rd$sys),
         compute = c(elapsed = cp$elapsed, user = cp$user, sys = cp$sys),
         write   = c(elapsed = wr$elapsed, user = wr$user, sys = wr$sys))
}

.wf_cov_multi <- function(paths, engine, threads, bp, write_bigwig) {
    flt <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
    worker <- if (identical(engine, "BamScale")) {
        force(threads); force(write_bigwig); force(flt)
        function(p) {
            ga  <- BamScale::bam_read(p, what = c("rname", "pos", "cigar", "strand"),
                                      as = "GAlignments", threads = threads, BPPARAM = NULL, param = flt)
            cvg <- GenomicAlignments::coverage(ga)
            if (isTRUE(write_bigwig)) { bw <- tempfile(fileext = ".bw"); rtracklayer::export.bw(cvg, bw); unlink(bw) }
            length(cvg)
        }
    } else {
        force(write_bigwig); force(flt)
        function(p) {
            ga  <- GenomicAlignments::readGAlignments(p, param = flt)
            cvg <- GenomicAlignments::coverage(ga)
            if (isTRUE(write_bigwig)) { bw <- tempfile(fileext = ".bw"); rtracklayer::export.bw(cvg, bw); unlink(bw) }
            length(cvg)
        }
    }
    master <- .wf_time_phase(
        if (is.null(bp)) lapply(paths, worker) else BiocParallel::bplapply(paths, worker, BPPARAM = bp)
    )
    list(total = c(elapsed = master$elapsed, user = master$user, sys = master$sys))
}

.wf_atac_single <- function(path, engine, threads) {
    FLAG <- .wf_atac_flag()
    reader <- if (identical(engine, "BamScale")) {
        # Pass a real ScanBamParam (a plain scanBamFlag integer via list(flag=)
        # is NOT honored by bam_read and would apply no filter). Read only
        # `isize` to mirror scanBam(what="isize") exactly.
        function() BamScale::bam_read(path, what = "isize", as = "data.frame",
                                      threads = threads, BPPARAM = NULL,
                                      param = Rsamtools::ScanBamParam(flag = FLAG))
    } else {
        function() Rsamtools::scanBam(path, param = Rsamtools::ScanBamParam(what = "isize", flag = FLAG))[[1]]
    }
    rd  <- .wf_time_phase(reader()); dat <- rd$value
    cp  <- .wf_time_phase(table(abs(dat$isize)))
    list(read    = c(elapsed = rd$elapsed, user = rd$user, sys = rd$sys),
         compute = c(elapsed = cp$elapsed, user = cp$user, sys = cp$sys),
         write   = c(elapsed = 0, user = 0, sys = 0))
}

.wf_atac_multi <- function(paths, engine, threads, bp) {
    FLAG <- .wf_atac_flag()
    worker <- if (identical(engine, "BamScale")) {
        force(threads); force(FLAG)
        function(p) {
            d <- BamScale::bam_read(p, what = "isize", as = "data.frame",
                                    threads = threads, BPPARAM = NULL,
                                    param = Rsamtools::ScanBamParam(flag = FLAG))
            length(table(abs(d$isize)))
        }
    } else {
        force(FLAG)
        function(p) {
            d <- Rsamtools::scanBam(p, param = Rsamtools::ScanBamParam(what = "isize", flag = FLAG))[[1]]
            length(table(abs(d$isize)))
        }
    }
    master <- .wf_time_phase(
        if (is.null(bp)) lapply(paths, worker) else BiocParallel::bplapply(paths, worker, BPPARAM = bp)
    )
    list(total = c(elapsed = master$elapsed, user = master$user, sys = master$sys))
}

# ---------------------------------------------------------------------------
# Correctness (equivalence anchors)
# ---------------------------------------------------------------------------

.wf_correctness_coverage <- function(path) {
    tryCatch({
        flt <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
        ga_std <- GenomicAlignments::readGAlignments(path, param = flt)
        ga_bs  <- BamScale::bam_read(path, what = c("rname", "pos", "cigar", "strand"),
                                     as = "GAlignments", threads = 1L, BPPARAM = NULL, param = flt)
        same_levels  <- identical(GenomeInfoDb::seqlevels(ga_std), GenomeInfoDb::seqlevels(ga_bs))
        same_lengths <- identical(GenomeInfoDb::seqlengths(ga_std), GenomeInfoDb::seqlengths(ga_bs))
        same_n       <- identical(length(ga_std), length(ga_bs))
        same_cov     <- identical(GenomicAlignments::coverage(ga_std), GenomicAlignments::coverage(ga_bs))
        data.frame(workload = "coverage", comparator = "GenomicAlignments::coverage",
                   passed = isTRUE(same_levels && same_lengths && same_n && same_cov),
                   detail = sprintf("levels=%s lengths=%s n=%s coverage=%s",
                                    same_levels, same_lengths, same_n, same_cov),
                   status = "ok", stringsAsFactors = FALSE)
    }, error = function(e) data.frame(workload = "coverage", comparator = "GenomicAlignments::coverage",
                                      passed = FALSE, detail = conditionMessage(e), status = "error",
                                      stringsAsFactors = FALSE))
}

.wf_correctness_atac <- function(path) {
    tryCatch({
        FLAG <- .wf_atac_flag()
        std <- Rsamtools::scanBam(path, param = Rsamtools::ScanBamParam(what = "isize", flag = FLAG))[[1]]
        bs  <- BamScale::bam_read(path, what = "isize", as = "data.frame",
                                  threads = 1L, BPPARAM = NULL,
                                  param = Rsamtools::ScanBamParam(flag = FLAG))
        t_std <- table(abs(std$isize)); t_bs <- table(abs(bs$isize))
        same <- identical(t_std, t_bs)
        data.frame(workload = "atacqc", comparator = "Rsamtools::scanBam(isize)",
                   passed = isTRUE(same),
                   detail = sprintf("fragment-size table identical=%s (std n=%d, bs n=%d)",
                                    same, sum(t_std), sum(t_bs)),
                   status = "ok", stringsAsFactors = FALSE)
    }, error = function(e) data.frame(workload = "atacqc", comparator = "Rsamtools::scanBam(isize)",
                                      passed = FALSE, detail = conditionMessage(e), status = "error",
                                      stringsAsFactors = FALSE))
}

# One-time anchor against the real ATACseqQC package (not on the timing path).
.wf_atacseqqc_anchor <- function(path) {
    if (!requireNamespace("ATACseqQC", quietly = TRUE)) {
        return(data.frame(check = "ATACseqQC::fragSizeDist", passed = NA,
                          detail = "ATACseqQC not installed; anchor skipped",
                          version = NA_character_, stringsAsFactors = FALSE))
    }
    tryCatch({
        FLAG <- .wf_atac_flag()
        ours <- table(abs(Rsamtools::scanBam(path, param = Rsamtools::ScanBamParam(what = "isize", flag = FLAG))[[1]]$isize))
        fd <- ATACseqQC::fragSizeDist(path, basename(path))
        pkg <- fd[[1]]
        # Align both as named integer vectors over the union of fragment lengths.
        keys <- union(names(ours), names(pkg))
        a <- setNames(as.integer(ours[keys]), keys); a[is.na(a)] <- 0L
        b <- setNames(as.integer(pkg[keys]),  keys); b[is.na(b)] <- 0L
        same <- identical(a, b)
        data.frame(check = "ATACseqQC::fragSizeDist", passed = isTRUE(same),
                   detail = sprintf("fragment-size table identical=%s (ours=%d, pkg=%d)", same, sum(a), sum(b)),
                   version = as.character(utils::packageVersion("ATACseqQC")), stringsAsFactors = FALSE)
    }, error = function(e) data.frame(check = "ATACseqQC::fragSizeDist", passed = FALSE,
                                      detail = conditionMessage(e),
                                      version = tryCatch(as.character(utils::packageVersion("ATACseqQC")), error = function(x) NA_character_),
                                      stringsAsFactors = FALSE))
}

# ---------------------------------------------------------------------------
# Case execution
# ---------------------------------------------------------------------------

.wf_run_case <- function(meta, run_once, iterations, warmup) {
    .wf_log(sprintf("[%s | %s | %s] req_t=%s workers=%s files=%d",
                    meta$scenario, meta$workload, meta$method,
                    meta$threads_requested, meta$bp_workers_effective, meta$n_files))
    agg <- .wf_bench_phases(run_once, iterations = iterations, warmup = warmup)
    row <- data.frame(
        scenario = meta$scenario, workload = meta$workload,
        method = meta$method, method_family = meta$method_family,
        comparison_track = if (!is.null(meta$comparison_track)) meta$comparison_track else "fair",
        threads_requested = meta$threads_requested, threads_effective = meta$threads_effective,
        bp_workers_requested = meta$bp_workers_requested, bp_workers_effective = meta$bp_workers_effective,
        # Total OS cores the arm actually used = OpenMP threads x BiocParallel
        # workers. Exposes the resource asymmetry the audit flagged (BamScale can
        # use more cores per file than the single-threaded standard reader).
        cores_used = meta$threads_effective * meta$bp_workers_effective,
        n_files = meta$n_files, n_records = meta$n_records, total_mb = meta$total_mb,
        iterations_planned = iterations,
        iterations_completed = if (!is.null(agg$iterations_completed)) agg$iterations_completed else 0L,
        status = agg$status,
        error_message = if (!is.null(agg$error_message)) agg$error_message else NA_character_,
        stringsAsFactors = FALSE
    )
    # Attach every timing/dispersion column the engine produced (read/compute/
    # write/total x median/min/max/sd/user/sys). Missing columns are NA-filled at
    # rbind time by .wf_rbind_fill, so error rows stay aligned.
    if (identical(agg$status, "ok")) {
        for (nm in names(agg)) if (grepl("_s$", nm)) row[[nm]] <- agg[[nm]]
    }
    total_s <- if (!is.null(row$total_s)) row$total_s else NA_real_
    read_s  <- if (!is.null(row$read_s))  row$read_s  else NA_real_
    row$total_mb_per_s      <- if (!is.na(total_s) && total_s > 0) row$total_mb  / total_s else NA_real_
    row$read_mb_per_s       <- if (!is.na(read_s)  && read_s  > 0) row$total_mb  / read_s  else NA_real_
    # Denominator is the unfiltered countBam total, not the filtered records each
    # workflow processes -> label as INPUT records/s (ratio-neutral, absolute only).
    row$input_records_per_s <- if (!is.na(total_s) && total_s > 0) row$n_records / total_s else NA_real_
    row
}

.wf_build_meta <- function(scenario, workload, method, method_family,
                           threads_req, threads_eff, workers_req, workers_eff,
                           n_files, n_records, total_mb, comparison_track = "fair") {
    list(scenario = scenario, workload = workload, method = method, method_family = method_family,
         threads_requested = threads_req, threads_effective = threads_eff,
         bp_workers_requested = workers_req, bp_workers_effective = workers_eff,
         n_files = n_files, n_records = n_records, total_mb = total_mb,
         comparison_track = comparison_track)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

args <- .wf_parse_args()

cfg <- list(
    threads       = .wf_int_vec(args$threads, c(1L, 4L, 12L, 24L, 48L)),
    workers       = .wf_int_vec(args$workers, c(1L, 2L, 4L, 8L, 12L)),
    n_files       = .wf_as_int(args$n_files, 12L),
    budget        = .wf_as_int(args$budget_threads, 48L),
    max_threads   = .wf_as_int(args$max_threads, 48L),
    iterations    = .wf_as_int(args$iterations, 5L),
    warmup        = .wf_as_bool(args$warmup, TRUE),
    bp_backend    = if (!is.null(args$bp_backend)) args$bp_backend else "snow",
    include_coverage = .wf_as_bool(args$include_coverage, TRUE),
    include_atacqc   = .wf_as_bool(args$include_atacqc, !is.null(args$atac_bam_dir)),
    include_single   = .wf_as_bool(args$include_single, TRUE),
    include_multi    = .wf_as_bool(args$include_multi, TRUE),
    write_bigwig  = .wf_as_bool(args$write_bigwig, TRUE),
    multi_write_bigwig = .wf_as_bool(args$multi_write_bigwig, .wf_as_bool(args$write_bigwig, TRUE)),
    allow_repeat  = .wf_as_bool(args$allow_repeat_files, TRUE),
    preflight     = .wf_as_bool(args$correctness_preflight, TRUE),
    atac_anchor   = .wf_as_bool(args$atac_anchor, FALSE),
    bam_dir       = args$bam_dir,
    atac_bam_dir  = args$atac_bam_dir,
    outdir        = if (!is.null(args$outdir)) args$outdir else "benchmark_results"
)

# Clamp the core budget to the machine so threads_each never oversubscribes.
cfg$budget <- max(1L, min(cfg$budget, cfg$max_threads, .wf_detect_cores()))

out_dir <- file.path(cfg$outdir, format(Sys.time(), "wf_run_%Y%m%d_%H%M%S"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
.wf_log("output dir: ", out_dir, " | core budget: ", cfg$budget)

# --- resolve BAM inputs --------------------------------------------------
cov_bams <- character()
if (cfg$include_coverage) {
    cov_bams <- if (!is.null(cfg$bam_dir)) .wf_list_bam_dir(cfg$bam_dir) else .wf_fetch_chipseqdbdata(max(cfg$n_files, 2L))
    cov_bams <- cov_bams[order(.wf_file_size_mb(cov_bams), decreasing = TRUE)]
    if (!length(cov_bams)) stop("No coverage BAMs resolved.", call. = FALSE)
}
atac_bams <- character()
if (cfg$include_atacqc) {
    if (is.null(cfg$atac_bam_dir)) stop("--include-atacqc requires --atac-bam-dir (real paired-end ATAC BAMs).", call. = FALSE)
    atac_bams <- .wf_list_bam_dir(cfg$atac_bam_dir)
    atac_bams <- atac_bams[order(.wf_file_size_mb(atac_bams), decreasing = TRUE)]
    if (!length(atac_bams)) stop("No ATAC BAMs found in --atac-bam-dir.", call. = FALSE)
    invisible(lapply(atac_bams, .wf_ensure_index))
}

# --- reference record counts + sizes (once per unique file) --------------
.wf_log("counting records for throughput reference ...")
all_bams <- unique(c(cov_bams, atac_bams))
rec_map <- setNames(vapply(all_bams, .wf_count_records, numeric(1)), all_bams)
mb_map  <- setNames(.wf_file_size_mb(all_bams), all_bams)

files_df <- data.frame(
    file = all_bams, basename = basename(all_bams),
    size_mb = round(unname(mb_map[all_bams]), 1),
    n_records = unname(rec_map[all_bams]),
    workflow = ifelse(all_bams %in% atac_bams, "atacqc", "coverage"),
    stringsAsFactors = FALSE
)

# --- correctness preflight ----------------------------------------------
correctness <- list()
if (cfg$preflight) {
    if (cfg$include_coverage) { .wf_log("preflight: coverage equivalence ..."); correctness[["cov"]] <- .wf_correctness_coverage(cov_bams[[1]]) }
    if (cfg$include_atacqc)   { .wf_log("preflight: atac fragment-size equivalence ..."); correctness[["atac"]] <- .wf_correctness_atac(atac_bams[[1]]) }
}
if (cfg$atac_anchor && cfg$include_atacqc) {
    .wf_log("ATACseqQC anchor ...")
    anchor <- .wf_atacseqqc_anchor(atac_bams[[1]])
    .wf_write_csv(anchor, file.path(out_dir, "atacseqqc_anchor.csv"))
}

# --- build + run cases ---------------------------------------------------
summary_rows <- list()
add <- function(row) summary_rows[[length(summary_rows) + 1L]] <<- row

run_single <- function(workload, bams) {
    single <- bams[[1]]
    n_rec <- unname(rec_map[single]); mb <- unname(mb_map[single])
    std_family <- if (identical(workload, "coverage")) "GenomicAlignments" else "Rsamtools"
    std_method <- if (identical(workload, "coverage")) "GenomicAlignments::readGAlignments -> coverage -> export.bw"
                  else "Rsamtools::scanBam(isize) -> fragment-size table"
    once_bs  <- function(t) if (identical(workload, "coverage"))
        function() .wf_cov_single(single, "BamScale", t, cfg$write_bigwig) else function() .wf_atac_single(single, "BamScale", t)
    # BamScale across the thread grid
    for (t in cfg$threads) {
        meta <- .wf_build_meta("single", workload, "BamScale", "BamScale", t, t, 1L, 1L, 1L, n_rec, mb, "single-file")
        add(.wf_run_case(meta, once_bs(t), cfg$iterations, cfg$warmup))
    }
    # Standard reader once (single-threaded, flat reference)
    once_std <- if (identical(workload, "coverage")) function() .wf_cov_single(single, "std", 1L, cfg$write_bigwig)
                else function() .wf_atac_single(single, "std", 1L)
    meta <- .wf_build_meta("single", workload, std_method, std_family, 1L, 1L, 1L, 1L, 1L, n_rec, mb, "single-file")
    add(.wf_run_case(meta, once_std, cfg$iterations, cfg$warmup))
}

run_multi <- function(workload, bams) {
    files <- .wf_select_n_files(bams, cfg$n_files, cfg$allow_repeat)
    n_rec <- sum(unname(rec_map[files]), na.rm = TRUE); mb <- sum(unname(mb_map[files]))
    std_family <- if (identical(workload, "coverage")) "GenomicAlignments" else "Rsamtools"
    std_method <- if (identical(workload, "coverage")) "GenomicAlignments + BiocParallel" else "Rsamtools::scanBam + BiocParallel"
    pkgs <- if (identical(workload, "coverage")) c("BamScale", "GenomicAlignments", "rtracklayer") else c("BamScale", "Rsamtools")
    once <- function(engine, threads, bp) if (identical(workload, "coverage"))
        function() .wf_cov_multi(files, engine, threads, bp, cfg$multi_write_bigwig)
        else function() .wf_atac_multi(files, engine, threads, bp)

    # One BamScale + one standard case at (W workers x threads_each / 1). The
    # PSOCK cluster is started and its reader namespaces pre-loaded OUTSIDE the
    # timed region, so cluster spin-up and per-worker library loading are not
    # charged to every iteration (audit: otherwise they dilute the ratio->1.0).
    run_pair <- function(w, threads_each, track) {
        bp <- .wf_make_bpparam(cfg$bp_backend, w)
        if (!is.null(bp)) {
            BiocParallel::bpstart(bp)
            try(BiocParallel::bplapply(seq_len(w), function(i, pk) {
                for (p in pk) suppressMessages(requireNamespace(p, quietly = TRUE)); TRUE
            }, pk = pkgs, BPPARAM = bp), silent = TRUE)
        }
        add(.wf_run_case(.wf_build_meta("multi", workload, "BamScale", "BamScale",
            threads_each, threads_each, w, w, length(files), n_rec, mb, track),
            once("BamScale", threads_each, bp), cfg$iterations, cfg$warmup))
        add(.wf_run_case(.wf_build_meta("multi", workload, std_method, std_family,
            1L, 1L, w, w, length(files), n_rec, mb, track),
            once("std", 1L, bp), cfg$iterations, cfg$warmup))
        if (!is.null(bp)) try(BiocParallel::bpstop(bp), silent = TRUE)
    }

    # Full-budget sweep (W>=2): BamScale spends the whole core budget
    # (floor(budget/W) threads x W workers ~= budget cores); standard uses W cores
    # (1 thread/worker). cores_used in the output exposes the asymmetry. W=1 is
    # excluded (it is the single-file case, not a multi-file worker point).
    workers_grid <- cfg$workers[cfg$workers >= 2L & cfg$workers <= length(files)]
    for (w in workers_grid) run_pair(w, max(1L, floor(cfg$budget / w)), "full-budget")

    # Matched-cores comparison: both arms use C = min(budget, files) cores as
    # C workers x 1 thread, so the head-to-head ratio is apples-to-apples.
    C <- min(cfg$budget, length(files))
    if (C > 1L) run_pair(C, 1L, "matched-cores")
}

for (wl in c(if (cfg$include_coverage) "coverage", if (cfg$include_atacqc) "atacqc")) {
    # Do not emit timing for a comparison the equivalence preflight failed -- the
    # arms would be measuring different record sets (audit finding).
    if (identical(wl, "atacqc") && cfg$preflight &&
        !is.null(correctness[["atac"]]) && !isTRUE(correctness[["atac"]]$passed)) {
        .wf_log("SKIPPING atacqc timing: correctness preflight FAILED (non-equivalent record sets).")
        next
    }
    bams <- if (identical(wl, "coverage")) cov_bams else atac_bams
    if (cfg$include_single) run_single(wl, bams)
    if (cfg$include_multi)  run_multi(wl, bams)
}

# --- write outputs -------------------------------------------------------
summary_df <- .wf_rbind_fill(summary_rows)
.wf_write_csv(summary_df, file.path(out_dir, "summary.csv"))
.wf_write_csv(summary_df, file.path(out_dir, "summary_qmd.csv"))
.wf_write_csv(files_df, file.path(out_dir, "files.csv"))
.wf_write_csv(files_df[, c("workflow", "n_records", "size_mb")], file.path(out_dir, "reference_counts.csv"))
if (length(correctness)) .wf_write_csv(do.call(rbind, correctness), file.path(out_dir, "correctness_workflow.csv"))
.wf_write_csv(.wf_host_info(), file.path(out_dir, "host_info.csv"))
.wf_write_csv(.wf_pkg_versions(), file.path(out_dir, "pkg_versions.csv"))

cfg_lines <- vapply(names(cfg), function(k) {
    v <- cfg[[k]]
    sprintf("%s: %s", k, paste(as.character(v), collapse = ","))
}, character(1))
writeLines(c(cfg_lines,
             paste0("cov_bams: ", paste(basename(cov_bams), collapse = ",")),
             paste0("atac_bams: ", paste(basename(atac_bams), collapse = ","))),
           file.path(out_dir, "config.txt"))
writeLines(capture.output(utils::sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

.wf_log("done. ", nrow(summary_df), " cases written to ", out_dir)
