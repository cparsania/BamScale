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
# New-API workloads: the four C++ aggregation endpoints as first-class arms.
# The BamScale side is a single fused call (read+compute+write inside C++), so
# it reports total-only; the standard side reuses the phase-decomposed recipes.
# ---------------------------------------------------------------------------

.wf_total_only <- function(expr_fn) {
    tm <- .wf_time_phase(expr_fn())
    list(total = c(elapsed = tm$elapsed, user = tm$user, sys = tm$sys))
}

.wf_mapq_single_std <- function(path) {
    rd <- .wf_time_phase(Rsamtools::scanBam(path, param = Rsamtools::ScanBamParam(what = "mapq"))[[1]])
    cp <- .wf_time_phase(table(rd$value$mapq, useNA = "ifany"))
    list(read    = c(elapsed = rd$elapsed, user = rd$user, sys = rd$sys),
         compute = c(elapsed = cp$elapsed, user = cp$user, sys = cp$sys),
         write   = c(elapsed = 0, user = 0, sys = 0))
}

# Registry: per API workload, the BamScale fused call, the standard phase recipe,
# the per-file worker closures for multi, and metadata. `t` = OpenMP threads.
.wf_api_registry <- function(cfg) {
    FLAG <- .wf_atac_flag()
    fragP <- Rsamtools::ScanBamParam(flag = FLAG)
    list(
        fastcov = list(
            include = cfg$include_fastcov,
            std_method = "GenomicAlignments::readGAlignments -> coverage",
            std_family = "GenomicAlignments",
            pkgs = c("BamScale", "GenomicAlignments"),
            bs_single = function(f, t) .wf_total_only(function()
                BamScale::bam_coverage(f, threads = t, BPPARAM = NULL)),
            std_single = function(f) .wf_cov_single(f, "std", 1L, write_bigwig = FALSE),
            bs_worker = function(t) { force(t); function(p)
                length(BamScale::bam_coverage(p, threads = t, BPPARAM = NULL)) },
            std_worker = function() function(p) {
                flt <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
                length(GenomicAlignments::coverage(GenomicAlignments::readGAlignments(p, param = flt)))
            }
        ),
        bigwig = list(
            include = cfg$include_bigwig,
            std_method = "GenomicAlignments::readGAlignments -> coverage -> export.bw",
            std_family = "GenomicAlignments",
            pkgs = c("BamScale", "GenomicAlignments", "rtracklayer"),
            bs_single = function(f, t) .wf_total_only(function() {
                bw <- tempfile(fileext = ".bw")
                on.exit(unlink(bw), add = TRUE)
                BamScale::bam_coverage_bigwig(f, bw, threads = t, parallel = TRUE, BPPARAM = NULL)
            }),
            std_single = function(f) .wf_cov_single(f, "std", 1L, write_bigwig = TRUE),
            bs_worker = function(t) { force(t); function(p) {
                bw <- tempfile(fileext = ".bw")
                on.exit(unlink(bw), add = TRUE)
                BamScale::bam_coverage_bigwig(p, bw, threads = t, parallel = TRUE, BPPARAM = NULL)
                1L
            } },
            std_worker = function() function(p) {
                flt <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
                cvg <- GenomicAlignments::coverage(GenomicAlignments::readGAlignments(p, param = flt))
                bw <- tempfile(fileext = ".bw")
                on.exit(unlink(bw), add = TRUE)
                rtracklayer::export.bw(cvg, bw)
                1L
            }
        ),
        fragsize = list(
            include = cfg$include_fragsize,
            std_method = "Rsamtools::scanBam(isize) -> fragment-size table",
            std_family = "Rsamtools",
            pkgs = c("BamScale", "Rsamtools"),
            bs_single = function(f, t) .wf_total_only(function()
                BamScale::fragment_sizes(f, param = fragP, threads = t, BPPARAM = NULL)),
            std_single = function(f) .wf_atac_single(f, "std", 1L),
            bs_worker = function(t) { force(t); force(fragP); function(p)
                nrow(BamScale::fragment_sizes(p, param = fragP, threads = t, BPPARAM = NULL)) },
            std_worker = function() { force(FLAG); function(p) {
                d <- Rsamtools::scanBam(p, param = Rsamtools::ScanBamParam(what = "isize", flag = FLAG))[[1]]
                length(table(abs(d$isize)))
            } }
        ),
        mapq = list(
            include = cfg$include_mapq,
            std_method = "Rsamtools::scanBam(mapq) -> table",
            std_family = "Rsamtools",
            pkgs = c("BamScale", "Rsamtools"),
            bs_single = function(f, t) .wf_total_only(function()
                BamScale::mapq_dist(f, threads = t, BPPARAM = NULL)),
            std_single = function(f) .wf_mapq_single_std(f),
            bs_worker = function(t) { force(t); function(p)
                nrow(BamScale::mapq_dist(p, threads = t, BPPARAM = NULL)) },
            std_worker = function() function(p) {
                d <- Rsamtools::scanBam(p, param = Rsamtools::ScanBamParam(what = "mapq"))[[1]]
                length(table(d$mapq, useNA = "ifany"))
            }
        )
    )
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

.wf_gate_row <- function(workload, comparator, passed, detail, status = "ok") {
    data.frame(workload = workload, comparator = comparator, passed = isTRUE(passed),
               detail = detail, status = status, stringsAsFactors = FALSE)
}

.wf_correctness_fastcov <- function(path) {
    tryCatch({
        flt <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
        std <- GenomicAlignments::coverage(GenomicAlignments::readGAlignments(path, param = flt))
        bs  <- BamScale::bam_coverage(path, threads = 2L, BPPARAM = NULL)
        same_levels  <- identical(names(std), names(bs))
        same_lengths <- identical(vapply(std, length, numeric(1)), vapply(bs, length, numeric(1)))
        same_cov     <- identical(bs, std)
        .wf_gate_row("fastcov", "GenomicAlignments::coverage",
                     same_levels && same_lengths && same_cov,
                     sprintf("levels=%s lengths=%s coverage_identical=%s",
                             same_levels, same_lengths, same_cov))
    }, error = function(e) .wf_gate_row("fastcov", "GenomicAlignments::coverage",
                                        FALSE, conditionMessage(e), "error"))
}

.wf_correctness_bigwig <- function(path) {
    tryCatch({
        flt <- Rsamtools::ScanBamParam(flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE))
        bw_bs  <- tempfile(fileext = ".bw"); bw_std <- tempfile(fileext = ".bw")
        on.exit(unlink(c(bw_bs, bw_std)), add = TRUE)
        BamScale::bam_coverage_bigwig(path, bw_bs, threads = 2L, parallel = TRUE, BPPARAM = NULL)
        cvg <- GenomicAlignments::coverage(GenomicAlignments::readGAlignments(path, param = flt))
        rtracklayer::export.bw(cvg, bw_std)
        # Both files pass through the same importer, so equal values => identical objects.
        imp_bs  <- rtracklayer::import(bw_bs,  as = "RleList")
        imp_std <- rtracklayer::import(bw_std, as = "RleList")
        common <- intersect(names(imp_bs), names(imp_std))
        auc <- function(rl) sum(vapply(rl, function(r)
            sum(as.numeric(S4Vectors::runValue(r)) * as.numeric(S4Vectors::runLength(r))), numeric(1)))
        same_vals <- identical(imp_bs[common], imp_std[common])
        same_auc  <- isTRUE(all.equal(auc(imp_bs[common]), auc(imp_std[common])))
        .wf_gate_row("bigwig", "rtracklayer::export.bw round-trip",
                     same_vals && same_auc,
                     sprintf("reimported_identical=%s auc_equal=%s contigs=%d",
                             same_vals, same_auc, length(common)))
    }, error = function(e) .wf_gate_row("bigwig", "rtracklayer::export.bw round-trip",
                                        FALSE, conditionMessage(e), "error"))
}

.wf_correctness_fragsize <- function(path) {
    tryCatch({
        FLAG <- .wf_atac_flag()
        std <- Rsamtools::scanBam(path, param = Rsamtools::ScanBamParam(what = "isize", flag = FLAG))[[1]]
        t_std <- table(abs(std$isize))
        fs <- BamScale::fragment_sizes(path, param = Rsamtools::ScanBamParam(flag = FLAG),
                                       threads = 2L, BPPARAM = NULL)
        keys <- union(names(t_std), as.character(fs$fragment_size))
        a <- setNames(as.integer(t_std[keys]), keys); a[is.na(a)] <- 0L
        b <- setNames(as.integer(fs$count[match(keys, as.character(fs$fragment_size))]), keys)
        b[is.na(b)] <- 0L
        same <- identical(a, b)
        .wf_gate_row("fragsize", "Rsamtools::scanBam(isize) table",
                     same, sprintf("fragment-size table identical=%s (std n=%d, bs n=%.0f)",
                                   same, sum(a), sum(fs$count)))
    }, error = function(e) .wf_gate_row("fragsize", "Rsamtools::scanBam(isize) table",
                                        FALSE, conditionMessage(e), "error"))
}

.wf_correctness_mapq <- function(path) {
    tryCatch({
        std <- Rsamtools::scanBam(path, param = Rsamtools::ScanBamParam(what = "mapq"))[[1]]
        t_std <- table(std$mapq, useNA = "ifany")
        # Rsamtools reports MAPQ 255 ("unavailable") as NA; mapq_dist keeps 255.
        nm <- names(t_std); nm[is.na(nm)] <- "255"
        names(t_std) <- nm
        md <- BamScale::mapq_dist(path, threads = 2L, BPPARAM = NULL)
        keys <- union(names(t_std), as.character(md$mapq))
        a <- setNames(as.integer(t_std[keys]), keys); a[is.na(a)] <- 0L
        b <- setNames(as.integer(md$count[match(keys, as.character(md$mapq))]), keys)
        b[is.na(b)] <- 0L
        same <- identical(a[order(as.integer(names(a)))], b[order(as.integer(names(b)))])
        .wf_gate_row("mapq", "Rsamtools::scanBam(mapq) table (NA<->255)",
                     same, sprintf("mapq table identical=%s (std n=%d, bs n=%.0f)",
                                   same, sum(a), sum(md$count)))
    }, error = function(e) .wf_gate_row("mapq", "Rsamtools::scanBam(mapq) table (NA<->255)",
                                        FALSE, conditionMessage(e), "error"))
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
    # New-API workloads (fused C++ aggregation endpoints); default off.
    include_fastcov  = .wf_as_bool(args$include_fastcov, FALSE),
    include_bigwig   = .wf_as_bool(args$include_bigwig, FALSE),
    include_fragsize = .wf_as_bool(args$include_fragsize, FALSE),
    include_mapq     = .wf_as_bool(args$include_mapq, FALSE),
    megadepth_bin    = args$megadepth_bin,
    megadepth_threads = .wf_int_vec(args$megadepth_threads, c(1L, 8L, 48L, 72L)),
    outdir        = if (!is.null(args$outdir)) args$outdir else "benchmark_results"
)

# Clamp the core budget to the machine so threads_each never oversubscribes.
cfg$budget <- max(1L, min(cfg$budget, cfg$max_threads, .wf_detect_cores()))

out_dir <- file.path(cfg$outdir, format(Sys.time(), "wf_run_%Y%m%d_%H%M%S"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
.wf_log("output dir: ", out_dir, " | core budget: ", cfg$budget)

cfg$any_api <- cfg$include_fastcov || cfg$include_bigwig ||
    cfg$include_fragsize || cfg$include_mapq

# --- resolve BAM inputs --------------------------------------------------
cov_bams <- character()
if (cfg$include_coverage) {
    cov_bams <- if (!is.null(cfg$bam_dir)) .wf_list_bam_dir(cfg$bam_dir) else .wf_fetch_chipseqdbdata(max(cfg$n_files, 2L))
    cov_bams <- cov_bams[order(.wf_file_size_mb(cov_bams), decreasing = TRUE)]
    if (!length(cov_bams)) stop("No coverage BAMs resolved.", call. = FALSE)
}
atac_bams <- character()
if (cfg$include_atacqc || cfg$any_api) {
    if (is.null(cfg$atac_bam_dir)) stop("--include-atacqc / API workloads require --atac-bam-dir (real paired-end ATAC BAMs).", call. = FALSE)
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
    if (cfg$include_fastcov)  { .wf_log("preflight: fastcov (bam_coverage) equivalence ..."); correctness[["fastcov"]] <- .wf_correctness_fastcov(atac_bams[[1]]) }
    if (cfg$include_bigwig)   { .wf_log("preflight: bigwig round-trip equivalence ..."); correctness[["bigwig"]] <- .wf_correctness_bigwig(atac_bams[[1]]) }
    if (cfg$include_fragsize) { .wf_log("preflight: fragsize (fragment_sizes) equivalence ..."); correctness[["fragsize"]] <- .wf_correctness_fragsize(atac_bams[[1]]) }
    if (cfg$include_mapq)     { .wf_log("preflight: mapq (mapq_dist) equivalence ..."); correctness[["mapq"]] <- .wf_correctness_mapq(atac_bams[[1]]) }
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

# --- new-API workload cases ----------------------------------------------

run_single_api <- function(wl, spec, bams) {
    single <- bams[[1]]
    n_rec <- unname(rec_map[single]); mb <- unname(mb_map[single])
    for (t in cfg$threads) {
        meta <- .wf_build_meta("single", wl, "BamScale", "BamScale", t, t, 1L, 1L, 1L, n_rec, mb, "single-file")
        add(.wf_run_case(meta, function() spec$bs_single(single, t), cfg$iterations, cfg$warmup))
    }
    meta <- .wf_build_meta("single", wl, spec$std_method, spec$std_family, 1L, 1L, 1L, 1L, 1L, n_rec, mb, "single-file")
    add(.wf_run_case(meta, function() spec$std_single(single), cfg$iterations, cfg$warmup))

    # Optional context arm: megadepth on the bigwig workload (its own tool, own
    # conventions -- report-only track, never a gated comparison).
    if (identical(wl, "bigwig") && !is.null(cfg$megadepth_bin) && file.exists(cfg$megadepth_bin)) {
        md_threads <- cfg$megadepth_threads[cfg$megadepth_threads <= cfg$max_threads]
        for (t in md_threads) {
            once_md <- local({ tt <- t; function() .wf_total_only(function() {
                pref <- tempfile()
                on.exit(unlink(paste0(pref, c(".all.bw", ".unique.bw", ".annotation.tsv"))), add = TRUE)
                system2(cfg$megadepth_bin,
                        c(shQuote(single), "--bigwig", "--double-count",
                          "--threads", tt, "--prefix", shQuote(pref), "--no-auc-stdout"),
                        stdout = FALSE, stderr = FALSE)
            }) })
            meta <- .wf_build_meta("single", wl, "megadepth --bigwig", "megadepth",
                                   t, t, 1L, 1L, 1L, n_rec, mb, "context")
            add(.wf_run_case(meta, once_md, cfg$iterations, cfg$warmup))
        }
    }
}

run_multi_api <- function(wl, spec, bams) {
    files <- .wf_select_n_files(bams, cfg$n_files, cfg$allow_repeat)
    n_rec <- sum(unname(rec_map[files]), na.rm = TRUE); mb <- sum(unname(mb_map[files]))

    run_pair <- function(w, threads_each, track) {
        bp <- .wf_make_bpparam(cfg$bp_backend, w)
        if (!is.null(bp)) {
            BiocParallel::bpstart(bp)
            try(BiocParallel::bplapply(seq_len(w), function(i, pk) {
                for (p in pk) suppressMessages(requireNamespace(p, quietly = TRUE)); TRUE
            }, pk = spec$pkgs, BPPARAM = bp), silent = TRUE)
        }
        bs_worker <- spec$bs_worker(threads_each)
        std_worker <- spec$std_worker()
        once <- function(worker) function() {
            master <- .wf_time_phase(
                if (is.null(bp)) lapply(files, worker) else BiocParallel::bplapply(files, worker, BPPARAM = bp)
            )
            list(total = c(elapsed = master$elapsed, user = master$user, sys = master$sys))
        }
        add(.wf_run_case(.wf_build_meta("multi", wl, "BamScale", "BamScale",
            threads_each, threads_each, w, w, length(files), n_rec, mb, track),
            once(bs_worker), cfg$iterations, cfg$warmup))
        add(.wf_run_case(.wf_build_meta("multi", wl, spec$std_method, spec$std_family,
            1L, 1L, w, w, length(files), n_rec, mb, track),
            once(std_worker), cfg$iterations, cfg$warmup))
        if (!is.null(bp)) try(BiocParallel::bpstop(bp), silent = TRUE)
    }

    workers_grid <- cfg$workers[cfg$workers >= 2L & cfg$workers <= length(files)]
    for (w in workers_grid) run_pair(w, max(1L, floor(cfg$budget / w)), "full-budget")
    C <- min(cfg$budget, length(files))
    if (C > 1L) run_pair(C, 1L, "matched-cores")
}

if (cfg$any_api) {
    api_reg <- .wf_api_registry(cfg)
    for (wl in names(api_reg)) {
        spec <- api_reg[[wl]]
        if (!isTRUE(spec$include)) next
        if (cfg$preflight && !is.null(correctness[[wl]]) && !isTRUE(correctness[[wl]]$passed)) {
            .wf_log(sprintf("SKIPPING %s timing: correctness preflight FAILED.", wl))
            next
        }
        if (cfg$include_single) run_single_api(wl, spec, atac_bams)
        if (cfg$include_multi)  run_multi_api(wl, spec, atac_bams)
    }
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
