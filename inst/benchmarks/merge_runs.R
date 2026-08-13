#!/usr/bin/env Rscript
# Merge a single-scenario run dir and a multi-scenario run dir (produced by the
# staged orchestrator) into one directory shaped like a normal harness run, so
# the manuscript qmd can consume it via a single param. Optionally filter to a
# workload subset (e.g. split a workflow run into coverage-only and atacqc-only
# dirs for the qmd's cov_dir / atac_dir params).
#   Usage: merge_runs.R --out=DIR --dirs=RUN1,RUN2[,...] [--workloads=coverage]
args <- list()
for (a in commandArgs(trailingOnly = TRUE)) {
    kv <- sub("^--", "", a)
    k <- gsub("-", "_", sub("=.*$", "", kv))
    args[[k]] <- sub("^[^=]*=", "", kv)
}
out  <- args$out;  if (is.null(out)) stop("--out required")
dirs <- strsplit(args$dirs, ",")[[1]]; if (!length(dirs)) stop("--dirs required")
wls  <- if (!is.null(args$workloads)) strsplit(args$workloads, ",")[[1]] else NULL

dir.create(out, showWarnings = FALSE, recursive = TRUE)
`%||%` <- function(a, b) if (is.null(a)) b else a

rd <- function(d, f) {
    p <- file.path(d, f)
    if (!file.exists(p)) return(NULL)
    tryCatch(read.csv(p, stringsAsFactors = FALSE), error = function(e) NULL)
}
rbind_fill <- function(lst) {
    lst <- Filter(Negate(is.null), lst)
    if (!length(lst)) return(NULL)
    cols <- unique(unlist(lapply(lst, names)))
    do.call(rbind, lapply(lst, function(d) {
        for (m in setdiff(cols, names(d))) d[[m]] <- NA
        d[, cols]
    }))
}

for (f in c("summary.csv", "summary_qmd.csv", "iterations.csv", "files.csv",
            "correctness_workflow.csv", "correctness_preflight.csv",
            "reference_counts.csv")) {
    m <- rbind_fill(lapply(dirs, rd, f))
    if (is.null(m)) next
    if (!is.null(wls) && "workload" %in% names(m)) m <- m[m$workload %in% wls, ]
    if (identical(f, "files.csv")) m <- m[!duplicated(m$file %||% m$basename), ]
    write.csv(m, file.path(out, f), row.names = FALSE)
}

# carry singletons from the first dir that has them
for (f in c("host_info.csv", "pkg_versions.csv", "atacseqqc_anchor.csv",
            "config.txt", "sessionInfo.txt")) {
    for (d in dirs) {
        p <- file.path(d, f)
        if (file.exists(p)) { file.copy(p, file.path(out, f), overwrite = TRUE); break }
    }
}
cat("merged", length(dirs), "run dirs into", out,
    if (!is.null(wls)) paste0(" (workloads: ", paste(wls, collapse = ","), ")") else "", "\n")
