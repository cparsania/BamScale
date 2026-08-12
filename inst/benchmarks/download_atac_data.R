#!/usr/bin/env Rscript
# download_atac_data.R
#
# Fetch a real, public paired-end ATAC-seq dataset for the atacqc workflow
# benchmark, index each BAM, and verify the reads are paired-end (required for
# the fragment-size distribution). There is no ATAC-seq data on this host, so
# this must be run once before `--include-atacqc=true`.
#
# Three ways to specify the data (in precedence order):
#   1. --urls=URL1,URL2,...             direct .bam download URLs
#   2. --experiment=ENCSRxxxxxxx        an ENCODE experiment; its alignment BAMs
#                                       are resolved via the ENCODE REST API
#   3. --files=ENCFFxxxxxxx,...         ENCODE file accessions (built into URLs)
#
# ENCODE reference: an experiment page lists files with output_type
# "alignments" and file_format "bam"; each file downloads from
#   https://www.encodeproject.org/files/<ACC>/@@download/<ACC>.bam
# Browse ATAC-seq experiments at https://www.encodeproject.org/search/?type=Experiment&assay_title=ATAC-seq
#
# Usage:
#   Rscript inst/benchmarks/download_atac_data.R \
#     --experiment=ENCSRxxxxxxx --assembly=GRCh38 --output-type=alignments \
#     --dest=/spectrum/GSCT/chiragp/018_BamScale/atac_bams --max-files=8

suppressWarnings(suppressMessages(ok <- requireNamespace("Rsamtools", quietly = TRUE)))
if (!ok) stop("Rsamtools is required.")

# download.file defaults to a 60s timeout, which truncates multi-GB BAMs. Raise
# it and use libcurl (follows the ENCODE -> signed-S3 redirect).
options(timeout = max(7200, getOption("timeout")))

parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
    out <- list()
    for (a in argv) {
        if (!grepl("^--", a)) next
        a <- sub("^--", "", a)
        k <- if (grepl("=", a)) sub("=.*$", "", a) else a
        v <- if (grepl("=", a)) sub("^[^=]*=", "", a) else "true"
        out[[gsub("-", "_", k)]] <- v
    }
    out
}
split_csv <- function(x) if (is.null(x)) character() else {
    v <- trimws(strsplit(x, ",", fixed = TRUE)[[1]]); v[nzchar(v)]
}

args <- parse_args()
dest       <- if (!is.null(args$dest)) args$dest else "atac_bams"
assembly   <- if (!is.null(args$assembly)) args$assembly else "GRCh38"
out_type   <- if (!is.null(args$output_type)) args$output_type else "alignments"
max_files  <- if (!is.null(args$max_files)) as.integer(args$max_files) else 8L
do_index   <- !(tolower(as.character(args$index)) %in% c("false", "0", "no"))

encode_file_url <- function(acc) sprintf("https://www.encodeproject.org/files/%s/@@download/%s.bam", acc, acc)

# --manifest=manifest.csv mode: rows with columns accession,url,md5sum (as
# written by select_encode_atac.R). Optional --role=single|multi filters rows.
# md5 is verified after download; a mismatch deletes the file and reports it.
manifest_df <- NULL
if (!is.null(args$manifest)) {
    manifest_df <- utils::read.csv(args$manifest, stringsAsFactors = FALSE)
    if (!is.null(args$role)) manifest_df <- manifest_df[manifest_df$role == args$role, ]
    if (!is.null(args$rows)) {
        idx <- suppressWarnings(as.integer(split_csv(args$rows)))
        manifest_df <- manifest_df[idx[!is.na(idx) & idx <= nrow(manifest_df)], ]
    }
}

resolve_urls <- function() {
    if (!is.null(manifest_df)) return(manifest_df$url)
    if (!is.null(args$urls)) return(split_csv(args$urls))
    if (!is.null(args$files)) return(vapply(split_csv(args$files), encode_file_url, character(1)))
    if (!is.null(args$experiment)) {
        if (!requireNamespace("jsonlite", quietly = TRUE)) {
            stop("jsonlite needed to resolve --experiment; either install it or pass --urls / --files.", call. = FALSE)
        }
        url <- sprintf("https://www.encodeproject.org/experiments/%s/?format=json", args$experiment)
        js <- jsonlite::fromJSON(url, simplifyVector = FALSE)
        files <- js$files
        keep <- Filter(function(f) {
            isTRUE(identical(f$file_format, "bam")) &&
                isTRUE(identical(f$output_type, out_type)) &&
                (is.null(f$assembly) || identical(f$assembly, assembly)) &&
                (is.null(f$status) || identical(f$status, "released"))
        }, files)
        if (!length(keep)) stop("No matching BAM files for experiment ", args$experiment,
                                " (output_type=", out_type, ", assembly=", assembly, ").", call. = FALSE)
        accs <- vapply(keep, function(f) f$accession, character(1))
        vapply(accs, encode_file_url, character(1))
    } else {
        stop("Specify one of --urls, --files, or --experiment. See header for ENCODE guidance.", call. = FALSE)
    }
}

urls <- head(unique(resolve_urls()), max_files)
if (!length(urls)) stop("No download URLs resolved.")
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("Resolved %d BAM URL(s); downloading to %s\n", length(urls), normalizePath(dest, mustWork = FALSE)))

is_paired <- function(bam) {
    fl <- tryCatch(
        Rsamtools::scanBam(Rsamtools::BamFile(bam, yieldSize = 2000L),
                           param = Rsamtools::ScanBamParam(what = "flag"))[[1]]$flag,
        error = function(e) integer())
    if (!length(fl)) return(NA)
    mean(bitwAnd(fl, 0x1L) != 0L) > 0.5
}

results <- list()
for (u in urls) {
    acc <- sub("\\.bam$", "", basename(u))
    dst <- file.path(dest, paste0(acc, ".bam"))
    if (!file.exists(dst)) {
        cat("  downloading ", basename(u), " ... ", sep = "")
        rc <- tryCatch({ utils::download.file(u, dst, mode = "wb", quiet = TRUE, method = "libcurl"); 0L },
                       error = function(e) { cat("FAILED (", conditionMessage(e), ")\n", sep = ""); 1L })
        if (rc != 0L) next
        cat("ok\n")
    } else {
        cat("  exists: ", basename(dst), "\n", sep = "")
    }
    # md5 verification against the manifest (delete on mismatch so a re-run retries)
    md5_ok <- NA
    if (!is.null(manifest_df)) {
        want <- manifest_df$md5sum[match(acc, manifest_df$accession)]
        if (length(want) == 1L && !is.na(want)) {
            got <- unname(tools::md5sum(dst))
            md5_ok <- identical(got, want)
            if (!isTRUE(md5_ok)) {
                cat("    MD5 MISMATCH (", got, " != ", want, ") -- deleting for retry\n", sep = "")
                unlink(dst)
                next
            }
        }
    }
    if (do_index) tryCatch(Rsamtools::indexBam(dst), error = function(e) cat("    index failed: ", conditionMessage(e), "\n", sep = ""))
    results[[dst]] <- list(size_mb = round(file.size(dst) / 1024^2, 1), paired = is_paired(dst), md5_ok = md5_ok)
}

cat("\nSummary:\n")
for (p in names(results)) {
    r <- results[[p]]
    cat(sprintf("  %-40s  %8.1f MB   paired=%s   md5_ok=%s\n",
                basename(p), r$size_mb, r$paired, r$md5_ok))
}
paired_all <- all(vapply(results, function(r) isTRUE(r$paired), logical(1)))
if (!paired_all) {
    cat("\nWARNING: some files are not paired-end. ATAC fragment-size QC requires paired-end reads.\n")
} else {
    cat(sprintf("\nAll %d files are paired-end. Use: --include-atacqc=true --atac-bam-dir=%s\n",
                length(results), normalizePath(dest, mustWork = FALSE)))
}
