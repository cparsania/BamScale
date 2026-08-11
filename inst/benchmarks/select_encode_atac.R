#!/usr/bin/env Rscript
# Select the FINAL manuscript benchmark dataset from the ENCODE portal:
# a homogeneous set of released GRCh38 paired-end ATAC-seq alignment BAMs.
#   - multi set : --n-multi files in [--multi-min-gb, --multi-max-gb], ONE file
#                 per experiment (replicates of the same experiment are near-
#                 duplicates), largest first, excluding --exclude accessions and
#                 the single file's experiment.
#   - single    : in [--single-min-gb, --single-max-gb], the file whose QC
#                 total_reads (AtacAlignmentQualityMetric) is closest ABOVE
#                 --single-target-reads. Picked by real read count, not size.
# Writes <out>/encode_atac_manifest.csv: the download input and, verbatim, the
# manuscript supplementary table (accession-citable, md5-verifiable).
#
# Usage:
#   Rscript select_encode_atac.R --out=/path/dataset_dir \
#     [--n-multi=24] [--multi-min-gb=2.5] [--multi-max-gb=4.5] \
#     [--single-min-gb=5] [--single-max-gb=9] [--single-target-reads=150000000] \
#     [--exclude=ENCFF646NWY,ENCFF415FEC]
suppressWarnings(suppressMessages({
    ok <- requireNamespace("jsonlite", quietly = TRUE)
    if (!ok) stop("jsonlite is required")
}))

parse_args <- function() {
    out <- list()
    for (a in commandArgs(trailingOnly = TRUE)) {
        if (!startsWith(a, "--")) next
        kv <- sub("^--", "", a)
        k <- gsub("-", "_", sub("=.*$", "", kv))
        v <- if (grepl("=", kv, fixed = TRUE)) sub("^[^=]*=", "", kv) else "TRUE"
        out[[k]] <- v
    }
    out
}
args <- parse_args()
gnum <- function(k, d) if (!is.null(args[[k]])) as.numeric(args[[k]]) else d
gchr <- function(k, d) if (!is.null(args[[k]])) args[[k]] else d

out_dir       <- gchr("out", stop("--out is required"))
n_multi       <- as.integer(gnum("n_multi", 24))
multi_min     <- gnum("multi_min_gb", 2.5) * 1e9
multi_max     <- gnum("multi_max_gb", 4.5) * 1e9
single_min    <- gnum("single_min_gb", 5) * 1e9
single_max    <- gnum("single_max_gb", 9) * 1e9
target_reads  <- gnum("single_target_reads", 1.5e8)
excl          <- strsplit(gchr("exclude", "ENCFF646NWY,ENCFF415FEC"), ",")[[1]]

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- query the ENCODE portal (File objects, embedded QC) --------------------
base <- paste0(
    "https://www.encodeproject.org/search/?type=File",
    "&file_format=bam&output_type=alignments&assembly=GRCh38&status=released",
    "&assay_term_name=ATAC-seq&mapped_run_type=paired-ended",
    "&format=json&limit=all&sort=-file_size",
    "&field=accession&field=file_size&field=md5sum&field=dataset",
    "&field=biosample_ontology.term_name&field=lab.title",
    "&field=quality_metrics&field=date_created"
)
message("Querying ENCODE portal ...")
js <- jsonlite::fromJSON(base, simplifyVector = FALSE)
hits <- js[["@graph"]]
message(sprintf("candidates: %d", length(hits)))
if (length(hits) < n_multi + 1L) stop("Too few candidates returned")

field <- function(h, ...) {
    v <- h
    for (k in c(...)) {
        if (is.null(v)) return(NA_character_)
        v <- v[[k]]
    }
    if (is.null(v)) NA_character_ else as.character(v)
}
qc_total_reads <- function(h) {
    qms <- h$quality_metrics
    if (is.null(qms)) return(NA_real_)
    for (qm in qms) {
        tr <- qm$total_reads
        if (!is.null(tr)) return(as.numeric(tr))
    }
    NA_real_
}

df <- do.call(rbind, lapply(hits, function(h) {
    data.frame(
        accession   = field(h, "accession"),
        experiment  = sub("^/experiments/([^/]+)/$", "\\1", field(h, "dataset")),
        biosample   = field(h, "biosample_ontology", "term_name"),
        lab         = field(h, "lab", "title"),
        file_size   = as.numeric(field(h, "file_size")),
        md5sum      = field(h, "md5sum"),
        total_reads = qc_total_reads(h),
        date_created = field(h, "date_created"),
        stringsAsFactors = FALSE
    )
}))
df <- df[!is.na(df$accession) & !is.na(df$file_size), ]
df <- df[!df$accession %in% excl, ]

# ---- single: read-count-closest-above-target within the size band -----------
sg <- df[df$file_size >= single_min & df$file_size <= single_max & !is.na(df$total_reads), ]
if (nrow(sg) == 0L) stop("No single-file candidates with QC total_reads in band")
above <- sg[sg$total_reads >= target_reads, ]
single <- if (nrow(above) > 0L) above[which.min(above$total_reads), ] else sg[which.max(sg$total_reads), ]
message(sprintf("single: %s (%.1f GB, %.0fM reads, %s, %s)",
                single$accession, single$file_size / 1e9, single$total_reads / 1e6,
                single$biosample, single$experiment))

# ---- multi: size band, one per experiment, largest first --------------------
mg <- df[df$file_size >= multi_min & df$file_size <= multi_max, ]
mg <- mg[mg$experiment != single$experiment, ]
mg <- mg[order(-mg$file_size), ]
mg <- mg[!duplicated(mg$experiment), ]
if (nrow(mg) < n_multi) stop(sprintf("Only %d multi candidates after dedupe", nrow(mg)))
multi <- mg[seq_len(n_multi), ]
message(sprintf("multi: %d files, %.1f GB total, %d distinct experiments, %d biosamples",
                nrow(multi), sum(multi$file_size) / 1e9,
                length(unique(multi$experiment)), length(unique(multi$biosample))))

manifest <- rbind(
    cbind(role = "single", single),
    cbind(role = "multi", multi)
)
manifest$url <- sprintf("https://www.encodeproject.org/files/%s/@@download/%s.bam",
                        manifest$accession, manifest$accession)
out_csv <- file.path(out_dir, "encode_atac_manifest.csv")
write.csv(manifest, out_csv, row.names = FALSE)
message("manifest: ", out_csv)
message(sprintf("TOTAL download: %.1f GB", sum(manifest$file_size) / 1e9))
