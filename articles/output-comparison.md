# BamScale vs Rsamtools Output Reproducibility

## Purpose

This vignette checks reproducibility by comparing **BamScale** and
**Rsamtools** outputs on the same BAM input.

- `step1` fields: `qname, flag, rname, pos, mapq, cigar`
- `seq+qual` fields: `qname, seq, qual`

To keep the comparison deterministic and fair, both tools are run in
single-thread mode for this check.

## Input BAM

``` r
bam <- ompBAM::example_BAM("Unsorted")

cat("BAM:", bam, "\n")
#> BAM: /home/runner/work/_temp/Library/ompBAM/extdata/THP1_ND_1.bam
cat("Exists:", file.exists(bam), "\n")
#> Exists: TRUE
cat("Size (MB):", round(file.size(bam) / 1024^2, 3), "\n")
#> Size (MB): 1.323
```

## Helper Functions

``` r
scanbam_to_df <- function(scan_rec, fields) {
  cols <- vector("list", length(fields))
  names(cols) <- fields

  lens <- integer(length(fields))

  for (i in seq_along(fields)) {
    nm <- fields[[i]]
    v <- scan_rec[[nm]]

    if (is.null(v)) {
      next
    }

    if (nm %in% c("seq", "qual", "rname")) {
      v <- as.character(v)
    }

    cols[[i]] <- v
    lens[[i]] <- length(v)
  }

  n <- if (any(lens > 0L)) max(lens) else 0L

  for (i in seq_along(fields)) {
    if (is.null(cols[[i]])) {
      cols[[i]] <- rep(NA, n)
    } else if (length(cols[[i]]) != n) {
      if (length(cols[[i]]) == 1L && n > 1L) {
        cols[[i]] <- rep(cols[[i]], n)
      } else {
        stop("Inconsistent field lengths in Rsamtools output for field: ", fields[[i]], call. = FALSE)
      }
    }
  }

  as.data.frame(cols, stringsAsFactors = FALSE)
}

normalize_df <- function(df) {
  out <- as.data.frame(df, stringsAsFactors = FALSE)

  for (nm in names(out)) {
    if (is.factor(out[[nm]])) out[[nm]] <- as.character(out[[nm]])
    if (is.integer(out[[nm]]) || is.numeric(out[[nm]])) {
      suppressWarnings(out[[nm]] <- as.integer(out[[nm]]))
    } else {
      out[[nm]] <- as.character(out[[nm]])
    }
  }

  out
}

compare_fields <- function(a, b, fields) {
  stopifnot(all(fields %in% names(a)), all(fields %in% names(b)))

  n <- min(nrow(a), nrow(b))
  if (n == 0L) {
    return(data.frame(
      field = fields,
      identical = NA,
      stringsAsFactors = FALSE
    ))
  }

  aa <- a[seq_len(n), fields, drop = FALSE]
  bb <- b[seq_len(n), fields, drop = FALSE]

  data.frame(
    field = fields,
    identical = vapply(fields, function(f) identical(aa[[f]], bb[[f]]), logical(1)),
    stringsAsFactors = FALSE
  )
}


preview_pair <- function(bs_df, rs_df, fields, n = 5L) {
  n <- min(as.integer(n), nrow(bs_df), nrow(rs_df))
  if (n <= 0L) {
    return(data.frame())
  }

  bs_preview <- bs_df[seq_len(n), fields, drop = FALSE]
  rs_preview <- rs_df[seq_len(n), fields, drop = FALSE]

  bs_preview$tool <- "BamScale"
  rs_preview$tool <- "Rsamtools"
  bs_preview$row_id <- seq_len(n)
  rs_preview$row_id <- seq_len(n)

  out <- rbind(bs_preview, rs_preview)
  out <- out[, c("row_id", "tool", fields), drop = FALSE]
  out
}

kable_fit <- function(x) {
  tbl <- knitr::kable(
    x,
    format = "html",
    table.attr = 'class="repro-table"'
  )
  knitr::asis_output(paste0('<div class="repro-scroll">', tbl, '</div>'))
}
```

## Step1 Output Comparison

``` r
fields_step1 <- c("qname", "flag", "rname", "pos", "mapq", "cigar")

bs_step1 <- BamScale::bam_read(
  file = bam,
  what = fields_step1,
  as = "data.frame",
  include_unmapped = TRUE,
  threads = 1L,
  auto_threads = FALSE
)

rs_step1 <- Rsamtools::scanBam(
  bam,
  param = Rsamtools::ScanBamParam(what = fields_step1)
)[[1]]

bs_step1_df <- normalize_df(bs_step1)
rs_step1_df <- normalize_df(scanbam_to_df(rs_step1, fields_step1))

step1_cmp <- compare_fields(bs_step1_df, rs_step1_df, fields_step1)

summary_step1 <- data.frame(
  rows_bamscale = nrow(bs_step1_df),
  rows_rsamtools = nrow(rs_step1_df),
  rows_compared = min(nrow(bs_step1_df), nrow(rs_step1_df)),
  all_fields_identical = all(step1_cmp$identical),
  stringsAsFactors = FALSE
)

kable_fit(summary_step1)
```

| rows_bamscale | rows_rsamtools | rows_compared | all_fields_identical |
|--------------:|---------------:|--------------:|:---------------------|
|         10000 |          10000 |         10000 | TRUE                 |

``` r
kable_fit(step1_cmp)
```

|       | field | identical |
|:------|:------|:----------|
| qname | qname | TRUE      |
| flag  | flag  | TRUE      |
| rname | rname | TRUE      |
| pos   | pos   | TRUE      |
| mapq  | mapq  | TRUE      |
| cigar | cigar | TRUE      |

``` r

cat("\nQuick output preview (first matched rows):\n")
#> 
#> Quick output preview (first matched rows):
kable_fit(preview_pair(bs_step1_df, rs_step1_df, fields_step1, n = 5L))
```

| row_id | tool      | qname                                     | flag | rname |      pos | mapq | cigar         |
|-------:|:----------|:------------------------------------------|-----:|:------|---------:|-----:|:--------------|
|      1 | BamScale  | ST-E00600:137:H77Y3CCXY:1:1101:6837:1309  |  163 | 19    |   572614 |  255 | 1S88M6798N61M |
|      2 | BamScale  | ST-E00600:137:H77Y3CCXY:1:1101:6837:1309  |   83 | 19    |   579499 |  255 | 1S148M1S      |
|      3 | BamScale  | ST-E00600:137:H77Y3CCXY:1:1101:10450:1309 |  163 | 6     | 44252112 |  255 | 60S86M4S      |
|      4 | BamScale  | ST-E00600:137:H77Y3CCXY:1:1101:10450:1309 |   83 | 6     | 44252124 |  255 | 1S146M3S      |
|      5 | BamScale  | ST-E00600:137:H77Y3CCXY:1:1101:15077:1309 |   99 | 12    | 46185884 |  255 | 1S149M        |
|      1 | Rsamtools | ST-E00600:137:H77Y3CCXY:1:1101:6837:1309  |  163 | 19    |   572614 |  255 | 1S88M6798N61M |
|      2 | Rsamtools | ST-E00600:137:H77Y3CCXY:1:1101:6837:1309  |   83 | 19    |   579499 |  255 | 1S148M1S      |
|      3 | Rsamtools | ST-E00600:137:H77Y3CCXY:1:1101:10450:1309 |  163 | 6     | 44252112 |  255 | 60S86M4S      |
|      4 | Rsamtools | ST-E00600:137:H77Y3CCXY:1:1101:10450:1309 |   83 | 6     | 44252124 |  255 | 1S146M3S      |
|      5 | Rsamtools | ST-E00600:137:H77Y3CCXY:1:1101:15077:1309 |   99 | 12    | 46185884 |  255 | 1S149M        |

## Seq + Qual Output Comparison

``` r
fields_seqqual <- c("qname", "seq", "qual")

bs_sq <- BamScale::bam_read(
  file = bam,
  what = fields_seqqual,
  as = "data.frame",
  seqqual_mode = "compatible",
  include_unmapped = TRUE,
  threads = 1L,
  auto_threads = FALSE
)

rs_sq <- Rsamtools::scanBam(
  bam,
  param = Rsamtools::ScanBamParam(what = fields_seqqual)
)[[1]]

bs_sq_df <- normalize_df(bs_sq)
rs_sq_df <- normalize_df(scanbam_to_df(rs_sq, fields_seqqual))

seqqual_cmp <- compare_fields(bs_sq_df, rs_sq_df, fields_seqqual)

summary_sq <- data.frame(
  rows_bamscale = nrow(bs_sq_df),
  rows_rsamtools = nrow(rs_sq_df),
  rows_compared = min(nrow(bs_sq_df), nrow(rs_sq_df)),
  all_fields_identical = all(seqqual_cmp$identical),
  stringsAsFactors = FALSE
)

kable_fit(summary_sq)
```

| rows_bamscale | rows_rsamtools | rows_compared | all_fields_identical |
|--------------:|---------------:|--------------:|:---------------------|
|         10000 |          10000 |         10000 | TRUE                 |

``` r
kable_fit(seqqual_cmp)
```

|       | field | identical |
|:------|:------|:----------|
| qname | qname | TRUE      |
| seq   | seq   | TRUE      |
| qual  | qual  | TRUE      |

``` r

cat("\nQuick output preview (first matched rows):\n")
#> 
#> Quick output preview (first matched rows):
kable_fit(preview_pair(
  bs_sq_df,
  rs_sq_df,
  fields_seqqual,
  n = 5L
))
```

| row_id | tool      | qname                                     | seq                                                                                                                                                    | qual                                                                                                                                                                |
|-------:|:----------|:------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|      1 | BamScale  | ST-E00600:137:H77Y3CCXY:1:1101:6837:1309  | NGACCGGCGAGGAATAGGAATCATGGCGGCTGCGCTGTTCGTGCTGCTGGGATTCGCGCTGCTGGGCACCCACGGAGCCTCCGGGGCTGCCGGCACAGTCTTCACTACCGTAGAAGACCTTGGCTCGAAGATACTCCTCACCTGCTCCTT | \#AA-AAFJFJJJJJJJJJAAJJJJJJJJJJJJJF\<-FJF-AAFJJJJFFJJFFF\<FFJ-\<\<AFAFF7-AFJFJJJJJ77FJAAJJ7AA777AJJAJFJ-A\<AA\<AAF7FFA-FAFA\<-77\<FJA–)77A7AA7–7AF\<)7777\<))–      |
|      2 | BamScale  | ST-E00600:137:H77Y3CCXY:1:1101:6837:1309  | TGCCGGCACAGTCTTCACTACCGTAGAAGACCTTGGCTCCAAGATACTCCTCACCTGCTCCTTGAATGACAGCGCCACAGAGGTCACAGGGCACCGCTGGCTGAAGGCGGGCGTGGTGCTGAAGGAGGACGCGCTGCCAGGCCAGAAAAN | -FAAA)-AJFF7JJFFAF-FJAF-\<A-7–JJJFA-F-A7-\<-F7-FFJFFJAFJFJJAJFJ\<FFJJJFAF7JJFFJFFJ7\<-\<AA-\<AF-JJA\<7\<777A7\<JA77-7JJJJJFJJJFJJJFFJJFJJJJJJJJJF7-JJJF-AFAAA#      |
|      3 | BamScale  | ST-E00600:137:H77Y3CCXY:1:1101:10450:1309 | NGGAGCTTCGAGGTGGTATATATGACCGAGCCCATTAATGAGTACTGTGTGCAGCAGCGTAAGCAATTTGATGGGAAGAGCGTGGTCTCAGTTAGCAAGGAGGGTCTGGAGCTGCCTGAGGATGAGGAGGAGCTGAAGAAGATGGAGGAG | \#AA-AFFJFJJAAFFJAJJJ\<AJJJFJJJJJA\<\<F7-F-\<-\<FFAFFJFAJJF\<\<FJ—-\<\<-7\<J-FF-AJJJF\<J\<JF-AFAA7-7AJJJF-7777A-7AA—777—7FJ–AFA7-7AAJFAA)7))7—A\<7\<–7—A7           |
|      4 | BamScale  | ST-E00600:137:H77Y3CCXY:1:1101:10450:1309 | TGGGAAGAGCCTGGTCTCAGTTACCAAGGAGGGTCTGGAGCTGCCTGAGGATGAGTAGGAGAAGAAGAAGATGGAAGAGAGCAAGGTAAAGTTTGAGAACCTCTGCAAGCTCATGAAAGAAATCTTAGATAAGAAGGTCGAGATGGTGAN | AJA-FFJJJFF\<7FF-AA\<7F\<AF\<-7FJFJJJJJFF-\<7-JFJJJJFFAJFJ\<A-JF-JJJAFFFFJJJJFJ7\<–F-A-\<7JF\<-FF\<\<—\<-77\<-A-FJJFJJFJFJJJJJJJJJJJ\<JJJJJJJJJJJJJJ\<7JJJF\<FFAAA# |
|      5 | BamScale  | ST-E00600:137:H77Y3CCXY:1:1101:15077:1309 | NGCATCTCTACCCCTACTGTCCAGTAGGTGGGATGTGGCTGGGCTGGACAGTCCAGATTATACGAACTGGCAACTCTGAACAAACACCCTCCCTGGAAACAATATTATATTTGATGGTTAGATTCTTTAGCAAACCTATTACATTATTCG | \#AAFFF-FJJJA-AFJJJJJJJJJJJJJJJJJJJJJJJJFJF7-7FJJJFJJJ\<\<F-FF-FJJFJFJ-F\<JA7-7-7AF-JJJFJJJFJJJ7FJJJ\<FJJJJJJJJJJJJFJJJJJ7FFJ-7-7-A-AFFF-F–7F7AJJJJJ-\<\<FJJ        |
|      1 | Rsamtools | ST-E00600:137:H77Y3CCXY:1:1101:6837:1309  | NGACCGGCGAGGAATAGGAATCATGGCGGCTGCGCTGTTCGTGCTGCTGGGATTCGCGCTGCTGGGCACCCACGGAGCCTCCGGGGCTGCCGGCACAGTCTTCACTACCGTAGAAGACCTTGGCTCGAAGATACTCCTCACCTGCTCCTT | \#AA-AAFJFJJJJJJJJJAAJJJJJJJJJJJJJF\<-FJF-AAFJJJJFFJJFFF\<FFJ-\<\<AFAFF7-AFJFJJJJJ77FJAAJJ7AA777AJJAJFJ-A\<AA\<AAF7FFA-FAFA\<-77\<FJA–)77A7AA7–7AF\<)7777\<))–      |
|      2 | Rsamtools | ST-E00600:137:H77Y3CCXY:1:1101:6837:1309  | TGCCGGCACAGTCTTCACTACCGTAGAAGACCTTGGCTCCAAGATACTCCTCACCTGCTCCTTGAATGACAGCGCCACAGAGGTCACAGGGCACCGCTGGCTGAAGGCGGGCGTGGTGCTGAAGGAGGACGCGCTGCCAGGCCAGAAAAN | -FAAA)-AJFF7JJFFAF-FJAF-\<A-7–JJJFA-F-A7-\<-F7-FFJFFJAFJFJJAJFJ\<FFJJJFAF7JJFFJFFJ7\<-\<AA-\<AF-JJA\<7\<777A7\<JA77-7JJJJJFJJJFJJJFFJJFJJJJJJJJJF7-JJJF-AFAAA#      |
|      3 | Rsamtools | ST-E00600:137:H77Y3CCXY:1:1101:10450:1309 | NGGAGCTTCGAGGTGGTATATATGACCGAGCCCATTAATGAGTACTGTGTGCAGCAGCGTAAGCAATTTGATGGGAAGAGCGTGGTCTCAGTTAGCAAGGAGGGTCTGGAGCTGCCTGAGGATGAGGAGGAGCTGAAGAAGATGGAGGAG | \#AA-AFFJFJJAAFFJAJJJ\<AJJJFJJJJJA\<\<F7-F-\<-\<FFAFFJFAJJF\<\<FJ—-\<\<-7\<J-FF-AJJJF\<J\<JF-AFAA7-7AJJJF-7777A-7AA—777—7FJ–AFA7-7AAJFAA)7))7—A\<7\<–7—A7           |
|      4 | Rsamtools | ST-E00600:137:H77Y3CCXY:1:1101:10450:1309 | TGGGAAGAGCCTGGTCTCAGTTACCAAGGAGGGTCTGGAGCTGCCTGAGGATGAGTAGGAGAAGAAGAAGATGGAAGAGAGCAAGGTAAAGTTTGAGAACCTCTGCAAGCTCATGAAAGAAATCTTAGATAAGAAGGTCGAGATGGTGAN | AJA-FFJJJFF\<7FF-AA\<7F\<AF\<-7FJFJJJJJFF-\<7-JFJJJJFFAJFJ\<A-JF-JJJAFFFFJJJJFJ7\<–F-A-\<7JF\<-FF\<\<—\<-77\<-A-FJJFJJFJFJJJJJJJJJJJ\<JJJJJJJJJJJJJJ\<7JJJF\<FFAAA# |
|      5 | Rsamtools | ST-E00600:137:H77Y3CCXY:1:1101:15077:1309 | NGCATCTCTACCCCTACTGTCCAGTAGGTGGGATGTGGCTGGGCTGGACAGTCCAGATTATACGAACTGGCAACTCTGAACAAACACCCTCCCTGGAAACAATATTATATTTGATGGTTAGATTCTTTAGCAAACCTATTACATTATTCG | \#AAFFF-FJJJA-AFJJJJJJJJJJJJJJJJJJJJJJJJFJF7-7FJJJFJJJ\<\<F-FF-FJJFJFJ-F\<JA7-7-7AF-JJJFJJJFJJJ7FJJJ\<FJJJJJJJJJJJJFJJJJJ7FFJ-7-7-A-AFFF-F–7F7AJJJJJ-\<\<FJJ        |

## Reproducibility Summary

``` r
repro_tbl <- rbind(
  data.frame(
    workload = "step1",
    rows_bamscale = nrow(bs_step1_df),
    rows_rsamtools = nrow(rs_step1_df),
    rows_compared = min(nrow(bs_step1_df), nrow(rs_step1_df)),
    all_fields_identical = all(step1_cmp$identical),
    stringsAsFactors = FALSE
  ),
  data.frame(
    workload = "seqqual",
    rows_bamscale = nrow(bs_sq_df),
    rows_rsamtools = nrow(rs_sq_df),
    rows_compared = min(nrow(bs_sq_df), nrow(rs_sq_df)),
    all_fields_identical = all(seqqual_cmp$identical),
    stringsAsFactors = FALSE
  )
)

kable_fit(repro_tbl)
```

| workload | rows_bamscale | rows_rsamtools | rows_compared | all_fields_identical |
|:---------|--------------:|---------------:|--------------:|:---------------------|
| step1    |         10000 |          10000 |         10000 | TRUE                 |
| seqqual  |         10000 |          10000 |         10000 | TRUE                 |

## Session Information

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] crayon_1.5.3         cli_3.6.5            knitr_1.51          
#>  [4] rlang_1.1.7          xfun_0.56            generics_0.1.4      
#>  [7] textshaping_1.0.4    jsonlite_2.0.0       S4Vectors_0.48.0    
#> [10] Biostrings_2.78.0    BamScale_0.1.0       htmltools_0.5.9     
#> [13] stats4_4.5.2         ragg_1.5.0           sass_0.4.10         
#> [16] rmarkdown_2.30       Seqinfo_1.0.0        evaluate_1.0.5      
#> [19] jquerylib_0.1.4      bitops_1.0-9         fastmap_1.2.0       
#> [22] IRanges_2.44.0       yaml_2.3.12          lifecycle_1.0.5     
#> [25] compiler_4.5.2       codetools_0.2-20     fs_1.6.6            
#> [28] Rcpp_1.1.1           XVector_0.50.0       BiocParallel_1.44.0 
#> [31] systemfonts_1.3.1    digest_0.6.39        R6_2.6.1            
#> [34] parallel_4.5.2       GenomicRanges_1.62.1 bslib_0.10.0        
#> [37] tools_4.5.2          ompBAM_1.14.0        Rsamtools_2.26.0    
#> [40] pkgdown_2.2.0        BiocGenerics_0.56.0  cachem_1.1.0        
#> [43] desc_1.4.3
```

## Notes

- Equality here is reported at the field level over rows compared in
  order.
- If row counts differ, inspect filtering assumptions (`flag`,
  `include_unmapped`, `which`, etc.) first.
- For performance benchmarking, keep this output check separate from
  timing runs.
