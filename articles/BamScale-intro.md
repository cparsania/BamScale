# Getting started with BamScale

## Introduction

`BamScale` provides multithreaded sequential BAM traversal for
R/Bioconductor workflows. The package is designed for users who already
rely on `Rsamtools`, `GenomicAlignments`, and `BiocParallel`, but need
faster BAM parsing without moving to a separate non-Bioconductor
workflow model.

The package targets three common classes of BAM access:

- metadata-oriented scans for filtering, fragment summaries, and quality
  control
- generation of alignment objects for downstream `GenomicAlignments`
  workflows
- sequence and quality extraction, including an optimized compact mode

The motivation for inclusion in Bioconductor is therefore
straightforward: `BamScale` accelerates a core data-access bottleneck
while preserving familiar Bioconductor inputs, filtering semantics, and
output types.

## Installation

At present, `BamScale` depends on the `ompBAM` engine and is intended
for use with an OpenMP-capable toolchain.

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("BamScale")
```

## Basic usage

``` r

library(BamScale)

bam <- ompBAM::example_BAM("Unsorted")
```

### Metadata-oriented BAM access

The most common BAM-access pattern is extraction of alignment metadata
such as read name, flag, reference name, position, mapping quality, and
CIGAR string.

``` r

x <- bam_read(
    file = bam,
    what = c("qname", "flag", "rname", "pos", "mapq", "cigar"),
    as = "data.frame",
    threads = 2
)

head(x)
#>                                       qname flag rname      pos mapq
#> 1  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309  163    19   572614  255
#> 2  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309   83    19   579499  255
#> 3 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309  163     6 44252112  255
#> 4 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309   83     6 44252124  255
#> 5 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309   99    12 46185884  255
#> 6 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309  147    12 46185968  255
#>           cigar
#> 1 1S88M6798N61M
#> 2      1S148M1S
#> 3      60S86M4S
#> 4      1S146M3S
#> 5        1S149M
#> 6      7S142M1S
```

### Filtering with `ScanBamParam`

`BamScale` accepts
[`Rsamtools::ScanBamParam`](https://rdrr.io/pkg/Rsamtools/man/ScanBamParam-class.html),
allowing existing filtering logic to be reused directly.

``` r

param <- Rsamtools::ScanBamParam(
    what = c("qname", "flag", "mapq"),
    mapqFilter = 20L,
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)
)

filtered <- bam_read(
    file = bam,
    param = param,
    as = "data.frame",
    threads = 2
)

head(filtered)
#>                                       qname flag mapq
#> 1  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309  163  255
#> 2  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309   83  255
#> 3 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309  163  255
#> 4 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309   83  255
#> 5 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309   99  255
#> 6 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309  147  255
```

### Alignment-object output

When downstream workflows expect `GenomicAlignments` objects,
[`bam_read()`](https://cparsania.github.io/BamScale/reference/bam_read.md)
can return those directly.

``` r

ga <- bam_read(
    file = bam,
    what = c("qname", "flag", "rname", "pos", "cigar", "strand"),
    as = "GAlignments",
    threads = 2
)

ga
#> GAlignments object with 10000 alignments and 2 metadata columns:
#>           seqnames strand              cigar    qwidth     start       end
#>              <Rle>  <Rle>        <character> <integer> <integer> <integer>
#>       [1]       19      +      1S88M6798N61M       150    572614    579560
#>       [2]       19      -           1S148M1S       150    579499    579646
#>       [3]        6      +           60S86M4S       150  44252112  44252197
#>       [4]        6      -           1S146M3S       150  44252124  44252269
#>       [5]       12      +             1S149M       150  46185884  46186032
#>       ...      ...    ...                ...       ...       ...       ...
#>    [9996]        2      -         16S122M12S       150 186680457 186680578
#>    [9997]       15      +           4S143M3S       150  90501177  90501319
#>    [9998]       15      -        26M2I56M66S       150  90501301  90501382
#>    [9999]        4      + 3M336N95M752N51M1S       150 121801488 121802724
#>   [10000]        4      -      1S117M121N32M       150 121802677 121802946
#>               width     njunc |                  qname      flag
#>           <integer> <integer> |            <character> <integer>
#>       [1]      6947         1 | ST-E00600:137:H77Y3C..       163
#>       [2]       148         0 | ST-E00600:137:H77Y3C..        83
#>       [3]        86         0 | ST-E00600:137:H77Y3C..       163
#>       [4]       146         0 | ST-E00600:137:H77Y3C..        83
#>       [5]       149         0 | ST-E00600:137:H77Y3C..        99
#>       ...       ...       ... .                    ...       ...
#>    [9996]       122         0 | ST-E00600:137:H77Y3C..        83
#>    [9997]       143         0 | ST-E00600:137:H77Y3C..       163
#>    [9998]        82         0 | ST-E00600:137:H77Y3C..        83
#>    [9999]      1237         2 | ST-E00600:137:H77Y3C..       163
#>   [10000]       270         1 | ST-E00600:137:H77Y3C..        83
#>   -------
#>   seqinfo: 25 sequences from an unspecified genome; no seqlengths
```

### Sequence and quality extraction

For `seq` and `qual`, `BamScale` supports both a
compatibility-preserving mode and an optimized compact mode.

``` r

seqqual_compatible <- bam_read(
    file = bam,
    what = c("qname", "seq", "qual"),
    as = "data.frame",
    seqqual_mode = "compatible",
    threads = 2
)

head(seqqual_compatible)
#>                                       qname
#> 1  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309
#> 2  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309
#> 3 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309
#> 4 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309
#> 5 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309
#> 6 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309
#>                                                                                                                                                      seq
#> 1 NGACCGGCGAGGAATAGGAATCATGGCGGCTGCGCTGTTCGTGCTGCTGGGATTCGCGCTGCTGGGCACCCACGGAGCCTCCGGGGCTGCCGGCACAGTCTTCACTACCGTAGAAGACCTTGGCTCGAAGATACTCCTCACCTGCTCCTT
#> 2 TGCCGGCACAGTCTTCACTACCGTAGAAGACCTTGGCTCCAAGATACTCCTCACCTGCTCCTTGAATGACAGCGCCACAGAGGTCACAGGGCACCGCTGGCTGAAGGCGGGCGTGGTGCTGAAGGAGGACGCGCTGCCAGGCCAGAAAAN
#> 3 NGGAGCTTCGAGGTGGTATATATGACCGAGCCCATTAATGAGTACTGTGTGCAGCAGCGTAAGCAATTTGATGGGAAGAGCGTGGTCTCAGTTAGCAAGGAGGGTCTGGAGCTGCCTGAGGATGAGGAGGAGCTGAAGAAGATGGAGGAG
#> 4 TGGGAAGAGCCTGGTCTCAGTTACCAAGGAGGGTCTGGAGCTGCCTGAGGATGAGTAGGAGAAGAAGAAGATGGAAGAGAGCAAGGTAAAGTTTGAGAACCTCTGCAAGCTCATGAAAGAAATCTTAGATAAGAAGGTCGAGATGGTGAN
#> 5 NGCATCTCTACCCCTACTGTCCAGTAGGTGGGATGTGGCTGGGCTGGACAGTCCAGATTATACGAACTGGCAACTCTGAACAAACACCCTCCCTGGAAACAATATTATATTTGATGGTTAGATTCTTTAGCAAACCTATTACATTATTCG
#> 6 AACAAAAACCCTCCCTAGCAACAATATTATATTTGATCGTTAGATTCTTTAGCAAACCTATTACCTTATTCGATGTCAGCTAACCTCTTGGTTTGATCATCTTTTCCAGCTGCTCTAGGGGCCTAACCACCTTGGATAACCTGGTCTCTN
#>                                                                                                                                                     qual
#> 1 #AA-AAFJFJJJJJJJJJAAJJJJJJJJJJJJJF<-FJF-AAFJJJJFFJJFFF<FFJ-<<AFAFF7-AFJFJJJJJ77FJAAJJ7AA777AJJAJFJ-A<AA<AAF7FFA-FAFA<-77<FJA--)77A7AA7--7AF<)7777<))--
#> 2 -FAAA)-AJFF7JJFFAF-FJAF-<A-7--JJJFA-F-A7-<-F7-FFJFFJAFJFJJAJFJ<FFJJJFAF7JJFFJFFJ7<-<AA-<AF-JJA<7<777A7<JA77-7JJJJJFJJJFJJJFFJJFJJJJJJJJJF7-JJJF-AFAAA#
#> 3 #AA-AFFJFJJAAFFJAJJJ<AJJJFJJJJJA<<F7-F-<-<FFAFFJFAJJF<<FJ----<<-7<J-FF-AJJJF<J<JF-AFAA7-7AJJJF-7777A-7AA---777---7FJ--AFA7-7AAJFAA)7))7---A<7<--7---A7
#> 4 AJA-FFJJJFF<7FF-AA<7F<AF<-7FJFJJJJJFF-<7-JFJJJJFFAJFJ<A-JF-JJJAFFFFJJJJFJ7<--F-A-<7JF<-FF<<---<-77<-A-FJJFJJFJFJJJJJJJJJJJ<JJJJJJJJJJJJJJ<7JJJF<FFAAA#
#> 5 #AAFFF-FJJJA-AFJJJJJJJJJJJJJJJJJJJJJJJJFJF7-7FJJJFJJJ<<F-FF-FJJFJFJ-F<JA7-7-7AF-JJJFJJJFJJJ7FJJJ<FJJJJJJJJJJJJFJJJJJ7FFJ-7-7-A-AFFF-F--7F7AJJJJJ-<<FJJ
#> 6 7A7-7--<A-)<)AA7-7-A<<7-FF<7<F<7-7A-F-JJAAFJJJJJJA<JF<A7JJJFAA-J<JJFJJA-7<JJJJFFJJAJJA7<<<7<FJJAJJJFJJFJJJA-JJFJJFJJJJJJJJFJJJJJJ<F-JJJJJFFJJFJJFFAAA#
```

Compact mode returns lower-level raw vectors for throughput-oriented
workflows. These values can be decoded back to ordinary strings
explicitly when needed.

``` r

seqqual_compact <- bam_read(
    file = bam,
    what = c("qname", "qwidth", "seq", "qual"),
    as = "data.frame",
    seqqual_mode = "compact",
    threads = 2
)

head(seqqual_compact)
#>                                       qname qwidth          seq         qual
#> 1  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309    150 f4, 12, .... 02, 20, ....
#> 2  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309    150 84, 22, .... 0c, 25, ....
#> 3 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309    150 f4, 41, .... 02, 20, ....
#> 4 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309    150 84, 44, .... 20, 29, ....
#> 5 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309    150 f4, 21, .... 02, 20, ....
#> 6 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309    150 11, 21, .... 16, 20, ....

seqqual_decoded <- decode_seqqual_compact(seqqual_compact)
head(seqqual_decoded)
#>                                       qname qwidth
#> 1  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309    150
#> 2  ST-E00600:137:H77Y3CCXY:1:1101:6837:1309    150
#> 3 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309    150
#> 4 ST-E00600:137:H77Y3CCXY:1:1101:10450:1309    150
#> 5 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309    150
#> 6 ST-E00600:137:H77Y3CCXY:1:1101:15077:1309    150
#>                                                                                                                                                      seq
#> 1 NGACCGGCGAGGAATAGGAATCATGGCGGCTGCGCTGTTCGTGCTGCTGGGATTCGCGCTGCTGGGCACCCACGGAGCCTCCGGGGCTGCCGGCACAGTCTTCACTACCGTAGAAGACCTTGGCTCGAAGATACTCCTCACCTGCTCCTT
#> 2 TGCCGGCACAGTCTTCACTACCGTAGAAGACCTTGGCTCCAAGATACTCCTCACCTGCTCCTTGAATGACAGCGCCACAGAGGTCACAGGGCACCGCTGGCTGAAGGCGGGCGTGGTGCTGAAGGAGGACGCGCTGCCAGGCCAGAAAAN
#> 3 NGGAGCTTCGAGGTGGTATATATGACCGAGCCCATTAATGAGTACTGTGTGCAGCAGCGTAAGCAATTTGATGGGAAGAGCGTGGTCTCAGTTAGCAAGGAGGGTCTGGAGCTGCCTGAGGATGAGGAGGAGCTGAAGAAGATGGAGGAG
#> 4 TGGGAAGAGCCTGGTCTCAGTTACCAAGGAGGGTCTGGAGCTGCCTGAGGATGAGTAGGAGAAGAAGAAGATGGAAGAGAGCAAGGTAAAGTTTGAGAACCTCTGCAAGCTCATGAAAGAAATCTTAGATAAGAAGGTCGAGATGGTGAN
#> 5 NGCATCTCTACCCCTACTGTCCAGTAGGTGGGATGTGGCTGGGCTGGACAGTCCAGATTATACGAACTGGCAACTCTGAACAAACACCCTCCCTGGAAACAATATTATATTTGATGGTTAGATTCTTTAGCAAACCTATTACATTATTCG
#> 6 AACAAAAACCCTCCCTAGCAACAATATTATATTTGATCGTTAGATTCTTTAGCAAACCTATTACCTTATTCGATGTCAGCTAACCTCTTGGTTTGATCATCTTTTCCAGCTGCTCTAGGGGCCTAACCACCTTGGATAACCTGGTCTCTN
#>                                                                                                                                                     qual
#> 1 #AA-AAFJFJJJJJJJJJAAJJJJJJJJJJJJJF<-FJF-AAFJJJJFFJJFFF<FFJ-<<AFAFF7-AFJFJJJJJ77FJAAJJ7AA777AJJAJFJ-A<AA<AAF7FFA-FAFA<-77<FJA--)77A7AA7--7AF<)7777<))--
#> 2 -FAAA)-AJFF7JJFFAF-FJAF-<A-7--JJJFA-F-A7-<-F7-FFJFFJAFJFJJAJFJ<FFJJJFAF7JJFFJFFJ7<-<AA-<AF-JJA<7<777A7<JA77-7JJJJJFJJJFJJJFFJJFJJJJJJJJJF7-JJJF-AFAAA#
#> 3 #AA-AFFJFJJAAFFJAJJJ<AJJJFJJJJJA<<F7-F-<-<FFAFFJFAJJF<<FJ----<<-7<J-FF-AJJJF<J<JF-AFAA7-7AJJJF-7777A-7AA---777---7FJ--AFA7-7AAJFAA)7))7---A<7<--7---A7
#> 4 AJA-FFJJJFF<7FF-AA<7F<AF<-7FJFJJJJJFF-<7-JFJJJJFFAJFJ<A-JF-JJJAFFFFJJJJFJ7<--F-A-<7JF<-FF<<---<-77<-A-FJJFJJFJFJJJJJJJJJJJ<JJJJJJJJJJJJJJ<7JJJF<FFAAA#
#> 5 #AAFFF-FJJJA-AFJJJJJJJJJJJJJJJJJJJJJJJJFJF7-7FJJJFJJJ<<F-FF-FJJFJFJ-F<JA7-7-7AF-JJJFJJJFJJJ7FJJJ<FJJJJJJJJJJJJFJJJJJ7FFJ-7-7-A-AFFF-F--7F7AJJJJJ-<<FJJ
#> 6 7A7-7--<A-)<)AA7-7-A<<7-FF<7<F<7-7A-F-JJAAFJJJJJJA<JF<A7JJJFAA-J<JJFJJA-7<JJJJFFJJAJJA7<<<7<FJJAJJJFJJFJJJA-JJFJJFJJJJJJJJFJJJJJJ<F-JJJJJFFJJFJJFFAAA#
```

### Fast count summaries

[`bam_count()`](https://cparsania.github.io/BamScale/reference/bam_count.md)
provides chromosome-level count summaries using the same BAM filtering
model.

``` r

counts <- bam_count(bam, threads = 2)
counts
#>    seqname seqlength count
#> 1        1 248956422  1308
#> 2       10 133797422   308
#> 3       11 135086622   600
#> 4       12 133275309   648
#> 5       13 114364328   162
#> 6       14 107043718   230
#> 7       15 101991189   294
#> 8       16  90338345   334
#> 9       17  83257441   570
#> 10      18  80373285    78
#> 11      19  58617616   600
#> 12       2 242193529   586
#> 13      20  64444167   168
#> 14      21  46709983    76
#> 15      22  50818468   200
#> 16       3 198295559   454
#> 17       4 190214555   252
#> 18       5 181538259   516
#> 19       6 170805979   824
#> 20       7 159345973   412
#> 21       8 145138636   452
#> 22       9 138394717   372
#> 23      MT     16569   288
#> 24       X 156040895   254
#> 25       Y  57227415    14
#> 26       *        NA     0
```

## Relationship to existing Bioconductor tools

`BamScale` is intended to complement, not replace, existing Bioconductor
packages.

- `Rsamtools` remains the canonical low-level BAM access layer in R and
  defines the filtering idioms that `BamScale` reuses through
  `ScanBamParam`.
- `GenomicAlignments` remains the standard package for alignment-centric
  downstream workflows, and `BamScale` supports direct generation of
  compatible alignment objects.
- `BiocParallel` remains the standard mechanism for file-level
  parallelism, and `BamScale` adds a second axis of parallelism through
  per-file OpenMP threads.

The main difference is therefore performance-oriented: `BamScale`
accelerates the traversal step itself while staying close to existing
Bioconductor usage patterns.

## Session information

``` r

sessionInfo()
#> R version 4.6.0 (2026-04-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
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
#> other attached packages:
#> [1] BamScale_0.99.9  BiocStyle_2.40.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] sass_0.4.10                 generics_0.1.4             
#>  [3] SparseArray_1.12.2          bitops_1.0-9               
#>  [5] lattice_0.22-9              digest_0.6.39              
#>  [7] evaluate_1.0.5              grid_4.6.0                 
#>  [9] bookdown_0.46               fastmap_1.2.0              
#> [11] jsonlite_2.0.0              Matrix_1.7-5               
#> [13] cigarillo_1.2.0             BiocManager_1.30.27        
#> [15] Biostrings_2.80.0           codetools_0.2-20           
#> [17] textshaping_1.0.5           jquerylib_0.1.4            
#> [19] abind_1.4-8                 cli_3.6.6                  
#> [21] rlang_1.2.0                 crayon_1.5.3               
#> [23] XVector_0.52.0              ompBAM_1.16.0              
#> [25] Biobase_2.72.0              cachem_1.1.0               
#> [27] DelayedArray_0.38.1         yaml_2.3.12                
#> [29] S4Arrays_1.12.0             tools_4.6.0                
#> [31] parallel_4.6.0              BiocParallel_1.46.0        
#> [33] Rsamtools_2.28.0            SummarizedExperiment_1.42.0
#> [35] BiocGenerics_0.58.0         R6_2.6.1                   
#> [37] matrixStats_1.5.0           stats4_4.6.0               
#> [39] lifecycle_1.0.5             Seqinfo_1.2.0              
#> [41] S4Vectors_0.50.0            fs_2.1.0                   
#> [43] IRanges_2.46.0              ragg_1.5.2                 
#> [45] desc_1.4.3                  pkgdown_2.2.0              
#> [47] bslib_0.10.0                Rcpp_1.1.1-1.1             
#> [49] systemfonts_1.3.2           xfun_0.57                  
#> [51] GenomicRanges_1.64.0        GenomicAlignments_1.48.0   
#> [53] MatrixGenerics_1.24.0       knitr_1.51                 
#> [55] htmltools_0.5.9             rmarkdown_2.31             
#> [57] compiler_4.6.0
```
