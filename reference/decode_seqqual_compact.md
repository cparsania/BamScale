# Decode compact `seq` and `qual` columns in BamScale output

Convenience wrapper for converting a compact
[`bam_read()`](https://cparsania.github.io/BamScale/reference/bam_read.md)
result back to ordinary sequence and quality strings.

## Usage

``` r
decode_seqqual_compact(
  x,
  seq_col = "seq",
  qual_col = "qual",
  qwidth_col = "qwidth"
)
```

## Arguments

- x:

  A `data.frame`,
  [`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html),
  or list-like object containing compact BamScale `seq` and/or `qual`
  columns.

- seq_col:

  Name of the compact sequence column.

- qual_col:

  Name of the compact quality column.

- qwidth_col:

  Name of the read-width column used to decode compact sequence bytes.

## Value

`x` with compact `seq` and/or `qual` columns replaced by decoded
character vectors. The input class is preserved.

## See also

[`decode_compact_seq()`](https://cparsania.github.io/BamScale/reference/decode_compact_seq.md),
[`decode_compact_qual()`](https://cparsania.github.io/BamScale/reference/decode_compact_qual.md),
[`bam_read()`](https://cparsania.github.io/BamScale/reference/bam_read.md)

## Examples

``` r
x <- data.frame(qwidth = 4L)
x$seq <- I(list(as.raw(c(0x12, 0x48))))
x$qual <- I(list(as.raw(c(0L, 1L, 2L, 3L))))
decode_seqqual_compact(x)
#>   qwidth  seq qual
#> 1      4 ACGT !"#$
```
