# Decode compact BamScale sequence output

Decodes `seq` values returned by
`bam_read(..., seqqual_mode = "compact")` back to ordinary character
strings.

## Usage

``` r
decode_compact_seq(seq, qwidth)
```

## Arguments

- seq:

  A list (or list-column) of `raw` vectors produced by compact BamScale
  sequence extraction.

- qwidth:

  Integer vector of read widths. This is required because compact
  sequence bytes use BAM's 4-bit packed encoding (two bases per byte).

## Value

A character vector containing decoded sequence strings.

## See also

[`decode_compact_qual()`](https://cparsania.github.io/BamScale/reference/decode_compact_qual.md),
[`decode_seqqual_compact()`](https://cparsania.github.io/BamScale/reference/decode_seqqual_compact.md),
[`bam_read()`](https://cparsania.github.io/BamScale/reference/bam_read.md)

## Examples

``` r
decode_compact_seq(
  seq = list(as.raw(c(0x12, 0x48))),
  qwidth = 4L
)
#> [1] "ACGT"
```
