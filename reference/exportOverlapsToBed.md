# Export Overlap Groups to BED Files

This function exports genomic overlap groups from
[`extractOverlaps()`](https://ckntav.github.io/gVenn/reference/extractOverlaps.md)
to BED format files, creating one BED file per overlap group.

## Usage

``` r
exportOverlapsToBed(
  grouped,
  output_dir = ".",
  output_prefix = "overlaps",
  with_date = TRUE,
  verbose = TRUE
)
```

## Arguments

- grouped:

  Genomic overlap groups from
  [`extractOverlaps()`](https://ckntav.github.io/gVenn/reference/extractOverlaps.md)
  (must be `GRangesList`).

- output_dir:

  A string specifying the output directory. Defaults to `"."`.

- output_prefix:

  A string specifying the filename prefix. Defaults to `"overlaps"`.

- with_date:

  Logical (default `TRUE`). Whether to prepend the current date to
  filenames.

- verbose:

  Logical. If `TRUE`, print messages. Default `TRUE`.

## Value

Invisibly returns a character vector of file paths created.

## Details

This function only works with genomic overlaps (i.e., when the input to
[`extractOverlaps()`](https://ckntav.github.io/gVenn/reference/extractOverlaps.md)
was a `GenomicOverlapResult` object, resulting in a `GRangesList`). It
does not work with set overlaps (character vectors). Each overlap group
will be saved as a separate BED file with the group identifier included
in the filename.
