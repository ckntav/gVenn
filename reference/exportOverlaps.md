# Export Overlap Groups to Excel

This function exports the output of
[`extractOverlaps()`](https://ckntav.github.io/gVenn/reference/extractOverlaps.md)
to an Excel file, creating one sheet per overlap group. Genomic overlaps
(`GRanges`) are converted to data frames before export.

## Usage

``` r
exportOverlaps(
  grouped,
  output_dir = ".",
  output_file = "overlap_groups",
  with_date = TRUE,
  verbose = TRUE
)
```

## Arguments

- grouped:

  Overlap groups from
  [`extractOverlaps()`](https://ckntav.github.io/gVenn/reference/extractOverlaps.md).

- output_dir:

  A string specifying the output directory. Defaults to `"."`.

- output_file:

  A string specifying the base filename (without extension). Defaults to
  `"overlap_groups"`.

- with_date:

  Logical (default `TRUE`). Whether to prepend the current date (from
  `today`) to the filename.

- verbose:

  Logical. If `TRUE`, print a message with the saved path. Default
  `TRUE`.

## Value

Overlap groups are saved to a Excel file on disk. Invisibly returns the
full path to the saved file.

## Examples

``` r
res <- computeOverlaps(list(A = letters[1:3], B = letters[2:4]))
grouped <- extractOverlaps(res)
exportOverlaps(grouped, output_dir = tempdir(), output_file = "overlap_groups")
#>  > Overlap groups saved in /tmp/RtmpfpCeIe/20251213_overlap_groups.xlsx
```
