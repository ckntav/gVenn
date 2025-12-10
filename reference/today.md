# Today's Date at Package Load Time

This variable stores the current date (in "yyyymmdd" format) at the time
the package is loaded. It is useful for reproducible filenames (e.g., in
[`saveViz()`](https://ckntav.github.io/gVenn/reference/saveViz.md)), and
is automatically set when the package is attached.

## Usage

``` r
today
```

## Format

A character string (e.g., "20250624").

## Examples

``` r
# Print the date stored at package load
library(gVenn)
today
#> [1] "20251210"

# Use it in a filename
paste0("venn_plot_", today, ".pdf")
#> [1] "venn_plot_20251210.pdf"
```
