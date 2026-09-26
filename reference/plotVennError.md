# Plot where a Venn diagram misrepresents the data

Redraws the Venn diagram produced by
[`plotVenn`](https://ckntav.github.io/gVenn/reference/plotVenn.md),
shaded by the signed error of each region, so that the regions the
diagram gets wrong can be read off the picture rather than inferred from
a single summary statistic. This is a thin wrapper around
[`error_plot`](https://jolars.github.io/eulerr/reference/error_plot.html).

## Usage

``` r
plotVennError(venn, ...)
```

## Arguments

- venn:

  A plot returned by
  [`plotVenn`](https://ckntav.github.io/gVenn/reference/plotVenn.md),
  carrying the `"euler_fit"` attribute that function attaches. Note that
  this is the
  [`plotVenn()`](https://ckntav.github.io/gVenn/reference/plotVenn.md)
  output, not the `GenomicOverlapResult` or `SetOverlapResult` that
  [`plotVenn()`](https://ckntav.github.io/gVenn/reference/plotVenn.md)
  itself takes.

- ...:

  Additional arguments passed to
  [`error_plot`](https://jolars.github.io/eulerr/reference/error_plot.html).

## Value

The plot returned by
[`error_plot`](https://jolars.github.io/eulerr/reference/error_plot.html),
a `gTree` that can be drawn, arranged with other grobs, or written out
with [`saveViz`](https://ckntav.github.io/gVenn/reference/saveViz.md).

## Details

Where `diagError` reports how far off the worst region is, this plot
shows which regions are off and in which direction, with each region's
`regionError` printed inside it. Use it when `diagError` exceeds 1e-6 or
when `attr(venn, "fit_diagnostics")$undrawnRegions` is not empty.

## See also

[`plotVenn`](https://ckntav.github.io/gVenn/reference/plotVenn.md) for
the diagram itself and its fit diagnostics,
[`plotUpSet`](https://ckntav.github.io/gVenn/reference/plotUpSet.md) for
a representation that does not approximate any count by an area.

## Examples

``` r
data(gene_list)
res_sets <- computeOverlaps(gene_list)

# Where does the diagram misrepresent the counts?
venn <- plotVenn(res_sets)
#> ✔ Venn diagError = 1.728e-14  (<= 1e-06)
#>   Access fit diagnostics with attr(<plotVenn output>, "fit_diagnostics")
plotVennError(venn)
```
