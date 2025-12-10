# Changelog

## gVenn 1.1.1

### New features

- Add `bg` parameter to
  [`saveViz()`](https://ckntav.github.io/gVenn/reference/saveViz.md) for
  controlling plot background color, including transparent backgrounds.
  Users can now save plots with `bg = "transparent"` for use in
  presentations or publications requiring transparent backgrounds.
- Add hex sticker logo created using the `hexSticker` R package.
- Update graphical abstract highlighting gVenn’s overlap visualization
  and extraction capabilities

### Documentation

- Add example to vignette demonstrating transparent background export
  using `bg = "transparent"` parameter in
  [`saveViz()`](https://ckntav.github.io/gVenn/reference/saveViz.md)

## gVenn 0.99.5

### Minor update

- Add `comb_col` parameter to
  [`plotUpSet()`](https://ckntav.github.io/gVenn/reference/plotUpSet.md)
  for customizing the color of combination matrix elements.

### Documentation

- Update UpSet plot example in the vignette with color customization of
  combination matrix elements using the `comb_col` parameter.

## gVenn 0.99.4

### New features

- Add
  [`exportOverlapsToBed()`](https://ckntav.github.io/gVenn/reference/exportOverlapsToBed.md)
  function to export genomic overlap groups to BED format files.

### Documentation

- Updated vignette to include information about
  [`exportOverlapsToBed()`](https://ckntav.github.io/gVenn/reference/exportOverlapsToBed.md)
  and guidance on choosing between
  [`exportOverlaps()`](https://ckntav.github.io/gVenn/reference/exportOverlaps.md)
  (Excel format) and
  [`exportOverlapsToBed()`](https://ckntav.github.io/gVenn/reference/exportOverlapsToBed.md)
  (BED format) based on data type and downstream needs.

## gVenn 0.99.3

### Minor updates

- Set default colors in
  [`plotVenn()`](https://ckntav.github.io/gVenn/reference/plotVenn.md).

## gVenn 0.99.2

### New features

- Add customization options for
  [`plotVenn()`](https://ckntav.github.io/gVenn/reference/plotVenn.md).

### Documentation

- Improved clarity in function documentation and examples.
- Enhanced vignette with additional customization examples for
  [`plotVenn()`](https://ckntav.github.io/gVenn/reference/plotVenn.md).

## gVenn 0.99.1

### Minor updates

- Package refinements and documentation improvements for Bioconductor
  submission.

## gVenn 0.99.0

### New features

- Initial release of the `gVenn` package.
- Introduced a workflow for overlap analysis:
  - [`computeOverlaps()`](https://ckntav.github.io/gVenn/reference/computeOverlaps.md)
    computes intersections across multiple sets of `GRanges` or gene
    lists, returning counts and membership categories.
  - [`extractOverlaps()`](https://ckntav.github.io/gVenn/reference/extractOverlaps.md)
    retrieves the actual elements (regions or genes) that belong to each
    overlap group for downstream analysis.
  - [`exportOverlaps()`](https://ckntav.github.io/gVenn/reference/exportOverlaps.md)
    exports overlap groups to an Excel file, creating one sheet per
    group and converting `GRanges` to data frames when needed.
- Added visualization functions:
  - [`plotVenn()`](https://ckntav.github.io/gVenn/reference/plotVenn.md)
    to draw proportional Venn diagrams based on overlaps between genomic
    regions (e.g., ChIP-seq peaks).
  - [`plotUpSet()`](https://ckntav.github.io/gVenn/reference/plotUpSet.md)
    to visualize complex overlaps with an UpSet plot.
- Added
  [`saveViz()`](https://ckntav.github.io/gVenn/reference/saveViz.md) to
  export visualizations to PDF, PNG, or SVG formats, with optional date
  tagging in filenames.
