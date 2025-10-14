# gVenn 0.99.3

## New features

- Add `exportOverlapsToBed()` function to export genomic overlap groups
to BED format files.

## Documentation

- Updated vignette to include information about `exportOverlapsToBed()` and
guidance on choosing between `exportOverlaps()` (Excel format) and
`exportOverlapsToBed()` (BED format) based on data type and downstream needs.

# gVenn 0.99.3

## Minor updates

- Set default colors in `plotVenn()`.

# gVenn 0.99.2

## New features

- Add customization options for `plotVenn()`.

## Documentation

- Improved clarity in function documentation and examples.
- Enhanced vignette with additional customization examples for `plotVenn()`.

# gVenn 0.99.1

## Minor updates

- Package refinements and documentation improvements for Bioconductor submission.

# gVenn 0.99.0

## New fFatures

- Initial release of the `gVenn` package.
- Introduced a workflow for overlap analysis:
  - `computeOverlaps()` computes intersections across multiple sets of 
    `GRanges` or gene lists, returning counts and membership categories.
  - `extractOverlaps()` retrieves the actual elements (regions or genes) that 
    belong to each overlap group for downstream analysis.
  - `exportOverlaps()` exports overlap groups to an Excel file, creating one 
    sheet per group and converting `GRanges` to data frames when needed.
- Added visualization functions:
  - `plotVenn()` to draw proportional Venn diagrams based on overlaps 
    between genomic regions (e.g., ChIP-seq peaks).
  - `plotUpSet()` to visualize complex overlaps with an UpSet plot.
- Added `saveViz()` to export visualizations to PDF, PNG, or SVG formats, 
  with optional date tagging in filenames.
