# gVenn 2.0.0

## New features

- Add a `mode` argument to `computeOverlaps()` controlling how genomic
  intervals are made non-redundant before they are classified. The default,
  `mode = "reduce"`, keeps the previous "reduce-then-classify" behavior.
  The new `mode = "disjoin"` reduces each set on its own, then partitions the
  union into non-overlapping segments with `GenomicRanges::disjoin()`, so that
  every segment is covered by exactly one combination of sets. Because the 
  resulting intervals are merged in one mode and disjoint in the other, the
  `reduced_regions` element of `GenomicOverlapResult` is renamed to `regions`.
- `plotVenn()` now attaches `eulerr`'s goodness-of-fit diagnostics (`stress`,
  `diagError`, `regionError`) to the returned plot as a `"fit_diagnostics"`
  attribute, retrievable with `attr(venn, "fit_diagnostics")`. A message
  reporting `stress` and `diagError` is printed by default; set the new
  `verbose = FALSE` argument to silence it.
- Add an `ignore.strand` argument to `computeOverlaps()`/`computeGenomicOverlaps()`,
  passed through to `GenomicRanges::reduce()`, `GenomicRanges::disjoin()`, and
  `IRanges::overlapsAny()`. Defaults to `FALSE` (previous, strand-aware
  behavior is unchanged); set to `TRUE` to disregard strand when merging or
  partitioning regions and when determining overlaps.
- `computeOverlaps()` now warns when two or more input genomic sets share no
  chromosome name at all, a common symptom of mismatched chromosome naming
  conventions (e.g. `"chr1"` vs `"1"`) or of comparing different genome
  assemblies, cases where overlaps would otherwise be silently and
  permanently empty. Genome-assembly conflicts on a *shared* chromosome name
  already error via `GenomicRanges::GRangesList()`, unchanged.

## Minor updates

- The right-hand annotation of `plotUpSet()` is now labelled according to the
  type of the input: `"Region size"` for a `GenomicOverlapResult` and
  `"Set size"` for a `SetOverlapResult`.

# gVenn 1.3.2

## Bug fixes

- Fix `plotVenn()` failing with `fills$fill must have length 1, n_sets, or
  n_subsets` when the data did not populate every region of the diagram.
  The default fill palette is now recycled to `length(fit$original.values)`
  (eulerr's `n_subsets`), so every region receives a color regardless of
  which combinations are populated.

# gVenn 1.3.1

## Bug fixes

- Fix `plotVenn()` failing on 2-set inputs. The default fill palette had
  a fixed length of 7, which violated eulerr's stricter validation
  (`fills$fill` must have length 1, `n_sets`, or `n_subsets`). The default
  is now recycled to match `n_sets`.

# gVenn 1.1.1

## New features

- Add `bg` parameter to `saveViz()` for controlling plot background color,
including transparent backgrounds. Users can now save plots with
`bg = "transparent"` for use in presentations or publications requiring
transparent backgrounds.
- Add hex sticker logo created using the `hexSticker` R package.
- Update graphical abstract highlighting gVenn's overlap visualization and 
extraction capabilities

## Documentation

- Add example to vignette demonstrating transparent background export using
`bg = "transparent"` parameter in `saveViz()`

# gVenn 0.99.5

## Minor update

- Add `comb_col` parameter to `plotUpSet()` for customizing the color of
combination matrix elements.

## Documentation

- Update UpSet plot example in the vignette with color customization of
combination matrix elements using the `comb_col` parameter.

# gVenn 0.99.4

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

## New features

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
