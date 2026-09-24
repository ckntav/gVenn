# gVenn 1.99.0

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
  `diagError`, `regionError`, `undrawnRegions`) to the returned plot as a
  `"fit_diagnostics"` attribute, and the fit itself as `"euler_fit"`. A message
  reports `diagError` and whether it falls above or below the 1e-6 threshold
  of Micallef and Rodgers (2014), who introduced the measure,
  and names the populated regions the diagram gives no area to at all, which
  `undrawnRegions` records. Set the new `verbose = FALSE` to silence it.
- Add `plotVennError()`, which redraws the diagram shaded by the signed error of
  each region, so the regions a diagram misrepresents can be read off the picture.
  A thin wrapper around  `eulerr::error_plot()`. It takes the plot returned
  by `plotVenn()` as input.
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
- `computeOverlaps()` labels overlap categories faster at large numbers of
  regions (~10x at 10^5-10^6 regions), by vectorizing the internal
  `defineCategories()` helper instead of looping row by row.

## Bug fixes

- Fix `plotVenn()` painting different regions of a four-set diagram in the same
  color. Since 1.3.2 the palette has been recycled to the number of regions, but
  it held seven colors against the fifteen regions of a four-set diagram, so
  eight of them repeated a color. Eight colors were appended. Two- and three-set
  diagrams draw on the unchanged first seven, so their output is identical.

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
