# Compute Overlaps Between Multiple Sets or Genomic Regions

`computeOverlaps()` is the main entry point for overlap analysis. It
accepts either genomic region objects (`GRanges`/`GRangesList`) or
ordinary sets (character/numeric vectors) and computes a binary overlap
matrix describing the presence or absence of each element across sets.

## Usage

``` r
computeOverlaps(x, mode = c("reduce", "disjoin"), ignore.strand = FALSE)
```

## Arguments

- x:

  Input sets. One of:

  - A `GRangesList` object.

  - A named list of `GRanges` objects.

  - A named list of atomic vectors (character, numeric, factor, etc.),
    all of the same type.

- mode:

  Character string selecting where the boundaries of the partition fall.
  One of:

  - `"reduce"` (default):
    [`GenomicRanges::reduce()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html)
    collapses the union of all intervals into a non-redundant collection
    of "reduced regions".

  - `"disjoin"`: each set's intervals are collapsed individually and the
    union is then cut at every set boundary with
    [`GenomicRanges::disjoin()`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html),
    yielding a larger number of smaller, position-exact "disjoint
    regions".

  See Examples for a chained configuration under both modes. Ignored
  (with a warning) for non-genomic inputs.

- ignore.strand:

  Logical, default `FALSE`. Controls whether strand is taken into
  account when regions are made non-redundant (`reduce()`/ `disjoin()`)
  and when overlaps against the input sets are determined
  (`overlapsAny()`). With the default `FALSE`, regions on opposite
  strands are never merged or considered overlapping. With `TRUE`,
  strand is disregarded throughout, and regions on opposite strands can
  be merged into a single, unstranded (`"*"`) region. Ignored (with a
  warning) for non-genomic inputs.

## Value

An S3 object encoding the overlap result whose class depends on the
input type:

- GenomicOverlapResult:

  Returned when the input is genomic (`GRangesList` or list of
  `GRanges`). A list with:

  - `regions`: A `GRanges` object containing the reduced regions
    (`mode = "reduce"`) or disjoint regions (`mode = "disjoin"`). Each
    region is annotated with an `intersect_category` column giving the
    binary code of its overlap group.

  - `overlap_matrix`: A logical matrix indicating whether each region
    overlaps each input set (rows = regions, columns = sets).

  - `mode`: The `mode` used to build the regions.

- SetOverlapResult:

  Returned when the input is a list of atomic vectors. A list with:

  - `unique_elements`: Character vector of all unique elements across
    the sets.

  - `overlap_matrix`: A logical matrix indicating whether each element
    is present in each set (rows = elements, columns = sets).

  - `intersect_category`: Character vector giving the binary code of the
    overlap group (e.g., `"110"`) of each element.

## Details

- When provided with genomic regions, the function builds a
  non-redundant set of intervals (see `mode`), then determines which
  original sets each region overlaps.

- When provided with ordinary sets (e.g., gene symbols), it collects all
  unique elements and records which sets contain them.

The resulting object encodes both the overlap matrix and compact
category labels (e.g., `"110"`) representing the overlap pattern of each
element. These results can be directly passed to visualization functions
such as
[`plotVenn()`](https://ckntav.github.io/gVenn/reference/plotVenn.md) or
[`plotUpSet()`](https://ckntav.github.io/gVenn/reference/plotUpSet.md).

Internally, `computeOverlaps()` dispatches to either
`computeGenomicOverlaps()` (for genomic inputs) or
`computeSetOverlaps()` (for ordinary sets). Users are encouraged to call
only `computeOverlaps()`.

### Chromosome names and genome assemblies

gVenn does not harmonize chromosome names or coordinate systems across
input sets: genomic sets are expected to follow the same naming
convention and to come from the same genome assembly. As a safeguard,
`computeOverlaps()` warns whenever two input sets share no chromosome
name at all, naming the offending pairs. Such a mismatch usually signals
incompatible naming conventions (e.g. `"chr1"` vs `"1"`) or different
assemblies, cases in which overlaps between those sets would otherwise
be silently and permanently empty rather than raising an error. Only
chromosomes actually present in each set are compared. Conflicting
assemblies declared on a *shared* chromosome name (e.g. `"chr1"` tagged
`hg38` in one set and `hg19` in another) already raise an error in
[`GenomicRanges::GRangesList()`](https://rdrr.io/pkg/GenomicRanges/man/GRangesList-class.html),
before any overlap is computed.

## See also

[`plotVenn`](https://ckntav.github.io/gVenn/reference/plotVenn.md),
[`plotUpSet`](https://ckntav.github.io/gVenn/reference/plotUpSet.md),
[`GRangesList`](https://rdrr.io/pkg/GenomicRanges/man/GRangesList-class.html),
[`reduce`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html),
[`disjoin`](https://rdrr.io/pkg/IRanges/man/inter-range-methods.html)

## Examples

``` r
# Example with gene sets (built-in dataset)
data(gene_list)
ov_sets <- computeOverlaps(gene_list)
head(ov_sets$overlap_matrix)
#>          random_genes_A random_genes_B random_genes_C
#> ALPP               TRUE          FALSE          FALSE
#> ACTG1P9            TRUE          FALSE          FALSE
#> AHSG               TRUE          FALSE          FALSE
#> ASIC2              TRUE          FALSE          FALSE
#> ACTG1P10           TRUE          FALSE          FALSE
#> ALAS1              TRUE          FALSE          FALSE
plotVenn(ov_sets)
#> ✔ Venn diagError = 8.693e-13  (<= 1e-06)
#>   Access fit diagnostics with attr(<plotVenn output>, "fit_diagnostics")


# Example with genomic regions (built-in dataset)
data(a549_chipseq_peaks)
ov_gr <- computeOverlaps(a549_chipseq_peaks)
head(ov_gr$overlap_matrix)
#>      MED1_Dex_chr7 BRD4_Dex_chr7 GR_Dex_chr7
#> [1,]         FALSE          TRUE       FALSE
#> [2,]         FALSE          TRUE       FALSE
#> [3,]         FALSE         FALSE        TRUE
#> [4,]          TRUE          TRUE        TRUE
#> [5,]         FALSE          TRUE       FALSE
#> [6,]         FALSE          TRUE       FALSE
plotVenn(ov_gr)
#> ✔ Venn diagError = 2.229e-12  (<= 1e-06)
#>   Access fit diagnostics with attr(<plotVenn output>, "fit_diagnostics")


# Chained overlaps: A-B and B-C overlap, but A and C do not
A <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
B <- GenomicRanges::GRanges("chr1", IRanges::IRanges(180, 300))
C <- GenomicRanges::GRanges("chr1", IRanges::IRanges(280, 400))

# "reduce" merges the chain into a single "111" region
computeOverlaps(list(A = A, B = B, C = C))$regions
#> GRanges object with 1 range and 1 metadata column:
#>       seqnames    ranges strand | intersect_category
#>          <Rle> <IRanges>  <Rle> |        <character>
#>   [1]     chr1   100-400      * |                111
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

# "disjoin" keeps A-B and B-C as separate two-way intersections
computeOverlaps(list(A = A, B = B, C = C), mode = "disjoin")$regions
#> GRanges object with 5 ranges and 1 metadata column:
#>       seqnames    ranges strand | intersect_category
#>          <Rle> <IRanges>  <Rle> |        <character>
#>   [1]     chr1   100-179      * |                100
#>   [2]     chr1   180-200      * |                110
#>   [3]     chr1   201-279      * |                010
#>   [4]     chr1   280-300      * |                011
#>   [5]     chr1   301-400      * |                001
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths

#  Overlapping regions on opposite strands: kept separate by default,
# merged when ignore.strand = TRUE
D <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200), strand = "+")
E <- GenomicRanges::GRanges("chr1", IRanges::IRanges(150, 250), strand = "-")
computeOverlaps(list(D = D, E = E))$regions
#> GRanges object with 2 ranges and 1 metadata column:
#>       seqnames    ranges strand | intersect_category
#>          <Rle> <IRanges>  <Rle> |        <character>
#>   [1]     chr1   100-200      + |                 10
#>   [2]     chr1   150-250      - |                 01
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
computeOverlaps(list(D = D, E = E), ignore.strand = TRUE)$regions
#> GRanges object with 1 range and 1 metadata column:
#>       seqnames    ranges strand | intersect_category
#>          <Rle> <IRanges>  <Rle> |        <character>
#>   [1]     chr1   100-250      * |                 11
#>   -------
#>   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
