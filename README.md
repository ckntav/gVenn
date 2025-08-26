
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gVenn

<!-- badges: start -->

<!-- badges: end -->

**Proportional Venn and UpSet diagrams for genomic regions and gene set
overlaps.**

**gVenn** stands for **gene/genomic Venn**.  
It provides tools to compute overlaps between sets of gene or genomic
regions and visualize them as Venn diagrams with areas proportional to
the number of overlapping elements. With seamless support for `GRanges`
and `GRangesList` objects, **gVenn** integrates naturally into
Bioconductor workflows such as ChIP-seq, ATAC-seq, or other
interval-based analyses, and produces clean, publication-ready figures.

## Installation

You can install the development version of gVenn from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("ckntav/gVenn")
```

## Quick start

This quick example demonstrates how to compute overlaps between ChIP-seq
peaks and visualize them with both a Venn diagram and an UpSet plot.

### Load example ChIP-seq data and compute overlaps

``` r
library(gVenn)

# Example dataset of ChIP-seq peaks (A549 cell line, 3 genomic regions)
data(a549_chipseq_peaks)
a549_chipseq_peaks
#> Loading required namespace: GenomicRanges
#> GRangesList object of length 3:
#> $MED1_Dex_chr7
#> GRanges object with 336 ranges and 2 metadata columns:
#>         seqnames              ranges strand |        name     score
#>            <Rle>           <IRanges>  <Rle> | <character> <numeric>
#>     [1]     chr7     1157024-1157513      * |        4997         0
#>     [2]     chr7     1520389-1521218      * |        4998         0
#>     [3]     chr7     1536927-1537642      * |        4999         0
#>     [4]     chr7     2309837-2310506      * |        5000         0
#>     [5]     chr7     3028013-3028396      * |        5001         0
#>     ...      ...                 ...    ... .         ...       ...
#>   [332]     chr7 158733134-158733544      * |        5328         0
#>   [333]     chr7 158818327-158819201      * |        5329         0
#>   [334]     chr7 158821150-158821448      * |        5330         0
#>   [335]     chr7 158863388-158864513      * |        5331         0
#>   [336]     chr7 159015348-159016094      * |        5332         0
#>   -------
#>   seqinfo: 24 sequences from an unspecified genome; no seqlengths
#> 
#> $BRD4_Dex_chr7
#> GRanges object with 604 ranges and 2 metadata columns:
#>         seqnames              ranges strand |        name     score
#>            <Rle>           <IRanges>  <Rle> | <character> <numeric>
#>     [1]     chr7       234690-235402      * |        9419         0
#>     [2]     chr7       538240-538633      * |        9420         0
#>     [3]     chr7     1156721-1157555      * |        9421         0
#>     [4]     chr7     1504294-1504733      * |        9422         0
#>     [5]     chr7     1506830-1507301      * |        9423         0
#>     ...      ...                 ...    ... .         ...       ...
#>   [600]     chr7 158829343-158830028      * |       10018         0
#>   [601]     chr7 158856251-158856723      * |       10019         0
#>   [602]     chr7 158863108-158864616      * |       10020         0
#>   [603]     chr7 159012435-159013222      * |       10021         0
#>   [604]     chr7 159015311-159016245      * |       10022         0
#>   -------
#>   seqinfo: 24 sequences from an unspecified genome; no seqlengths
#> 
#> $GR_Dex_chr7
#> GRanges object with 450 ranges and 2 metadata columns:
#>         seqnames              ranges strand |        name     score
#>            <Rle>           <IRanges>  <Rle> | <character> <numeric>
#>     [1]     chr7       729847-730122      * |        6571         0
#>     [2]     chr7     1156806-1157495      * |        6572         0
#>     [3]     chr7     1520508-1521044      * |        6573         0
#>     [4]     chr7     2309959-2310483      * |        6574         0
#>     [5]     chr7     2860620-2860960      * |        6575         0
#>     ...      ...                 ...    ... .         ...       ...
#>   [446]     chr7 158733144-158733534      * |        7016         0
#>   [447]     chr7 158818350-158819168      * |        7017         0
#>   [448]     chr7 158821076-158821582      * |        7018         0
#>   [449]     chr7 158863549-158864364      * |        7019         0
#>   [450]     chr7 159015407-159016007      * |        7020         0
#>   -------
#>   seqinfo: 24 sequences from an unspecified genome; no seqlengths

# Compute overlaps
ov <- computeOverlaps(a549_chipseq_peaks)
```

### Visualize

``` r
# Draw Venn diagram
plotVenn(ov)
```

<img src="man/figures/README-example_venn-1.png" width="100%" />

``` r
# Draw UpSet plot (useful for larger numbers of sets)
plotUpSet(ov)
```

<img src="man/figures/README-example_upset-1.png" width="100%" />

### Extract elements per overlap group

``` r
groups <- extractOverlaps(ov)
groups
#> GRangesList object of length 7:
#> $group_010
#> GRanges object with 267 ranges and 1 metadata column:
#>         seqnames              ranges strand | intersect_category
#>            <Rle>           <IRanges>  <Rle> |        <character>
#>     [1]     chr7       234690-235402      * |                010
#>     [2]     chr7       538240-538633      * |                010
#>     [3]     chr7     1504294-1504733      * |                010
#>     [4]     chr7     1506830-1507301      * |                010
#>     [5]     chr7     1513353-1513690      * |                010
#>     ...      ...                 ...    ... .                ...
#>   [263]     chr7 155618941-155619523      * |                010
#>   [264]     chr7 155644241-155644737      * |                010
#>   [265]     chr7 158829343-158830028      * |                010
#>   [266]     chr7 158856251-158856723      * |                010
#>   [267]     chr7 159012435-159013222      * |                010
#>   -------
#>   seqinfo: 24 sequences from an unspecified genome; no seqlengths
#> 
#> ...
#> <6 more elements>
```

## Contributing

Pull requests are welcome.
