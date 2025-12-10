# Extract Overlap Groups from Genomic or Set Overlap Results

This function extracts subsets of intersecting elements grouped by their
overlap category (e.g., "110"). For genomic overlaps, it returns a
`GRangesList`; for set overlaps, it returns a named list of character
vectors.

## Usage

``` r
extractOverlaps(overlap_object)
```

## Arguments

- overlap_object:

  A `GenomicOverlapsResult` or `SetOverlapsResult` object.

## Value

A named list of grouped intersecting elements:

- If input is a `GenomicOverlapsResult`, a `GRangesList` split by
  `intersect_category`.

- If input is a `SetOverlapsResult`, a named `list` of character vectors
  grouped by `intersect_category`.

## Examples

``` r
# Example with gene sets (built-in dataset)
data(gene_list)
res_sets <- computeOverlaps(gene_list)
group_gene <- extractOverlaps(res_sets)
group_gene
#> $group_001
#>  [1] "ACTN1"  "ALDOA"  "CRYBG1" "AK4"    "ACYP2"  "AFP"    "ACTN2"  "AMPD3" 
#>  [9] "ACADS"  "ACY1"   "ADH1B"  "ACTBP7" "ADH6"   "ACTBP4" "ANXA2"  "AGXT"  
#> [17] "TLE5"  
#> 
#> $group_010
#>  [1] "AFM"       "ADPRH"     "AIF1"      "ACVR2A"    "ACTA1"     "PLIN2"    
#>  [7] "ALDH1A3"   "ALK"       "ACTG2"     "ADCY8"     "ABR"       "ADCYAP1"  
#> [13] "ADRA2A"    "ADORA2B"   "ADRA1B"    "ANXA5"     "ACLY"      "ANK3"     
#> [19] "ADCY6"     "ACR"       "ADAM10"    "AARS1"     "ACVR2B"    "ACO2"     
#> [25] "ADH1C"     "PARP1"     "ADD2"      "ADARB1"    "ADCY2"     "ANXA1"    
#> [31] "ALCAM"     "ADORA2A"   "AMFR"      "AMBN"      "NAT2"      "A1BG"     
#> [37] "ACADL"     "ADH5"      "ACTG1P8"   "AADAC"     "ACOX1"     "ALDH9A1"  
#> [43] "ANXA2P2"   "ADORA2BP1" "ADRB3"    
#> 
#> $group_100
#>  [1] "ALPP"     "ACTG1P9"  "AHSG"     "ASIC2"    "ACTG1P10" "ALAS1"   
#>  [7] "AKT2"     "PARP1P1"  "ABCD1"    "SLC25A6"  "AAMP"     "ADCP1"   
#> [13] "ACADVL"   "ACTG1"    "ANGPT2"   "AGTR1"    "ACACB"    "ACTBP9"  
#> [19] "ALDH1B1"  "ADAR"     "ABCD2"    "AMHR2"    "ABCB7"    "ABCA1"   
#> [25] "PARP4"    "ACTG1P1"  "JAG1"     "ACTA2"    "ADH7"     "AP1B1"   
#> [31] "ACVR1"    "ACTN4"    "A2MP1"    "ABCA4"    "ALAD"     "ADRA1A"  
#> [37] "ADCY5"    "ALDOB"    "AP2B1"    "AMELY"    "ABL1"     "ACTC1"   
#> [43] "AK2"      "ALOX12B"  "ACTN3"    "AIC"      "ALB"      "NATP"    
#> [49] "ANG"      "AHR"      "ABCA2"    "ALPL"     "ANXA2P1"  "AMELX"   
#> [55] "AHCY"     "PARP1P2"  "ALOX5"    "AMPD1"    "AFA"      "ACADSB"  
#> [61] "AIH3"     "ACAN"     "AGA"      "AMY1C"    "ADSS2"    "ALDH2"   
#> [67] "ALOX15B" 
#> 
#> $group_011
#>  [1] "AMY1A"  "ACRV1"  "ALAS2"  "ABCA3"  "ALPG"   "AMY2A"  "ADH4"   "ADARB2"
#>  [9] "NR0B1"  "AMYP1"  "A2M"    "AGTR2"  "ALOX15" "ACP1"   "ADH1A"  "AF8T"  
#> 
#> $group_101
#> [1] "NAT1"   "ACAA1"  "AGT"    "AMD1P2"
#> 
#> $group_110
#>  [1] "ADRB1"    "ABAT"     "ALOX5AP"  "ADD1"     "ACVR1B"   "AANAT"   
#>  [7] "ADSL"     "ADCY7"    "ALX3"     "ALOX12P1" "ANGPT1"   "ACTG1P6" 
#> [13] "ADAM8"    "ACHE"     "ANCR"     "ACP3"     "ACP5"     "APLNR"   
#> [19] "ACTBP8"   "ADCY1"    "ADA"     
#> 
#> $group_111
#>  [1] "ACP2"      "ALDH3A1"   "ACTB"      "ACACA"     "ASIC1"     "SLC25A5"  
#>  [7] "ACTL6A"    "AMY2B"     "AMH"       "AMPH"      "ADK"       "ALDH3A2"  
#> [13] "ACTG1P3"   "ACO1"      "ACTG1P7"   "ALPI"      "ANXA4"     "AGL"      
#> [19] "ADRB2"     "ABCF1"     "ABO"       "AMD1"      "ALS3"      "ALOX12"   
#> [25] "AMBP"      "AMPD2"     "ALDH1A1"   "AFG3L1P"   "ADFN"      "ADCYAP1R1"
#> [31] "ADD3"      "ALOX12P2"  "BIN1"     
#> 

# Example with genomic regions (built-in dataset)
data(a549_chipseq_peaks)
res_genomic <- computeOverlaps(a549_chipseq_peaks)
group_genomic <- extractOverlaps(res_genomic)
group_genomic
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
