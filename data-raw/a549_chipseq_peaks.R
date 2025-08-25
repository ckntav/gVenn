## code to prepare `a549_chipseq_peaks` dataset goes here

# Load required packages
library(rtracklayer)   # to import BED files as GRanges
library(GenomicRanges) # to construct and manipulate GRanges/GRangesList
library(usethis)       # to save processed datasets into data/ for the package

# Locate the gzipped BED files inside inst/extdata/ of the installed package
bed_med1 <- system.file("extdata", "A549_MED1_Dex.stdchr.bed.gz", package = "gVenn")
bed_brd4 <- system.file("extdata", "A549_BRD4_Dex.stdchr.bed.gz", package = "gVenn")
bed_gr   <- system.file("extdata", "A549_GR_Dex.stdchr.bed.gz", package = "gVenn")

# Helper function: keep only peaks on chr8 (or "8" depending on genome build)
keep_chr8 <- function(gr) {
    gr[seqnames(gr) %in% c("chr8", "8")]
}

# Import full consensus peak sets
med1_full <- rtracklayer::import(bed_med1)
brd4_full <- rtracklayer::import(bed_brd4)
gr_full   <- rtracklayer::import(bed_gr)

# Subset each peak set to chr8 only (keeps the dataset small and fast for examples)
med1_chr8 <- keep_chr8(med1_full)
brd4_chr8 <- keep_chr8(brd4_full)
gr_chr8 <- keep_chr8(gr_full)

# Combine into a GRangesList with descriptive names
a549_chipseq_peaks <- GRangesList("MED1_Dex_chr8" = med1_chr8,
                                  "BRD4_Dex_chr8" = brd4_chr8,
                                  "GR_Dex_chr8" = gr_chr8)

# Save the dataset into data/a549_chipseq_peaks.rda
#   - overwrite = TRUE allows updating
#   - compress = "xz" ensures small package size (Bioconductor friendly)
usethis::use_data(a549_chipseq_peaks, overwrite = TRUE, compress = "xz")
