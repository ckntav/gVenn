## code to prepare `a549_chipseq_peaks` dataset goes here

# Load required packages
library(rtracklayer)   # to import BED files as GRanges
library(GenomicRanges) # to construct and manipulate GRanges/GRangesList
library(usethis)       # to save processed datasets into data/ for the package

# Locate the gzipped BED files inside inst/extdata/ of the installed package
bed_med1 <- system.file("extdata", "A549_MED1_Dex.stdchr.bed.gz", package = "gVenn")
bed_brd4 <- system.file("extdata", "A549_BRD4_Dex.stdchr.bed.gz", package = "gVenn")
bed_gr   <- system.file("extdata", "A549_GR_Dex.stdchr.bed.gz", package = "gVenn")

# Helper function: keep only peaks on chr7 (or "8" depending on genome build)
keep_chr7 <- function(gr) {
    gr[seqnames(gr) %in% c("chr7", "7")]
}

# Import full consensus peak sets
med1_full <- rtracklayer::import(bed_med1)
brd4_full <- rtracklayer::import(bed_brd4)
gr_full   <- rtracklayer::import(bed_gr)

# Subset each peak set to chr7 only (keeps the dataset small and fast for examples)
med1_chr7 <- keep_chr7(med1_full)
brd4_chr7 <- keep_chr7(brd4_full)
gr_chr7 <- keep_chr7(gr_full)

# Combine into a GRangesList with descriptive names
a549_chipseq_peaks <- GRangesList("MED1_Dex_chr7" = med1_chr7,
                                  "BRD4_Dex_chr7" = brd4_chr7,
                                  "GR_Dex_chr7" = gr_chr7)

# Save the dataset into data/a549_chipseq_peaks.rda
#   - overwrite = TRUE allows updating
#   - compress = "xz" ensures small package size (Bioconductor friendly)
usethis::use_data(a549_chipseq_peaks, overwrite = TRUE, compress = "xz")
