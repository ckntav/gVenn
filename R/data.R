#' A549 ChIP-seq Consensus Peak Subsets (Dex, chr7)
#'
#' Example consensus peak subsets for MED1, BRD4, and GR after dexamethasone
#' treatment in A549 cells. Each set has been restricted to peaks on
#' \code{chr7} to keep the dataset small and suitable for examples and tests.
#'
#' @format A \code{GRangesList} with 3 named elements:
#' \describe{
#'   \item{MED1_Dex_chr7}{Consensus MED1 peaks (chr7 subset).}
#'   \item{BRD4_Dex_chr7}{Consensus BRD4 peaks (chr7 subset).}
#'   \item{GR_Dex_chr7}{Consensus GR peaks (chr7 subset).}
#' }
#'
#' @details
#' The original full consensus peak sets are available as gzipped BED files in
#' \code{inst/extdata/}:
#' \itemize{
#'   \item \code{A549_MED1_Dex.stdchr.bed.gz}
#'   \item \code{A549_BRD4_Dex.stdchr.bed.gz}
#'   \item \code{A549_GR_Dex.stdchr.bed.gz}
#' }
#' These are not trimmed, but for package efficiency the dataset here
#' (\code{a549_chipseq_peaks}) only includes the chr7 subsets.
#'
#' @source Internal consensus peak sets processed in A549 cells after
#' dexamethasone stimulation.
#'
#' @references
#' Tav C, Fournier É, Fournier M, Khadangi F, Baguette A, Côté MC, Silveira MAD,
#' Bérubé-Simard F-A, Bourque G, Droit A, Bilodeau S (2023).
#' "Glucocorticoid stimulation induces regionalized gene responses within
#' topologically associating domains." \emph{Frontiers in Genetics}.
#' \doi{10.3389/fgene.2023.1237092}
#'
#' @examples
#' # Load dataset
#' data(a549_chipseq_peaks)
#' a549_chipseq_peaks
#'
#' # Compute overlaps and plot
#' ov <- computeOverlaps(a549_chipseq_peaks)
#' plotVenn(ov)
"a549_chipseq_peaks"

#' Example Gene Lists with Overlaps
#'
#' A synthetic dataset of three gene lists, created from the first 250 human
#' gene symbols in \pkg{org.Hs.eg.db}.
#'
#' @format A named \code{list} of length 3. Each element is a character vector
#' of gene symbols:
#' \describe{
#'   \item{random_genes_A}{125 gene symbols.}
#'   \item{random_genes_B}{115 gene symbols.}
#'   \item{random_genes_C}{70 gene symbols.}
#' }
#'
#' @source Generated from \pkg{org.Hs.eg.db} (keys of type \code{SYMBOL}),
#' using a reproducible random seed.
#'
#' @examples
#' data(gene_list)
#'
#' # Inspect the list
#' str(gene_list)
#'
#' # Compute overlaps and plot
#' ov <- computeOverlaps(gene_list)
#' plotVenn(ov)
"gene_list"
