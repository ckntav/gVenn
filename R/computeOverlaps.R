#' Encode Overlap Patterns into Category Strings
#'
#' Internal helper that converts a logical or numeric matrix of presence/absence
#' values into compact category codes. Each row of the matrix is collapsed into
#' a character string (e.g., \code{"110"}, \code{"101"}), representing which sets
#' are present (1) or absent (0) for that element.
#'
#' @param data A logical or numeric matrix where rows represent elements
#'     (e.g., genomic regions or genes) and columns represent sets. Entries must
#'     be interpretable as 0/1 (e.g., \code{TRUE}/\code{FALSE}, \code{1}/\code{0}).
#'
#' @return A character vector of category strings, one per row of \code{data}.
#'
#' @details This function is used internally by
#'     \code{\link{computeGenomicOverlaps}} and
#'     \code{\link{computeSetOverlaps}} to label elements with their overlap
#'     pattern.
#'
#' @examples
#' m <- matrix(c(TRUE, FALSE, TRUE,
#'               TRUE, TRUE, FALSE),
#'             nrow = 2, byrow = TRUE)
#' define_categories(m)
#' # Returns: c("101", "110")
#'
#' @keywords internal
#' @noRd
define_categories <- function(data) {
    categories <- apply(data, 1, function(row) {
        paste0(as.integer(row), collapse = "")
    })
    return(categories)
}

#' Compute Genomic Overlaps Across GRanges Sets
#'
#' This function computes overlaps across multiple genomic region sets provided
#' as a `GRangesList` or a list of `GRanges` objects.
#' It reduces all regions into a unified, non-redundant set and determines which
#'  original sets each region overlaps.
#' This facilitates the analysis and visualization of genomic intersection
#' patterns (e.g., using Venn or UpSet plots).
#'
#' @param genomic_regions A `GRangesList` or a named list of `GRanges` objects.
#'   Each element should represent a genomic region set (e.g., ChIP-seq peaks,
#'   annotated genes, etc.).
#'
#' @return An object of class `GenomicOverlapsResult`, which is a list with the
#' following components:
#' \describe{
#'   \item{reduced_regions}{A `GRanges` object containing the reduced (merged)
#'   genomic intervals across all sets.
#'   Each region is annotated with an `intersect_category` string representing
#'   the overlap pattern (e.g., `"110"`).}
#'   \item{overlap_matrix}{A logical matrix indicating which reduced regions
#'   overlap with which input sets.
#'   Rows correspond to reduced regions; columns correspond to the input sets.}
#' }
#'
#' @details Internally, the function uses `GenomicRanges::reduce()` to merge
#' overlapping or adjacent regions across all sets.
#'   It then determines overlaps between each reduced region and the original
#'   input sets using `IRanges::overlapsAny()`.
#'   The resulting matrix can be used to generate set diagrams or for further
#'   statistical analysis.
#'
#' @seealso \code{\link[GenomicRanges]{GRangesList}},
#' \code{\link[GenomicRanges]{reduce}}, \code{\link[IRanges]{overlapsAny}},
#'   \code{\link{plotVenn}}, \code{\link{plotUpSet}}
#'
#' @examples
#' gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 500),
#'                               width = 100))
#' gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(150, 700),
#'                               width = 100))
#' gr3 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(900),
#'                               width = 100))
#'
#' peak_sets <- list(H3K27ac = gr1, MED1 = gr2, BRD4 = gr3)
#' overlap_result <- computeGenomicOverlaps(peak_sets)
#'
#' head(overlap_result$overlap_matrix)
#' GenomicRanges::mcols(overlap_result$reduced_regions)$intersect_category
#'
#' @export
computeGenomicOverlaps <- function(genomic_regions) {
    if (inherits(genomic_regions, "list")) {
        genomic_regions <- GenomicRanges::GRangesList(genomic_regions)
    } else if (!inherits(genomic_regions, "GRangesList")) {
        stop("Input must be a list of GRanges or a GRangesList.")
    }

    reduced_regions <- GenomicRanges::reduce(unlist(genomic_regions))
    overlap_matrix <- matrix(FALSE,
                             nrow = length(reduced_regions),
                             ncol = length(genomic_regions))

    for (i in seq_along(genomic_regions)) {
        overlap_matrix[, i] <- IRanges::overlapsAny(reduced_regions,
                                                    genomic_regions[[i]])
    }

    colnames(overlap_matrix) <- names(genomic_regions)

    intersect_category <- define_categories(overlap_matrix)
    GenomicRanges::mcols(reduced_regions)$intersect_category <- intersect_category

    res <- list(
        reduced_regions = reduced_regions,
        overlap_matrix = overlap_matrix
    )
    class(res) <- "GenomicOverlapResult"
    return(res)
}

#' Compute Overlaps Between Named Sets
#'
#' This function computes overlaps across a list of character vectors (e.g., gene symbols, transcript IDs, region names),
#' returning a binary matrix of presence/absence and overlap categories per element.
#'
#' @param named_sets A named list of character vectors, where each vector contains identifiers (e.g., gene symbols) belonging to a set.
#'
#' @return An object of class `SetOverlapsResult`, a list with the following components:
#' \describe{
#'   \item{unique_elements}{A character vector of all unique elements across the input sets.}
#'   \item{overlap_matrix}{A logical matrix indicating for each element (rows) whether it is present in each set (columns).}
#'   \item{intersect_category}{A character vector encoding the pattern of overlaps per element (e.g., "110", "101").}
#' }
#'
#' @examples
#' gene_sets <- list(
#'   TF1_targets = c("TP53", "BRCA1", "MYC"),
#'   TF2_targets = c("MYC", "ESR1"),
#'   TF3_targets = c("TP53", "GATA3")
#' )
#'
#' res <- computeSetOverlaps(gene_sets)
#' head(res$overlap_matrix)
#' table(res$intersect_category)
#'
#' # Can be passed to plotVenn() or plotUpSet()
#' plotVenn(res)
#'
#' @export
computeSetOverlaps <- function(named_sets) {
    stopifnot(is.list(named_sets),
              all(vapply(named_sets, is.character, logical(1))))

    all_elements <- unique(unlist(named_sets))
    overlap_matrix <- matrix(FALSE,
                             nrow = length(all_elements),
                             ncol = length(named_sets))
    rownames(overlap_matrix) <- all_elements
    colnames(overlap_matrix) <- names(named_sets)

    for (i in seq_along(named_sets)) {
        overlap_matrix[all_elements %in% named_sets[[i]], i] <- TRUE
    }

    intersect_category <- define_categories(overlap_matrix)

    res <- list(
        unique_elements = all_elements,
        overlap_matrix = overlap_matrix,
        intersect_category = intersect_category
    )
    class(res) <- "SetOverlapResult"
    return(res)
}

#' Compute Overlaps Between Multiple Sets or Genomic Regions
#'
#' `computeOverlaps()` is the main entry point for overlap analysis. It accepts
#' either genomic region objects (`GRanges`/`GRangesList`) or ordinary sets
#' (character/numeric vectors) and computes a binary overlap matrix describing
#' the presence or absence of each element across sets.
#'
#' - When provided with genomic regions, the function merges all intervals into
#'   a non-redundant set (`reduce()`), then determines which original sets each
#'   region overlaps.
#' - When provided with ordinary sets (e.g., gene symbols), it collects all
#'   unique elements and records which sets contain them.
#'
#' The resulting object encodes both the overlap matrix and compact category
#' labels (e.g., `"110"`) representing the overlap pattern of each element.
#' These results can be directly passed to visualization functions such as
#' `plotVenn()` or `plotUpSet()`.
#'
#' @param x Input sets. One of:
#'   \itemize{
#'     \item A `GRangesList` object.
#'     \item A named list of `GRanges` objects.
#'     \item A named list of atomic vectors (character, numeric, factor, etc.),
#'       all of the same type.
#'   }
#'
#' @return
#' An S3 object encoding the overlap result whose class depends on the input type:
#'
#' \describe{
#'   \item{GenomicOverlapResult}{Returned when the input is genomic
#'       (`GRangesList` or list of `GRanges`). A list with:
#'       \itemize{
#'         \item \code{reduced_regions}: A `GRanges` object containing the
#'             merged (non-redundant) intervals. Each region is annotated with
#'             an \code{intersect_category} column.
#'         \item \code{overlap_matrix}: A logical matrix indicating whether each
#'             reduced region overlaps each input set (rows = regions,
#'             columns = sets).
#'       }}
#'   \item{SetOverlapResult}{Returned when the input is a list of atomic
#'       vectors. A list with:
#'       \itemize{
#'         \item \code{unique_elements}: Character vector of all unique elements
#'             across the sets.
#'         \item \code{overlap_matrix}: A logical matrix indicating whether each
#'             element is present in each set (rows = elements, columns = sets).
#'         \item \code{intersect_category}: Character vector of category codes
#'             (e.g., `"110"`) for each element.
#'       }}
#' }
#'
#' @details
#' Internally, `computeOverlaps()` dispatches to either
#' `computeGenomicOverlaps()` (for genomic inputs) or
#' `computeSetOverlaps()` (for ordinary sets). Users are encouraged to call
#' only `computeOverlaps()`.
#'
#' @examples
#' # Example with simple sets
#' sets <- list(A = letters[1:4], B = letters[3:6])
#' ov1 <- computeOverlaps(sets)
#' head(ov1$overlap_matrix)
#'
#' # Example with genomic regions
#' if (requireNamespace("GenomicRanges", quietly = TRUE)) {
#'     gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 50), width = 20))
#'     gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(15, 90), width = 20))
#'     ov2 <- computeOverlaps(list(A = gr1, B = gr2))
#'     head(ov2$overlap_matrix)
#' }
#'
#' @seealso \code{\link{plotVenn}}, \code{\link{plotUpSet}},
#'     \code{\link[GenomicRanges]{GRangesList}}, \code{\link[GenomicRanges]{reduce}}
#'
#' @export
computeOverlaps <- function(x) {
    if (missing(x) || is.null(x)) {
        stop("'x' must be provided.", call. = FALSE)
    }

    # ---- direct GRangesList -------------------------------------------------
    if (methods::is(x, "GRangesList")) {
        if (is.null(names(x))) names(x) <- paste0("set", seq_along(x))
        return(computeGenomicOverlaps(x))
    }

    # ---- list inputs --------------------------------------------------------
    if (!is.list(x)) {
        stop("'x' must be a GRangesList, a list of GRanges, or a list of atomic vectors.",
             call. = FALSE)
    }
    if (length(x) == 0L) {
        stop("'x' is an empty list.", call. = FALSE)
    }
    if (is.null(names(x))) {
        names(x) <- paste0("set", seq_along(x))
    }

    is_gr <- vapply(x, function(e) methods::is(e, "GRanges"), logical(1))
    if (all(is_gr)) {
        # list of GRanges -> genomic
        return(computeGenomicOverlaps(x))
    }

    # atomic vectors (genes/ids). Allow numeric/factor but coerce to character.
    is_atomic_vec <- vapply(x, function(e) is.atomic(e) && !is.list(e), logical(1))
    if (all(is_atomic_vec)) {
        x_chr <- lapply(x, function(e) as.character(e))
        return(computeSetOverlaps(x_chr))
    }

    # mixed or unsupported
    bad_idx <- which(!(is_gr | is_atomic_vec))[1]
    bad_cls <- paste(class(x[[bad_idx]]), collapse = "/")
    stop("All elements of 'x' must be the same type: either all GRanges or all ",
         "atomic vectors. Found unsupported element of class: ", bad_cls, ".",
         call. = FALSE)
}
