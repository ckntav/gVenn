#' Encode Overlap Patterns into Category Strings
#'
#' Internal helper that converts a logical or numeric matrix of presence/absence
#' values into compact category codes. Each row of the matrix is collapsed into
#' a character string (e.g., \code{"110"}, \code{"101"}), representing which
#' sets are present (1) or absent (0) for that element.
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
#' defineCategories(m)
#' # Returns: c("101", "110")
#'
#' @keywords internal
#' @noRd
defineCategories <- function(data) {
    categories <- apply(data, 1, function(row) {
        paste0(as.integer(row), collapse = "")
    })
    return(categories)
}

#' Compute Genomic Overlaps Across GRanges Sets
#'
#' This function computes overlaps across multiple genomic region sets provided
#' as a `GRangesList` or a list of `GRanges` objects.
#' It builds a unified, non-redundant set of intervals and determines which
#'  original sets each interval overlaps.
#' This facilitates the analysis and visualization of genomic intersection
#' patterns (e.g., using Venn or UpSet plots).
#'
#' @param genomic_regions A `GRangesList` or a named list of `GRanges` objects.
#'   Each element should represent a genomic region set (e.g., ChIP-seq peaks,
#'   annotated genes, etc.).
#' @param mode Character string, either `"reduce"` (default) or `"disjoin"`.
#'   See the `mode` documentation of \code{\link{computeOverlaps}}.
#'
#' @return An object of class `GenomicOverlapsResult`, which is a list with the
#' following components:
#' \describe{
#'   \item{regions}{A `GRanges` object containing the merged
#'   (`mode = "reduce"`) or disjoint (`mode = "disjoin"`) genomic intervals
#'   across all sets.
#'   Each region is annotated with an `intersect_category` string representing
#'   the overlap pattern (e.g., `"110"`).}
#'   \item{overlap_matrix}{A logical matrix indicating which regions
#'   overlap with which input sets.
#'   Rows correspond to regions; columns correspond to the input sets.}
#'   \item{mode}{The `mode` used to build the regions.}
#' }
#'
#' @details With `mode = "reduce"`, the function uses `GenomicRanges::reduce()`
#'   to merge overlapping or adjacent regions across all sets. With
#'   `mode = "disjoin"`, each set is first reduced on its own (to drop
#'   within-set redundancy), then `GenomicRanges::disjoin()` partitions the
#'   union into non-overlapping segments delimited by every set boundary.
#'   In both cases, overlaps between the resulting regions and the original
#'   input sets are determined with `IRanges::overlapsAny()`.
#'
#' @seealso \code{\link[GenomicRanges]{GRangesList}},
#' \code{\link[GenomicRanges]{reduce}}, \code{\link[GenomicRanges]{disjoin}},
#' \code{\link[IRanges]{overlapsAny}},
#'   \code{\link{plotVenn}}, \code{\link{plotUpSet}}
#'
#' @examples
#' library(gVenn)
#'
#' # Example dataset of A549 ChIP-seq peaks (3 sets)
#' data(a549_chipseq_peaks)
#'
#' ov <- computeGenomicOverlaps(a549_chipseq_peaks)
#'
#' # Inspect the overlap matrix
#' head(ov$overlap_matrix)
#'
#' # Check the intersection category assigned to each region
#' GenomicRanges::mcols(ov$regions)$intersect_category
#'
#' # Visualize with a Venn diagram
#' plotVenn(ov)
#'
#' @keywords internal
#' @noRd
computeGenomicOverlaps <- function(genomic_regions, mode = c("reduce", "disjoin")) {
    mode <- match.arg(mode)

    if (inherits(genomic_regions, "list")) {
        genomic_regions <- GenomicRanges::GRangesList(genomic_regions)
    } else if (!inherits(genomic_regions, "GRangesList")) {
        stop("Input must be a list of GRanges or a GRangesList.")
    }

    if (mode == "reduce") {
        regions <- GenomicRanges::reduce(unlist(genomic_regions))
    } else {
        # Reduce each set on its own first, so that intervals that are
        # redundant *within* a set do not introduce spurious breakpoints,
        # then split the union at every remaining set boundary.
        per_set <- GenomicRanges::reduce(genomic_regions)
        regions <- GenomicRanges::disjoin(unlist(per_set))
    }

    overlap_matrix <- matrix(FALSE,
                             nrow = length(regions),
                             ncol = length(genomic_regions))

    for (i in seq_along(genomic_regions)) {
        overlap_matrix[, i] <- IRanges::overlapsAny(regions,
                                                    genomic_regions[[i]])
    }

    colnames(overlap_matrix) <- names(genomic_regions)

    intersect_category <- defineCategories(overlap_matrix)
    GenomicRanges::mcols(regions)$intersect_category <- intersect_category

    res <- list(
        regions = regions,
        overlap_matrix = overlap_matrix,
        mode = mode
    )
    class(res) <- "GenomicOverlapResult"
    return(res)
}

#' Compute Overlaps Between Named Sets
#'
#' This function computes overlaps across a list of character vectors
#' (e.g., gene symbols, transcript IDs, region names),
#' returning a binary matrix of presence/absence and overlap categories per
#' element.
#'
#' @param named_sets A named list of character vectors, where each vector
#' contains identifiers (e.g., gene symbols) belonging to a set.
#'
#' @return An object of class `SetOverlapsResult`, a list with the following
#' components:
#' \describe{
#'   \item{unique_elements}{A character vector of all unique elements across
#'   the input sets.}
#'   \item{overlap_matrix}{A logical matrix indicating for each element (rows)
#'   whether it is present in each set (columns).}
#'   \item{intersect_category}{A character vector encoding the pattern of
#'   overlaps per element (e.g., "110", "101").}
#' }
#'
#' @examples
#' library(gVenn)
#'
#' # Example gene lists dataset (3 sets with overlaps)
#' data(gene_list)
#'
#' res <- computeSetOverlaps(gene_list)
#'
#' # Inspect the overlap matrix
#' head(res$overlap_matrix)
#'
#' # Summarize overlap categories
#' table(res$intersect_category)
#'
#' # Visualize with a Venn diagram
#' plotVenn(res)
#'
#' @keywords internal
#' @noRd
computeSetOverlaps <- function(named_sets) {
    stopifnot(is.list(named_sets), all(vapply(named_sets, is.character, logical(1))))

    all_elements <- unique(unlist(named_sets))
    overlap_matrix <- matrix(FALSE,
                             nrow = length(all_elements),
                             ncol = length(named_sets))
    rownames(overlap_matrix) <- all_elements
    colnames(overlap_matrix) <- names(named_sets)

    for (i in seq_along(named_sets)) {
        overlap_matrix[all_elements %in% named_sets[[i]], i] <- TRUE
    }

    intersect_category <- defineCategories(overlap_matrix)

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
#' - When provided with genomic regions, the function builds a non-redundant
#'   set of intervals (see `mode`), then determines which original sets each
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
#' @param mode Character string controlling how genomic intervals are made
#'   non-redundant before they are classified. One of:
#'   \itemize{
#'     \item `"reduce"` (default): all intervals from all sets are merged with
#'       `GenomicRanges::reduce()`, and each merged region is classified by the
#'       sets it overlaps. Region counts then correspond to merged loci.
#'     \item `"disjoin"`: each set is reduced on its own, and the union is then
#'       partitioned into non-overlapping segments with
#'       `GenomicRanges::disjoin()`. Every segment is covered by exactly one
#'       combination of sets, so a category such as `"111"` is reported only
#'       for positions genuinely shared by all three sets.
#'   }
#'   Ignored (with a warning) for non-genomic inputs.
#'
#' @return
#' An S3 object encoding the overlap result whose class depends on the input
#' type:
#'
#' \describe{
#'   \item{GenomicOverlapResult}{Returned when the input is genomic
#'       (`GRangesList` or list of `GRanges`). A list with:
#'       \itemize{
#'         \item \code{regions}: A `GRanges` object containing the
#'             non-redundant intervals, merged when `mode = "reduce"` and
#'             disjoint when `mode = "disjoin"`. Each region is annotated with
#'             an \code{intersect_category} column.
#'         \item \code{overlap_matrix}: A logical matrix indicating whether each
#'             region overlaps each input set (rows = regions,
#'             columns = sets).
#'         \item \code{mode}: The `mode` used to build the regions.
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
#' ## Choosing between `"reduce"` and `"disjoin"`
#'
#' The two modes answer different questions, and they differ whenever
#' overlaps are chained: an interval of A overlapping an interval of B, which
#' in turn overlaps an interval of C, without A and C sharing any position.
#'
#' `mode = "reduce"` merges such a chain into a single connected region and
#' labels it `"111"`, reporting a three-way intersection even though no base
#' pair is common to the three sets. This is the historical behavior, and it
#' is the natural one when the sets describe the same underlying features
#' (e.g., peaks called on replicates of the same experiment) and the question
#' is "which loci are shared?".
#'
#' `mode = "disjoin"` splits the chain at every set boundary, yielding one
#' `"110"` segment, one `"011"` segment, and the set-specific remainders. Use
#' it when categories must reflect genuinely shared genomic positions, which
#' is what an intersection is usually taken to mean. Note that the resulting
#' counts are counts of segments, not of input intervals, so a single long
#' interval may contribute to several categories.
#'
#' @examples
#' # Example with gene sets (built-in dataset)
#' data(gene_list)
#' ov_sets <- computeOverlaps(gene_list)
#' head(ov_sets$overlap_matrix)
#' plotVenn(ov_sets)
#'
#' # Example with genomic regions (built-in dataset)
#' data(a549_chipseq_peaks)
#' ov_gr <- computeOverlaps(a549_chipseq_peaks)
#' head(ov_gr$overlap_matrix)
#' plotVenn(ov_gr)
#'
#' # Chained overlaps: A-B and B-C overlap, but A and C do not
#' A <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
#' B <- GenomicRanges::GRanges("chr1", IRanges::IRanges(180, 300))
#' C <- GenomicRanges::GRanges("chr1", IRanges::IRanges(280, 400))
#'
#' # "reduce" merges the chain into a single "111" region
#' computeOverlaps(list(A = A, B = B, C = C))$regions
#'
#' # "disjoin" keeps A-B and B-C as separate two-way intersections
#' computeOverlaps(list(A = A, B = B, C = C), mode = "disjoin")$regions
#'
#' @seealso \code{\link{plotVenn}}, \code{\link{plotUpSet}},
#'     \code{\link[GenomicRanges]{GRangesList}},
#'     \code{\link[GenomicRanges]{reduce}},
#'     \code{\link[GenomicRanges]{disjoin}}
#'
#' @export
computeOverlaps <- function(x, mode = c("reduce", "disjoin")) {
    mode <- match.arg(mode)

    if (missing(x) || is.null(x)) {
        stop("'x' must be provided.", call. = FALSE)
    }

    # ---- direct GRangesList -------------------------------------------------
    if (methods::is(x, "GRangesList")) {
        if (is.null(names(x))) names(x) <- paste0("set", seq_along(x))
        return(computeGenomicOverlaps(x, mode = mode))
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
        return(computeGenomicOverlaps(x, mode = mode))
    }

    # atomic vectors (genes/ids). Allow numeric/factor but coerce to character.
    is_atomic_vec <- vapply(x, function(e) is.atomic(e) && !is.list(e), logical(1))
    if (all(is_atomic_vec)) {
        if (mode != "reduce") {
            warning("'mode' applies to genomic inputs only; ignored for sets ",
                    "of identifiers.", call. = FALSE)
        }
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
