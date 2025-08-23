define_categories <- function(data) {
    categories <- apply(data, 1, function(row) {
        paste0(as.integer(row), collapse = "")
    })
    return(categories)
}

#' Compute Genomic Overlaps Across GRanges Sets
#'
#' Internal engine used by \code{\link{compute_overlaps}} when the input
#' consists of \code{GRanges} or a \code{GRangesList}. Not intended for direct
#' user calls.
#'
#' @param genomic_regions A \code{GRangesList} or a named list of \code{GRanges}
#'     objects. Each element represents one genomic region set (e.g., ChIP-seq
#'     peaks).
#'
#' @return An object of class \code{GenomicOverlapResult}, a list with:
#' \describe{
#'   \item{reduced_regions}{\code{GRanges} of the merged intervals across all
#'     sets, annotated with \code{intersect_category} (e.g., \code{"110"}).}
#'   \item{overlap_matrix}{Logical matrix; rows are reduced regions, columns are
#'     input sets.}
#' }
#'
#' @details Uses \code{GenomicRanges::reduce()} to merge regions, then
#' \code{IRanges::overlapsAny()} to build the overlap matrix. Prefer calling
#' \code{\link{compute_overlaps}} instead.
#'
#' @seealso \code{\link{compute_overlaps}}
#'
#' @keywords internal
#' @noRd
compute_genomic_overlaps <- function(genomic_regions) {
    if (inherits(genomic_regions, "list")) {
        genomic_regions <- GenomicRanges::GRangesList(genomic_regions)
    } else if (!inherits(genomic_regions, "GRangesList")) {
        stop("Input must be a list of GRanges or a GRangesList.")
    }

    reduced_regions <- GenomicRanges::reduce(unlist(genomic_regions))
    overlap_matrix <- matrix(
        FALSE,
        nrow = length(reduced_regions),
        ncol = length(genomic_regions)
    )

    for (i in seq_along(genomic_regions)) {
        overlap_matrix[, i] <- IRanges::overlapsAny(
            reduced_regions,
            genomic_regions[[i]]
        )
    }

    colnames(overlap_matrix) <- names(genomic_regions)

    intersect_category <- define_categories(overlap_matrix)
    GenomicRanges::mcols(reduced_regions)$intersect_category <- intersect_category

    res <- list(
        reduced_regions = reduced_regions,
        overlap_matrix  = overlap_matrix
    )
    class(res) <- "GenomicOverlapResult"
    return(res)
}

#' Compute Overlaps Between Named Atomic Sets
#'
#' Internal engine used by \code{\link{compute_overlaps}} when the input is a
#' named list of atomic vectors (e.g., gene symbols). Not intended for direct
#' user calls.
#'
#' @param named_sets A named list of atomic vectors (character, numeric, factor,
#'     etc.). Each vector contains the members of one set.
#'
#' @return An object of class \code{SetOverlapResult}, a list with:
#' \describe{
#'   \item{unique_elements}{Character vector of all unique elements.}
#'   \item{overlap_matrix}{Logical matrix; rows are elements, columns are sets.}
#'   \item{intersect_category}{Character vector encoding the overlap pattern per
#'     element (e.g., \code{"110"}, \code{"101"}).}
#' }
#'
#' @details Elements are coerced to \code{character} before computing the
#' presence/absence matrix. Prefer calling \code{\link{compute_overlaps}}
#' instead.
#'
#' @seealso \code{\link{compute_overlaps}}
#'
#' @keywords internal
#' @noRd
compute_set_overlaps <- function(named_sets) {
    stopifnot(
        is.list(named_sets),
        all(vapply(named_sets, is.atomic, logical(1))),
        !is.null(names(named_sets)) || length(named_sets) == 0L
    )

    all_elements <- unique(unlist(named_sets))
    overlap_matrix <- matrix(
        FALSE,
        nrow = length(all_elements),
        ncol = length(named_sets)
    )
    rownames(overlap_matrix) <- all_elements
    colnames(overlap_matrix) <- names(named_sets)

    for (i in seq_along(named_sets)) {
        overlap_matrix[all_elements %in% as.character(named_sets[[i]]), i] <- TRUE
    }

    intersect_category <- define_categories(overlap_matrix)

    res <- list(
        unique_elements    = all_elements,
        overlap_matrix     = overlap_matrix,
        intersect_category = intersect_category
    )
    class(res) <- "SetOverlapResult"
    return(res)
}

#' Dispatch overlap computation based on input type
#'
#' Wrapper that inspects the input and calls either
#' \code{\link{compute_genomic_overlaps}} (for GRanges/GRangesList)
#' or \code{\link{compute_set_overlaps}} (for atomic vectors).
#'
#' @param x A \code{GRangesList}, a named list of \code{GRanges}, or a named
#'     list of atomic vectors (e.g., gene symbols). All elements must be of the
#'     same kind.
#'
#' @return The object returned by the chosen backend:
#'     \itemize{
#'         \item \code{GenomicOverlapResult} from \code{compute_genomic_overlaps()}
#'         \item \code{SetOverlapResult} from \code{compute_set_overlaps()}
#'     }
#'
#' @examples
#' # Sets
#' sets <- list(A = letters[1:4], B = letters[3:6])
#' ov1 <- compute_overlaps(sets)
#'
#' # Genomic
#' if (requireNamespace("GenomicRanges", quietly = TRUE)) {
#'     gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 50), width = 20))
#'     gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(15, 90), width = 20))
#'     ov2 <- compute_overlaps(list(A = gr1, B = gr2))
#' }
#'
#' @seealso \code{\link[GenomicRanges]{GRangesList}},
#'     \code{\link[GenomicRanges]{reduce}},
#'     \code{\link[IRanges]{overlapsAny}},
#'     \code{\link{plotVenn}},
#'     \code{\link{plotUpset}}
#'
#' @export
compute_overlaps <- function(x) {
    if (missing(x) || is.null(x)) {
        stop("'x' must be provided.", call. = FALSE)
    }

    # ---- direct GRangesList -------------------------------------------------
    if (methods::is(x, "GRangesList")) {
        if (is.null(names(x))) names(x) <- paste0("set", seq_along(x))
        return(compute_genomic_overlaps(x))
    }

    # ---- list inputs --------------------------------------------------------
    if (!is.list(x)) {
        stop(
            "'x' must be a GRangesList, a list of GRanges, or a list of atomic vectors.",
            call. = FALSE
        )
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
        return(compute_genomic_overlaps(x))
    }

    # atomic vectors (genes/ids). Allow numeric/factor but coerce to character.
    is_atomic_vec <- vapply(x, function(e) is.atomic(e) && !is.list(e), logical(1))
    if (all(is_atomic_vec)) {
        x_chr <- lapply(x, function(e) as.character(e))
        return(compute_set_overlaps(x_chr))
    }

    # mixed or unsupported
    bad_idx <- which(!(is_gr | is_atomic_vec))[1]
    bad_cls <- paste(class(x[[bad_idx]]), collapse = "/")
    stop(
        "All elements of 'x' must be the same type: either all GRanges or all ",
        "atomic vectors. Found unsupported element of class: ", bad_cls, ".",
        call. = FALSE
    )
}
