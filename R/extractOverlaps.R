#' Extract Overlap Groups from Genomic or Set Overlap Results
#'
#' This function extracts the elements of each overlap group, identified by
#' its binary code (e.g., `"110"`). For genomic overlaps, the elements are the
#' reduced or disjoint regions built by \code{\link{computeOverlaps}}, not the
#' input intervals.
#' For genomic overlaps, it returns a `GRangesList`; for set overlaps, it
#' returns a named list of character vectors.
#'
#' @param overlap_object A `GenomicOverlapResult` or `SetOverlapResult`
#' object.
#'
#' @return A named list with one entry per overlap group, named by its binary
#' code prefixed with `"group_"` (e.g., `"group_110"`) and ordered by the
#' number of sets involved:
#' - If input is a `GenomicOverlapResult`, a `GRangesList` split
#' by `intersect_category`.
#' - If input is a `SetOverlapResult`, a named `list` of character vectors
#' grouped by `intersect_category`.
#'
#' @export
#'
#' @examples
#' # Example with gene sets (built-in dataset)
#' data(gene_list)
#' res_sets <- computeOverlaps(gene_list)
#' group_gene <- extractOverlaps(res_sets)
#' group_gene
#'
#' # Example with genomic regions (built-in dataset)
#' data(a549_chipseq_peaks)
#' res_genomic <- computeOverlaps(a549_chipseq_peaks)
#' group_genomic <- extractOverlaps(res_genomic)
#' group_genomic
extractOverlaps <- function(overlap_object) {
    if (inherits(overlap_object, "GenomicOverlapResult")) {
        regions <- overlap_object[["regions"]]
        categories <- unique(regions$intersect_category)
        sorted_categories <- categories[order(stringr::str_count(categories, "1"))]

        intersect_regions_list <- GenomicRanges::GRangesList()
        for (category in sorted_categories) {
            intersect_regions_list[[paste0("group_", category)]] <- regions[regions$intersect_category == category]
        }

        return(intersect_regions_list)

    } else if (inherits(overlap_object, "SetOverlapResult")) {
        elements <- overlap_object[["unique_elements"]]
        categories <- overlap_object[["intersect_category"]]
        sorted_categories <- sort(unique(categories), decreasing = FALSE)
        sorted_categories <- sorted_categories[order(stringr::str_count(sorted_categories, "1"))]

        grouped <- split(elements, categories)
        grouped <- grouped[sorted_categories]  # reorder

        names(grouped) <- paste0("group_", names(grouped))
        return(grouped)

    } else {
        stop("Input must be a GenomicOverlapResult or SetOverlapResult object.")
    }
}
