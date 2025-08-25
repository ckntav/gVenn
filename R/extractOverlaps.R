#' Count the Number of "1"s in Strings
#'
#' This function counts how many times the character `"1"` appears
#' in each element of a character vector. It is mainly used internally
#' when working with binary category codes (e.g., `"101"`, `"110"`)
#' that represent overlap patterns across sets.
#'
#' @param x A character vector where each element is a string
#'   potentially containing `"1"` characters.
#'
#' @return An integer vector of the same length as `x`,
#'   giving the count of `"1"` in each string.
#'
#' @examples
#' countOnes(c("101", "110", "000", "111"))
#' # Returns: 2 2 0 3
#'
#' @seealso [stringr::str_count()]
#'
#' @keywords internal
countOnes <- function(x) {
    stringr::str_count(x, "1")
}

#' Extract Overlap Groups from Genomic or Set Overlap Results
#'
#' This function extracts subsets of intersecting elements grouped by their overlap category (e.g., "110").
#' For genomic overlaps, it returns a `GRangesList`; for set overlaps, it returns a named list of character vectors.
#'
#' @param overlap_object A `GenomicOverlapsResult` or `SetOverlapsResult` object.
#'
#' @return A named list of grouped intersecting elements:
#' - If input is a `GenomicOverlapsResult`, a `GRangesList` split by `intersect_category`.
#' - If input is a `SetOverlapsResult`, a named `list` of character vectors grouped by `intersect_category`.
#'
#' @export
#'
#' @examples
#' # For GenomicOverlapsResult
#' gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 500), width = 100))
#' gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(150, 700), width = 100))
#' gr3 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(900), width = 100))
#' peak_sets <- list(H3K27ac = gr1, MED1 = gr2, BRD4 = gr3)
#' res1 <- computeOverlaps(peak_sets)
#' extractOverlaps(res1)
#'
#' # For SetOverlapsResult
#' gene_sets <- list(
#'   TF1 = c("TP53", "BRCA1", "MYC", "GATA3", "FOXA1"),
#'   TF2 = c("MYC", "ESR1", "GATA3", "FOXA1", "AR"),
#'   TF3 = c("TP53", "GATA3", "AR", "NANOG", "SOX2")
#' )
#' res2 <- computeOverlaps(gene_sets)
#' extractOverlaps(res2)
extractOverlaps <- function(overlap_object) {
    if (inherits(overlap_object, "GenomicOverlapResult")) {
        reduced_regions <- overlap_object[["reduced_regions"]]
        categories <- unique(reduced_regions$intersect_category)
        sorted_categories <- categories[order(countOnes(categories))]

        intersect_regions_list <- GenomicRanges::GRangesList()
        for (category in sorted_categories) {
            intersect_regions_list[[paste0("group_", category)]] <- reduced_regions[reduced_regions$intersect_category == category]
        }

        return(intersect_regions_list)

    } else if (inherits(overlap_object, "SetOverlapResult")) {
        elements <- overlap_object[["unique_elements"]]
        categories <- overlap_object[["intersect_category"]]
        sorted_categories <- sort(unique(categories), decreasing = FALSE)
        sorted_categories <- sorted_categories[order(countOnes(sorted_categories))]

        grouped <- split(elements, categories)
        grouped <- grouped[sorted_categories]  # reorder

        names(grouped) <- paste0("group_", names(grouped))
        return(grouped)

    } else {
        stop("Input must be a GenomicOverlapResult or SetOverlapResult object.")
    }
}

#' Export Overlap Groups to Excel
#'
#' This function exports the output of `extractOverlaps()` to an Excel file,
#' creating one sheet per overlap group. Genomic overlaps (`GRanges`) are converted
#' to data frames before export.
#'
#' @param grouped Overlap groups from `extractOverlaps()`.
#' @param output_dir A string specifying the output directory. Defaults to `"."`.
#' @param output_file A string specifying the base filename (without extension). Defaults to `"overlap_groups"`.
#' @param with_date Logical (default `TRUE`). Whether to prepend the current date (from `today`) to the filename.
#' @param verbose Logical. If `TRUE`, print a message with the saved path. Default `TRUE`.
#'
#' @return Overlap groups are saved to a Excel file on disk. Invisibly returns the full path to the saved file.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'     res <- computeOverlaps(list(A = letters[1:3], B = letters[2:4]))
#'     grouped <- extractOverlaps(res)
#'     export_overlap_groups(grouped, "overlap_groups.xlsx")
#'     }
exportOverlaps <- function(grouped,
                           output_dir = ".",
                           output_file = "overlap_groups",
                           with_date = TRUE,
                           verbose = TRUE) {
    stopifnot(requireNamespace("writexl", quietly = TRUE))

    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    if (isTRUE(with_date)) {
        output_file <- paste0(today, "_", output_file)
    }

    filepath <- file.path(output_dir, paste0(output_file, ".xlsx"))

    # ensure sheet names are valid (Excel ≤ 31 chars, unique)
    sheet_names <- substr(names(grouped), 1, 31)
    sheet_names <- make.unique(sheet_names)

    # convert GRanges to data.frame if needed
    to_write <- lapply(grouped, function(x) {
        if (methods::is(x, "GRanges")) {
            as.data.frame(x)
        } else {
            data.frame(element = x)
        }
    })
    names(to_write) <- sheet_names

    writexl::write_xlsx(to_write, path = filepath)
    if (isTRUE(verbose)) {
        message(" > Overlap groups saved in ", filepath)
    }
    invisible(filepath)
}
