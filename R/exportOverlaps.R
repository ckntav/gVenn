#' Export Overlap Groups to Excel
#'
#' This function exports the output of `extractOverlaps()` to an Excel file,
#' creating one sheet per overlap group. Genomic overlaps (`GRanges`) are
#' converted to data frames before export.
#'
#' @param grouped Overlap groups from `extractOverlaps()`.
#' @param output_dir A string specifying the output directory. Defaults to `"."`.
#' @param output_file A string specifying the base filename (without extension).
#' Defaults to `"overlap_groups"`.
#' @param with_date Logical (default `TRUE`). Whether to prepend the current
#' date (from `today`) to the filename.
#' @param verbose Logical. If `TRUE`, print a message with the saved path.
#' Default `TRUE`.
#'
#' @return Overlap groups are saved to a Excel file on disk. Invisibly returns
#' the full path to the saved file.
#'
#' @export
#'
#' @examples
#' res <- computeOverlaps(list(A = letters[1:3], B = letters[2:4]))
#' grouped <- extractOverlaps(res)
#' exportOverlaps(grouped, output_dir = tempdir(), output_file = "overlap_groups")
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

#' Export Overlap Groups to BED Files
#'
#' This function exports genomic overlap groups from `extractOverlaps()` to
#' BED format files, creating one BED file per overlap group.
#'
#' @details
#' This function only works with genomic overlaps (i.e., when the input to
#' `extractOverlaps()` was a `GenomicOverlapResult` object, resulting in a
#' `GRangesList`). It does not work with set overlaps (character vectors).
#' Each overlap group will be saved as a separate BED file with the group
#' identifier included in the filename.
#'
#' @param grouped Genomic overlap groups from `extractOverlaps()`
#'   (must be `GRangesList`).
#' @param output_dir A string specifying the output directory. Defaults to `"."`.
#' @param output_prefix A string specifying the filename prefix.
#'   Defaults to `"overlaps"`.
#' @param with_date Logical (default `TRUE`). Whether to prepend the current
#'   date to filenames.
#' @param verbose Logical. If `TRUE`, print messages. Default `TRUE`.
#'
#' @return Invisibly returns a character vector of file paths created.
#'
#' @export
exportOverlapsToBed <- function(grouped,
                                output_dir = ".",
                                output_prefix = "overlaps",
                                with_date = TRUE,
                                verbose = TRUE) {
    # Check if data is genomic
    if (!all(sapply(grouped, function(x) methods::is(x, "GRanges")))) {
        stop("BED format export only works with genomic overlaps (GRanges objects).")
    }

    stopifnot(requireNamespace("rtracklayer", quietly = TRUE))

    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    if (isTRUE(with_date)) {
        output_prefix <- paste0(today, "_", output_prefix)
    }

    filepaths <- character(length(grouped))
    for (i in seq_along(grouped)) {
        group_name <- names(grouped)[i]
        bed_file <- file.path(output_dir, paste0(output_prefix, "_", group_name, ".bed"))
        rtracklayer::export.bed(grouped[[i]], con = bed_file)
        filepaths[i] <- bed_file
    }

    if (isTRUE(verbose)) {
        message(" > ", length(filepaths), " BED files saved in ", output_dir)
    }

    invisible(filepaths)
}

