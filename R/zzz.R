#' Today's Date at Package Load Time
#'
#' This variable stores the current date (in "yyyymmdd" format) at the time the
#' package is loaded. It is useful for reproducible filenames (e.g., in
#' `saveViz()`), and is automatically set when the package is attached.
#'
#' @format A character string (e.g., "20250624").
#' @export
#'
#' @examples
#' # Print the date stored at package load
#' library(gVenn)
#' today
#'
#' # Use it in a filename
#' paste0("venn_plot_", today, ".pdf")
today <- NULL

get_today <- function() {
    gsub("-", "", as.character(lubridate::today()))
}

.onLoad <- function(libname, pkgname) {
    # Assign to the namespace (for internal access)
    ns <- asNamespace(pkgname)
    assign("today", get_today(), envir = ns)

    # Assign to the package environment (for user access)
    assign("today", get_today(), envir = parent.env(environment()))
}
