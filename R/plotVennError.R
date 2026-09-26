#' Plot where a Venn diagram misrepresents the data
#'
#' Redraws the Venn diagram produced by \code{\link{plotVenn}}, shaded by the
#' signed error of each region, so that the regions the diagram gets wrong can
#' be read off the picture rather than inferred from a single summary
#' statistic. This is a thin wrapper around \code{\link[eulerr]{error_plot}}.
#'
#' Where `diagError` reports how far off the worst region is, this plot shows
#' which regions are off and in which direction, with each region's
#' `regionError` printed inside it. Use it when `diagError` exceeds 1e-6.
#'
#' @param venn A plot returned by \code{\link{plotVenn}}, carrying the
#'   `"euler_fit"` attribute that function attaches. Note that this is the
#'   `plotVenn()` output, not the `GenomicOverlapResult` or
#'   `SetOverlapResult` that `plotVenn()` itself takes.
#' @param ... Additional arguments passed to \code{\link[eulerr]{error_plot}}.
#'
#' @return The plot returned by \code{\link[eulerr]{error_plot}}, a `gTree`
#'   that can be drawn, arranged with other grobs, or written out with
#'   \code{\link{saveViz}}.
#' @export
#'
#' @seealso \code{\link{plotVenn}} for the diagram itself and its fit
#'   diagnostics, \code{\link{plotUpSet}} for a representation that does not
#'   approximate any count by an area.
#'
#' @examples
#' data(gene_list)
#' res_sets <- computeOverlaps(gene_list)
#'
#' # Where does the diagram misrepresent the counts?
#' venn <- plotVenn(res_sets)
#' plotVennError(venn)
plotVennError <- function(venn, ...) {

    if (inherits(venn, "GenomicOverlapResult") ||
        inherits(venn, "SetOverlapResult")) {
        stop("`venn` must be the plot returned by plotVenn(), not the overlap ",
             "object plotVenn() takes. plotVennError() reuses the fit stored ",
             "on that plot rather than computing a new one, so that it shades ",
             "the diagram you were shown. Use:\n",
             "    venn <- plotVenn(x)\n",
             "    plotVennError(venn)")
    }

    fit <- attr(venn, "euler_fit")

    if (!inherits(fit, "euler")) {
        stop("`venn` does not carry the eulerr fit plotVennError() needs. It ",
             "must be the plot returned by plotVenn(), with its \"euler_fit\" ",
             "attribute intact. Operations that rebuild or strip attributes ",
             "can drop it, in which case call plotVenn() again and pass its ",
             "return value straight to plotVennError().")
    }

    eulerr::error_plot(fit, ...)
}
