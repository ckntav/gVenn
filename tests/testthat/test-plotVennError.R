# Minimal overlap objects, mirroring the helpers in test-plotVenn.R
make_genomic_overlap <- function() {
    om <- matrix(
        c(
            1, 0, 0,
            1, 1, 0,
            0, 1, 1,
            1, 0, 1,
            0, 0, 1
        ),
        ncol = 3, byrow = TRUE
    )
    colnames(om) <- c("H3K27ac", "MED1", "BRD4")
    obj <- list(overlap_matrix = om)
    class(obj) <- "GenomicOverlapResult"
    obj
}

test_that("plotVenn() attaches the eulerr fit itself", {
    skip_if_not_installed("eulerr")
    obj <- make_genomic_overlap()
    p <- suppressMessages(plotVenn(obj))
    fit <- attr(p, "euler_fit")
    expect_true(inherits(fit, "euler"))
    # the diagnostics must describe this very fit, not a separate one
    diag <- attr(p, "fit_diagnostics")
    expect_identical(diag$diagError, fit$diagError)
    expect_identical(diag$stress, fit$stress)
})

test_that("plotVennError() returns a grid grob", {
    skip_if_not_installed("eulerr")
    obj <- make_genomic_overlap()
    p <- suppressMessages(plotVenn(obj))
    ep <- plotVennError(p)
    expect_true(inherits(ep, "eulergram"))
    expect_true(inherits(ep, "grob"))
})

test_that("plotVennError() reuses plotVenn()'s fit instead of refitting", {
    skip_if_not_installed("eulerr")
    obj <- make_genomic_overlap()
    p <- suppressMessages(plotVenn(obj))

    # eulerr::euler() seeds its optimiser from R's RNG, so a refit can place
    # the diagram differently. If plotVennError() refitted, advancing the RNG
    # between calls could change the picture; reusing the stored fit cannot.
    set.seed(1)
    a <- plotVennError(p)
    set.seed(99)
    runif(10)
    b <- plotVennError(p)

    # grid stamps each grob with an auto-incrementing name, so compare the
    # drawn coordinates rather than the trees as a whole
    geom <- function(g) {
        panel <- g$children[["canvas.grob"]]$children[["diagram.grob.1"]]
        shapes <- panel$children[grepl("^(fills\\.grob|edges\\.grob)", 
                                       names(panel$children))]
        lapply(shapes, function(x) list(x = as.numeric(x$x), y = as.numeric(x$y)))
    }

    ga <- geom(a)
    # guard against the comparison below passing on two empty lists
    expect_gt(length(ga), 0)
    expect_true(all(vapply(ga, function(s) length(s$x) > 0, logical(1))))

    expect_equal(ga, geom(b))
})

test_that("plotVennError() rejects the overlap object plotVenn() takes", {
    skip_if_not_installed("eulerr")
    expect_error(
        plotVennError(make_genomic_overlap()),
        "must be the plot returned by plotVenn"
    )
})

test_that("plotVennError() rejects input without the stored fit", {
    skip_if_not_installed("eulerr")
    p <- suppressMessages(plotVenn(make_genomic_overlap()))
    attr(p, "euler_fit") <- NULL
    expect_error(plotVennError(p), "does not carry the eulerr fit")
})
