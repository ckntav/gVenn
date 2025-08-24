test_that("plotVenn() rejects unsupported input classes", {
    bad <- list(overlap_matrix = matrix(TRUE, nrow = 2, ncol = 2))
    expect_error(
        plotVenn(bad),
        "GenomicOverlapResult or SetOverlapResult"
    )
})

# Helpers to fabricate minimal overlap objects without relying on computeOverlaps()
make_set_overlap <- function() {
    # Two sets, four elements => simple logical matrix works for eulerr
    om <- matrix(
        c(
            TRUE,  FALSE,
            TRUE,  TRUE,
            FALSE, TRUE,
            FALSE, FALSE
        ),
        ncol = 2, byrow = TRUE
    )
    colnames(om) <- c("A", "B")
    obj <- list(overlap_matrix = om)
    class(obj) <- "SetOverlapResult"
    obj
}

make_genomic_overlap <- function() {
    # Three sets; each row is an element; columns = set membership
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

test_that("plotVenn() returns a grid grob for SetOverlapResult", {
    skip_if_not_installed("eulerr")
    obj <- make_set_overlap()
    p <- plotVenn(obj)
    expect_true(inherits(p, "eulergram"))     # class defined by eulerr
    expect_true(inherits(p, "grob"))          # grid object
})

test_that("plotVenn() returns a grid grob for GenomicOverlapResult", {
    skip_if_not_installed("eulerr")
    obj <- make_genomic_overlap()
    p <- plotVenn(obj)
    expect_true(inherits(p, "eulergram"))
    expect_true(inherits(p, "grob"))
})

test_that("plotVenn() accepts legend = FALSE and custom fills", {
    obj <- make_genomic_overlap()
    # supply exactly three fills to match three sets
    custom_fills <- c("#2B70AB", "#FFB027", "#3EA742")
    expect_no_error(plotVenn(obj, legend = FALSE, fill = custom_fills))
})

test_that("plotVenn() works with labels = TRUE", {
    obj <- make_set_overlap()
    expect_no_error(plotVenn(obj, labels = TRUE))
})
