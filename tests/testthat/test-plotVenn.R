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

test_that("plotVenn() attaches eulerr fit diagnostics", {
    skip_if_not_installed("eulerr")
    obj <- make_genomic_overlap()
    p <- suppressMessages(plotVenn(obj))
    diag <- attr(p, "fit_diagnostics")
    expect_type(diag, "list")
    expect_named(diag, c("stress", "diagError", "regionError", "undrawnRegions"))
    expect_true(is.numeric(diag$stress) && length(diag$stress) == 1)
    expect_true(is.numeric(diag$diagError) && length(diag$diagError) == 1)
    expect_type(diag$undrawnRegions, "character")
})

test_that("plotVenn() reports regions the diagram gives no area to", {
    skip_if_not_installed("eulerr")
    # Large exclusive sets with small intersections: eulerr erases several
    # regions here, yet diagError stays below a thousandth, so a reader
    # watching the statistic alone would not notice the omission
    counts <- c(A = 19689, B = 17691, C = 36466, D = 7210,
                "A&B" = 222, "A&C" = 293, "A&D" = 264, "B&C" = 36,
                "B&D" = 164, "C&D" = 43, "A&B&C" = 1, "A&B&D" = 30,
                "A&C&D" = 30, "B&C&D" = 26, "A&B&C&D" = 1)
    ids <- paste0("id", seq_len(sum(counts)))
    nms <- unique(unlist(strsplit(names(counts), "&", fixed = TRUE)))
    sets <- stats::setNames(lapply(nms, function(x) character(0)), nms)
    i <- 1L
    for (r in names(counts)) {
        members <- strsplit(r, "&", fixed = TRUE)[[1]]
        block <- ids[i:(i + counts[[r]] - 1L)]
        for (s in members) sets[[s]] <- c(sets[[s]], block)
        i <- i + counts[[r]]
    }
    obj <- computeOverlaps(sets)

    set.seed(1111)
    p <- suppressMessages(plotVenn(obj))
    diag <- attr(p, "fit_diagnostics")

    expect_gt(length(diag$undrawnRegions), 0)
    expect_lt(diag$diagError, 0.001)   # small, yet regions have been erased

    # and the message names them rather than reporting diagError alone
    set.seed(1111)
    expect_message(plotVenn(obj), "regions are given no area at all")
    set.seed(1111)
    expect_message(plotVenn(obj), "> 1e-06")
})

test_that("plotVenn() messages fit diagnostics when verbose = TRUE (default)", {
    obj <- make_set_overlap()
    expect_message(plotVenn(obj), "Venn diagError")
})

test_that("plotVenn() silences fit diagnostics message when verbose = FALSE", {
    obj <- make_set_overlap()
    expect_no_message(plotVenn(obj, verbose = FALSE))
})

test_that("plotVenn() fits two sets with circles, exactly and reproducibly", {
    skip_if_not_installed("eulerr")
    # Two circles represent any two-set data exactly. Fitting ellipses here adds
    # two free parameters that buy nothing and open local minima: seed 7 used to
    # land in one, reporting diagError = 2.2e-01 and erasing A&B entirely, while
    # neighbouring seeds fitted the same data to 1e-12.
    two <- list(A = paste0("g", 1:500), B = paste0("g", 300:900))
    obj <- computeOverlaps(two)

    diag_errors <- vapply(1:12, function(s) {
        set.seed(s)
        attr(suppressMessages(plotVenn(obj)), "fit_diagnostics")$diagError
    }, numeric(1))
    expect_true(all(diag_errors <= 1e-6))

    set.seed(7)
    diag <- attr(suppressMessages(plotVenn(obj)), "fit_diagnostics")
    expect_length(diag$undrawnRegions, 0)
})
