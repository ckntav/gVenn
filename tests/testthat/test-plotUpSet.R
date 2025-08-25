test_that("plotUpSet() works with SetOverlapResult", {
    skip_on_cran()
    skip_if_not_installed("ComplexHeatmap")

    # tiny, deterministic toy data
    gene_sets <- list(
        TF1 = c("A", "B", "C"),
        TF2 = c("B", "C", "D"),
        TF3 = c("C", "E")
    )

    res_sets <- computeOverlaps(gene_sets)
    expect_true(inherits(res_sets, "SetOverlapResult"))

    upset <- plotUpSet(res_sets)
    # ComplexHeatmap can return a Heatmap or HeatmapList depending on options
    expect_true(inherits(upset, "Heatmap") || inherits(upset, "HeatmapList"))

    # If it is a single Heatmap, check that annotations were attached
    if (methods::is(upset, "Heatmap")) {
        # these are S4 slots on Heatmap objects
        expect_true(methods::is(upset@top_annotation, "HeatmapAnnotation"))
        expect_true(methods::is(upset@right_annotation, "HeatmapAnnotation"))
    }

    # Sanity-check the combination logic on the same overlap matrix
    # (this doesn't reach *inside* the plot object, but confirms inputs behave)
    combMat <- ComplexHeatmap::make_comb_mat(res_sets[["overlap_matrix"]])
    # Known intersections for the toy example
    # sets: TF1, TF2, TF3  (order as produced by make_comb_mat)
    # C is in all three -> one triple
    # B in TF1 & TF2 only
    # (A only TF1), (D only TF2), (E only TF3)
    cs <- ComplexHeatmap::comb_size(combMat)
    # total elements across unique combos should be 5
    expect_equal(sum(cs), 5L)
    # must contain at least one triple intersection
    expect_true(any(ComplexHeatmap::comb_degree(combMat) == 3))
})

test_that("plotUpSet() works with GenomicOverlapResult", {
    skip_on_cran()
    skip_if_not_installed("ComplexHeatmap")
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")

    gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 500, 900), width = 100))
    gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(150, 520, 1000), width = 100))
    gr3 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(300, 550, 1100), width = 100))
    peak_sets <- list(H3K27ac = gr1, MED1 = gr2, BRD4 = gr3)

    res_genomic <- computeOverlaps(peak_sets)
    expect_true(inherits(res_genomic, "GenomicOverlapResult"))

    upset <- plotUpSet(res_genomic)
    expect_true(inherits(upset, "Heatmap") || inherits(upset, "HeatmapList"))

    if (methods::is(upset, "Heatmap")) {
        expect_true(methods::is(upset@top_annotation, "HeatmapAnnotation"))
        expect_true(methods::is(upset@right_annotation, "HeatmapAnnotation"))
    }

    # Input matrix should be consumable by ComplexHeatmap::make_comb_mat
    combMat <- ComplexHeatmap::make_comb_mat(res_genomic[["overlap_matrix"]])
    # At least one intersection must exist
    expect_true(sum(ComplexHeatmap::comb_size(combMat)) >= length(unique(IRanges::start(IRanges::reduce(IRanges::IRanges(c(100, 150, 300, 500, 520, 550, 900, 1000, 1100), width = 100))))))
})

test_that("plotUpSet() errors clearly on wrong input", {
    expect_error(plotUpSet(list(a = 1)), "GenomicOverlapResult|SetOverlapResult")
})
