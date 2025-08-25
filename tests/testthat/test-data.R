test_that("a549_chipseq_peaks dataset is valid", {
    data("a549_chipseq_peaks", package = "gVenn")
    expect_s4_class(a549_chipseq_peaks, "GRangesList")
    expect_true(length(a549_chipseq_peaks) == 3)
})

test_that("gene_list dataset is valid", {
    data("gene_list", package = "gVenn")
    expect_type(gene_list, "list")
    expect_true(length(gene_list) == 3)
    expect_true(all(nzchar(gene_list))) # no empty strings
})
