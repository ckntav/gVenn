test_that("exportOverlapsToBed creates BED files for each group", {
    skip_if_not_installed("rtracklayer")

    temp_dir <- tempdir()
    gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 150))
    gr2 <- GenomicRanges::GRanges("chr2", IRanges::IRanges(200, 250))

    grouped <- list(group_10 = gr1, group_11 = gr2)

    filepaths <- exportOverlapsToBed(
        grouped,
        output_dir = temp_dir,
        with_date = FALSE,
        verbose = FALSE
    )

    expect_length(filepaths, 2)
    expect_true(all(file.exists(filepaths)))
    expect_match(basename(filepaths[1]), "group_10\\.bed$")
    expect_match(basename(filepaths[2]), "group_11\\.bed$")

    unlink(filepaths)
})

test_that("exportOverlapsToBed throws error for non-GRanges input", {
    skip_if_not_installed("rtracklayer")

    temp_dir <- tempdir()
    grouped <- list(group_10 = c("gene1", "gene2"))

    expect_error(
        exportOverlapsToBed(grouped, output_dir = temp_dir, verbose = FALSE),
        "BED format export only works with genomic overlaps"
    )
})

test_that("exportOverlapsToBed handles mixed input (should fail)", {
    skip_if_not_installed("rtracklayer")

    temp_dir <- tempdir()
    gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 150))

    grouped <- list(group_10 = gr1, group_11 = c("gene1", "gene2"))

    expect_error(
        exportOverlapsToBed(grouped, output_dir = temp_dir, verbose = FALSE),
        "BED format export only works with genomic overlaps"
    )
})

test_that("exportOverlapsToBed includes date prefix when requested", {
    skip_if_not_installed("rtracklayer")

    temp_dir <- tempdir()
    gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 150))
    grouped <- list(group_10 = gr1)

    filepaths <- exportOverlapsToBed(
        grouped,
        output_dir = temp_dir,
        output_prefix = "test",
        with_date = TRUE,
        verbose = FALSE
    )

    expect_match(basename(filepaths[1]), ".*_test_group_10\\.bed$")
    unlink(filepaths)
})

test_that("exportOverlapsToBed creates output directory", {
    skip_if_not_installed("rtracklayer")

    temp_dir <- file.path(tempdir(), "new_bed_dir")
    on.exit(unlink(temp_dir, recursive = TRUE))

    gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 150))
    grouped <- list(group_10 = gr1)

    expect_false(dir.exists(temp_dir))
    exportOverlapsToBed(grouped, output_dir = temp_dir,
                        with_date = FALSE, verbose = FALSE)
    expect_true(dir.exists(temp_dir))
})

test_that("exportOverlapsToBed verbose messaging works", {
    skip_if_not_installed("rtracklayer")

    temp_dir <- tempdir()
    gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 150))
    grouped <- list(group_10 = gr1)

    expect_message(
        exportOverlapsToBed(grouped, output_dir = temp_dir,
                            with_date = FALSE, verbose = TRUE),
        "BED files saved in"
    )

    expect_silent(
        exportOverlapsToBed(grouped, output_dir = temp_dir,
                            with_date = FALSE, verbose = FALSE)
    )
})

test_that("exportOverlapsToBed BED files are valid", {
    skip_if_not_installed("rtracklayer")

    temp_dir <- tempdir()
    gr1 <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr2"),
        ranges = IRanges::IRanges(start = c(100, 200), width = 50)
    )
    grouped <- list(group_10 = gr1)

    filepaths <- exportOverlapsToBed(grouped, output_dir = temp_dir,
                                     with_date = FALSE, verbose = FALSE)

    # Read back and verify
    reimported <- rtracklayer::import.bed(filepaths[1])
    expect_s4_class(reimported, "GRanges")
    expect_length(reimported, 2)

    unlink(filepaths)
})
