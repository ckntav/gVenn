test_that("extractOverlaps() groups and orders SetOverlapResult correctly", {
    # Build a minimal fake SetOverlapResult with the fields extractOverlaps() expects
    elements <- c("a", "b", "c", "d", "e")
    categories <- c("100", "110", "010", "011", "001")
    # order inside extractOverlaps(): sort(unique(categories)) then by count of '1'
    # sort(unique) -> "001","010","011","100","110"
    # ones count -> 1,1,2,1,2  => "001","010","100","011","110"

    set_obj <- structure(
        list(unique_elements = elements,
             intersect_category = categories),
        class = "SetOverlapResult"
    )

    grouped <- extractOverlaps(set_obj)
    expect_type(grouped, "list")
    expect_true(all(startsWith(names(grouped), "group_")))
    expect_identical(
        names(grouped),
        c("group_001", "group_010", "group_100", "group_011", "group_110")
    )

    # Check contents by category
    expect_identical(grouped[["group_100"]], "a")
    expect_identical(sort(grouped[["group_011"]]), sort(c("d")))  # only "d" has 011
    expect_setequal(unlist(grouped), elements)
})

test_that("extractOverlaps() groups and orders GenomicOverlapResult correctly", {
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")

    # Minimal GRanges with an 'intersect_category' mcol
    gr <- GenomicRanges::GRanges(
        seqnames = rep("chr1", 3),
        ranges = IRanges::IRanges(c(1, 101, 201), width = 50)
    )
    gr$intersect_category <- c("001", "110", "010")
    # unique(categories) preserved input order: "001","110","010"
    # ordered by ones count => "001"(1), "010"(1), "110"(2)

    genomic_obj <- structure(
        list(reduced_regions = gr),
        class = "GenomicOverlapResult"
    )

    grouped <- extractOverlaps(genomic_obj)
    expect_s4_class(grouped, "GRangesList")
    expect_identical(names(grouped), c("group_001", "group_010", "group_110"))

    # Each group contains only the expected ranges
    expect_equal(length(grouped[["group_001"]]), 1L)
    expect_equal(length(grouped[["group_010"]]), 1L)
    expect_equal(length(grouped[["group_110"]]), 1L)
    expect_identical(grouped[["group_001"]]$intersect_category, "001")
    expect_identical(grouped[["group_010"]]$intersect_category, "010")
    expect_identical(grouped[["group_110"]]$intersect_category, "110")
})

test_that("extractOverlaps() errors on wrong input class", {
    expect_error(extractOverlaps(list(foo = 1)), "GenomicOverlapResult|SetOverlapResult")
})

test_that("exportOverlaps() writes an .xlsx and returns the path (Set lists)", {
    skip_on_cran()
    skip_if_not_installed("writexl")

    tmpdir <- withr::local_tempdir()
    grouped <- list(
        group_001 = c("x", "y"),
        group_010 = c("z"),
        group_110 = c("u", "v", "w")
    )

    path <- exportOverlaps(grouped,
                           output_dir = tmpdir,
                           output_file = "groups_set",
                           with_date = FALSE,
                           verbose = FALSE)
    expect_true(file.exists(path))
    expect_match(basename(path), "^groups_set\\.xlsx$")
})

test_that("exportOverlaps() emits a message when verbose = TRUE", {
    skip_on_cran()
    skip_if_not_installed("writexl")

    tmpdir <- withr::local_tempdir()
    grouped <- list(group_001 = letters[1:2])

    expect_message(
        exportOverlaps(grouped,
                       output_dir = tmpdir,
                       output_file = "msg_groups",
                       with_date = FALSE,
                       verbose = TRUE),
        "Overlap groups saved"
    )
})
