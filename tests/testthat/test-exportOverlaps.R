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
