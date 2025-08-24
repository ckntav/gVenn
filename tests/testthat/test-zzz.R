test_that("get_today returns YYYYMMDD string", {
    val <- get_today()
    expect_type(val, "character")
    expect_equal(nchar(val), 8L)
    expect_match(val, "^[0-9]{8}$")
})

test_that("exported `today` looks like YYYYMMDD at runtime", {
    pkgname <- "gVenn"
    exported <- getExportedValue(pkgname, "today")

    expect_type(exported, "character")
    expect_equal(length(exported), 1L)
    expect_equal(nchar(exported), 8L)
    expect_match(exported, "^[0-9]{8}$")
})
