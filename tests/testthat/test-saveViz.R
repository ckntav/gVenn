test_that("saveViz() saves a ggplot to PDF/PNG/SVG and returns the filepath", {
    skip_on_cran()
    skip_if_not_installed("ggplot2")
    tmpdir <- withr::local_tempdir()

    p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()

    # PDF
    pdf_path <- saveViz(p,
                        output_dir = tmpdir,
                        output_file = "p_scatter",
                        format = "pdf",
                        with_date = FALSE,
                        verbose = FALSE)
    expect_true(file.exists(pdf_path))
    expect_match(pdf_path, "\\.pdf$")

    # PNG
    png_path <- saveViz(p,
                        output_dir = tmpdir,
                        output_file = "p_scatter_png",
                        format = "png",
                        with_date = FALSE,
                        height = 3,
                        verbose = FALSE)
    expect_true(file.exists(png_path))
    expect_match(png_path, "\\.png$")

    # SVG
    svg_path <- saveViz(p,
                        output_dir = tmpdir,
                        output_file = "p_scatter_svg",
                        format = "svg",
                        with_date = FALSE,
                        verbose = FALSE)
    expect_true(file.exists(svg_path))
    expect_match(svg_path, "\\.svg$")
})

test_that("saveViz() creates output directory recursively", {
    skip_on_cran()
    skip_if_not_installed("ggplot2")
    base <- withr::local_tempdir()
    nested_dir <- file.path(base, "a", "very", "deep", "folder")

    p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()

    out_path <- saveViz(p,
                        output_dir = nested_dir,
                        output_file = "deep_plot",
                        format = "pdf",
                        with_date = FALSE,
                        verbose = FALSE)

    expect_true(dir.exists(nested_dir))
    expect_true(file.exists(out_path))
})

test_that("saveViz() prints a message when verbose = TRUE, silent when FALSE", {
    skip_on_cran()
    skip_if_not_installed("ggplot2")
    tmpdir <- withr::local_tempdir()
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()

    expect_message(
        saveViz(p, output_dir = tmpdir, output_file = "msg_plot", format = "pdf",
                with_date = FALSE, verbose = TRUE),
        "Visualization \\(pdf\\) saved"
    )

    expect_silent(
        saveViz(p, output_dir = tmpdir, output_file = "silent_plot", format = "pdf",
                with_date = FALSE, verbose = FALSE)
    )
})

test_that("saveViz() errors on unsupported format", {
    skip_on_cran()
    skip_if_not_installed("ggplot2")
    tmpdir <- withr::local_tempdir()
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(mpg, wt)) + ggplot2::geom_point()

    expect_error(
        saveViz(p, output_dir = tmpdir, output_file = "badfmt", format = "jpeg",
                with_date = FALSE, verbose = FALSE),
        "should be one of"
    )
})
