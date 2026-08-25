test_that("defineCategories encodes rows to 0/1 strings", {
    m <- matrix(c(TRUE, FALSE, TRUE,
                  TRUE,  TRUE,  FALSE),
                nrow = 2, byrow = TRUE)
    out <- defineCategories(m)
    expect_type(out, "character")
    expect_identical(out, c("101", "110"))

    m2 <- matrix(c(1, 0, 1,
                   0, 0, 0),
                 nrow = 2, byrow = TRUE)
    expect_identical(defineCategories(m2), c("101", "000"))

    # edge: single row, single column
    expect_identical(defineCategories(matrix(TRUE, nrow = 1)), "1")
})

test_that("computeSetOverlaps works for simple named sets", {
    sets <- list(
        A = c("a", "b", "c"),
        B = c("b", "c", "d"),
        C = c("c", "e")
    )
    res <- computeSetOverlaps(sets)

    # class + structure
    expect_s3_class(res, "SetOverlapResult")
    expect_named(res, c("unique_elements", "overlap_matrix", "intersect_category"))

    # unique elements are union of all
    expect_setequal(res$unique_elements, c("a", "b", "c", "d", "e"))

    # matrix shape and names
    expect_true(is.matrix(res$overlap_matrix))
    expect_type(res$overlap_matrix, "logical")
    expect_equal(nrow(res$overlap_matrix), length(res$unique_elements))
    expect_equal(colnames(res$overlap_matrix), c("A", "B", "C"))
    expect_setequal(rownames(res$overlap_matrix), res$unique_elements)

    # categories length matches rows; spot-check a few
    expect_length(res$intersect_category, nrow(res$overlap_matrix))
    m <- res$overlap_matrix[rownames(res$overlap_matrix) %in% c("a", "b", "c", "d", "e"), , drop = FALSE]
    cats <- res$intersect_category[match(rownames(m), rownames(res$overlap_matrix))]

    # A only -> 100
    expect_identical(unname(cats[rownames(m) == "a"]), "100")
    # A∩B -> 110
    expect_identical(unname(cats[rownames(m) == "b"]), "110")
    # A∩B∩C -> 111
    expect_identical(unname(cats[rownames(m) == "c"]), "111")
    # B only -> 010
    expect_identical(unname(cats[rownames(m) == "d"]), "010")
    # C only -> 001
    expect_identical(unname(cats[rownames(m) == "e"]), "001")
})

test_that("computeSetOverlaps input validation", {
    # not a list
    expect_error(computeSetOverlaps("foo"), class = "simpleError")

    # list but non-character element
    bad <- list(A = letters[1:3], B = 1:3)
    expect_error(computeSetOverlaps(bad), class = "simpleError")
})

test_that("computeGenomicOverlaps merges and annotates correctly", {
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")

    gr1 <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(c(100, 300), width = 100) # [100-199], [300-399]
    )
    gr2 <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(c(150, 700), width = 100) # [150-249] overlaps gr1[1]; [700-799] no overlap
    )
    gr3 <- GenomicRanges::GRanges(
        seqnames = "chr1",
        ranges = IRanges::IRanges(900, width = 100)         # [900-999] singleton
    )

    res <- computeGenomicOverlaps(list(A = gr1, B = gr2, C = gr3))

    # class + structure
    expect_s3_class(res, "GenomicOverlapResult")
    expect_named(res, c("regions", "overlap_matrix", "mode"))
    expect_identical(res$mode, "reduce")

    # reduced regions should split into merged chunks:
    # merged of [100-199] U [150-249] -> [100-249]
    # [300-399]
    # [700-799]
    # [900-999]
    rr <- res$regions
    expect_true(inherits(rr, "GRanges"))
    expect_equal(length(rr), 4L)

    # overlap matrix shape and names
    om <- res$overlap_matrix
    expect_true(is.matrix(om))
    expect_type(om, "logical")
    expect_equal(nrow(om), length(rr))
    expect_equal(colnames(om), c("A", "B", "C"))

    # categories exist in mcols
    cats <- GenomicRanges::mcols(rr)$intersect_category
    expect_type(cats, "character")
    expect_length(cats, length(rr))

    # Spot-check expected patterns per region (order follows reduce/unlist)
    # rr[1] ~ [100-249] -> overlaps A and B -> "110"
    # rr[2] ~ [300-399] -> overlaps A only  -> "100"
    # rr[3] ~ [700-799] -> overlaps B only  -> "010"
    # rr[4] ~ [900-999] -> overlaps C only  -> "001"
    expect_identical(cats, c("110", "100", "010", "001"))

    # Matrix rows should match categories
    expect_identical(defineCategories(om), cats)
})

test_that("computeGenomicOverlaps accepts GRangesList and list<GRanges>", {
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")

    grA <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
    grB <- GenomicRanges::GRanges("chr1", IRanges::IRanges(20, 42))

    # list<GRanges>
    res1 <- computeGenomicOverlaps(list(A = grA, B = grB))
    expect_s3_class(res1, "GenomicOverlapResult")
    expect_equal(colnames(res1$overlap_matrix), c("A", "B"))

    # GRangesList
    grl <- GenomicRanges::GRangesList(A = grA, B = grB)
    res2 <- computeGenomicOverlaps(grl)
    expect_s3_class(res2, "GenomicOverlapResult")
    expect_equal(colnames(res2$overlap_matrix), c("A", "B"))
})

test_that("computeOverlaps dispatches correctly", {
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")

    # sets (atomic vectors)
    sets <- list(A = letters[1:3], B = letters[2:4])
    r1 <- computeOverlaps(sets)
    expect_s3_class(r1, "SetOverlapResult")

    # list<GRanges>
    gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 50), width = 10))
    gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(5, 90), width = 10))
    r2 <- computeOverlaps(list(A = gr1, B = gr2))
    expect_s3_class(r2, "GenomicOverlapResult")

    # GRangesList
    grl <- GenomicRanges::GRangesList(A = gr1, B = gr2)
    r3 <- computeOverlaps(grl)
    expect_s3_class(r3, "GenomicOverlapResult")
})

test_that("computeOverlaps auto-names unnamed inputs", {
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")

    # unnamed sets
    sets <- list(letters[1:3], letters[2:4])
    r1 <- computeOverlaps(sets)
    expect_equal(colnames(r1$overlap_matrix), c("set1", "set2"))

    # unnamed GRanges list
    gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
    gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(20, 42))
    r2 <- computeOverlaps(list(gr1, gr2))
    expect_equal(colnames(r2$overlap_matrix), c("set1", "set2"))
})

test_that("computeOverlaps errors on bad/mixed inputs", {
    # missing or NULL
    expect_error(computeOverlaps(), class = "simpleError")
    expect_error(computeOverlaps(NULL), class = "simpleError")

    # not a list or GRangesList
    expect_error(computeOverlaps("foo"), class = "simpleError")

    # empty list
    expect_error(computeOverlaps(list()), class = "simpleError")

    # mixed list: GRanges + character
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 10))
    bad <- list(A = gr, B = letters[1:3])
    expect_error(computeOverlaps(bad), class = "simpleError")
})

test_that("mode = 'disjoin' splits chained overlaps into pairwise ones", {
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")

    # A overlaps B, B overlaps C, but A and C share no position.
    A <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
    B <- GenomicRanges::GRanges("chr1", IRanges::IRanges(180, 300))
    C <- GenomicRanges::GRanges("chr1", IRanges::IRanges(280, 400))
    sets <- list(A = A, B = B, C = C)

    # "reduce" merges the chain into one region called a 3-way intersection
    red <- computeOverlaps(sets, mode = "reduce")
    expect_identical(red$mode, "reduce")
    expect_equal(length(red$regions), 1L)
    expect_identical(GenomicRanges::mcols(red$regions)$intersect_category,
                     "111")

    # "disjoin" keeps A∩B and B∩C separate, and reports no 3-way intersection
    dis <- computeOverlaps(sets, mode = "disjoin")
    expect_s3_class(dis, "GenomicOverlapResult")
    expect_identical(dis$mode, "disjoin")

    rr <- dis$regions
    cats <- GenomicRanges::mcols(rr)$intersect_category
    expect_identical(cats, c("100", "110", "010", "011", "001"))
    expect_false("111" %in% cats)
    expect_identical(defineCategories(dis$overlap_matrix), cats)
    expect_equal(colnames(dis$overlap_matrix), c("A", "B", "C"))

    # segments are non-overlapping and cover the same span as the reduced ones
    expect_equal(length(GenomicRanges::reduce(rr)), 1L)
    expect_equal(sum(GenomicRanges::width(rr)),
                 sum(GenomicRanges::width(red$regions)))
    expect_identical(
        as.character(GenomicRanges::start(rr)),
        as.character(c(100, 180, 201, 280, 301))
    )
    expect_identical(
        as.character(GenomicRanges::end(rr)),
        as.character(c(179, 200, 279, 300, 400))
    )
})

test_that("mode = 'disjoin' ignores within-set redundancy", {
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")

    # Two overlapping intervals within A must not be split into extra segments
    A <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(100, 150), c(200, 250)))
    B <- GenomicRanges::GRanges("chr1", IRanges::IRanges(500, 600))

    res <- computeOverlaps(list(A = A, B = B), mode = "disjoin")
    rr <- res$regions
    expect_equal(length(rr), 2L)
    expect_identical(GenomicRanges::mcols(rr)$intersect_category, c("10", "01"))
    expect_identical(as.character(GenomicRanges::end(rr)[1]), "250")
})

test_that("mode = 'disjoin' separates bookended intervals from different sets", {
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")

    # Adjacent but not overlapping: reduce() merges them, disjoin() does not
    A <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
    B <- GenomicRanges::GRanges("chr1", IRanges::IRanges(201, 300))

    red <- computeOverlaps(list(A = A, B = B), mode = "reduce")
    expect_identical(GenomicRanges::mcols(red$regions)$intersect_category,
                     "11")

    dis <- computeOverlaps(list(A = A, B = B), mode = "disjoin")
    expect_identical(GenomicRanges::mcols(dis$regions)$intersect_category,
                     c("10", "01"))
})

test_that("mode is validated and ignored for non-genomic inputs", {
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")

    gr1 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(100, 200))
    gr2 <- GenomicRanges::GRanges("chr1", IRanges::IRanges(180, 300))

    # unknown mode
    expect_error(computeOverlaps(list(A = gr1, B = gr2), mode = "merge"),
                 class = "simpleError")

    # GRangesList honours mode too
    grl <- GenomicRanges::GRangesList(A = gr1, B = gr2)
    expect_identical(computeOverlaps(grl, mode = "disjoin")$mode, "disjoin")

    # sets of identifiers: mode is irrelevant, warn if explicitly set
    sets <- list(A = letters[1:3], B = letters[2:4])
    expect_silent(computeOverlaps(sets))
    expect_warning(computeOverlaps(sets, mode = "disjoin"), "genomic inputs only")
})

test_that("category strings match overlap_matrix for both engines", {
    skip_if_not_installed("GenomicRanges")
    skip_if_not_installed("IRanges")

    # sets
    sets <- list(A = c("x", "y"), B = c("y"))
    sres <- computeSetOverlaps(sets)
    expect_identical(defineCategories(sres$overlap_matrix), sres$intersect_category)

    # genomic
    grA <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(1, 100), width = 20))
    grB <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(10, 300), width = 20))
    gres <- computeGenomicOverlaps(list(A = grA, B = grB))
    cats <- GenomicRanges::mcols(gres$regions)$intersect_category
    expect_identical(defineCategories(gres$overlap_matrix), cats)
})
