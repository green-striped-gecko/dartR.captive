# Tests of the reviewed gl.report.kin.classes
# (function-review/reports/dartR.captive/gl.report.kin.classes.md).

skip_without_testset2 <- function() {
  skip_if_not(exists("testset2.gl") && exists("testset2.gs"),
              "testset2.gl/gs need dartR.data >= 1.2.5")
}

test_that("SNP structure, read-only, class counts", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  x0 <- x
  kin <- gl.kin(x, verbose = 0)
  expect_length(capture.output(r <- gl.report.kin.classes(x, kin = kin,
                                                          verbose = 0)), 0)
  expect_identical(x, x0)
  expect_named(r, c("pairs", "conflicts"))
  expect_named(r$pairs, c("id1", "id2", "kinship", "kinship.adj", "class"))
  expect_equal(as.vector(table(r$pairs$class)), c(127, 30, 1000, 3233))
  expect_equal(names(table(r$pairs$class)),
               c("full-sib", "parent-offspring", "second-degree",
                 "third-degree"))
  baseline <- stats::median(kin[row(kin) != col(kin)])
  expect_equal(r$pairs$kinship.adj, r$pairs$kinship - baseline)
  expect_equal(nrow(r$conflicts), 19)
  expect_true(all(r$conflicts$class == "second-degree"))
  ab <- r$pairs[r$pairs$id1 == "CB_AB_01" & r$pairs$id2 == "CB_AB_02", ]
  expect_equal(ab$class, "full-sib")
  ra <- gl.report.kin.classes(x, kin = kin, all.pairs = TRUE, verbose = 0)
  expect_equal(nrow(ra$pairs), nInd(x) * (nInd(x) - 1) / 2)
})

test_that("gl.grm input is halved to the same classes", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  g <- gl.grm(x, plotheatmap = FALSE, verbose = 0)
  r <- gl.report.kin.classes(x, kin = g, verbose = 0)
  expect_equal(as.vector(table(r$pairs$class)), c(127, 30, 1000, 3233))
})

test_that("SilicoDArT keeps first-degree", {
  skip_without_testset2()
  r <- gl.report.kin.classes(testset2.gs, verbose = 0)
  expect_equal(as.vector(table(r$pairs$class)), c(192, 944, 3325))
  expect_equal(nrow(r$conflicts), 39)
})

test_that("NA kinship affects only its own pairs, with a warning", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  kin <- gl.kin(x, verbose = 0)
  kin[1, 2] <- kin[2, 1] <- NA
  expect_output(r <- gl.report.kin.classes(x, kin = kin, verbose = 1),
                "1 pairs have missing kinship")
  expect_false(anyNA(r$pairs$class))
  expect_gt(nrow(r$pairs), 4000)
  ra <- gl.report.kin.classes(x, kin = kin, all.pairs = TRUE, verbose = 0)
  expect_equal(sum(is.na(ra$pairs$class)), 1)
  expect_equal(nrow(ra$pairs), nInd(x) * (nInd(x) - 1) / 2)
})

test_that("fewer than three individuals error; fewer than ten warn", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  expect_error(gl.report.kin.classes(testset2.gl[1:2, ], verbose = 0),
               "three individuals")
  expect_output(gl.report.kin.classes(testset2.gl[1:5, ], verbose = 1),
                "median baseline is unreliable")
})

test_that("population-structure warning at verbose >= 1 only", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  kin <- gl.kin(x, verbose = 0)
  expect_output(gl.report.kin.classes(x, kin = kin, verbose = 1),
                "31 populations in x")
  cb <- x[pop(x) == "EmmacCaptBred", ]
  out <- capture.output(gl.report.kin.classes(cb, verbose = 1))
  expect_false(any(grepl("populations in x", out)))
})

test_that("low call rate warns; filtering resolves most conflicts", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  kin <- gl.kin(x, verbose = 0)
  expect_output(gl.report.kin.classes(x, kin = kin, verbose = 1),
                "24 individuals have call rate below 0.8")
  xf <- gl.filter.callrate(x, method = "loc", threshold = 0.95, verbose = 0)
  r <- gl.report.kin.classes(xf, kin = gl.kin(xf, verbose = 0), verbose = 0)
  expect_equal(nrow(r$conflicts), 5)
})

test_that("conflicts is a zero-row data frame when there are none", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  kin <- gl.kin(x, verbose = 0)
  x@other$ind.metrics$sire <- NULL
  r <- gl.report.kin.classes(x, kin = kin, verbose = 0)
  expect_s3_class(r$conflicts, "data.frame")
  expect_equal(nrow(r$conflicts), 0)
  expect_named(r$conflicts,
               c("offspring", "parent", "kinship", "kinship.adj", "class"))
})

test_that("input errors", {
  skip_without_testset2()
  x <- testset2.gl
  expect_error(gl.report.kin.classes(x, oh.thresh = 1, verbose = 0),
               "oh.thresh")
  expect_error(gl.report.kin.classes(x[1, ], verbose = 0), "three")
})
