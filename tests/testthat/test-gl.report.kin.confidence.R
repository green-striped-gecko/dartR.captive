# Tests of the reviewed gl.report.kin.confidence
# (function-review/reports/dartR.captive/gl.report.kin.confidence.md).

skip_without_testset2 <- function() {
  skip_if_not(exists("testset2.gl") && exists("testset2.gs"),
              "testset2.gl/gs need dartR.data >= 1.2.5")
}

test_that("SNP structure, read-only, seeded values", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  x0 <- x
  kin <- gl.kin(x, verbose = 0)
  set.seed(1)
  expect_length(capture.output(
    r <- gl.report.kin.confidence(x, kin = kin, nboots = 50, verbose = 0)
  ), 0)
  expect_identical(x, x0)
  expect_named(r, c("se", "ci.width", "summary"))
  expect_equal(dim(r$se), c(nInd(x), nInd(x)))
  expect_equal(dimnames(r$se), list(indNames(x), indNames(x)))
  expect_named(r$summary, c("mean.se", "median.se", "q90.se"))
  # G/2, as gl.kin returns (was G/2 - mean(diag(G) - 1): median 0.0258)
  expect_equal(r$summary$median.se, 0.02202817, tolerance = 1e-6)
  expect_equal(r$se["CB_AB_01", "CB_AB_02"], 0.01981219, tolerance = 1e-6)
})

test_that("seeded runs are reproducible", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  kin <- gl.kin(x, verbose = 0)
  set.seed(3)
  r1 <- gl.report.kin.confidence(x, kin = kin, nboots = 5, verbose = 0)
  set.seed(3)
  r2 <- gl.report.kin.confidence(x, kin = kin, nboots = 5, verbose = 0)
  expect_identical(r1, r2)
})

test_that("SilicoDArT", {
  skip_without_testset2()
  set.seed(1)
  r <- gl.report.kin.confidence(testset2.gs, nboots = 30, verbose = 0)
  # gl.kin's dominant estimator (was band correlation / 2: median 0.0182)
  expect_equal(r$summary$median.se, 0.01786112, tolerance = 1e-6)
})

test_that("the bootstrapped estimators reproduce gl.kin on all loci", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  kin <- gl.kin(x, verbose = 0)
  m <- as.matrix(x)
  p <- colMeans(m, na.rm = TRUE) / 2
  m <- m[, !is.na(p)]
  p <- p[!is.na(p)]
  Z <- sweep(m, 2, 2 * p)
  Z[is.na(Z)] <- 0
  G <- tcrossprod(Z) / (2 * sum(p * (1 - p)))
  d <- G / 2 - unclass(kin)[, ]
  off <- row(d) != col(d)
  expect_lt(max(abs(d[off])), 1e-4)
  expect_lt(max(abs(diag(d))), 0.01)
})

test_that("warnings: low call rate, estimator mismatch, NA kin", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  kin <- gl.kin(x, verbose = 0)
  expect_output(gl.report.kin.confidence(x, kin = kin, nboots = 2,
                                         verbose = 1),
                "24 individuals have call rate below 0.8")
  k2 <- kin
  attr(k2, "method") <- "emibd9"
  expect_output(gl.report.kin.confidence(x, kin = k2, nboots = 2,
                                         verbose = 1),
                "method 'emibd9'")
  k3 <- kin
  k3[1, 2] <- k3[2, 1] <- NA
  out <- capture.output(gl.report.kin.confidence(x, kin = k3, nboots = 3,
                                                 verbose = 3))
  expect_false(any(grepl("NaN", out)))
})

test_that("input errors", {
  skip_without_testset2()
  x <- testset2.gl
  expect_error(gl.report.kin.confidence(x, nboots = 1, verbose = 0), "nboots")
  expect_error(gl.report.kin.confidence(x, conf = 1, verbose = 0), "conf")
})
