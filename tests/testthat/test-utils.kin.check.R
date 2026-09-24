# Tests of the reference-kinship guard
# (function-review/reports/dartR.captive/kinship-suite-rowsum-audit.md).
# Kinship is centred on the individuals it is estimated on, so mean kinship
# and whole-set gene diversity need kinship from a wider reference.

skip_without_testset2 <- function() {
  skip_if_not(exists("testset2.gl") && exists("testset2.gs"),
              "testset2.gl/gs need dartR.data >= 1.2.5")
}
captive <- function(x) {
  gl.keep.pop(x, pop.list = "EmmacCaptBred", verbose = 0)
}

test_that("grm kinship rows sum to 0 over the estimation set", {
  skip_if_not_installed("rrBLUP")
  x <- gl.filter.allna(platypus.gl, verbose = 0)
  kin <- gl.kin(x, verbose = 0)
  expect_lt(max(abs(rowMeans(kin))), 1e-12)
  expect_identical(attr(kin, "ref.ids"), indNames(x))
})

test_that("a superset kin is restricted to x, keeping its attributes", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  full <- testset2.gl
  cb <- captive(full)
  K <- gl.kin(full, verbose = 0)
  k <- utils.kin.check(cb, K, need.reference = TRUE)
  expect_identical(dimnames(k), list(indNames(cb), indNames(cb)))
  expect_equal(as.vector(k), as.vector(K[indNames(cb), indNames(cb)]))
  expect_equal(attr(k, "scale"), "kinship")
  expect_equal(attr(k, "method"), "grm")
  expect_identical(attr(k, "ref.ids"), indNames(full))
  # rows no longer average 0 over the managed group
  expect_gt(min(rowMeans(k)), 0.03)
})

test_that("self-referenced kin is rejected when a reference is needed", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  cb <- captive(testset2.gl)
  msg <- "wider reference"
  expect_error(utils.kin.check(cb, NULL, need.reference = TRUE), msg)
  k.self <- gl.kin(cb, verbose = 0)
  expect_error(utils.kin.check(cb, k.self, need.reference = TRUE), msg)
  # attributes lost by subsetting: caught by the zero row means
  k.bare <- k.self[indNames(cb), indNames(cb)]
  attributes(k.bare) <- attributes(k.bare)[c("dim", "dimnames")]
  expect_error(utils.kin.check(cb, k.bare, need.reference = TRUE), msg)
  # without need.reference the old behaviour stands
  expect_silent(k <- utils.kin.check(cb, NULL))
  expect_identical(rownames(k), indNames(cb))
})

test_that("missing individuals are an error", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  K <- gl.kin(testset2.gl, verbose = 0)
  expect_error(utils.kin.check(testset2.gl, K[1:10, 1:10]),
               "include every indNames")
})

test_that("management functions stop without reference kin, run with it", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  pdf(NULL)
  on.exit(grDevices::dev.off())
  full <- testset2.gl
  cb <- captive(full)
  K <- gl.kin(full, verbose = 0)
  msg <- "wider reference"
  expect_error(gl.report.ind.remove(cb, verbose = 0), msg)
  expect_error(gl.report.repro.targets(cb, verbose = 0), msg)
  expect_error(gl.select.pairs(cb, verbose = 0), msg)
  expect_error(gl.report.mate.suitability(cb, verbose = 0), msg)
  expect_error(gl.report.gd.projection(cb, ne = 20, plot.display = FALSE,
                                       verbose = 0), msg)
  # superset kin gives the same result as subsetting by hand
  a <- gl.report.repro.targets(cb, kin = K, verbose = 0)
  b <- gl.report.repro.targets(cb, kin = K[indNames(cb), indNames(cb)],
                               verbose = 0)
  expect_equal(a, b)
  expect_equal(attr(gl.select.pairs(cb, kin = K, verbose = 0), "gd.start"),
               1 - mean(K[indNames(cb), indNames(cb)]))
  p <- gl.report.gd.projection(cb, kin = K, ne = 20, plot.display = FALSE,
                               verbose = 0)
  expect_equal(p$summary$gd.now, 1 - mean(K[indNames(cb), indNames(cb)]))
  # gd.now supplied directly needs no kinship
  expect_equal(gl.report.gd.projection(gd.now = 0.3, ne = 25, years = 50,
                                       plot.display = FALSE,
                                       verbose = 0)$summary$gd.now, 0.3)
})

test_that("kin.sets meanMK is the within-set block mean", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  K <- gl.kin(x, verbose = 0)
  s <- gl.report.kin.sets(x, kin = K, verbose = 0)$sets
  cb <- indNames(x)[pop(x) == "EmmacCaptBred"]
  expect_equal(s$meanMK[s$pop == "EmmacCaptBred"], mean(K[cb, cb]))
  expect_equal(s$meanMK, 1 - s$GD.w)
})

test_that("gl.report.kinship warns for one self-referenced population", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  pdf(NULL)
  on.exit(grDevices::dev.off())
  cb <- captive(testset2.gl)
  out <- capture.output(gl.report.kinship(cb, verbose = 1))
  expect_true(any(grepl("one population", out)))
  out <- capture.output(gl.report.kinship(cb, kin = gl.kin(testset2.gl,
                                                           verbose = 0),
                                          verbose = 1))
  expect_false(any(grepl("one population", out)))
})
