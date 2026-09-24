# Tests of the reviewed gl.report.gd.projection
# (function-review/reports/dartR.captive/gl.report.gd.projection.md).

skip_without_testset2 <- function() {
  skip_if_not(exists("testset2.gl") && exists("testset2.gs"),
              "testset2.gl/gs need dartR.data >= 1.2.5")
}

test_that("genlight input, structure, read-only, formulas", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  kin <- gl.kin(testset2.gl, verbose = 0)
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  cb0 <- cb
  expect_length(capture.output(
    p <- gl.report.gd.projection(cb, kin = kin, ne = 50, years = 100,
                                 plot.display = FALSE, verbose = 0)), 0)
  expect_identical(cb, cb0)
  expect_named(p, c("projection", "summary"))
  expect_named(p$projection, c("year", "gd", "prop.retained"))
  expect_named(p$summary, c("gd.now", "ne", "years.to.target",
                            "ne.required", "years.to.target.source",
                            "ne.required.source"))
  expect_equal(nrow(p$projection), 101)
  expect_equal(p$summary$gd.now,
               1 - mean(kin[indNames(cb), indNames(cb)]))
  expect_equal(p$summary$gd.now, 0.9337007, tolerance = 1e-6)
  lam <- 1 - 1 / 100
  expect_equal(p$projection$prop.retained, lam^(0:100))
  expect_equal(p$summary$years.to.target, log(0.9) / log(lam))
  expect_equal(p$summary$ne.required, 1 / (2 * (1 - 0.9^(1 / 100))))
  # source-relative: gene diversity itself falls to 0.9
  g <- p$summary$gd.now
  expect_equal(p$summary$years.to.target.source, log(0.9 / g) / log(lam))
  expect_equal(p$summary$years.to.target.source, 3.66, tolerance = 1e-3)
  expect_equal(p$summary$ne.required.source,
               1 / (2 * (1 - (0.9 / g)^(1 / 100))))
  expect_equal(p$summary$ne.required.source, 1360.4, tolerance = 1e-4)
})

test_that("source-relative outputs are NA when the goal is already missed", {
  expect_output(p <- gl.report.gd.projection(gd.now = 0.85, ne = 50,
                                             plot.display = FALSE,
                                             verbose = 1),
                "already at or below gd.target")
  expect_true(is.na(p$summary$years.to.target.source))
  expect_true(is.na(p$summary$ne.required.source))
  expect_false(is.na(p$summary$years.to.target))
})

test_that("retention outputs do not depend on gd.now", {
  a <- gl.report.gd.projection(gd.now = 0.5, ne = 50, plot.display = FALSE,
                               verbose = 0)
  b <- gl.report.gd.projection(gd.now = 0.95, ne = 50, plot.display = FALSE,
                               verbose = 0)
  expect_equal(a$projection$prop.retained, b$projection$prop.retained)
  expect_equal(a$summary$years.to.target, b$summary$years.to.target)
  expect_equal(a$summary$ne.required, b$summary$ne.required)
})

test_that("gen.length scales time", {
  p <- gl.report.gd.projection(gd.now = 0.9, ne = 25, gen.length = 5,
                               years = 50, plot.display = FALSE, verbose = 0)
  lam <- 1 - 1 / 50
  expect_equal(p$projection$gd, 0.9 * lam^((0:50) / 5))
  expect_equal(p$summary$years.to.target, 5 * log(0.9) / log(lam))
  expect_equal(p$summary$ne.required, 1 / (2 * (1 - 0.9^(5 / 50))))
})

test_that("NA kinship is ignored with a warning; low call rate warns", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  kin <- gl.kin(testset2.gl, verbose = 0)
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  ids <- indNames(cb)
  kin[ids[1], ids[2]] <- kin[ids[2], ids[1]] <- NA
  # before the review: "missing value where TRUE/FALSE needed"
  out <- capture.output(p <- gl.report.gd.projection(cb, kin = kin, ne = 50,
                                                     plot.display = FALSE,
                                                     verbose = 1))
  expect_true(any(grepl("1 missing kinship values", out)))
  expect_true(any(grepl("24 individuals have call rate below 0.8", out)))
  sub <- kin[ids, ids]
  expect_equal(p$summary$gd.now, 1 - mean(sub, na.rm = TRUE))
})

test_that("input errors and plot file", {
  skip_without_testset2()
  expect_error(gl.report.gd.projection(ne = 50, verbose = 0), "Supply either")
  expect_error(gl.report.gd.projection(gd.now = 0.9, ne = 0.5, verbose = 0),
               "ne must be")
  expect_error(gl.report.gd.projection(gd.now = 1.2, ne = 50, verbose = 0),
               "gd.now")
  expect_error(gl.report.gd.projection(gd.now = 0.9, ne = 50, gd.target = 1,
                                       verbose = 0), "gd.target")
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  expect_error(gl.report.gd.projection(cb, ne = 50, verbose = 0),
               "wider reference")
  # scalar arguments are checked before any kinship work
  expect_error(gl.report.gd.projection(cb, verbose = 0), "ne must be")
  d <- tempfile()
  dir.create(d)
  gl.report.gd.projection(gd.now = 0.9, ne = 25, plot.display = FALSE,
                          plot.file = "gp", plot.dir = d, verbose = 0)
  expect_true(file.exists(file.path(d, "gp.RDS")))
})
