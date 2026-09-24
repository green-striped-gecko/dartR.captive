# Tests of the reviewed gl.report.kin.groups
# (function-review/reports/dartR.captive/gl.report.kin.groups.md).

skip_without_testset2 <- function() {
  skip_if_not(exists("testset2.gl") && exists("testset2.gs"),
              "testset2.gl/gs need dartR.data >= 1.2.5")
}

test_that("reference kinship, structure, read-only, values", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  kin <- gl.kin(testset2.gl, verbose = 0)
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  cb0 <- cb
  expect_length(capture.output(
    r <- gl.report.kin.groups(cb, group.col = "cohort",
                              kin = kin[indNames(cb), indNames(cb)],
                              verbose = 0)), 0)
  expect_identical(cb, cb0)
  expect_named(r, c("groups", "kin.groups", "gd"))
  expect_named(r$groups, c("group", "n", "MK", "meanF"))
  expect_equal(r$groups$n, c(5, 4, 5, 5, 5))
  expect_equal(r$groups$MK,
               c(0.0912, 0.0622, 0.0472, 0.0782, 0.0519), tolerance = 1e-3)
  expect_equal(r$groups$meanF,
               c(-0.2725, -0.2514, -0.3365, -0.3821, -0.3834),
               tolerance = 1e-3)
  sub <- kin[indNames(cb), indNames(cb)]
  expect_equal(r$gd, 1 - mean(sub))
  expect_equal(r$kin.groups["F1_AB", "F1_AE"], 0.0823, tolerance = 1e-3)
  # a larger reference matrix is restricted to x (utils.kin.check, #113)
  rl <- gl.report.kin.groups(cb, group.col = "cohort", kin = kin,
                             verbose = 0)
  expect_equal(rl$groups, r$groups)
})

test_that("self-referenced kinship is an error", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  # before the review: MK = 0 for every group and GD = 1
  expect_error(gl.report.kin.groups(testset2.gl, verbose = 0),
               "wider reference")
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  expect_error(gl.report.kin.groups(cb, group.col = "cohort",
                                    kin = gl.kin(cb, verbose = 0),
                                    verbose = 0), "wider reference")
})

test_that("NA kinship is ignored in the block means, with a warning", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  kin <- gl.kin(testset2.gl, verbose = 0)
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  ids <- indNames(cb)
  kin[ids[1], ids[2]] <- kin[ids[2], ids[1]] <- NA
  expect_output(r <- gl.report.kin.groups(cb, group.col = "cohort",
                                          kin = kin, verbose = 1),
                "1 pairs have missing kinship")
  expect_false(anyNA(r$groups$MK))
  expect_false(is.na(r$gd))
})

test_that("low call rate warns", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  kin <- gl.kin(testset2.gl, verbose = 0)
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  expect_output(gl.report.kin.groups(cb, group.col = "cohort", kin = kin,
                                     verbose = 1),
                "24 individuals have call rate below 0.8")
})

test_that("SilicoDArT meanF is 0 by construction", {
  skip_without_testset2()
  gs <- testset2.gs
  kin <- gl.kin(gs, verbose = 0)
  sub <- gs[pop(gs) %in% levels(pop(gs))[1:3], ]
  r <- gl.report.kin.groups(sub, kin = kin, verbose = 0)
  expect_true(all(r$groups$meanF == 0))
})

test_that("group column handling", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  kin <- gl.kin(testset2.gl, verbose = 0)
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  expect_error(gl.report.kin.groups(cb, group.col = "nope", kin = kin,
                                    verbose = 0), "not found")
  cb@other$ind.metrics$cohort[1:2] <- NA
  expect_output(r <- gl.report.kin.groups(cb, group.col = "cohort",
                                          kin = kin, verbose = 1),
                "2 individual\\(s\\) with missing group value dropped")
  expect_equal(sum(r$groups$n), 22)
})
