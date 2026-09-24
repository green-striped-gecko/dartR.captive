# Tests of the reviewed gl.report.repro.targets
# (function-review/reports/dartR.captive/gl.report.repro.targets.md).

skip_without_testset2 <- function() {
  skip_if_not(exists("testset2.gl") && exists("testset2.gs"),
              "testset2.gl/gs need dartR.data >= 1.2.5")
}

cb.setup <- function() {
  kin <- gl.kin(testset2.gl, verbose = 0)
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  list(cb = cb, kin = kin[indNames(cb), indNames(cb)])
}

test_that("structure, read-only, allocation", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  s <- cb.setup()
  cb0 <- s$cb
  tg <- gl.report.repro.targets(s$cb, kin = s$kin, verbose = 0)
  expect_identical(s$cb, cb0)
  expect_named(tg, c("id", "pop", "sex", "MK", "target"))
  expect_equal(nrow(tg), 24)
  expect_equal(as.vector(tapply(tg$target, tg$sex, sum)), c(24, 24))
  expect_equal(tg$MK, unname(rowMeans(s$kin)[tg$id]), ignore_attr = TRUE)
  tt <- setNames(tg$target, tg$id)
  expect_equal(unname(tt[c("CB_CD_01", "CB_AB_01", "CB_AB_02", "CB_X_05")]),
               c(4, 0, 0, 2))
  # expected mean kinship of the offspring generation: 0.0582 against
  # 0.0667 for equal contributions
  # each parent passes on half of each offspring genome: c = target / 48
  c <- tg$target / 48
  k <- s$kin[tg$id, tg$id]
  expect_equal(as.numeric(t(c) %*% k %*% c), 0.05815443, tolerance = 1e-5)
  tg12 <- gl.report.repro.targets(s$cb, kin = s$kin, n.target = 12,
                                  verbose = 0)
  expect_equal(as.vector(tapply(tg12$target, tg12$sex, sum)), c(12, 12))
})

test_that("n.target handling, SilicoDArT", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  s <- cb.setup()
  # an invalid n.target is an error (was reset to nInd(x); change 4)
  for (bad in list(0, -3, NA_real_, "12", c(10, 12))) {
    expect_error(gl.report.repro.targets(s$cb, kin = s$kin,
                                         n.target = bad, verbose = 0),
                 "n.target must be")
  }
  expect_output(tg <- gl.report.repro.targets(s$cb, kin = s$kin,
                                              n.target = 12.7, verbose = 1),
                "rounded down to 12")
  expect_equal(sum(tg$target), 24)
  kg <- gl.kin(testset2.gs, verbose = 0)
  cbs <- gl.keep.pop(testset2.gs, pop.list = "EmmacCaptBred", verbose = 0)
  expect_equal(sum(gl.report.repro.targets(cbs, kin = kg,
                                           verbose = 0)$target), 48)
})

test_that("NA sex, NA kinship, low call rate, self-reference", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  s <- cb.setup()
  # NA sex counts as unknown (was 'missing value'; change 3)
  cb2 <- s$cb
  cb2@other$ind.metrics$sex[1] <- NA
  expect_output(tg <- gl.report.repro.targets(cb2, kin = s$kin, verbose = 1),
                "1 individuals of unknown sex")
  expect_equal(nrow(tg), 23)
  expect_false(indNames(cb2)[1] %in% tg$id)
  expect_equal(as.vector(tapply(tg$target, tg$sex, sum)), c(24, 24))
  # NA kinship is ignored in the mean kinship (was 'missing value'; change 2)
  k2 <- s$kin
  k2[1, 2] <- k2[2, 1] <- NA
  expect_output(tg <- gl.report.repro.targets(s$cb, kin = k2, verbose = 1),
                "1 missing \\(NA\\) kinship")
  expect_false(anyNA(tg$MK))
  expect_equal(as.vector(tapply(tg$target, tg$sex, sum)), c(24, 24))
  # low call rate warns (change 1)
  expect_output(gl.report.repro.targets(s$cb, kin = s$kin, verbose = 1),
                "24 individuals have call rate below 0.8")
  expect_error(gl.report.repro.targets(s$cb, verbose = 0), "wider")
})
