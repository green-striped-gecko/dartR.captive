# Tests of the reviewed gl.report.mate.suitability
# (function-review/reports/dartR.captive/gl.report.mate.suitability.md).

skip_without_testset2 <- function() {
  skip_if_not(exists("testset2.gl") && exists("testset2.gs"),
              "testset2.gl/gs need dartR.data >= 1.2.5")
}

cb.setup <- function() {
  kin <- gl.kin(testset2.gl, verbose = 0)
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  list(cb = cb, kin = kin[indNames(cb), indNames(cb)])
}

test_that("structure, read-only, component values", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  s <- cb.setup()
  cb0 <- s$cb
  r <- gl.report.mate.suitability(s$cb, kin = s$kin, verbose = 0)
  expect_identical(s$cb, cb0)
  expect_named(r, c("msi", "f.off", "dgd", "mkdiff", "completeness",
                    "rank.dgd", "rank.mkdiff", "rank.f", "msi.base",
                    "settings"))
  expect_equal(dim(r$msi), c(11, 13))
  expect_equal(as.vector(table(factor(r$msi,
                                      levels = c(1:6, "NoWay")))),
               c(5, 3, 17, 50, 12, 25, 31))
  expect_equal(r$settings$f.noway, 0.125)
  expect_equal(r$settings$f.breakpoint, 0.0600394, tolerance = 1e-5)
  expect_equal(r$f.off["CB_AB_01", "CB_AB_02"], 0.2137832, tolerance = 1e-5)
  expect_equal(r$msi["CB_AB_01", "CB_CD_02"], "6")
  # dgd equals the closed form for one virtual offspring
  k <- s$kin
  n <- nrow(k)
  m <- rownames(r$dgd)
  f <- colnames(r$dgd)
  cf <- (sum(k) / n^2) -
    (sum(k) + outer(rowSums(k)[m], rowSums(k)[f], "+") +
       0.5 * (1 + k[m, f])) / (n + 1)^2
  expect_equal(r$dgd, cf)
  expect_equal(r$mkdiff,
               abs(outer(rowMeans(k)[m], rowMeans(k)[f], "-")),
               ignore_attr = TRUE)
})

test_that("unknown.breaks, invalid f.noway, SilicoDArT", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  s <- cb.setup()
  r <- gl.report.mate.suitability(s$cb, kin = s$kin,
                                  unknown.breaks = c(0.0625, 0.125, 0.25, 0.5),
                                  verbose = 0)
  expect_equal(as.vector(table(r$msi)), c(47, 65, 31))
  # an out-of-range f.noway is an error (was reset to 0.125; change 4)
  expect_error(gl.report.mate.suitability(s$cb, kin = s$kin, f.noway = 5,
                                          verbose = 0),
               "f.noway must be")
  kg <- gl.kin(testset2.gs, verbose = 0)
  cbs <- gl.keep.pop(testset2.gs, pop.list = "EmmacCaptBred", verbose = 0)
  rs <- gl.report.mate.suitability(cbs, kin = kg[indNames(cbs), indNames(cbs)],
                                   verbose = 0)
  expect_equal(sum(rs$msi == "NoWay"), 34)
})

test_that("NA sex, NA kinship, low call rate, self-reference", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  s <- cb.setup()
  r0 <- gl.report.mate.suitability(s$cb, kin = s$kin, verbose = 0)
  # NA sex counts as unknown (was 'subscript out of bounds'; change 3)
  cb2 <- s$cb
  cb2@other$ind.metrics$sex[1] <- NA
  expect_output(r <- gl.report.mate.suitability(cb2, kin = s$kin,
                                                verbose = 1),
                "1 individuals of unknown sex")
  expect_equal(sum(dim(r$msi)), sum(dim(r0$msi)) - 1)
  # NA kinship is ignored in the means; the NA pairing gets MSI NA
  # (was 'missing value where TRUE/FALSE needed'; change 2)
  m <- rownames(r0$msi)[1]
  f <- colnames(r0$msi)[1]
  k2 <- s$kin
  k2[m, f] <- k2[f, m] <- NA
  expect_output(r <- gl.report.mate.suitability(s$cb, kin = k2, verbose = 1),
                "1 missing \\(NA\\) kinship")
  expect_true(is.na(r$msi[m, f]))
  expect_equal(sum(is.na(r$msi)), 1)
  expect_false(anyNA(r$dgd))
  # one NA barely moves the other pairings' dgd
  expect_equal(r$dgd[-1, -1], r0$dgd[-1, -1], tolerance = 1e-3)
  # low call rate warns (change 1)
  expect_output(gl.report.mate.suitability(s$cb, kin = s$kin, verbose = 1),
                "24 individuals have call rate below 0.8")
  expect_error(gl.report.mate.suitability(s$cb, verbose = 0),
               "wider")
})
