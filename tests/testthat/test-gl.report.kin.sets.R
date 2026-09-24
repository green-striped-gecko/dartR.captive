# Tests of the reviewed gl.report.kin.sets
# (function-review/reports/dartR.captive/gl.report.kin.sets.md).

skip_without_testset2 <- function() {
  skip_if_not(exists("testset2.gl") && exists("testset2.gs"),
              "testset2.gl/gs need dartR.data >= 1.2.5")
}

test_that("structure, read-only, values match the formulas", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  x0 <- x
  kin <- gl.kin(x, verbose = 0)
  expect_length(capture.output(
    r <- gl.report.kin.sets(x, kin = kin, verbose = 0)), 0)
  expect_identical(x, x0)
  expect_named(r, c("sets", "mkb", "fst"))
  expect_named(r$sets, c("pop", "n", "meanMK", "GD.w", "meanF"))
  expect_equal(nrow(r$sets), nPop(x))
  cb <- r$sets[r$sets$pop == "EmmacCaptBred", ]
  expect_equal(cb$GD.w, 0.9337007, tolerance = 1e-6)
  expect_equal(cb$meanF, -0.3284069, tolerance = 1e-6)
  b <- split(indNames(x), pop(x))
  s <- b[["EmmacCaptBred"]]
  t <- b[["EmmacMaclGeor"]]
  fb <- mean(kin[s, t])
  ft <- mean(kin[c(s, t), c(s, t)])
  expect_equal(r$mkb["EmmacCaptBred", "EmmacMaclGeor"], fb)
  expect_equal(r$fst["EmmacCaptBred", "EmmacMaclGeor"],
               1 - (1 - ft) / (1 - fb))
  expect_true(all(diag(r$fst) == 0))
  expect_equal(unname(diag(r$mkb)), r$sets$meanMK)
})

test_that("random halves of one population give Fst > 0", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  kin <- gl.kin(testset2.gl, verbose = 0)
  y <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  set.seed(1)
  p <- as.character(pop(y))
  p[sample(24, 12)] <- "half2"
  pop(y) <- p
  r <- gl.report.kin.sets(y, kin = kin, verbose = 0)
  expect_equal(r$fst[1, 2], 0.009798068, tolerance = 1e-6)
})

test_that("one population needs reference kinship", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  one <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  # before the review: GD.w = 1 and meanMK = 0 by construction
  expect_error(gl.report.kin.sets(one, verbose = 0), "wider reference")
  r <- gl.report.kin.sets(one, kin = gl.kin(testset2.gl, verbose = 0),
                          verbose = 0)
  expect_equal(r$sets$GD.w, 0.9337007, tolerance = 1e-6)
})

test_that("NA kinship is ignored with a warning; low call rate warns", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  kin <- gl.kin(x, verbose = 0)
  s <- indNames(x)[pop(x) == "EmmacCaptBred"]
  kin[s[1], s[2]] <- kin[s[2], s[1]] <- NA
  # before the review: GD.w NA for the colony and 60 Fst values NA
  out <- capture.output(r <- gl.report.kin.sets(x, kin = kin, verbose = 1))
  expect_true(any(grepl("1 missing kinship values", out)))
  expect_true(any(grepl("24 individuals have call rate below 0.8", out)))
  expect_false(anyNA(r$sets$GD.w))
  expect_false(anyNA(r$fst))
  expect_equal(r$sets$GD.w[r$sets$pop == "EmmacCaptBred"],
               1 - mean(kin[s, s], na.rm = TRUE))
})

test_that("SilicoDArT meanF is 0", {
  skip_without_testset2()
  r <- gl.report.kin.sets(testset2.gs, verbose = 0)
  expect_true(all(r$sets$meanF == 0))
})
