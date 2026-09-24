# Tests of the reviewed gl.report.kinship
# (function-review/reports/dartR.captive/gl.report.kinship.md).

skip_without_testset2 <- function() {
  skip_if_not(exists("testset2.gl") && exists("testset2.gs"),
              "testset2.gl/gs need dartR.data >= 1.2.5")
}

test_that("structure, read-only, silent at verbose = 0", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  h <- length(x@other$history)
  expect_length(capture.output(r <- gl.report.kinship(x, verbose = 0)), 0)
  expect_named(r, c("ind", "pop"))
  expect_named(r$ind, c("id", "pop", "sex", "MK", "MKrank", "F"))
  expect_named(r$pop, c("pop", "n", "meanMK", "GD", "FGE", "meanF"))
  expect_equal(nrow(r$ind), nInd(x))
  expect_equal(nrow(r$pop), nPop(x) + 1)
  expect_equal(r$pop$pop[1], "overall")
  expect_equal(length(x@other$history), h)
})

test_that("MK is within population; ranks within population and sex", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  kin <- gl.kin(x, verbose = 0)
  r <- gl.report.kinship(x, kin = kin, verbose = 0)
  cb <- indNames(x)[pop(x) == "EmmacCaptBred"]
  i <- match(cb, r$ind$id)
  expect_equal(r$ind$MK[i], unname(rowMeans(kin[cb, cb])))
  # the dataset-wide row means are 0 by construction (why MK is per
  # population); the old MKrank ranked these rounding residuals
  expect_lt(max(abs(rowMeans(kin))), 1e-12)
  expect_equal(r$ind$F, unname(2 * diag(kin) - 1))
  # ranks 1..n within each population x sex group
  grp <- paste(r$ind$pop, r$ind$sex)
  ok <- tapply(r$ind$MKrank, grp, function(v) all(sort(v) == seq_along(v)))
  expect_true(all(ok))
  # rank 1 male of the captive group is its lowest within-group MK
  m <- r$ind[i, ][r$ind$sex[i] == "Male", ]
  expect_equal(m$id[m$MKrank == 1], m$id[which.min(m$MK)])
  # the ind and pop tables agree
  expect_equal(mean(r$ind$MK[i]),
               r$pop$meanMK[r$pop$pop == "EmmacCaptBred"])
})

test_that("overall row: GD and FGE NA; population rows unchanged", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  r <- gl.report.kinship(testset2.gl, verbose = 0)
  expect_true(is.na(r$pop$GD[1]) && is.na(r$pop$FGE[1]))
  expect_lt(abs(r$pop$meanMK[1]), 1e-12)
  expect_false(anyNA(r$pop$FGE[-1]))
  cbp <- r$pop[r$pop$pop == "EmmacCaptBred", ]
  expect_equal(cbp$GD, 0.9336865, tolerance = 1e-6)
  expect_equal(cbp$FGE, 7.539946, tolerance = 1e-6)
  expect_equal(cbp$FGE, 1 / (2 * cbp$meanMK))
})

test_that("dominant data: no spurious FGE warning, F is 0", {
  skip_without_testset2()
  out <- capture.output(r <- gl.report.kinship(testset2.gs, verbose = 1))
  expect_true(is.na(r$pop$FGE[1]))
  expect_false(any(grepl("FGE reported as NA", out)))
  expect_true(all(r$ind$F == 0))
})

test_that("NA, empty and 'Unknown' sexes form one rank group", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  x@other$ind.metrics$sex[1:3] <- NA
  r <- gl.report.kinship(x, verbose = 0)
  unk <- is.na(r$ind$sex) | r$ind$sex == "Unknown"
  for (p in unique(r$ind$pop[unk])) {
    v <- r$ind$MKrank[unk & r$ind$pop == p]
    expect_equal(sort(v), seq_along(v))
  }
})

test_that("missing kinship values are ignored with a warning", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  kin <- gl.kin(x, verbose = 0)
  kin[1, 2] <- kin[2, 1] <- NA
  out <- capture.output(r <- gl.report.kinship(x, kin = kin, verbose = 1))
  expect_false(anyNA(r$ind$MK))
  expect_false(is.na(r$pop$meanMK[1]))
  expect_true(any(grepl("1 pair\\(s\\).*missing kinship", out)))
})
