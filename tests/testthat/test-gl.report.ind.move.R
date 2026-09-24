# Tests of the reviewed gl.report.ind.move
# (function-review/reports/dartR.captive/gl.report.ind.move.md).

skip_without_testset2 <- function() {
  skip_if_not(exists("testset2.gl") && exists("testset2.gs"),
              "testset2.gl/gs need dartR.data >= 1.2.5")
}

pl <- c("EmmacCaptBred", "EmmacMaclGeor", "EmmacBurnBara")

test_that("structure, read-only, values match the formula", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  kin <- gl.kin(testset2.gl, verbose = 0)
  ex <- gl.keep.pop(testset2.gl, pop.list = pl, verbose = 0)
  ex0 <- ex
  expect_length(capture.output(
    r <- gl.report.ind.move(ex, kin = kin, verbose = 0)), 0)
  expect_identical(ex, ex0)
  expect_named(r, c("id", "from", "to", "dgd.source", "dgd.dest", "net"))
  expect_equal(nrow(r), nInd(ex) * 2)
  expect_false(is.unsorted(rev(r$net)))
  expect_equal(r$net, r$dgd.source + r$dgd.dest)
  k <- kin[indNames(ex), indNames(ex)]
  b <- split(indNames(ex), pop(ex))
  gd <- function(ids) 1 - mean(k[ids, ids])
  for (i in 1:10) {
    s <- b[[r$from[i]]]
    t <- b[[r$to[i]]]
    expect_equal(r$dgd.source[i], gd(setdiff(s, r$id[i])) - gd(s))
    expect_equal(r$dgd.dest[i], gd(c(t, r$id[i])) - gd(t))
  }
  expect_equal(sum(r$from[1:10] == "EmmacCaptBred"), 9)
})

test_that("NA kinship is ignored with a warning; low call rate warns", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  kin <- gl.kin(testset2.gl, verbose = 0)
  ex <- gl.keep.pop(testset2.gl, pop.list = pl, verbose = 0)
  cb <- indNames(ex)[pop(ex) == "EmmacCaptBred"]
  kin[cb[1], cb[2]] <- kin[cb[2], cb[1]] <- NA
  # before the review 70 of 92 moves (all into or out of EmmacCaptBred)
  # were NA
  out <- capture.output(r <- gl.report.ind.move(ex, kin = kin, verbose = 1))
  expect_true(any(grepl("1 pairs have missing kinship", out)))
  expect_true(any(grepl("24 individuals have call rate below 0.8", out)))
  expect_false(anyNA(r$net))
})

test_that("SilicoDArT, singletons and input errors", {
  skip_without_testset2()
  ex.gs <- gl.keep.pop(testset2.gs, pop.list = pl, verbose = 0)
  r <- gl.report.ind.move(ex.gs, verbose = 0)
  expect_equal(nrow(r), nInd(ex.gs) * 2)
  one <- gl.keep.pop(testset2.gs, pop.list = "EmmacCaptBred", verbose = 0)
  expect_error(gl.report.ind.move(one, verbose = 0), "two populations")
})
