# Tests of the reviewed gl.report.ind.remove
# (function-review/reports/dartR.captive/gl.report.ind.remove.md).

skip_without_testset2 <- function() {
  skip_if_not(exists("testset2.gl") && exists("testset2.gs"),
              "testset2.gl/gs need dartR.data >= 1.2.5")
}

test_that("structure, read-only, values, greedy set", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  kin <- gl.kin(testset2.gl, verbose = 0)
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  cb0 <- cb
  expect_length(capture.output(
    r <- gl.report.ind.remove(cb, kin = kin, verbose = 0)), 0)
  expect_identical(cb, cb0)
  expect_named(r, c("ranking", "removal.set"))
  expect_named(r$ranking, c("id", "pop", "MK", "dGD"))
  expect_named(r$removal.set, c("step", "id", "gd.after"))
  k <- kin[indNames(cb), indNames(cb)]
  g <- function(ids) 1 - mean(k[ids, ids])
  d <- vapply(r$ranking$id, function(i)
    g(setdiff(indNames(cb), i)) - g(indNames(cb)), numeric(1))
  expect_equal(r$ranking$dGD, unname(d))
  expect_equal(r$ranking$MK, unname(rowMeans(k)[r$ranking$id]))
  expect_equal(r$removal.set$id,
               c("CB_AB_01", "CB_AB_02", "CB_X_03", "CB_AB_05", "CB_X_04",
                 "CB_X_01", "CB_AE_03"))
  expect_equal(r$removal.set$gd.after[7], 0.9419, tolerance = 1e-4)
  expect_equal(r$removal.set$gd.after,
               vapply(seq_len(7), function(s)
                 g(setdiff(indNames(cb), r$removal.set$id[1:s])),
                 numeric(1)))
})

test_that("self-referenced kinship is an error (#113)", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  expect_error(gl.report.ind.remove(cb, verbose = 0), "wider reference")
})

test_that("NA kinship is ignored with a warning; greedy matches brute force", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  kin <- gl.kin(testset2.gl, verbose = 0)
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  ids <- indNames(cb)
  kin[ids[1], ids[2]] <- kin[ids[2], ids[1]] <- NA
  # before the review: "missing value where TRUE/FALSE needed"
  out <- capture.output(r <- gl.report.ind.remove(cb, kin = kin, verbose = 1))
  expect_true(any(grepl("1 missing kinship values", out)))
  expect_true(any(grepl("24 individuals have call rate below 0.8", out)))
  expect_false(anyNA(r$ranking$dGD))
  # the running-sum greedy equals a brute-force greedy on the NA matrix
  k <- kin[ids, ids]
  g <- function(keep) 1 - mean(k[keep, keep], na.rm = TRUE)
  rem <- ids
  cur <- g(rem)
  bf <- character(0)
  repeat {
    gains <- vapply(rem, function(i) g(setdiff(rem, i)), numeric(1))
    if (max(gains) <= cur) break
    bf <- c(bf, rem[which.max(gains)])
    cur <- max(gains)
    rem <- setdiff(rem, bf)
  }
  expect_equal(r$removal.set$id, bf)
  expect_equal(r$removal.set$gd.after[length(bf)], cur)
})

test_that("invalid n.best is an error; n.best caps", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  kin <- gl.kin(testset2.gl, verbose = 0)
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  # before the review n.best = 0 removed without limit (7 animals)
  expect_error(gl.report.ind.remove(cb, kin = kin, n.best = 0, verbose = 0),
               "n.best")
  expect_error(gl.report.ind.remove(cb, kin = kin, n.best = "2", verbose = 0),
               "n.best")
  r2 <- gl.report.ind.remove(cb, kin = kin, n.best = 2, verbose = 0)
  expect_equal(nrow(r2$removal.set), 2)
})

test_that("SilicoDArT runs", {
  skip_without_testset2()
  cb.gs <- gl.keep.pop(testset2.gs, pop.list = "EmmacCaptBred", verbose = 0)
  r <- gl.report.ind.remove(cb.gs, kin = gl.kin(testset2.gs, verbose = 0),
                            n.best = 3, verbose = 0)
  expect_equal(nrow(r$ranking), nInd(cb.gs))
  expect_lte(nrow(r$removal.set), 3)
})
