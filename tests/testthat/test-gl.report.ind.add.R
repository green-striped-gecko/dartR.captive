# Tests of the reviewed gl.report.ind.add
# (function-review/reports/dartR.captive/gl.report.ind.add.md).

skip_without_testset2 <- function() {
  skip_if_not(exists("testset2.gl") && exists("testset2.gs"),
              "testset2.gl/gs need dartR.data >= 1.2.5")
}

test_that("structure, read-only, values match the formula", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  x0 <- x
  kin <- gl.kin(x, verbose = 0)
  expect_length(capture.output(
    r <- gl.report.ind.add(x, candidates = "EmmacMaclGeor",
                           target.pop = "EmmacCaptBred", kin = kin,
                           verbose = 0)), 0)
  expect_identical(x, x0)
  expect_named(r, c("id", "from", "dgd", "rank"))
  expect_equal(nrow(r), 11)
  expect_equal(r$rank, 1:11)
  expect_false(is.unsorted(rev(r$dgd)))
  expect_equal(r$id[1], "UC_01060")
  t.ids <- indNames(x)[pop(x) == "EmmacCaptBred"]
  f <- function(c) mean(kin[t.ids, t.ids]) -
    mean(kin[c(t.ids, c), c(t.ids, c)])
  expect_equal(r$dgd, unname(vapply(r$id, f, numeric(1))))
})

test_that("recorded founders rank last among wild candidates", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  kin <- gl.kin(x, verbose = 0)
  wild <- indNames(x)[pop(x) != "EmmacCaptBred"]
  r <- gl.report.ind.add(x, candidates = wild, target.pop = "EmmacCaptBred",
                         kin = kin, verbose = 0)
  expect_equal(tail(r$id, 2), c("AA000307", "AA019158"))
})

test_that("NA kinship is ignored with a warning", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  kin <- gl.kin(x, verbose = 0)
  t.ids <- indNames(x)[pop(x) == "EmmacCaptBred"]
  kin[t.ids[1], t.ids[2]] <- kin[t.ids[2], t.ids[1]] <- NA
  # before the review every dgd was NA
  expect_output(r <- gl.report.ind.add(x, candidates = "EmmacMaclGeor",
                                       target.pop = "EmmacCaptBred",
                                       kin = kin, verbose = 1),
                "1 pairs have missing kinship")
  expect_false(anyNA(r$dgd))
  # the shared engine keeps NA propagation by default
  expect_true(is.na(utils.kin.dgd(kin[t.ids, t.ids])))
})

test_that("duplicate candidates are evaluated once", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  kin <- gl.kin(x, verbose = 0)
  r <- gl.report.ind.add(x, candidates = c("AA010915", "AA010915"),
                         target.pop = "EmmacCaptBred", kin = kin, verbose = 0)
  expect_equal(nrow(r), 1)
})

test_that("low call rate warns", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  kin <- gl.kin(x, verbose = 0)
  expect_output(gl.report.ind.add(x, candidates = "EmmacMaclGeor",
                                  target.pop = "EmmacCaptBred", kin = kin,
                                  verbose = 1),
                "24 target or candidate individuals have call rate below 0.8")
})

test_that("SilicoDArT and input errors", {
  skip_without_testset2()
  r <- gl.report.ind.add(testset2.gs, candidates = "EmmacMaclGeor",
                         target.pop = "EmmacCaptBred", verbose = 0)
  expect_equal(nrow(r), 11)
  x <- testset2.gs
  expect_error(gl.report.ind.add(x, candidates = "EmmacMaclGeor",
                                 target.pop = "nope", verbose = 0),
               "target.pop")
  expect_error(gl.report.ind.add(x, candidates = "nobody",
                                 target.pop = "EmmacCaptBred", verbose = 0),
               "not found")
  expect_error(gl.report.ind.add(x, candidates = "EmmacCaptBred",
                                 target.pop = "EmmacCaptBred", verbose = 0),
               "already belong")
})
