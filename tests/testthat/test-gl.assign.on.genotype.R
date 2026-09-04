test_that("gl.assign.on.genotype returns a dartR genlight with the unknown plus qualifying populations (baseline snapshot)", {
  sub <- testset.gl[1:30, 1:200]

  # Baseline captured 2026-09-04 against dartR.captive commit ad61194.
  # This snapshots CURRENT behaviour, not asserted-correct behaviour.
  res <- gl.assign.on.genotype(sub, unknown = "UC_00126", nmin = 5, verbose = 0)

  expect_s4_class(res, "dartR")
  expect_equal(nInd(res), 11)
  expect_equal(nLoc(res), 10)
  expect_setequal(as.character(pop(res)), c(rep("EmmacBurdMist", 10), "unknown"))
  expect_true(all(ploidy(res) == 2))
  expect_true("UC_00126" %in% indNames(res))
})

test_that("gl.assign.on.genotype errors when unknown is not in the genlight object", {
  sub <- testset.gl[1:30, 1:200]

  expect_error(
    gl.assign.on.genotype(sub, unknown = "not_a_real_individual", verbose = 0),
    "not present in the dataset"
  )
})

test_that("gl.assign.on.genotype errors when pop(x) is NULL", {
  sub <- testset.gl[1:30, 1:200]
  pop(sub) <- NULL

  expect_error(
    gl.assign.on.genotype(sub, unknown = "UC_00126", verbose = 0),
    "Population assignments"
  )
})

test_that("gl.assign.on.genotype errors on duplicate individual names", {
  sub <- testset.gl[1:30, 1:200]
  indNames(sub)[2] <- indNames(sub)[1]

  expect_error(
    gl.assign.on.genotype(sub, unknown = indNames(sub)[1], verbose = 0),
    "Duplicate individual names"
  )
})

test_that("gl.assign.on.genotype returns the unknown individual only, with a warning, when no population meets aic.threshold (F1 fix)", {
  # Approved fix: matches the documented @return contract. Previously
  # crashed inside gl.keep.pop with an opaque "no populations listed to
  # keep!" error.
  sub <- testset.gl[1:30, 1:200]

  # dartR warnings here are cat(warn(...)) console output, not R condition
  # objects, so this is checked via captured output rather than expect_warning.
  out <- capture.output(
    res <- gl.assign.on.genotype(sub, unknown = "UC_00126", nmin = 5, aic.threshold = 1, verbose = 1)
  )
  expect_true(any(grepl("No population met the AIC weight threshold", out)))
  expect_equal(nInd(res), 1)
  expect_equal(indNames(res), "UC_00126")
})

test_that("gl.assign.on.genotype caps n.best at the number of candidate populations instead of crashing (F2 fix)", {
  # Approved fix: n.best = 5 with only 2 candidate populations available
  # now retains both (capped) instead of crashing on NA indices produced
  # by unchecked 1:n.best indexing.
  sub <- testset.gl[1:30, 1:200]

  res <- gl.assign.on.genotype(sub, unknown = "UC_00126", nmin = 5, n.best = 5, verbose = 0)

  expect_equal(nInd(res), 20)
  expect_setequal(popNames(res), c("EmmacBurdMist", "EmmacCoopEulb", "unknown"))
})

test_that("gl.assign.on.genotype rejects SilicoDArT (ploidy = 1) input with a clear error (F3 fix)", {
  # Approved fix: the genotype-likelihood model is diploid-HWE-only and is
  # not valid for presence/absence data. Previously ran to completion and
  # returned a genlight object with no indication the statistics were
  # invalid for this data type.
  sub <- testset.gs[1:30, 1:200]

  expect_error(
    gl.assign.on.genotype(sub, unknown = indNames(sub)[2], nmin = 5, verbose = 0),
    "not valid for"
  )
})
