test_that("gl.assign.grm returns a named numeric vector, one entry per candidate population", {
  sub <- testset.gl[1:30, 1:200]
  res <- gl.assign.grm(sub, unknown = "UC_00126", verbose = 0)

  expect_type(res, "double")
  expect_null(dim(res))
  # UC_00126's own population (EmmacMaclGeor) has only itself as a member,
  # so it is silently absent from the result once removed.
  expect_false("EmmacMaclGeor" %in% names(res))
  expect_setequal(
    names(res),
    c("EmmacBurdMist", "EmmacMDBCond", "EmmacMDBBowm", "EmmacMDBForb",
      "EmmacMDBMaci", "EmmacBurnBara", "EmmacMDBSanf", "EmmacCoopEulb")
  )
})

test_that("gl.assign.grm output is numerically stable across repeat runs (baseline snapshot)", {
  sub <- testset.gl[1:30, 1:200]
  res <- gl.assign.grm(sub, unknown = "UC_00126", verbose = 0)

  # Baseline captured 2026-09-01 against dartR.captive commit 571d6a6.
  # This snapshots CURRENT behaviour, not asserted-correct behaviour.
  expect_equal(unname(res["EmmacBurdMist"]), 0.24418249, tolerance = 1e-6)
  expect_equal(unname(res["EmmacCoopEulb"]), -0.44421730, tolerance = 1e-6)
  expect_equal(unname(res[1]), max(res), tolerance = 1e-6) # sorted decreasing
})

test_that("gl.assign.grm errors when unknown is not in the genlight object", {
  sub <- testset.gl[1:30, 1:200]

  expect_error(
    gl.assign.grm(sub, unknown = "not_a_real_individual", verbose = 0),
    "not in the genlight object"
  )
})

test_that("gl.assign.grm errors on SilicoDArT input via the inherited gl.grm gate", {
  sub <- testset.gs[1:10, 1:100]

  expect_error(
    gl.assign.grm(sub, unknown = indNames(sub)[1], verbose = 0),
    "SilicoDArT"
  )
})

test_that("gl.assign.grm errors clearly when pop(x) is NULL (F1 fix)", {
  # Approved fix: gl.assign.grm now guards is.null(pop(x)) before doing any
  # work, matching gl.assign.pa, instead of failing deep inside order()
  # with an opaque "argument 1 is not a vector" error.
  sub <- testset.gl[1:30, 1:200]
  pop(sub) <- NULL

  expect_error(
    gl.assign.grm(sub, unknown = "AA019237", verbose = 0),
    "Population assignments"
  )
})

test_that("gl.assign.grm errors on duplicate individual names (F4 fix)", {
  # Approved fix: gl.assign.grm now guards duplicated(indNames(x)) before
  # doing any work, matching gl.assign.pa, instead of silently folding two
  # different individuals into population 'unknown'.
  sub <- testset.gl[1:30, 1:200]
  indNames(sub)[2] <- indNames(sub)[1]

  expect_error(
    gl.assign.grm(sub, unknown = indNames(sub)[1], verbose = 0),
    "Duplicate individual names"
  )
})
