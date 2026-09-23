# Tests of the reviewed gl.run.colony
# (function-review/reports/dartR.captive/gl.run.colony.md). Needs the COLONY2 binary;
# set DARTR_COLONY_PATH to its folder (default ~/programs).

colony_dir <- function() {
  Sys.getenv("DARTR_COLONY_PATH", path.expand("~/programs"))
}
skip_without_colony <- function() {
  exe <- switch(Sys.info()[["sysname"]],
                Darwin = "colony2s.out",
                Linux = "colony2s.ifort.out",
                Windows = "Colony2p.exe")
  skip_if_not(file.exists(file.path(colony_dir(), exe)), "COLONY not found")
}
run_fixture <- function() {
  x <- testset.gl[1:20, 1:30]
  x@other$ind.metrics$offspring <- rep(c("yes", "no", "no"), c(10, 5, 5))
  x@other$ind.metrics$father <- rep(c("no", "yes", "no"), c(10, 5, 5))
  x@other$ind.metrics$mother <- rep(c("no", "no", "yes"), c(10, 5, 5))
  x
}

test_that("change 4/5: output lands in outpath; returns files and BestConfig", {
  skip_without_colony()
  wd <- withr::local_tempdir()
  withr::local_dir(wd)
  out <- withr::local_tempdir()
  res <- gl.run.colony(run_fixture(), colony.path = colony_dir(),
                       outpath = out, seed = 1, length.run = 1, verbose = 0)
  expect_named(res, c("files", "best.config"))
  expect_true(file.exists(file.path(out, "my_project.BestConfig")))
  expect_false(file.exists(file.path(wd, "my_project.BestConfig")))
  expect_true(all(normalizePath(dirname(res$files)) == normalizePath(out)))
  expect_equal(colnames(res$best.config),
               c("OffspringID", "FatherID", "MotherID", "CloneIndex",
                 "ClusterIndex"))
  expect_equal(nrow(res$best.config), 10)
  expect_equal(normalizePath(getwd()), normalizePath(wd))
})

test_that("change 1: default outpath = NULL runs", {
  skip_without_colony()
  withr::local_dir(withr::local_tempdir())
  res <- gl.run.colony(run_fixture(), colony.path = colony_dir(), seed = 1,
                       length.run = 1, verbose = 0)
  expect_equal(nrow(res$best.config), 10)
})

test_that("change 3: outpath with a space runs", {
  skip_without_colony()
  withr::local_dir(withr::local_tempdir())
  out <- file.path(withr::local_tempdir(), "my out")
  dir.create(out)
  res <- gl.run.colony(run_fixture(), colony.path = colony_dir(),
                       outpath = out, seed = 1, length.run = 1, verbose = 0)
  expect_equal(nrow(res$best.config), 10)
})

test_that("change 2: stale output files do not hide a failed run", {
  skip_without_colony()
  withr::local_dir(withr::local_tempdir())
  out <- withr::local_tempdir()
  file.create(file.path(out, "my_project.BestConfig"))
  Sys.setFileTime(file.path(out, "my_project.BestConfig"),
                  Sys.time() - 3600)
  # an empty offspring genotype block makes COLONY stop on a data error
  x <- run_fixture()
  local_mocked_bindings(
    gl2colony = function(...) {
      f <- file.path(out, "colony2.dat")
      writeLines("broken", f)
      f
    }
  )
  expect_error(
    gl.run.colony(x, colony.path = colony_dir(), outpath = out, seed = 1,
                  verbose = 0),
    "COLONY"
  )
})

test_that("change 6: missing executable names the expected file", {
  withr::local_dir(withr::local_tempdir())
  expect_error(
    gl.run.colony(run_fixture(), colony.path = tempdir(),
                  outpath = tempdir(), seed = 1, verbose = 0),
    "COLONY executable not found"
  )
})

test_that("change 6: nothing prints at verbose = 0", {
  skip_without_colony()
  withr::local_dir(withr::local_tempdir())
  out <- capture.output(
    gl.run.colony(run_fixture(), colony.path = colony_dir(),
                  outpath = withr::local_tempdir(), seed = 1,
                  length.run = 1, verbose = 0)
  )
  expect_length(out, 0)
})

test_that("A1: long and spaced names are restored in best.config", {
  skip_without_colony()
  withr::local_dir(withr::local_tempdir())
  x <- run_fixture()
  nm <- indNames(x)
  # offspring 1-2 share their first 20 characters; father 11 has a space
  indNames(x)[1] <- paste0(strrep("A", 20), "_offspring_one")
  indNames(x)[2] <- paste0(strrep("A", 20), "_offspring_two")
  indNames(x)[11] <- "father eleven"
  out <- withr::local_tempdir()
  res <- gl.run.colony(x, colony.path = colony_dir(), outpath = out,
                       seed = 1, length.run = 1, verbose = 0)
  expect_setequal(res$best.config$OffspringID, indNames(x)[1:10])
  expect_true(all(res$best.config$FatherID %in% c(indNames(x)[11:15]) |
                    startsWith(res$best.config$FatherID, "*")))
  map <- read.csv(file.path(out, "my_project.IDmap.csv"))
  expect_equal(map$name, indNames(x))
  expect_true(file.path(out, "my_project.IDmap.csv") %in% res$files)
})

test_that("A1: short names are used as-is and no map is written", {
  skip_without_colony()
  withr::local_dir(withr::local_tempdir())
  out <- withr::local_tempdir()
  res <- gl.run.colony(run_fixture(), colony.path = colony_dir(),
                       outpath = out, seed = 1, length.run = 1, verbose = 0)
  expect_setequal(res$best.config$OffspringID, indNames(run_fixture())[1:10])
  expect_false(file.exists(file.path(out, "my_project.IDmap.csv")))
})
