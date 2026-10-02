# Tests for gl.parentage.power(). Run inside the package (report()/warn() and the
# utils.* helpers resolve from the dartR.base namespace there). For a standalone
# run, bind those helpers into the environment first (see dev notes).

test_that("gl.parentage.power returns the expected structure", {
  gl <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                    verbose = 0)
  out <- gl.parentage.power(gl, n.offspring = 20, plot.display = FALSE,
                            verbose = 0)
  expect_type(out, "list")
  expect_named(out, c("power", "settings", "plot"))
  expect_s3_class(out$power, "data.frame")
  expect_true(all(c("n.loci", "true.retained", "correct.unique", "false.rate")
                  %in% colnames(out$power)))
  expect_s3_class(out$plot, "ggplot")
})

test_that("power metrics are valid proportions", {
  gl <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                    verbose = 0)
  out <- gl.parentage.power(gl, n.offspring = 20, plot.display = FALSE,
                            verbose = 0)
  expect_true(all(out$power$true.retained >= 0 & out$power$true.retained <= 1))
  expect_true(all(out$power$correct.unique >= 0 & out$power$correct.unique <= 1))
  expect_true(all(out$power$false.rate >= 0 & out$power$false.rate <= 1))
})

test_that("the true parent is always retained under zero error", {
  gl <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                    verbose = 0)
  out <- gl.parentage.power(gl, n.offspring = 30, error.rate = 0,
                            missing.rate = 0, plot.display = FALSE, verbose = 0)
  # With no error and a zero mismatch tolerance the true parent/pair can never
  # be excluded by Mendelian incompatibility.
  expect_true(all(out$power$true.retained == 1))
})

test_that("both single-parent and pair modes run", {
  gl <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                    verbose = 0)
  op <- gl.parentage.power(gl, n.offspring = 15, pairs = TRUE,
                           plot.display = FALSE, verbose = 0)
  os <- gl.parentage.power(gl, n.offspring = 15, pairs = FALSE,
                           plot.display = FALSE, verbose = 0)
  expect_true(op$settings$pairs)
  expect_false(os$settings$pairs)
})

test_that("engine and settings are recorded", {
  gl <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                    verbose = 0)
  out <- gl.parentage.power(gl, n.offspring = 10, plot.display = FALSE,
                            verbose = 0)
  expect_equal(out$settings$engine, "exclusion")
  expect_equal(out$settings$n.offspring, 10)
})

test_that("fewer than two candidates is a fatal error", {
  gl <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                    verbose = 0)
  expect_error(gl.parentage.power(gl[1, ], n.offspring = 5, verbose = 0))
})
