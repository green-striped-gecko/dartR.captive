test_that("gl.relatedness plumbs through to relatedness_cpp", {
  skip_if_not_installed("dartR.coancestry")
  skip_if_not_installed("dartR.data")
  gl  <- dartR.data::platypus.gl[1:40, 1:1000]
  est <- c("wang", "lynchrd")
  got <- gl.relatedness(gl, estimators = est, plot.out = FALSE, verbose = 0)
  snp <- as.matrix(gl); storage.mode(snp) <- "integer"
  ref <- dartR.coancestry::relatedness_cpp(snp, estimators = est,
            max_dyads = as.integer(choose(40, 2)))
  expect_equal(got$dyads$wang,    ref$dyads$wang)
  expect_equal(got$dyads$lynchrd, ref$dyads$lynchrd)
})

test_that("gl.relatedness rejects an unknown estimator", {
  skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:10, 1:100]
  expect_error(gl.relatedness(gl, estimators = "bogus", verbose = 0),
               regexp = "Unknown estimator")
})
