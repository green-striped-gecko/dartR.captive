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

test_that("gl.relatedness returns symmetric per-estimator matrices keyed by name", {
  skip_if_not_installed("dartR.coancestry"); skip_if_not_installed("dartR.data")
  gl  <- dartR.data::platypus.gl[1:40, 1:1000]
  got <- gl.relatedness(gl, estimators = "wang", plot.out = FALSE, verbose = 0)
  M <- got$wang
  expect_equal(dim(M), c(40L, 40L))
  expect_identical(rownames(M), indNames(gl))
  expect_identical(colnames(M), indNames(gl))
  expect_true(all(is.na(diag(M))))
  expect_equal(M[lower.tri(M)], t(M)[lower.tri(M)])          # symmetric
  d <- got$dyads
  i <- match(d$ind1[1], indNames(gl)); j <- match(d$ind2[1], indNames(gl))
  expect_equal(M[i, j], d$wang[1])                            # matrix matches long table
  expect_true(all(d$ind1 %in% indNames(gl)))                 # ind cols are names now
})

test_that("gl.relatedness default estimators and optional frames", {
  skip_if_not_installed("dartR.coancestry"); skip_if_not_installed("dartR.data")
  gl  <- dartR.data::platypus.gl[1:40, 1:1000]
  got <- gl.relatedness(gl, plot.out = FALSE, verbose = 0)
  expect_true(all(c("wang","lynchli","lynchrd","ritland","quellergt","loiselle")
                  %in% names(got)))
  expect_true("dyads" %in% names(got))
  expect_null(got$delta19); expect_null(got$trio_delta); expect_null(got$inbreeding)
})

test_that("gl.relatedness pre-filters and produces CIs when bootstrapping", {
  skip_if_not_installed("dartR.coancestry"); skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:40, 1:1000]   # contains 405 monomorphic/all-NA loci
  expect_warning(
    got <- gl.relatedness(gl, estimators = "wang", n.boots = 20,
                          plot.out = FALSE, verbose = 1),
    regexp = "pre-filter|monomorphic")
  expect_true(all(c("wang_lo", "wang_hi") %in% names(got$dyads)))
})
