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
  # dartR convention: the pre-filter reports through cat(warn(...)) gated on
  # verbose, it never signals an R warning condition.
  expect_no_warning(
    got <- gl.relatedness(gl, estimators = "wang", n.boots = 20,
                          plot.out = FALSE, verbose = 0))
  expect_true(all(c("wang_lo", "wang_hi") %in% names(got$dyads)))
})

test_that("gl.relatedness pre-filter reports only at verbose >= 1", {
  skip_if_not_installed("dartR.coancestry"); skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:40, 1:1000]
  expect_output(
    gl.relatedness(gl, estimators = "wang", n.boots = 5,
                   plot.out = FALSE, verbose = 1),
    regexp = "pre-filter")
  quiet <- capture.output(
    gl.relatedness(gl, estimators = "wang", n.boots = 5,
                   plot.out = FALSE, verbose = 0))
  expect_false(any(grepl("pre-filter", quiet)))
})

test_that("gl.relatedness plotting is optional; bad plot.stat errors", {
  skip_if_not_installed("dartR.coancestry"); skip_if_not_installed("dartR.data")
  gl  <- dartR.data::platypus.gl[1:40, 1:1000]
  res <- gl.relatedness(gl, estimators = c("wang", "lynchrd"),
                        plot.out = FALSE, verbose = 0)
  expect_true(is.list(res) && "wang" %in% names(res))
  expect_error(
    gl.relatedness(gl, estimators = "wang", plot.stat = "trioml",
                   plot.out = TRUE, verbose = 0),
    regexp = "plot.stat")
})

test_that("gl.relatedness populates and name-maps optional frames", {
  skip_if_not_installed("dartR.coancestry"); skip_if_not_installed("dartR.data")
  gl  <- dartR.data::platypus.gl[1:20, 1:300]
  got <- gl.relatedness(gl, estimators = c("dyadml", "inbreeding"),
                        allow.inbreeding = TRUE, plot.out = FALSE, verbose = 0)
  expect_true("dyadml" %in% names(got))
  expect_false(is.null(got$delta19))
  expect_false(is.null(got$inbreeding))
  expect_true(all(got$delta19$ind1 %in% indNames(gl)))   # mapped to names
  expect_true(all(got$inbreeding$ind %in% indNames(gl)))
  expect_true(all(c("LH", "LR", "dyadml_F") %in% names(got$inbreeding)))
})

test_that("gl.relatedness with only 'inbreeding' returns no dyad matrices", {
  skip_if_not_installed("dartR.coancestry"); skip_if_not_installed("dartR.data")
  gl  <- dartR.data::platypus.gl[1:20, 1:300]
  got <- gl.relatedness(gl, estimators = "inbreeding", plot.out = FALSE, verbose = 0)
  expect_false(is.null(got$inbreeding))
  expect_true(all(!c("wang","lynchli","lynchrd","ritland","quellergt",
                     "loiselle","dyadml","trioml") %in% names(got)))
})

# Tests added by the function review
# (function-review/reports/dartR.captive/gl.relatedness.md).

test_that("agreement with related::coancestry (Coancestry Fortran)", {
  skip_if_not_installed("dartR.coancestry"); skip_if_not_installed("dartR.data")
  skip_if_not_installed("related")
  # related reads and writes files in the working directory
  withr::local_dir(withr::local_tempdir())
  # named through a variable: related is not a declared dependency
  coancestry <- getExportedValue("related", "coancestry")
  gl <- gl.filter.monomorphs(dartR.data::platypus.gl[1:30, ], verbose = 0)
  gl <- gl.filter.callrate(gl, method = "loc", threshold = 1, verbose = 0)
  gl <- gl.filter.monomorphs(gl, verbose = 0)
  est <- c("wang", "lynchli", "lynchrd", "ritland", "quellergt")
  got <- gl.relatedness(gl, estimators = est, plot.out = FALSE, verbose = 0)
  m <- as.matrix(gl)
  gd <- data.frame(id = indNames(gl))
  for (j in seq_len(ncol(m))) {
    gd[[paste0("L", j, "a")]] <- ifelse(m[, j] == 2, 2, 1)
    gd[[paste0("L", j, "b")]] <- ifelse(m[, j] == 0, 1, 2)
  }
  invisible(capture.output(cc <- suppressWarnings(
    coancestry(gd, lynchli = 1, lynchrd = 1, quellergt = 1,
                        ritland = 1, wang = 1))))
  co <- cc$relatedness
  key <- function(a, b) paste(pmin(a, b), pmax(a, b))
  ix <- match(key(got$dyads$ind1, got$dyads$ind2),
              key(co$ind1.id, co$ind2.id))
  # related rounds its output to 4 decimals
  for (e in c("lynchrd", "ritland", "quellergt")) {
    expect_lt(max(abs(got$dyads[[e]] - co[[e]][ix])), 1e-4)
  }
  # wang and lynchli sit above related by an amount that shrinks with n
  # (documented in Details)
  expect_equal(mean(got$dyads$wang - co$wang[ix]), 0.0469, tolerance = 0.01)
  expect_equal(mean(got$dyads$lynchli - co$lynchli[ix]), 0.0183,
               tolerance = 0.01)
})

test_that("output contract: invisible, tagged, NA for pairs without data", {
  skip_if_not_installed("dartR.coancestry"); skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:12, 1:300]
  vis <- withVisible(gl.relatedness(gl, estimators = "wang",
                                    plot.out = FALSE, verbose = 0))
  expect_false(vis$visible)                                   # change 7
  expect_equal(attr(vis$value$wang, "scale"), "relatedness")  # change 4
  # the kinship series halves the tagged matrix
  expect_equal(as.vector(utils.kin.as.kinship(vis$value$wang)),
               as.vector(vis$value$wang) / 2)
  # an individual with no calls: NA in every estimator (was 0 / NaN; change 2)
  gx <- gl
  gx@gen[[1]] <- new("SNPbin", rep(NA_integer_, nLoc(gx)), ploidy = 2L)
  expect_output(r <- gl.relatedness(gx, estimators = c("wang", "dyadml"),
                                    plot.out = FALSE, verbose = 1),
                "11 pairs share no called locus")
  expect_true(all(is.na(r$wang[1, ])))
  expect_true(all(is.na(r$dyadml[1, ])))
  expect_false(anyNA(r$wang[2:12, 2:12][upper.tri(diag(11))]))
  expect_equal(sum(is.na(r$dyads$wang)), 11)
  # with n.boots > 0 the dropped individual keeps an NA row (was 11 x 11)
  rb <- gl.relatedness(gx, estimators = "wang", n.boots = 3,
                       plot.out = FALSE, verbose = 0)
  expect_equal(dim(rb$wang), c(12L, 12L))
  expect_identical(rownames(rb$wang), indNames(gx))
  expect_true(all(is.na(rb$wang[1, ])))
})

test_that("arguments are validated before the engine runs (change 5)", {
  skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:10, 1:100]
  bad <- list(list(n.boots = NA), list(n.boots = -1), list(num.trios = 0),
              list(n.threads = 0), list(rng.seed = NA),
              list(allow.inbreeding = NA), list(plot.out = "yes"),
              list(plot.stat = "trioml"))
  for (b in bad) {
    expect_error(do.call(gl.relatedness,
                         c(list(gl, estimators = "wang", verbose = 0), b)),
                 paste(names(b), "must be"))
  }
})

test_that("the heatmap is saved without being displayed (change 6)", {
  skip_if_not_installed("dartR.coancestry"); skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:10, 1:200]
  d <- withr::local_tempdir()
  n.dev <- length(grDevices::dev.list())
  out <- capture.output(
    gl.relatedness(gl, estimators = "wang", plot.out = FALSE,
                   plot.file = "h", plot.dir = d, verbose = 0))
  expect_length(list.files(d), 1)
  expect_equal(length(grDevices::dev.list()), n.dev)
  expect_false(any(grepl("gl.colors", out)))
})
