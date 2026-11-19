# Slow file: skipped on CRAN to keep the check time under 10 minutes; runs
# with devtools::test() and on GitHub Actions.
skip_on_cran()

# File-backed objects: gl.grm computes rrBLUP::A.mat's default relationship
# matrix in blocks of loci (dartR.base::utils.fbm.amat / utils.fbm.grm)
# without decoding the genotype matrix; it agrees with the gen-backed result
# (rrBLUP::A.mat) up to the rounding of the summed cross-product. Other A.mat
# options convert the object explicitly, with a warning, and then give the
# gen-backed result exactly.

x <- gl.compliance.check(platypus.gl[, 1:300], verbose = 0)
f <- suppressWarnings(gl.gen2fbm(x, verbose = 0))

grm <- function(obj, ...) {
  out <- utils::capture.output(
    g <- suppressWarnings(gl.grm(obj, plotheatmap = FALSE, verbose = 3, ...)))
  list(g = g, out = out)
}

test_that("gl.grm: file-backed equals gen-backed (default A.mat)", {
  skip_if_not_installed("bigstatsr")
  skip_if_not_installed("rrBLUP")
  old <- options(bigstatsr.block.sizeGB = getOption("bigstatsr.block.sizeGB"))
  on.exit(options(old), add = TRUE)
  for (bs in c(1, 1e-6)) {                   # 1e-6: one locus per block
    options(bigstatsr.block.sizeGB = bs)
    for (args in list(list(), list(min.MAF = 0.05))) {
      a <- do.call(grm, c(list(x), args))
      b <- do.call(grm, c(list(f), args))
      expect_equal(b$g, a$g, tolerance = 1e-12)
      expect_identical(attributes(b$g), attributes(a$g))
      expect_identical(b$out, a$out)
    }
  }
})

test_that("gl.grm never decodes a file-backed object for default A.mat", {
  skip_if_not_installed("bigstatsr")
  skip_if_not_installed("rrBLUP")
  local_mocked_bindings(.fbm_as_matrix = function(x) {
    stop("the whole FBM was decoded")
  }, .package = "dartR.base")
  expect_no_error(grm(f))
})

test_that("gl.grm converts explicitly for other A.mat options", {
  skip_if_not_installed("bigstatsr")
  skip_if_not_installed("rrBLUP")
  a <- grm(x, shrink = TRUE)
  b <- grm(f, shrink = TRUE)
  expect_identical(b$g, a$g)
  expect_true(any(grepl("whole genotype matrix", b$out)))
  expect_false(any(grepl("whole genotype matrix", a$out)))
})

# Exact cases: counts read in blocks of loci, identical results
test_that("gl.report.parent.offspring and gl.report.errorrate: identical", {
  skip_if_not_installed("bigstatsr")
  pdf(NULL); on.exit(dev.off(), add = TRUE)
  old <- options(bigstatsr.block.sizeGB = getOption("bigstatsr.block.sizeGB"))
  on.exit(options(old), add = TRUE)
  y <- gl.compliance.check(platypus.gl[1:40, 1:400], verbose = 0)
  m <- as.matrix(y)
  m[2, 1:150] <- NA                          # a poorly scored individual
  y@gen <- new("genlight", m, ploidy = ploidy(y))@gen
  fy <- suppressWarnings(gl.gen2fbm(y, verbose = 0))
  reps <- rep(1:20, 2)                       # 20 "replicate" pairs
  calls <- list(
    po = function(z) gl.report.parent.offspring(z,
                                                verbose = 3),
    er = function(z) gl.report.errorrate(z, replicates = reps,
                                         plot.display = FALSE, verbose = 3)
  )
  for (bs in c(1, 1e-6)) {
    options(bigstatsr.block.sizeGB = bs)
    for (nm in names(calls)) {
      g.out <- utils::capture.output(g <- suppressWarnings(calls[[nm]](y)))
      f.out <- utils::capture.output(f <- suppressWarnings(calls[[nm]](fy)))
      if (is.list(g) && !is.data.frame(g)) {
        g[vapply(g, inherits, logical(1), "ggplot")] <- NULL
        f[vapply(f, inherits, logical(1), "ggplot")] <- NULL
      }
      expect_identical(f, g, label = paste(nm, bs))
      expect_identical(f.out, g.out, label = paste(nm, bs, "output"))
    }
  }
  local_mocked_bindings(.fbm_as_matrix = function(x) {
    stop("the whole FBM was decoded")
  }, .package = "dartR.base")
  for (nm in names(calls)) {
    expect_no_error(utils::capture.output(suppressWarnings(calls[[nm]](fy))),
                    message = nm)
  }
})

# Full-matrix case: explicit decode with a memory guard
test_that("gl.parentage.power decodes explicitly, within the memory limit", {
  skip_if_not_installed("bigstatsr")
  gl <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                    verbose = 0)
  fgl <- suppressWarnings(gl.gen2fbm(gl, verbose = 0))
  run <- function(z, verbose = 0) {
    set.seed(1)
    gl.parentage.power(z, n.offspring = 20, plot.display = FALSE,
                       verbose = verbose)
  }
  g <- run(gl)
  f <- run(fgl)
  g$plot <- f$plot <- NULL
  expect_identical(f, g)
  out <- utils::capture.output(run(fgl, verbose = 2))
  expect_true(any(grepl("whole genotype matrix", out)))
  old <- options(dartR_fbm_max_decode_gb = 1e-9)  # below the test size
  on.exit(options(old), add = TRUE)
  expect_error(run(fgl), "dartR_fbm_max_decode_gb")
  expect_no_error(run(gl))                   # gen-backed: no limit
})
