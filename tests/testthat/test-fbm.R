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
