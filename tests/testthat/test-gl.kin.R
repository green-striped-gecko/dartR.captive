# Tests of the reviewed gl.kin
# (function-review/reports/dartR.captive/gl.kin.md). The EMIBD9 test needs
# the binary; set DARTR_EMIBD9_PATH to its folder
# (default ~/programs/emibd9-v1.0).

emibd9_dir <- function() {
  Sys.getenv("DARTR_EMIBD9_PATH", path.expand("~/programs/emibd9-v1.0"))
}
skip_without_emibd9 <- function() {
  exe <- if (Sys.info()[["sysname"]] == "Windows") "EM_IBD_P.exe" else
    "EM_IBD_P"
  skip_if_not(file.exists(file.path(emibd9_dir(), exe)), "EMIBD9 not found")
}
skip_without_testset2 <- function() {
  skip_if_not(exists("testset2.gl") && exists("testset2.gs"),
              "testset2.gl/gs need dartR.data >= 1.2.5")
}
# captive-bred offspring of testset2.gl paired with their sampled sire
sire_pairs <- function(x) {
  im <- x@other$ind.metrics
  ok <- im$sire %in% indNames(x) & im$id %in% indNames(x)
  cbind(as.character(im$id[ok]), as.character(im$sire[ok]))
}

test_that("grm: contract shape, attributes, input untouched", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  x <- testset2.gl
  h <- length(x@other$history)
  kin <- gl.kin(x, verbose = 0)
  expect_true(is.matrix(kin) && is.numeric(kin))
  expect_identical(dimnames(kin), list(indNames(x), indNames(x)))
  expect_equal(attr(kin, "method"), "grm")
  expect_equal(attr(kin, "datatype"), "SNP")
  expect_equal(attr(kin, "nLoc"), nLoc(x))
  expect_equal(attr(kin, "scale"), "kinship")
  expect_equal(length(x@other$history), h)
})

test_that("grm: kinship is G / 2, diagonal included", {
  skip_if_not_installed("rrBLUP")
  x <- gl.filter.allna(platypus.gl, verbose = 0)
  kin <- gl.kin(x, verbose = 0)
  G <- gl.grm(x, plotheatmap = FALSE, verbose = 0)
  expect_equal(as.vector(kin), as.vector(G) / 2)
})

test_that("grm: known parent-offspring pairs have kinship near 0.25", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  # before the fix: mean 0.097, none above the first-degree break 0.1875
  x <- gl.filter.callrate(testset2.gl, threshold = 0.95, verbose = 0)
  kin <- gl.kin(x, verbose = 0)
  po <- sire_pairs(x)
  expect_equal(nrow(po), 24)
  expect_equal(mean(kin[po]), 0.25, tolerance = 0.03)
  expect_true(all(kin[po] > 0.1875))
  # unrelated pairs centre on 0, so a virtual offspring is not "outbred"
  expect_lt(abs(median(kin[upper.tri(kin)])), 0.02)
})

test_that("dominant: fixed diagonal, tagged, warns on pairs above 0.5", {
  skip_without_testset2()
  x <- testset2.gs
  out <- capture.output(kin <- gl.kin(x, verbose = 1))
  expect_equal(attr(kin, "method"), "dominant")
  expect_equal(attr(kin, "scale"), "kinship")
  expect_true(all(diag(kin) == 0.5))
  expect_equal(mean(kin[upper.tri(kin)]), -0.002251, tolerance = 1e-4)
  expect_true(any(grepl("16 pair\\(s\\).*above 0.5", out)))
  expect_length(capture.output(gl.kin(x, verbose = 0)), 0)
})

test_that("method validation errors", {
  skip_without_testset2()
  expect_error(gl.kin(testset2.gl, method = "dominant", verbose = 0),
               "presence/absence")
  expect_error(gl.kin(testset2.gs, method = "grm", verbose = 0),
               "requires SNP")
  expect_error(gl.kin(testset2.gl, method = "foo", verbose = 0),
               "must be one of")
})

test_that("verbose = 3 prints mean pairwise kinship, not a constant GD", {
  skip_if_not_installed("rrBLUP")
  x <- gl.filter.allna(platypus.gl, verbose = 0)
  expect_length(capture.output(gl.kin(x, verbose = 0)), 0)
  out <- capture.output(kin <- gl.kin(x, verbose = 3))
  expect_false(any(grepl("Gene diversity", out)))
  line <- grep("Mean pairwise kinship", out, value = TRUE)
  expect_length(line, 1)
  expect_match(line, format(round(mean(kin[upper.tri(kin)]), 4)),
               fixed = TRUE)
})

test_that("emibd9: returns $rel unscaled, found via emibd9.path", {
  skip_without_emibd9()
  skip_without_testset2()
  x <- gl.keep.ind(testset2.gl,
                   indNames(testset2.gl)[pop(testset2.gl) == "EmmacCaptBred"],
                   verbose = 0)
  x <- gl.filter.monomorphs(gl.filter.callrate(x, threshold = 0.9,
                                               verbose = 0), verbose = 0)
  # the working directory is not the binary folder
  wd <- setwd(tempdir())
  on.exit(setwd(wd))
  res <- gl.run.EMIBD9(x, emibd9.path = emibd9_dir(), plot.out = FALSE,
                       verbose = 0)
  kin <- gl.kin(x, method = "emibd9", emibd9.path = emibd9_dir(),
                verbose = 0)
  expect_equal(as.vector(kin), as.vector(res$rel))
  expect_true(all(diag(kin) == 0.5))
  expect_equal(attr(kin, "scale"), "kinship")
  expect_identical(dimnames(kin), list(indNames(x), indNames(x)))
})
