# Tests of the scale-tag contract: gl.grm tags "relatedness",
# gl.run.EMIBD9 and gl.kin tag "kinship", and utils.kin.as.kinship /
# utils.kin.check convert by the tag.

toy <- function(scale = NULL) {
  m <- matrix(c(1, 0.5, 0.5, 1), 2, dimnames = list(c("a", "b"), c("a", "b")))
  attr(m, "scale") <- scale
  m
}

test_that("relatedness is halved, kinship and untagged pass unchanged", {
  k <- utils.kin.as.kinship(toy("relatedness"))
  expect_equal(as.vector(k), c(0.5, 0.25, 0.25, 0.5))
  expect_equal(attr(k, "scale"), "kinship")
  expect_equal(as.vector(utils.kin.as.kinship(toy("kinship"))),
               as.vector(toy()))
  k <- utils.kin.as.kinship(toy())
  expect_equal(as.vector(k), as.vector(toy()))
  expect_equal(attr(k, "scale"), "kinship")
  expect_identical(dimnames(k), dimnames(toy()))
})

test_that("unknown scale is a fatal error", {
  expect_error(utils.kin.as.kinship(toy("distance")), "unknown matrix scale")
})

test_that("conversion message only at verbose >= 2", {
  expect_length(
    capture.output(k <- utils.kin.as.kinship(toy("relatedness"))), 0)
  out <- capture.output(k <- utils.kin.as.kinship(toy("relatedness"),
                                             verbose = 2))
  expect_true(any(grepl("halved to kinship", out)))
})

test_that("gl.grm tags relatedness; gl.kin output is kinship", {
  skip_if_not_installed("rrBLUP")
  x <- gl.filter.allna(platypus.gl, verbose = 0)
  G <- gl.grm(x, plotheatmap = FALSE, verbose = 0)
  expect_equal(attr(G, "scale"), "relatedness")
  kin <- gl.kin(x, verbose = 0)
  expect_equal(attr(kin, "scale"), "kinship")
  expect_equal(as.vector(kin), as.vector(G) / 2)
})

test_that("utils.kin.check turns gl.grm output into gl.kin output", {
  skip_if_not_installed("rrBLUP")
  # before: gl.grm output passed as kin was used as kinship (2x too large)
  x <- gl.filter.allna(platypus.gl, verbose = 0)
  G <- gl.grm(x, plotheatmap = FALSE, verbose = 0)
  expect_equal(as.vector(utils.kin.check(x, G)),
               as.vector(gl.kin(x, verbose = 0)))
  attr(G, "scale") <- "distance"
  expect_error(utils.kin.check(x, G), "unknown matrix scale")
})
