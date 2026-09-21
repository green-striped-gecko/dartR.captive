test_that("gl.grm returns a square matrix sized to nInd, no plot", {
  sub <- testset.gl[1:10, 1:100]
  G <- gl.grm(sub, plotheatmap = FALSE, verbose = 0)

  expect_equal(dim(G), c(10, 10))
  expect_equal(rownames(G), indNames(sub))
  expect_true(isSymmetric(unname(G)))
})

test_that("gl.grm output is numerically stable across repeat runs (baseline snapshot)", {
  sub <- testset.gl[1:10, 1:100]
  G <- gl.grm(sub, plotheatmap = FALSE, verbose = 0)

  # Baseline captured 2026-08-26 against dartR.captive commit bdb623d.
  # This snapshots CURRENT behaviour, not asserted-correct behaviour.
  expect_equal(unname(G[1, 1]), 0.10859585, tolerance = 1e-6)
  expect_equal(unname(G[2, 2]), 2.10321355, tolerance = 1e-6)
  expect_equal(unname(G[1, 2]), -0.17634953, tolerance = 1e-6)
  expect_equal(range(diag(G)), c(0.07693525, 4.211809), tolerance = 1e-5)
})

test_that("gl.grm errors on SilicoDArT data instead of returning a meaningless matrix (F1 fix)", {
  # Approved fix: gl.grm now gates on datatype and refuses SilicoDArT
  # (presence/absence, ploidy 1) input, since the additive-relationship
  # algorithm and the documented 1..2 diagonal range only hold for SNP
  # dosage data.
  sub <- testset.gs[1:10, 1:100]

  expect_error(
    gl.grm(sub, plotheatmap = FALSE, verbose = 0),
    "SilicoDArT"
  )
})

test_that("gl.grm no longer errors when plotheatmap = FALSE and plot.file is set (F2 fix)", {
  # Approved fix: the plot-save block now only runs inside
  # `if (plotheatmap == TRUE)`, so `p3` is always defined where it's used.
  # plotheatmap = FALSE + plot.file set now computes and returns the matrix
  # and warns (at verbose >= 1) that nothing was saved.
  sub <- testset.gl[1:10, 1:100]

  msg <- capture.output(
    G <- gl.grm(sub, plotheatmap = FALSE, plot.file = "grm_test",
                plot.dir = tempdir(), verbose = 1)
  )
  expect_true(any(grepl("plotheatmap = FALSE", msg)))
  expect_equal(dim(G), c(10, 10))
  expect_false(file.exists(file.path(tempdir(), "grm_test.RDS")))
})

test_that("gl.grm assigns a default population when none is set", {
  sub <- testset.gl[1:10, 1:100]
  pop(sub) <- NULL

  G <- gl.grm(sub, plotheatmap = FALSE, verbose = 0)
  expect_equal(dim(G), c(10, 10))
})
