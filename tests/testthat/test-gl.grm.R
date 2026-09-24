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

  # Snapshot of CURRENT behaviour, not asserted-correct behaviour.
  # Values re-anchored 2026-09-24 to testset.gl from dartR.data >= 1.2.4
  # (1.2.4 changed genotypes in this subset; the 2026-08-26 values came
  # from an older testset.gl). The fingerprint separates a change in the
  # test data from a change in gl.grm.
  m <- as.matrix(sub)
  expect_equal(c(sum(m, na.rm = TRUE), sum(is.na(m))), c(513, 111),
               label = "testset.gl[1:10, 1:100] fingerprint (testset.gl changed?)")
  expect_equal(unname(G[1, 1]), 0.02240585, tolerance = 1e-6)
  expect_equal(unname(G[2, 2]), 0.96305974, tolerance = 1e-6)
  expect_equal(unname(G[1, 2]), -0.03638501, tolerance = 1e-6)
  expect_equal(range(diag(G)), c(0.0158735, 0.9630597), tolerance = 1e-5)
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

test_that("single-copy loci are kept, so the matrix is platform independent", {
  # value from x86 Linux/Windows; arm64 macOS gave 0.943133 before the
  # tolerance because mean() put some single-copy loci below 1/(2n)
  x <- platypus.gl[1:12, 1:200]
  G <- gl.grm(x, plotheatmap = FALSE, verbose = 0)
  expect_equal(round(unname(G["T27", c("T27", "T35", "SDS4", "SDS12")]), 6),
               c(0.985509, -0.173415, -0.098468, -0.090907))
  # a user-supplied min.MAF is passed through unchanged
  G2 <- gl.grm(x, plotheatmap = FALSE, verbose = 0, min.MAF = 0.1)
  expect_equal(unclass(G2)[, ], rrBLUP::A.mat(as.matrix(x) - 1,
                                               min.MAF = 0.1)[, ],
               ignore_attr = TRUE)
})
