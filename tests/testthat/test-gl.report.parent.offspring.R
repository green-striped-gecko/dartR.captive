# Tests of the reviewed gl.report.parent.offspring
# (function-review/reports/dartR.captive/gl.report.parent.offspring.md).

test_that("platypus.gl: same pairs as before review, proportions, p (1, 3)", {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  res <- gl.report.parent.offspring(platypus.gl, verbose = 0)
  expect_equal(nrow(res), 10)
  expect_equal(names(res), c("Outlier", "ind1", "ind2", "n.loci", "prop",
                             "zscore", "p"))
  expect_equal(sort(unique(res$Outlier)), c(0, 1, 3, 4))
  expect_setequal(res$ind1[res$Outlier == 0], c("T5", "T42"))
  expect_equal(res$prop, res$Outlier / res$n.loci)
  # p is the standard-normal lower tail of the standardised proportion
  expect_equal(res$p, pnorm(res$zscore))
  expect_equal(res$p[1], 1.06002e-04, tolerance = 1e-4)
})

test_that("counts match an independent matrix computation (4)", {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- gl.filter.rdepth(gl.filter.reproducibility(platypus.gl, threshold = 1,
                                                  verbose = 0),
                        lower = 12, verbose = 0)
  m <- as.matrix(x)
  h0 <- (m == 0) * 1
  h2 <- (m == 2) * 1
  h0[is.na(h0)] <- 0
  h2[is.na(h2)] <- 0
  cnt <- h0 %*% t(h2) + h2 %*% t(h0)
  res <- gl.report.parent.offspring(platypus.gl, verbose = 0)
  expect_equal(res$Outlier, cnt[cbind(res$ind1, res$ind2)])
})

test_that("range below 1.5 flags more pairs (2)", {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  n <- vapply(c(0.5, 1, 1.5, 3), function(r) {
    nrow(gl.report.parent.offspring(platypus.gl, range = r, verbose = 0))
  }, numeric(1))
  expect_equal(n, c(218, 55, 10, 0))
})

test_that("missing data no longer makes an individual everyone's relative (3)", {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- platypus.gl
  g <- as.matrix(x)
  set.seed(1)
  g["T10", sample(nLoc(x), 800)] <- NA
  xx <- new("genlight", g, ind.names = indNames(x), loc.names = locNames(x),
            ploidy = 2)
  xx@other <- x@other
  res <- gl.report.parent.offspring(xx, verbose = 0)
  # 44 of 54 flagged pairs involved T10 before the review
  expect_lt(sum(res$ind1 == "T10" | res$ind2 == "T10"), 10)
})

test_that("no pairs returns an empty data frame with all columns", {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  res <- gl.report.parent.offspring(testset.gl[1:10, 1:100], verbose = 0)
  expect_equal(nrow(res), 0)
  expect_equal(names(res), c("Outlier", "ind1", "ind2", "n.loci", "prop",
                             "zscore", "p"))
})

test_that("SilicoDArT errors (5); nothing printed at verbose 0 (6)", {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_error(gl.report.parent.offspring(testset.gs[1:10, 1:50],
                                          verbose = 0),
               "Only SNP data")
  out <- capture.output(
    res <- gl.report.parent.offspring(testset.gl[1:10, 1:100], verbose = 0)
  )
  expect_length(out, 0)
})
