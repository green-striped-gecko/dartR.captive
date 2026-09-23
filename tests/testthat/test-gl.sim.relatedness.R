# Tests of the reviewed gl.sim.relatedness
# (function-review/reports/dartR.captive/gl.sim.relatedness.md).
# Needs the EMIBD9 binary; set DARTR_EMIBD9_PATH to its folder
# (default ~/programs/emibd9-v1.0).

emibd9_dir <- function() {
  Sys.getenv("DARTR_EMIBD9_PATH", path.expand("~/programs/emibd9-v1.0"))
}
skip_without_emibd9 <- function() {
  exe <- if (Sys.info()[["sysname"]] == "Windows") "EM_IBD_P.exe" else
    "EM_IBD_P"
  skip_if_not(file.exists(file.path(emibd9_dir(), exe)), "EMIBD9 not found")
}
# an unrelated population without missing data, independent of the
# dartR.data version of testset.gl
sim_population <- function() {
  set.seed(10)
  p <- gl.filter.monomorphs(gl.filter.allna(platypus.gl, verbose = 0),
                            verbose = 0)
  dartR.sim::gl.sim.ind(p, n = 40, popname = "sim")
}

test_that("returns values, mean, quantile interval, CI of mean, plot (2, 3)", {
  skip_without_emibd9()
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- sim_population()
  set.seed(1)
  out <- capture.output(
    res <- gl.sim.relatedness(x, rel = "half.sib", nboots = 6,
                              emibd9.path = emibd9_dir(),
                              plot.out = FALSE, verbose = 0)
  )
  expect_named(res, c("values", "mean", "interval", "ci.mean", "plot"))
  expect_length(res$values, 6)
  expect_equal(res$mean, mean(res$values))
  expect_equal(unname(res$interval),
               unname(quantile(res$values, c(0.025, 0.975))))
  expect_s3_class(res$plot, "ggplot")
  # nothing printed at verbose = 0, including EMIBD9 and gl.sim.offspring
  expect_length(out, 0)
})

test_that("full sibs are two offspring of the same parents (1)", {
  skip_without_emibd9()
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- sim_population()
  set.seed(1)
  res <- gl.sim.relatedness(x, rel = "full.sib", nboots = 6,
                            emibd9.path = emibd9_dir(),
                            plot.out = FALSE, verbose = 0)
  # expected kinship 0.25; parent-offspring averaging gave sd ~0.008
  expect_gt(res$mean, 0.2)
  expect_lt(res$mean, 0.3)
  expect_gt(sd(res$values), 0.01)
})

test_that("first cousins sit near 0.0625 (6)", {
  skip_without_emibd9()
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- sim_population()
  set.seed(1)
  res <- gl.sim.relatedness(x, rel = "first.cousin", nboots = 6,
                            emibd9.path = emibd9_dir(),
                            plot.out = FALSE, verbose = 0)
  expect_gt(res$mean, 0.02)
  expect_lt(res$mean, 0.11)
})

test_that("rel is validated and SilicoDArT errors (6)", {
  expect_error(gl.sim.relatedness(testset.gl[1:5, 1:10], rel = "fullsib",
                                  verbose = 0), "should be one of")
  expect_error(gl.sim.relatedness(testset.gs[1:5, 1:10], verbose = 0),
               "Only SNP data")
})
