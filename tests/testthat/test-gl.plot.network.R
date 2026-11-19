# Tests of the reviewed gl.plot.network
# (function-review/reports/dartR.captive/gl.plot.network.md).

pn_fixture <- function() {
  t <- gl.keep.ind(platypus.gl, ind.list = indNames(platypus.gl)[1:20],
                   verbose = 0)
  list(x = t, G = gl.grm(t, plotheatmap = FALSE, verbose = 0))
}

test_that("similarity input keeps the top alpha proportion of pairs (2, 7)", {
  skip_if_not_installed("igraph")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  f <- pn_fixture()
  net <- gl.plot.network(f$G, f$x, alpha = 0.05, verbose = 0)
  expect_s3_class(net, "igraph")
  g <- f$G[upper.tri(f$G)]
  expect_equal(igraph::ecount(net), sum(g >= quantile(g, 0.95)))
  expect_true(all(igraph::E(net)$weight >= quantile(g, 0.95)))
})

test_that("distance input keeps the closest pairs (2)", {
  skip_if_not_installed("igraph")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  f <- pn_fixture()
  D <- gl.dist.ind(f$x, method = "euclidean", plot.display = FALSE,
                   verbose = 0)
  d <- as.matrix(D)[upper.tri(as.matrix(D))]
  net <- gl.plot.network(D, f$x, alpha = 0.05, type = "distance",
                         verbose = 0)
  expect_true(all(igraph::E(net)$weight <= quantile(d, 0.05)))
  expect_true(min(d) %in% igraph::E(net)$weight)
})

test_that("x = NULL works; quiet at verbose 0 (1, 6)", {
  skip_if_not_installed("igraph")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  f <- pn_fixture()
  out <- capture.output(net <- gl.plot.network(f$G, verbose = 0))
  expect_s3_class(net, "igraph")
  expect_length(out, 0)
})

test_that("mismatched names and an invalid method error (5, 6)", {
  skip_if_not_installed("igraph")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  f <- pn_fixture()
  x2 <- f$x
  indNames(x2)[1] <- "zzz"
  expect_error(gl.plot.network(f$G, x2, verbose = 0), "individual names")
  expect_error(gl.plot.network(f$G, f$x, method = "xx", verbose = 0),
               "should be one of")
})
