test_that("gl.grm.network returns a plot and a labelled matrix (baseline)", {
  sub <- testset.gl[1:4, 1:50]
  nm <- indNames(sub)
  G <- diag(1, 4)
  dimnames(G) <- list(nm, nm)
  G[nm[1], nm[2]] <- G[nm[2], nm[1]] <- 0.6  # kinship 0.3
  G[nm[3], nm[4]] <- G[nm[4], nm[3]] <- 0.4  # kinship 0.2

  res <- gl.grm.network(G, sub, verbose = 0)

  expect_type(res, "list")
  expect_s3_class(res[[1]], "ggplot")
  expect_equal(dim(res[[2]]), c(4, 4))
  expect_equal(rownames(res[[2]]), nm)
  # self-loop rows are hard-coded to kinship 0 (current behaviour)
  expect_equal(unname(diag(res[[2]])), c(0, 0, 0, 0))
})

test_that("gl.grm.network F1 fix: print() no longer crashes with a single qualifying edge", {
  # Approved fix: the manual `breaks`/`labels` computation fed to
  # scale_colour_gradientn() is removed; the scale now derives its own
  # breaks, which handles a single-point data range correctly.
  sub <- testset.gl[1:3, 1:50]
  nm <- indNames(sub)
  G <- diag(1, 3)
  dimnames(G) <- list(nm, nm)
  G[nm[1], nm[2]] <- G[nm[2], nm[1]] <- 0.26  # only pair above default threshold

  res <- gl.grm.network(G, sub, verbose = 0)
  expect_s3_class(res[[1]], "ggplot")
  expect_silent(ggplot2::ggplot_build(res[[1]]))
})

test_that("gl.grm.network F2 fix: categorise=TRUE no longer crashes with a low threshold", {
  # Approved fix (option 2a): the undocumented 'First Cousins' bucket
  # (0.038-0.1) is dropped, so categorisation always yields at most the 3
  # documented buckets, matching the 3-color `color.categories` default.
  # The formerly-crashing pair (kinship 0.06) now gets no category (NA).
  sub <- testset.gl[1:5, 1:50]
  nm <- indNames(sub)
  G <- diag(1, 5)
  dimnames(G) <- list(nm, nm)
  G[nm[1], nm[2]] <- G[nm[2], nm[1]] <- 0.35  # "Same Individual"
  G[nm[1], nm[3]] <- G[nm[3], nm[1]] <- 0.25  # "Full Siblings..."
  G[nm[1], nm[4]] <- G[nm[4], nm[1]] <- 0.15  # "Half Siblings"
  G[nm[1], nm[5]] <- G[nm[5], nm[1]] <- 0.06  # no longer categorised

  res <- gl.grm.network(G, sub, categorise = TRUE, kinship.threshold = 0.05,
                        verbose = 0)
  expect_s3_class(res[[1]], "ggplot")
  expect_silent(ggplot2::ggplot_build(res[[1]]))
})

test_that("gl.grm.network F3 fix: an individual in 2+ above-threshold pairs is no longer duplicated", {
  # Approved fix: links_plot is deduplicated to one row per individual
  # (keeping the strongest relationship) before the plotcord merge.
  sub <- testset.gl[1:4, 1:50]
  nm <- indNames(sub)
  G <- diag(1, 4)
  dimnames(G) <- list(nm, nm)
  G[nm[1], nm[2]] <- G[nm[2], nm[1]] <- 0.6  # edge 1-2
  G[nm[1], nm[3]] <- G[nm[3], nm[1]] <- 0.4  # edge 1-3: individual 1 now in 2 pairs

  res <- gl.grm.network(G, sub, kinship.threshold = 0.125, verbose = 0)
  built <- ggplot2::ggplot_build(res[[1]])
  point_rows <- nrow(built$data[[2]])

  expect_equal(point_rows, nInd(sub))
})
