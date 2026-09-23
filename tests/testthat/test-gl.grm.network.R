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

# ---- Second-pass review (2026-09-23): approved changes 1-7 ----

test_that("r2 change 1: no pair above threshold plots nodes only, no error", {
  sub <- testset.gl[1:5, 1:50]
  nm <- indNames(sub)
  G <- diag(1, 5)
  dimnames(G) <- list(nm, nm)
  res <- gl.grm.network(G, sub, verbose = 0)
  expect_s3_class(res$plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(res$plot))
})

test_that("r2 change 2: categorise colours follow the documented order", {
  sub <- testset.gl[1:5, 1:50]
  nm <- indNames(sub)
  G <- diag(1, 5)
  dimnames(G) <- list(nm, nm)
  G[nm[1], nm[2]] <- G[nm[2], nm[1]] <- 0.7  # kinship 0.35
  G[nm[3], nm[4]] <- G[nm[4], nm[3]] <- 0.5  # kinship 0.25
  G[nm[5], nm[1]] <- G[nm[1], nm[5]] <- 0.3  # kinship 0.15
  cols_of <- function(G) {
    b <- ggplot2::ggplot_build(
      gl.grm.network(G, sub, categorise = TRUE, verbose = 0)$plot
    )
    setNames(b$data[[1]]$colour, b$plot$layers[[1]]$data$cat)
  }
  cols <- cols_of(G)
  expect_equal(unname(cols["Same Individual"]), "#E63E94")
  expect_equal(unname(cols["Full Siblings\nParent-Offspring"]), "#E5D44C")
  expect_equal(unname(cols["Half Siblings"]), "#3ED2E6")
  # dropping the full-sib pair must not shift the other colours
  G[nm[3], nm[4]] <- G[nm[4], nm[3]] <- 0
  cols2 <- cols_of(G)
  expect_equal(unname(cols2["Same Individual"]), "#E63E94")
  expect_equal(unname(cols2["Half Siblings"]), "#3ED2E6")
})

test_that("r2 change 3: default path returns kinship = G / 2", {
  sub <- testset.gl[1:3, 1:50]
  nm <- sort(indNames(sub))
  G <- diag(1.2, 3)
  dimnames(G) <- list(nm, nm)
  G[nm[2], nm[1]] <- G[nm[1], nm[2]] <- 0.4
  m <- gl.grm.network(G, sub, verbose = 0)$kinship
  expect_equal(m[nm[2], nm[1]], 0.2)
  # a relatedness of 0.2 (kinship 0.1) is below the default threshold
  G[nm[2], nm[1]] <- G[nm[1], nm[2]] <- 0.2
  b <- ggplot2::ggplot_build(gl.grm.network(G, sub, verbose = 0)$plot)
  expect_equal(nrow(b$data[[1]]), 0)
})

test_that("r2 change 4: G that does not match x errors clearly", {
  sub <- testset.gl[1:5, 1:50]
  nm <- indNames(sub)
  G <- diag(1, 4)
  dimnames(G) <- list(nm[1:4], nm[1:4])
  expect_error(gl.grm.network(G, sub, verbose = 0), "individual names of x")
  G <- diag(1, 5)
  expect_error(gl.grm.network(G, sub, verbose = 0), "individual names of x")
})

test_that("r2 change 5: return is a named list; matrix lower-triangular", {
  sub <- testset.gl[1:3, 1:50]
  nm <- indNames(sub)
  G <- diag(1.2, 3)
  dimnames(G) <- list(nm, nm)
  G[nm[1], nm[2]] <- G[nm[2], nm[1]] <- 0.4
  res <- gl.grm.network(G, sub, verbose = 0)
  expect_named(res, c("plot", "kinship"))
  expect_true(all(is.na(res$kinship[upper.tri(res$kinship)])))
  expect_equal(unname(diag(res$kinship)), c(0, 0, 0))
})

test_that("r2 change 6: invalid method errors", {
  sub <- testset.gl[1:3, 1:50]
  nm <- indNames(sub)
  G <- diag(1, 3)
  dimnames(G) <- list(nm, nm)
  expect_error(gl.grm.network(G, sub, method = "xx", verbose = 0))
})

test_that("G tagged as kinship (gl.run.EMIBD9 $rel) is not halved", {
  sub <- testset.gl[1:3, 1:50]
  nm <- sort(indNames(sub))
  G <- diag(0.5, 3)
  dimnames(G) <- list(nm, nm)
  G[nm[2], nm[1]] <- G[nm[1], nm[2]] <- 0.25  # full-sib kinship
  attr(G, "scale") <- "kinship"
  m <- gl.grm.network(G, sub, verbose = 0)$kinship
  expect_equal(m[nm[2], nm[1]], 0.25)
  # untagged, the same matrix is read as relatedness and halved
  attr(G, "scale") <- NULL
  m <- gl.grm.network(G, sub, verbose = 0)$kinship
  expect_equal(m[nm[2], nm[1]], 0.125)
  # standardise (Goudet et al. 2018): pairs 0.25, 0, 0, mean 1/12, so
  # (0.25 - 1/12) / (1 - 1/12) = 2/11; tagged and untagged scale alike
  attr(G, "scale") <- "kinship"
  m <- gl.grm.network(G, sub, standardise = TRUE, verbose = 0)$kinship
  expect_equal(m[nm[2], nm[1]], 2 / 11)
  expect_equal(m[nm[3], nm[1]], (0 - 1 / 12) / (1 - 1 / 12))
  attr(G, "scale") <- NULL
  m <- gl.grm.network(2 * G, sub, standardise = TRUE, verbose = 0)$kinship
  expect_equal(m[nm[2], nm[1]], 2 / 11)
})

test_that("standardise = TRUE recovers pedigree kinship on filtered data", {
  skip_if_not_installed("rrBLUP")
  skip_if_not_installed("igraph")
  skip_if_not(exists("testset2.gl"), "testset2.gl needs dartR.data >= 1.2.5")
  # before: kinship - mean inbreeding (0.158) gave parent-offspring 0.097
  x <- gl.filter.callrate(testset2.gl, threshold = 0.95, verbose = 0)
  im <- x@other$ind.metrics
  ok <- im$sire %in% indNames(x)
  po <- cbind(as.character(im$id[ok]), as.character(im$sire[ok]))
  G <- gl.grm(x, plotheatmap = FALSE, verbose = 0)
  pdf(NULL)
  on.exit(grDevices::dev.off())
  m <- gl.grm.network(G, x, standardise = TRUE, verbose = 0)$kinship
  # kinship holds the lower triangle only
  v <- pmax(m[po], m[po[, 2:1]], na.rm = TRUE)
  expect_equal(mean(v), 0.25, tolerance = 0.03)
  expect_true(all(v > 0.1875))
})
