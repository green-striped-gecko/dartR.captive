# Characterization test for utils.kin.dominant(). Captured on commit
# 586d93e, then updated only for the changes approved in
# function-review/reports/dartR.captive/utils.kin.md. A failing expectation
# here means behaviour changed; it does not mean the old behaviour was
# correct.

test_that("equals a pairwise loop over shared informative loci", {
  gs <- testset.gs[1:40, ]
  m <- as.matrix(gs)
  mu <- colMeans(m, na.rm = TRUE)
  p <- !is.na(mu) & mu > 0 & mu < 1
  m <- m[, p]
  mu <- mu[p]
  naive <- matrix(NA_real_, nrow(m), nrow(m))
  for (i in seq_len(nrow(m))) {
    for (j in seq_len(nrow(m))) {
      s <- !is.na(m[i, ]) & !is.na(m[j, ])
      if (any(s)) {
        naive[i, j] <- mean((m[i, s] - mu[s]) * (m[j, s] - mu[s]) /
                              (mu[s] * (1 - mu[s])))
      }
    }
  }
  naive <- naive / 2
  diag(naive) <- 0.5
  k <- utils.kin.dominant(gs, verbose = 0)
  expect_equal(unname(k), naive)
  expect_identical(rownames(k), indNames(gs))
  expect_identical(colnames(k), indNames(gs))
})

test_that("rejects SNP data", {
  expect_error(utils.kin.dominant(testset.gl[1:5, ], verbose = 0),
               "SilicoDArT")
})

test_that("flag messages inside gl.kin", {
  # change 2: the internal helper previously printed its own
  # Starting/Completed lines inside gl.kin's at verbose >= 1
  out <- capture.output(invisible(gl.kin(testset.gs[1:15, ], verbose = 1)))
  expect_false(any(grepl("utils.kin.dominant", out)))
  expect_true(any(grepl("Completed: gl.kin", out)))
})
