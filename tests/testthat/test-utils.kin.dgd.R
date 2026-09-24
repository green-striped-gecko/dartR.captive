# Characterization test for utils.kin.dgd(). Captured on commit 586d93e,
# then updated only for the changes approved in
# function-review/reports/dartR.captive/utils.kin.md.
# A failing expectation here means behaviour changed; it does not mean the
# old behaviour was correct.

kin_toy <- function() {
  set.seed(3)
  n <- 6
  A <- matrix(stats::runif(n * n, 0, 0.1), n)
  A <- (A + t(A)) / 2
  diag(A) <- 0.5 + stats::runif(n, 0, 0.05)
  dimnames(A) <- list(letters[1:n], letters[1:n])
  A
}

# Independent construction: append each virtual offspring explicitly
brute_gd <- function(A, pairs) {
  K <- A
  for (r in seq_len(nrow(pairs))) {
    a <- pairs[r, 1]
    b <- pairs[r, 2]
    v <- (K[a, ] + K[b, ]) / 2
    ids <- c(rownames(K), paste0("offspring_", r))
    K <- rbind(cbind(K, v), c(v, 0.5 * (1 + K[a, b])))
    dimnames(K) <- list(ids, ids)
  }
  1 - mean(K)
}

test_that("GD of the unmodified matrix", {
  A <- kin_toy()
  expect_equal(utils.kin.dgd(A), 1 - mean(A))
})

test_that("removals", {
  A <- kin_toy()
  expect_equal(utils.kin.dgd(A, drop = c("a", "f")),
               1 - mean(A[-c(1, 6), -c(1, 6)]))
  expect_error(utils.kin.dgd(A, drop = "z"), "not present")
  expect_error(utils.kin.dgd(A, drop = letters[1:6]), "empty")
})

test_that("virtual offspring, chained and selfed", {
  A <- kin_toy()
  pairs <- rbind(c("a", "b"), c("c", "d"),
                 c("offspring_1", "offspring_2"), c("e", "e"))
  expect_equal(utils.kin.dgd(A, add.pairs = pairs), brute_gd(A, pairs))
  expect_equal(utils.kin.dgd(A, add.pairs = data.frame(p1 = "a", p2 = "b")),
               brute_gd(A, rbind(c("a", "b"))))
  # a dropped individual cannot be a parent
  expect_error(utils.kin.dgd(A, drop = "a", add.pairs = rbind(c("a", "b"))),
               "not found")
})

test_that("missing kinship", {
  A <- kin_toy()
  A[1, 2] <- A[2, 1] <- NA
  expect_true(is.na(utils.kin.dgd(A)))
  expect_equal(utils.kin.dgd(A, na.rm = TRUE), 1 - mean(A, na.rm = TRUE))
})

test_that("zero-row add.pairs", {
  A <- kin_toy()
  # change 1: an empty table of pairs previously errored; it adds nothing
  expect_equal(utils.kin.dgd(A, add.pairs = matrix(character(0), ncol = 2)),
               1 - mean(A))
})
