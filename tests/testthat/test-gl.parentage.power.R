# Tests for gl.parentage.power(). Run inside the package (report()/warn() and the
# utils.* helpers resolve from the dartR.base namespace there). For a standalone
# run, bind those helpers into the environment first (see dev notes).

test_that("gl.parentage.power returns the expected structure", {
  gl <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                    verbose = 0)
  out <- gl.parentage.power(gl, n.offspring = 20, plot.display = FALSE,
                            verbose = 0)
  expect_type(out, "list")
  expect_named(out, c("power", "settings", "plot"))
  expect_s3_class(out$power, "data.frame")
  expect_true(all(c("n.loci", "true.retained", "correct.unique", "false.rate")
                  %in% colnames(out$power)))
  expect_s3_class(out$plot, "ggplot")
})

test_that("power metrics are valid proportions", {
  gl <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                    verbose = 0)
  out <- gl.parentage.power(gl, n.offspring = 20, plot.display = FALSE,
                            verbose = 0)
  expect_true(all(out$power$true.retained >= 0 & out$power$true.retained <= 1))
  expect_true(all(out$power$correct.unique >= 0 & out$power$correct.unique <= 1))
  expect_true(all(out$power$false.rate >= 0 & out$power$false.rate <= 1))
})

test_that("the true parent is always retained under zero error", {
  gl <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                    verbose = 0)
  out <- gl.parentage.power(gl, n.offspring = 30, error.rate = 0,
                            missing.rate = 0, plot.display = FALSE, verbose = 0)
  # With no error and a zero mismatch tolerance the true parent/pair can never
  # be excluded by Mendelian incompatibility.
  expect_true(all(out$power$true.retained == 1))
})

test_that("both single-parent and pair modes run", {
  gl <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                    verbose = 0)
  op <- gl.parentage.power(gl, n.offspring = 15, pairs = TRUE,
                           plot.display = FALSE, verbose = 0)
  os <- gl.parentage.power(gl, n.offspring = 15, pairs = FALSE,
                           plot.display = FALSE, verbose = 0)
  expect_true(op$settings$pairs)
  expect_false(os$settings$pairs)
})

test_that("engine and settings are recorded", {
  gl <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                    verbose = 0)
  out <- gl.parentage.power(gl, n.offspring = 10, plot.display = FALSE,
                            verbose = 0)
  expect_equal(out$settings$engine, "exclusion")
  expect_equal(out$settings$n.offspring, 10)
})

test_that("fewer than two candidates is a fatal error", {
  gl <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                    verbose = 0)
  expect_error(gl.parentage.power(gl[1, ], n.offspring = 5, verbose = 0))
})

# Regression tests for the review of PR #147 (F1-F4, D1) ----------------------

pp_fixture <- function(m) {
  g <- new("genlight", m, ploidy = 2)
  pop(g) <- rep("p", nrow(m))
  g
}

test_that("pair uniqueness checks every alternative pair (F1)", {
  # A and C are identical, so B x C explains an A x B offspring as well as the
  # true pair does. Only A x C offspring are uniquely assigned: 1/3, where
  # sampling alternatives with replacement gave about 1/2.
  cand <- rbind(A = rep(0, 10), B = rep(2, 10), C = rep(0, 10))
  score <- function(off, tp) utils.parentage.score.pair(off, cand, tp, 0)
  ab <- score(rep(1, 10), c(1, 2))
  ac <- score(rep(0, 10), c(1, 3))
  bc <- score(rep(1, 10), c(2, 3))
  expect_equal(c(ab$unique, ac$unique, bc$unique), c(FALSE, TRUE, FALSE))
  expect_equal(c(ab$false.rate, ac$false.rate, bc$false.rate), c(0.5, 0, 0.5))

  colnames(cand) <- paste0("L", 1:10)
  set.seed(1)
  out <- gl.parentage.power(pp_fixture(cand), n.offspring = 600,
                            plot.display = FALSE, verbose = 0)
  expect_equal(out$power$true.retained, 1)
  expect_lt(abs(out$power$correct.unique - 1 / 3), 0.06)
  # every non-unique offspring has exactly one compatible alternative of two
  expect_equal(out$power$false.rate, (1 - out$power$correct.unique) / 2)
})

test_that("the all-pairs mismatch matrix matches the trio test (F1)", {
  set.seed(42)
  cand <- matrix(sample(c(0, 1, 2, NA), 12 * 40, replace = TRUE,
                        prob = c(0.3, 0.3, 0.3, 0.1)), nrow = 12)
  for (k in 1:5) {
    off <- sample(c(0, 1, 2, NA), 40, replace = TRUE,
                  prob = c(0.3, 0.3, 0.3, 0.1))
    mm <- utils.parentage.pair.mismatch(off, cand)
    brute <- outer(1:12, 1:12, Vectorize(function(a, b)
      utils.parentage.trio.mismatch(off, cand[a, ], cand[b, ])))
    expect_equal(mm[upper.tri(mm)], brute[upper.tri(brute)])
    # a pair never has fewer trio mismatches than either member alone
    single <- utils.parentage.single.mismatch(off, cand)
    expect_true(all((brute >= outer(single, single, pmax))[upper.tri(brute)]))
  }
})

test_that("a parent outside the mismatch tolerance is never assigned (F2)", {
  cand <- rbind(A = c(0, 0, 1, 1, 1), B = c(1, 1, 0, 0, 0))
  off <- rep(2, 5)
  expect_equal(unname(utils.parentage.single.mismatch(off, cand)), c(2, 3))
  for (tp in list(c(1, 2), c(2, 1))) {
    sc <- utils.parentage.score.single(off, cand, target = tp[1],
                                      other = tp[2], max.mismatch = 1)
    expect_false(sc$retained)
    expect_false(sc$unique)
  }
  # the sole best candidate, but outside the tolerance
  cand2 <- rbind(A = c(0, 0, 1, 1, 1), C = rep(0, 5), B = c(1, 1, 0, 0, 0))
  sc <- utils.parentage.score.single(off, cand2, target = 1, other = 3,
                                    max.mismatch = 1)
  expect_false(sc$unique)
})

test_that("correct.unique never exceeds true.retained under error (F2)", {
  x <- gl.keep.pop(platypus.gl, pop.list = popNames(platypus.gl)[1],
                   verbose = 0)
  set.seed(7)
  for (pr in c(TRUE, FALSE)) {
    out <- gl.parentage.power(x, n.offspring = 40, pairs = pr,
                              error.rate = 0.05, missing.rate = 0.05,
                              max.mismatch = 2, plot.display = FALSE,
                              verbose = 0)
    expect_true(all(out$power$correct.unique <= out$power$true.retained))
  }
})

test_that("single mode scores the sampled parent with the other unsampled (D1)", {
  # error-free offspring of A x B: both true parents have zero mismatches, so
  # with both among the candidates they would always tie
  cand <- rbind(A = c(0, 0, 1, 2), B = c(0, 1, 0, 2), C = c(2, 1, 1, 1))
  off <- c(0, 0, 0, 2)
  expect_equal(unname(utils.parentage.single.mismatch(off, cand)), c(0, 0, 1))
  sc <- utils.parentage.score.single(off, cand, target = 1, other = 2,
                                    max.mismatch = 0)
  expect_true(sc$unique)
  expect_equal(sc$false.rate, 0)

  # error-free power is no longer stuck at zero
  x <- gl.keep.pop(platypus.gl, pop.list = popNames(platypus.gl)[1],
                   verbose = 0)
  set.seed(3)
  out <- gl.parentage.power(x, n.offspring = 40, pairs = FALSE,
                            plot.display = FALSE, verbose = 0)
  expect_true(all(out$power$true.retained == 1))
  expect_gt(max(out$power$correct.unique), 0.5)
})

test_that("COLONY runs once per replicate and averages them (F3)", {
  x <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                   verbose = 0)[1:2, ]
  cand <- indNames(x)
  calls <- list()
  # with two candidates every true pair is (cand[1], cand[2])
  assign <- list(cand, c("#1", "#2"), c(cand[1], "#1"))
  local_mocked_bindings(gl.run.colony = function(x, outpath, ...) {
    calls[[length(calls) + 1]] <<- list(
      outpath = outpath,
      off = as.matrix(x)[grepl("^off", indNames(x)), , drop = FALSE])
    k <- length(calls)
    offs <- indNames(x)[grepl("^off", indNames(x))]
    list(best.config = data.frame(OffspringID = offs,
                                  FatherID = assign[[k]][1],
                                  MotherID = assign[[k]][2]))
  })
  set.seed(11)
  out <- gl.parentage.power(x, n.offspring = 6, engine = "colony", n.rep = 3,
                            plot.display = FALSE, verbose = 0)
  expect_length(calls, 3)
  expect_equal(basename(vapply(calls, `[[`, "", "outpath")),
               c("rep1", "rep2", "rep3"))
  expect_false(identical(calls[[1]]$off, calls[[2]]$off))
  expect_equal(out$settings$n.rep, 3)
  expect_equal(nrow(out$power), 1)
  expect_equal(out$power$true.retained, 2 / 3)
  expect_equal(out$power$correct.unique, 1 / 3)
  expect_equal(out$power$false.rate, 0)

  calls <- list()
  gl.parentage.power(x, n.offspring = 6, engine = "colony", n.rep = 1,
                     plot.display = FALSE, verbose = 0)
  expect_length(calls, 1)
})

test_that("COLONY offspring matrix keeps its shape with one locus (F4)", {
  x <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
                   verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  seen <- NULL
  local_mocked_bindings(gl.run.colony = function(x, ...) {
    seen <<- x
    offs <- indNames(x)[grepl("^off", indNames(x))]
    list(best.config = data.frame(OffspringID = offs, FatherID = "#1",
                                  MotherID = "#2"))
  })
  for (case in list(c(loci = 1, n.off = 2), c(loci = 5, n.off = 1),
                    c(loci = 1, n.off = 1))) {
    xs <- x[, seq_len(case[["loci"]])]
    gl.parentage.power(xs, n.offspring = case[["n.off"]], engine = "colony",
                       plot.display = FALSE, verbose = 0)
    expect_equal(nLoc(seen), case[["loci"]])
    expect_equal(nInd(seen), nInd(x) + case[["n.off"]])
    expect_equal(tail(indNames(seen), case[["n.off"]]),
                 paste0("off", seq_len(case[["n.off"]])))
  }
})
