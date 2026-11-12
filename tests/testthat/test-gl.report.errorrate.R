# Tests of gl.report.errorrate

# n.pairs replicate pairs at nloc SNPs; each allele copy of each sample is
# miscalled with probability eps (the model of the function). dropout: a
# true heterozygote is called homozygous with this probability instead.
sim_reps <- function(n.pairs = 20, nloc = 2000, eps = 0.02, dropout = 0,
                     seed = 1) {
  set.seed(seed)
  q <- runif(nloc, 0.05, 0.95)
  call <- function(a1, a2) {
    if (eps > 0) {
      a1 <- ifelse(runif(length(a1)) < eps, 1 - a1, a1)
      a2 <- ifelse(runif(length(a2)) < eps, 1 - a2, a2)
    }
    g <- a1 + a2
    if (dropout > 0) {
      d <- g == 1 & runif(length(g)) < dropout
      g[d] <- sample(c(0, 2), sum(d), replace = TRUE)
    }
    g
  }
  m <- matrix(NA_real_, 2 * n.pairs, nloc)
  for (p in seq_len(n.pairs)) {
    a1 <- rbinom(nloc, 1, q)
    a2 <- rbinom(nloc, 1, q)
    m[2 * p - 1, ] <- call(a1, a2)
    m[2 * p, ] <- call(a1, a2)
  }
  ids <- paste0("S", rep(seq_len(n.pairs), each = 2), "_rep", 1:2)
  gl <- new("genlight", m, ind.names = ids,
            loc.names = paste0("snp", seq_len(nloc)), ploidy = 2)
  list(gl = gl, reps = sub("_rep[12]$", "", ids))
}

run_er <- function(...) {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  gl.report.errorrate(..., plot.display = FALSE, verbose = 0)
}

test_that("the pooled estimate recovers the simulated error rate", {
  s <- sim_reps(eps = 0.02)
  res <- run_er(s$gl, replicates = s$reps)
  po <- res$pooled
  expect_equal(po$epsilon, 0.02, tolerance = 0.1)        # within 10 %
  expect_true(po$lower < 0.02 && po$upper > 0.02)
  expect_true(po$lower < po$epsilon && po$epsilon < po$upper)
  expect_equal(po$n_pairs, 20)
  expect_equal(po$n_comparisons, 20 * 2000)
  # per pair: estimates scatter around the truth, mostly covered by the CI
  pp <- res$per.pair
  expect_true(mean(pp$lower <= 0.02 & pp$upper >= 0.02) >= 0.8)
  expect_equal(po$mean_pairs, mean(pp$epsilon))
  expect_equal(po$median_pairs, median(pp$epsilon))
  # the model is the one simulated, so it fits
  expect_gt(res$fit$p_value, 0.001)
  # true genotype frequencies, estimated jointly
  expect_equal(sum(res$genotype.freq), 1)
  expect_named(res$genotype.freq, c("0", "1", "2"))
})

test_that("identical replicates give an error rate of 0", {
  s <- sim_reps(n.pairs = 5, nloc = 500, eps = 0)
  res <- run_er(s$gl, replicates = s$reps)
  expect_true(all(res$per.pair$discordant == 0))
  expect_true(all(res$per.pair$epsilon == 0))
  expect_true(all(res$per.pair$lower == 0))
  expect_true(all(res$per.pair$upper > 0))
  expect_equal(res$pooled$epsilon, 0)
  expect_equal(res$pooled$lower, 0)
  # a smaller alpha gives a wider interval
  wide <- run_er(s$gl, replicates = s$reps, alpha = 0.01)
  expect_gt(wide$pooled$upper, res$pooled$upper)
})

test_that("per-pair counts and the observed and expected tables add up", {
  s <- sim_reps(n.pairs = 6, nloc = 800, eps = 0.03)
  m <- as.matrix(s$gl)
  m[1, 1:50] <- NA            # missing in one sample only
  m[2, 41:60] <- NA           # 10 loci missing in both
  gl <- new("genlight", m, ind.names = indNames(s$gl),
            loc.names = locNames(s$gl), ploidy = 2)
  res <- run_er(gl, replicates = s$reps)
  pp <- res$per.pair
  expect_named(pp, c("ind1", "ind2", "group", "loci", "discordant",
                     "discordance_rate", "missing_discordance", "epsilon",
                     "lower", "upper"))
  both <- !is.na(m[1, ]) & !is.na(m[2, ])
  expect_equal(pp$loci[1], sum(both))
  expect_equal(pp$discordant[1], sum(m[1, both] != m[2, both]))
  expect_equal(pp$discordance_rate, pp$discordant / pp$loci)
  # loci 1-40 and 51-60 called in exactly one sample, of 790 called in at
  # least one
  expect_equal(pp$missing_discordance[1], 50 / 790)
  expect_true(all(pp$missing_discordance[-1] == 0))
  # pooled tables: symmetric, one count per SNP comparison
  expect_equal(unname(res$observed), unname(t(res$observed)))
  expect_equal(sum(res$observed), sum(pp$loci))
  expect_equal(sum(res$expected), sum(res$observed))
  expect_equal(sum(diag(res$observed)), sum(pp$loci - pp$discordant))
  expect_identical(dimnames(res$observed),
                   list(sample1 = c("0", "1", "2"), sample2 = c("0", "1", "2")))
  expect_equal(res$fit$df, 2)
})

test_that("pairs can be given by name or by index, replicates by group", {
  s <- sim_reps(n.pairs = 4, nloc = 600, eps = 0.02)
  ids <- indNames(s$gl)
  by_grp <- run_er(s$gl, replicates = s$reps)
  by_idx <- run_er(s$gl, pairs = data.frame(c(1, 3, 5, 7), c(2, 4, 6, 8)))
  by_nam <- run_er(s$gl, pairs = data.frame(ids[c(1, 3, 5, 7)],
                                            ids[c(2, 4, 6, 8)]))
  by_fac <- run_er(s$gl, pairs = data.frame(factor(ids[c(1, 3, 5, 7)]),
                                            factor(ids[c(2, 4, 6, 8)])))
  expect_equal(by_idx$per.pair, by_nam$per.pair)
  expect_equal(by_fac$per.pair, by_nam$per.pair)
  expect_equal(by_idx$pooled, by_grp$pooled)
  expect_equal(by_grp$per.pair$group, paste0("S", 1:4))
  expect_true(all(is.na(by_idx$per.pair$group)))
  expect_equal(by_grp$per.pair$ind1, ids[c(1, 3, 5, 7)])

  # a group of three gives three pairs; NA and single members are ignored
  reps <- s$reps
  reps[3:4] <- "S1"           # S1 now has four members: 6 pairs
  reps[5] <- NA
  res <- run_er(s$gl, replicates = reps)
  expect_equal(sum(res$per.pair$group == "S1"), choose(4, 2))
  expect_false("S3" %in% res$per.pair$group)
  expect_equal(nrow(res$per.pair), choose(4, 2) + 1)
})

test_that("allelic dropout makes the symmetric model fit poorly", {
  s <- sim_reps(n.pairs = 20, nloc = 2000, eps = 0.005, dropout = 0.1)
  res <- run_er(s$gl, replicates = s$reps)
  expect_lt(res$fit$p_value, 0.001)
  # dropout is absorbed into a larger epsilon ...
  expect_gt(res$pooled$epsilon, 0.01)
  # ... and shows as excess 0 vs 2 discordance (a heterozygote that drops to
  # a different homozygote in each replicate)
  expect_gt(res$observed["0", "2"], 2 * res$expected["0", "2"])
})

test_that("pairs with too few loci are dropped, with a warning", {
  s <- sim_reps(n.pairs = 3, nloc = 300, eps = 0.02)
  m <- as.matrix(s$gl)
  m[1, 1:250] <- NA
  gl <- new("genlight", m, ind.names = indNames(s$gl),
            loc.names = locNames(s$gl), ploidy = 2)
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_output(
    res <- gl.report.errorrate(gl, replicates = s$reps, plot.display = FALSE,
                               verbose = 1),
    "1 pairs have fewer than 100")
  expect_true(is.na(res$per.pair$epsilon[1]))
  expect_equal(res$pooled$n_pairs, 2)
  expect_error(run_er(gl, replicates = s$reps, min.loci = 1000),
               "No pair has")
})

test_that("invalid arguments give clear errors", {
  s <- sim_reps(n.pairs = 2, nloc = 200)
  gl <- s$gl
  expect_error(run_er(gl), "either pairs or replicates")
  expect_error(run_er(gl, pairs = data.frame(1, 2), replicates = s$reps),
               "either pairs or replicates")
  expect_error(run_er(gl, pairs = data.frame(1, 1)), "same individual twice")
  expect_error(run_er(gl, pairs = data.frame(1, 9)), "between 1 and 4")
  expect_error(run_er(gl, pairs = data.frame("S1_rep1", "nobody")),
               "not found in x: nobody")
  expect_error(run_er(gl, pairs = data.frame(1:2)), "two columns")
  expect_error(run_er(gl, replicates = s$reps[-1]), "one entry per individual")
  expect_error(run_er(gl, replicates = indNames(gl)), "two or more members")
  expect_error(run_er(gl, replicates = s$reps, alpha = 1), "alpha")
  expect_error(run_er(gl, replicates = s$reps, min.loci = 0), "min.loci")
  expect_error(run_er(gl, replicates = s$reps, plot.colors = "red"),
               "two colours")
})

test_that("the plot is returned and saved; the function works with do.call", {
  s <- sim_reps(n.pairs = 3, nloc = 300)
  dir <- file.path(tempdir(), "test_errorrate_plots")
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE))
  res <- run_er(s$gl, replicates = s$reps, plot.file = "errorrate",
                plot.dir = dir)
  expect_s3_class(res$plot, "ggplot")
  expect_true(any(grepl("errorrate", list.files(dir))))
  pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_output(
    res2 <- do.call(gl.report.errorrate,
                    list(s$gl, replicates = s$reps, plot.display = FALSE,
                         verbose = 1)),
    "Completed: gl.report.errorrate")
  expect_equal(res2$pooled, res$pooled)
})

test_that("the example runs on bandicoot.gl", {
  res <- run_er(bandicoot.gl, pairs = data.frame(1, 2))
  expect_equal(nrow(res$per.pair), 1)
  expect_true(res$pooled$epsilon >= 0 && res$pooled$epsilon < 0.5)
})
