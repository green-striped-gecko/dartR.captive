# Tests of gl.test.offspring

# three candidate fathers, three candidate mothers and four offspring of
# father 1 with mothers 1 and 2
make_trios <- function(nloc = 1000) {
  set.seed(42)
  x <- platypus.gl[1:6, 1:nloc]
  dads <- x[1:3, ]
  mums <- x[4:6, ]
  kids <- gl.sim.offspring(fathers = dads[1, ], mothers = mums[1:2, ],
                           noffpermother = 2, verbose = 0)
  list(x = x, dads = dads, mums = mums, kids = kids)
}

run_test <- function(...) {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  gl.test.offspring(..., plot.display = FALSE, verbose = 0)
}

test_that("candidate mode returns one row per trio with the documented columns", {
  d <- make_trios()
  res <- run_test(p1 = d$dads, p2 = d$mums, off = d$kids)
  expect_named(res, c("results", "mismatch.loci", "plot", "heatmap",
                      "heatmap.info", "allele.freq", "settings", "plot.info"))
  expect_equal(nrow(res$results), 3 * 3 * nInd(d$kids))
  expect_named(res$results, c("test_id", "offspring", "parent1", "parent2",
                              "loci", "observed_mismatches",
                              "expected_mismatches", "prediction_lower",
                              "prediction_upper", "expected_mismatch_rate",
                              "p_value", "LOD", "reject_H0", "delta_LOD"))
  expect_equal(res$settings$mode, "candidate")
  expect_equal(res$settings$n.tests, nrow(res$results))
  expect_s3_class(res$plot, "ggplot")
  expect_s3_class(res$heatmap, "ggplot")
  expect_length(res$allele.freq, nLoc(d$kids))
  expect_identical(names(res$allele.freq), locNames(d$kids))
})

test_that("the true parents are not excluded and have the highest LOD", {
  d <- make_trios()
  res <- run_test(p1 = d$dads, p2 = d$mums, off = d$kids)$results
  true_mum <- setNames(rep(indNames(d$mums)[1:2], each = 2), indNames(d$kids))
  is_true <- res$parent1 == indNames(d$dads)[1] &
    res$parent2 == true_mum[res$offspring]
  expect_equal(sum(is_true), nInd(d$kids))
  expect_false(any(res$reject_H0[is_true]))
  expect_true(all(res$delta_LOD[is_true] == 0))
  expect_true(all(res$LOD[is_true] > 0))
  # every unrelated parent pair is excluded
  unrelated <- res$parent1 != indNames(d$dads)[1] &
    !(res$parent2 %in% indNames(d$mums)[1:2])
  expect_true(all(res$reject_H0[unrelated]))
  # delta_LOD is 0 for the best trio of each offspring and negative otherwise
  expect_true(all(tapply(res$delta_LOD, res$offspring, max) == 0))
  expect_true(all(res$delta_LOD <= 0))
})

test_that("p-values, prediction intervals and expected counts are consistent", {
  d <- make_trios()
  r <- run_test(p1 = d$dads, p2 = d$mums, off = d$kids)$results
  expect_true(all(r$p_value >= 0 & r$p_value <= 1))
  expect_true(all(r$prediction_lower <= r$prediction_upper))
  expect_true(all(r$expected_mismatches >= r$prediction_lower - 1 &
                    r$expected_mismatches <= r$prediction_upper + 1))
  expect_equal(r$expected_mismatch_rate, r$expected_mismatches / r$loci)
  expect_identical(r$reject_H0, r$p_value < 0.05)
  # zero observed mismatches can never be rejected
  expect_true(all(r$p_value[r$observed_mismatches == 0] == 1))
})

test_that("mismatch.loci lists the mismatching SNPs of each trio", {
  d <- make_trios()
  res <- run_test(p1 = d$dads, p2 = d$mums, off = d$kids)
  expect_length(res$mismatch.loci, nrow(res$results))
  expect_identical(names(res$mismatch.loci), paste0("test_", res$results$test_id))
  expect_equal(unname(lengths(res$mismatch.loci)), res$results$observed_mismatches)
  expect_true(all(unlist(res$mismatch.loci) %in% locNames(d$kids)))
  # a listed SNP is Mendelian-incompatible: e.g. both parents 0, offspring 2
  g <- function(obj, ind, loc) as.matrix(obj)[ind, loc]
  r1 <- which(res$results$observed_mismatches > 0)[1]
  loc <- res$mismatch.loci[[r1]][1]
  row <- res$results[r1, ]
  ga <- g(d$dads, row$parent1, loc); gb <- g(d$mums, row$parent2, loc)
  go <- g(d$kids, row$offspring, loc)
  possible <- unique(c(outer(if (ga == 0) 0 else if (ga == 2) 1 else 0:1,
                             if (gb == 0) 0 else if (gb == 2) 1 else 0:1, "+")))
  expect_false(go %in% possible)

  res2 <- run_test(p1 = d$dads, p2 = d$mums, off = d$kids, return.loci = FALSE)
  expect_null(res2$mismatch.loci)
  expect_identical(res2$results, res$results)
})

test_that("swapping p1 and p2 gives the same tests", {
  d <- make_trios()
  a <- run_test(p1 = d$dads, p2 = d$mums, off = d$kids)$results
  b <- run_test(p1 = d$mums, p2 = d$dads, off = d$kids)$results
  key <- function(r) paste(r$offspring, pmin(r$parent1, r$parent2),
                           pmax(r$parent1, r$parent2))
  b <- b[match(key(a), key(b)), ]
  expect_equal(a$LOD, b$LOD)
  expect_equal(a$p_value, b$p_value)
  expect_equal(a$observed_mismatches, b$observed_mismatches)
})

test_that("individuals in both p1 and p2 are tested once per parent pair", {
  d <- make_trios()
  parents <- rbind(d$dads, d$mums)
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_message(
    res <- gl.test.offspring(p1 = parents, p2 = parents, off = d$kids,
                             plot.display = FALSE, verbose = 2),
    "swapped parents")
  r <- res$results
  expect_equal(nrow(r), choose(nInd(parents), 2) * nInd(d$kids))
  expect_false(any(duplicated(paste(r$offspring, pmin(r$parent1, r$parent2),
                                    pmax(r$parent1, r$parent2)))))
  # same individual with different genotypes in p1 and p2
  p2 <- parents
  p2@gen[[1]] <- d$kids@gen[[1]]
  expect_error(run_test(p1 = parents, p2 = p2, off = d$kids),
               "different genotypes")
})

test_that("exploratory mode tests every individual against all pairs of the others", {
  d <- make_trios(500)
  x <- rbind(d$x[c(1, 4, 5)], d$kids[1:2])
  res <- run_test(x = x)
  n <- nInd(x)
  expect_equal(nrow(res$results), n * (n - 1) * (n - 2) / 2)
  expect_equal(res$settings$mode, "exploratory")
  r <- res$results
  expect_true(all(r$offspring != r$parent1 & r$offspring != r$parent2 &
                    r$parent1 != r$parent2))
  expect_error(run_test(x = d$x[1:2]), ">= 3 individuals")
  expect_error(run_test(x = x, max.tests = 10), "Too many trios")
})

test_that("missing genotypes reduce the loci per trio; min.loci gives NA", {
  d <- make_trios(500)
  kids <- d$kids
  m <- as.matrix(kids)
  m[1, 1:100] <- NA
  kids_na <- new("genlight", m, ind.names = indNames(kids),
                 loc.names = locNames(kids), ploidy = 2)
  r <- run_test(p1 = d$dads, p2 = d$mums, off = kids_na)$results
  first <- r$offspring == indNames(kids)[1]
  expect_true(all(r$loci[first] <= 400))
  expect_true(all(r$loci[!first] > 400))

  res <- run_test(p1 = d$dads, p2 = d$mums, off = kids, min.loci = 1000)
  expect_true(all(is.na(res$results$p_value)))
  expect_true(all(is.na(res$results$LOD)))
  expect_null(res$plot)
  expect_null(res$heatmap)
})

test_that("allele.freq can be one value, a vector in locus order or named", {
  d <- make_trios(300)
  args <- list(p1 = d$dads, p2 = d$mums, off = d$kids)
  est <- do.call(run_test, args)
  expect_match(est$settings$allele.freq.source, "estimated")
  af <- est$allele.freq
  pos <- do.call(run_test, c(args, list(allele.freq = unname(af))))
  named <- do.call(run_test, c(args, list(allele.freq = rev(af))))
  expect_equal(pos$results, est$results)
  expect_equal(named$results, est$results)
  expect_equal(pos$settings$allele.freq.source, "user-supplied")
  one <- do.call(run_test, c(args, list(allele.freq = 0.3)))
  expect_true(all(one$allele.freq == 0.3))
  expect_error(do.call(run_test, c(args, list(allele.freq = c(0.1, 0.2)))),
               "expected 1 or")
  expect_error(do.call(run_test, c(args, list(allele.freq = 1.5))),
               "in \\[0,1\\]")
  expect_error(do.call(run_test, c(args, list(allele.freq = NA_real_))),
               "finite numeric")
})

test_that("invalid arguments give clear errors", {
  d <- make_trios(200)
  expect_error(run_test(p1 = d$dads, p2 = d$mums, off = d$kids, epsilon = 0.5),
               "epsilon")
  expect_error(run_test(p1 = d$dads, p2 = d$mums, off = d$kids, alpha = 0),
               "alpha")
  expect_error(run_test(p1 = d$dads, p2 = d$mums, off = d$kids, min.loci = 0),
               "min.loci")
  expect_error(run_test(x = d$x, p1 = d$dads), "not both")
  expect_error(run_test(p1 = d$dads, off = d$kids), "all three")
  expect_error(run_test(p1 = d$dads, p2 = d$mums, off = d$kids[, 1:100]),
               "same number of SNPs")
  expect_error(run_test(p1 = d$dads, p2 = d$mums, off = d$kids, max.tests = 5),
               "Too many candidate trios")
})

test_that("plot.file saves the candidate plot and the heatmap", {
  d <- make_trios(300)
  dir <- file.path(tempdir(), "test_offspring_plots")
  on.exit(unlink(dir, recursive = TRUE))
  run_test(p1 = d$dads, p2 = d$mums, off = d$kids,
           plot.file = "parentage", plot.dir = dir)
  expect_setequal(list.files(dir),
                  c("parentage.png", "parentage_heatmap.png"))
})

test_that(".pb_pmf matches poibin::dpoibin", {
  skip_if_not_installed("poibin")
  set.seed(1)
  pp <- rbind(runif(300, 0, 0.05), c(runif(290, 0, 0.02), runif(10, 0.2, 0.9)),
              rep(0, 300))
  pmf <- dartR.captive:::.pb_pmf(pp)
  for (r in 1:2) {
    ref <- poibin::dpoibin(seq_len(nrow(pmf)) - 1, pp[r, ])
    expect_equal(pmf[, r], ref, tolerance = 1e-10)
  }
  # all-zero probabilities: all mass at 0
  expect_equal(pmf[1, 3], 1)
  expect_equal(sum(pmf[-1, 3]), 0, tolerance = 1e-12)
})
