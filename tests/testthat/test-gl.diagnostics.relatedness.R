# Tests of the reviewed gl.diagnostics.relatedness and its helpers
# (function-review/reports/dartR.captive/gl.diagnostics.relatedness.md).
# Tests that run the relatedness estimators need the non-CRAN engine
# 'dartR.coancestry' (used by gl.relatedness).

hand_pedigree <- function() {
  # F1-F4 founders; O1, O2 full sibs; O3 unrelated to them; O4 = O1 x O2
  data.frame(id = c("F1", "F2", "F3", "F4", "O1", "O2", "O3", "O4"),
             dad = c(NA, NA, NA, NA, "F1", "F1", "F3", "O1"),
             mom = c(NA, NA, NA, NA, "F2", "F2", "F4", "O2"))
}

test_that("missing parents are not shared parents (1)", {
  cl <- as.data.frame(CleanupExtractParents(hand_pedigree()[1:7, ]))
  expect_equal(sum(cl$relationship == "half_sibs"), 0)
  expect_equal(sum(cl$relationship == "half_first_cousins"), 0)
  expect_equal(sum(cl$relationship == "full_sibs"), 1)
})

test_that("pedigree kinship matches kinship2, including inbreeding (2)", {
  ped <- hand_pedigree()
  K <- pedigreeKinship(ped)
  expect_equal(K["O1", "O2"], 0.25)
  expect_equal(K["O1", "O3"], 0)
  expect_equal(K["O4", "O4"], 0.625)
  skip_if_not_installed("kinship2")
  kp <- kinship2::kinship(id = ped$id, dadid = ped$dad, momid = ped$mom)
  expect_equal(unname(K[ped$id, ped$id]), unname(kp[ped$id, ped$id]))
})

test_that("RMSE is the root mean square error against rel (4)", {
  df <- data.frame(RelDegree = "full_sibs", rel = c(0.25, 0.25),
                   wang = c(0.2, 0.35))
  out <- calcRMSE(list(df), "wang")[[1]]
  expect_equal(out["wang", "full_sibs"], sqrt(mean(c(-0.05, 0.1)^2)))
  expect_true(is.na(out["wang", "half_sibs"]))
})

test_that("all pairs are kept (3)", {
  skip_if_not_installed("dartR.coancestry")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- gl.filter.allna(testset.gl[1:15, ], verbose = 0)
  capture.output(
    res <- gl.diagnostics.relatedness(x, which_tests = "wang", verbose = 0)
  )
  expect_equal(nrow(res@MergedDf[[1]]) / 2, choose(15, 2))
})

test_that("attached pedigree: one row per pair, classes and rel (1, 2, 3)", {
  skip_if_not_installed("dartR.coancestry")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- gl.filter.allna(testset.gl[1:15, ], verbose = 0)
  ids <- indNames(x)
  x@other$ind.metrics$id <- ids
  x@other$ind.metrics$dad <- c(0, 0, 0, 0, ids[1], ids[1], ids[3], rep(0, 8))
  x@other$ind.metrics$mom <- c(0, 0, 0, 0, ids[2], ids[2], ids[4], rep(0, 8))
  capture.output(
    res <- gl.diagnostics.relatedness(x, which_tests = "wang",
                                      includedPed = TRUE, rmseOut = TRUE,
                                      verbose = 0)
  )
  m <- res@MergedDf[[1]]
  expect_equal(colnames(m), c("ind1", "ind2", "RelDegree", "rel", "wang",
                              "rrBLUP"))
  expect_equal(nrow(m), choose(15, 2))
  expect_equal(as.vector(table(m$RelDegree)[c("full_sibs",
                                               "parent_offspring",
                                               "unrelated")]),
               c(1, 6, 98))
  expect_true(all(m$rel[m$RelDegree == "unrelated"] == 0))
})

test_that("simulation runs with the default variable files (5)", {
  skip_if_not_installed("dartR.coancestry")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- gl.filter.allna(testset.gl[1:30, ], verbose = 0)
  set.seed(1)
  capture.output(
    res <- gl.diagnostics.relatedness(x, which_tests = "wang",
                                      run_sim = TRUE, rmseOut = TRUE, Ne = 50,
                                      verbose = 0)
  )
  m <- res@MergedDf[[1]]
  # analysisUnit = "generation": all pairs within each stored generation
  gen <- table(res@SimOutput[[1]]@other$ind.metrics$generation)
  expect_equal(nrow(m), sum(choose(gen, 2)))
  expect_false(anyDuplicated(m[, c("ind1", "ind2")]) > 0)
  expect_false(anyNA(m$generation))
  expect_true("unrelated" %in% m$RelDegree)
})

test_that("analysisUnit = 'pooled' estimates all pairs across generations", {
  skip_if_not_installed("dartR.coancestry")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- gl.filter.allna(testset.gl[1:20, 1:200], verbose = 0)
  set.seed(1)
  capture.output(
    res <- gl.diagnostics.relatedness(x, which_tests = "wang",
                                      run_sim = TRUE, Ne = 30,
                                      analysisUnit = "pooled", verbose = 0)
  )
  m <- res@MergedDf[[1]]
  expect_equal(nrow(m), choose(nInd(res@SimOutput[[1]]), 2))
  expect_true("parent_offspring" %in% m$RelDegree)
  expect_error(gl.diagnostics.relatedness(x, analysisUnit = "bad",
                                          verbose = 0))
})

test_that("SilicoDArT input errors (6)", {
  skip_if_not_installed("dartR.coancestry")
  expect_error(gl.diagnostics.relatedness(testset.gs[1:10, 1:50],
                                          verbose = 0),
               "Only SNP data")
})

test_that("estimates are gl.relatedness halved to kinship", {
  skip_if_not_installed("dartR.coancestry")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- gl.filter.allna(testset.gl[1:15, ], verbose = 0)
  res <- gl.diagnostics.relatedness(x, which_tests = c("wang", "loiselle"),
                                    verbose = 0)
  m <- tidyr::pivot_wider(res@MergedDf[[1]], names_from = "variable",
                          values_from = "value")
  r <- gl.relatedness(x, estimators = c("wang", "loiselle"),
                      plot.out = FALSE, verbose = 0)
  ij <- cbind(as.character(m$ind1), as.character(m$ind2))
  expect_equal(m$wang, r$wang[ij] / 2)
  expect_equal(m$loiselle, r$loiselle[ij] / 2)
})

test_that("unknown estimators in which_tests error", {
  skip_if_not_installed("dartR.coancestry")
  expect_error(
    gl.diagnostics.relatedness(testset.gl[1:10, ], which_tests = "foo",
                               verbose = 0),
    "Unknown estimator"
  )
})

test_that("relatives across generations get their own class, not unrelated", {
  # G1 x G2 -> P1, P2 (full sibs); P1 x M1 -> C1; P2 x M2 -> C2;
  # C1 x M3 -> D1; H1 = G1 x M4 (half sib of P1, P2)
  ped <- data.frame(
    id  = c("G1", "G2", "M1", "M2", "M3", "M4", "P1", "P2", "H1", "C1",
            "C2", "D1"),
    dad = c(NA, NA, NA, NA, NA, NA, "G1", "G1", "G1", "P1", "P2", "C1"),
    mom = c(NA, NA, NA, NA, NA, NA, "G2", "G2", "M4", "M1", "M2", "M3"))
  ids <- ped$id
  pairs <- t(combn(ids, 2))
  est <- data.frame(ind1 = pairs[, 1], ind2 = pairs[, 2], variable = "wang",
                    value = 0)
  m <- mergePedigreeTruth(est, ped)
  cls <- function(a, b) m$RelDegree[m$ind1 == min(a, b) & m$ind2 == max(a, b)]
  kin <- function(a, b) m$rel[m$ind1 == min(a, b) & m$ind2 == max(a, b)]
  expect_equal(cls("G1", "C1"), "grandparent_grandchild")
  expect_equal(kin("G1", "C1"), 0.125)
  expect_equal(cls("P2", "C1"), "avuncular")
  expect_equal(kin("P2", "C1"), 0.125)
  expect_equal(cls("H1", "C1"), "half_avuncular")
  expect_equal(cls("G2", "D1"), "great_grandparent_grandchild")
  expect_equal(cls("C1", "C2"), "full_first_cousins")
  expect_equal(cls("P2", "D1"), "other_relatives")   # great-avuncular
  expect_equal(cls("M1", "M2"), "unrelated")
  expect_true(all(m$rel[m$RelDegree == "unrelated"] == 0))
  expect_true(all(m$rel[m$RelDegree == "other_relatives"] > 0))
})

test_that("default variable files mirror x", {
  skip_if_not_installed("dartR.coancestry")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- gl.filter.allna(testset.gl[1:30, 1:200], verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  pop(x) <- rep(c("A", "B"), c(16, 14))
  set.seed(3)
  capture.output(
    res <- gl.diagnostics.relatedness(x, which_tests = "wang", run_sim = TRUE,
                                      Ne = c(30, 20),
                                      verbose = 0)
  )
  sim <- res@SimOutput[[1]]
  expect_equal(nLoc(sim), nLoc(x))
  # every stored generation has the sample sizes of x
  gen <- sim@other$ind.metrics$generation
  expect_true(all(table(gen, pop(sim))[, "A"] == 16))
  expect_true(all(table(gen, pop(sim))[, "B"] == 14))
  m <- res@MergedDf[[1]]
  expect_true("half_sibs" %in% m$RelDegree)
  expect_true(all(m$rel[m$RelDegree == "unrelated"] == 0))
})

test_that("a missing variable file takes real_freq from the one given", {
  skip_if_not_installed("dartR.coancestry")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- gl.filter.allna(testset.gl[1:20, 1:100], verbose = 0)
  capture.output(
    res <- gl.diagnostics.relatedness(
      x, which_tests = "wang", run_sim = TRUE,
      sim_variables = system.file("extdata", "sim_variables.csv",
                                  package = "dartR.sim"),
      verbose = 0)
  )
  # shipped sim file: real_freq = FALSE, one population of 50 in every
  # stored generation (3, plus the founders with dartR.sim >= store_founders)
  gen <- res@SimOutput[[1]]@other$ind.metrics$generation
  expect_true(all(table(gen) == 50))
})

test_that("founder inbreeding enters the pedigree kinship", {
  ped <- data.frame(id = c("A", "B", "C", "D"), dad = c(NA, NA, "A", "A"),
                    mom = c(NA, NA, "B", "B"), F = c(0.2, 0.1, NA, NA))
  K <- pedigreeKinship(ped)
  expect_equal(K["A", "A"], 0.6)
  expect_equal(K["C", "D"], (2 + 0.2 + 0.1) / 8)
  expect_equal(K["A", "C"], (0.6 + 0) / 2)
  expect_equal(pedigreeKinship(ped[, 1:3])["C", "D"], 0.25)
})

test_that("bias is the mean of estimate minus rel by class", {
  df <- data.frame(RelDegree = "full_sibs", rel = c(0.25, 0.25),
                   wang = c(0.2, 0.35))
  out <- calcBias(list(df), "wang")[[1]]
  expect_equal(out["wang", "full_sibs"], 0.025)
  expect_true(is.na(out["wang", "half_sibs"]))
})

test_that("copyMissing copies the missing pattern of a same-population donor", {
  x <- gl.filter.allna(testset.gl[1:10, 1:50], verbose = 0)
  pop(x) <- rep(c("A", "B"), each = 5)
  full <- gl.impute(x, method = "frequency", verbose = 0)
  sim <- rbind(full, full)
  indNames(sim) <- paste0("s", 1:20)
  pop(sim) <- rep(c("A", "B"), each = 5, times = 2)
  set.seed(1)
  out <- copyMissing(sim, x)
  mx <- is.na(as.matrix(x)); mo <- is.na(as.matrix(out))
  pats.A <- apply(mx[1:5, ], 1, paste, collapse = "")
  expect_true(all(apply(mo[pop(out) == "A", ], 1, paste, collapse = "") %in%
                  pats.A))
  expect_equal(as.matrix(out)[!mo], as.matrix(sim)[!mo])
})

test_that("simulation keeps parents, copies missing data and reports bias", {
  skip_if_not_installed("dartR.coancestry")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- gl.filter.allna(testset.gl[1:30, 1:200], verbose = 0)
  x <- gl.filter.monomorphs(x, verbose = 0)
  set.seed(5)
  capture.output(
    res <- gl.diagnostics.relatedness(x, which_tests = c("wang", "lynchrd"),
                                      run_sim = TRUE, biasOut = TRUE, Ne = 50,
                                      verbose = 0)
  )
  sim <- res@SimOutput[[1]]
  im <- sim@other$ind.metrics
  expect_equal(rownames(im), indNames(sim))
  expect_true(all(c("generation") %in% colnames(im)))
  expect_false(is.null(sim@other$sim.vars))
  expect_gt(mean(is.na(as.matrix(sim))), 0)
  expect_false(is.null(res@corOutList@biasPlot))
})

test_that("ind.metrics from generations bind by name with plain row names", {
  a <- data.frame(sex = "m", phenotype = "c", pat = NA, mat = NA,
                  F_founder = 0.1, row.names = "f1")
  b <- data.frame(sex = "f", phenotype = "c", pat = "f1", mat = "f2",
                  row.names = "o1")
  out <- bindIndMetrics(list(generation_0 = a, generation_1 = b))
  expect_equal(rownames(out), c("f1", "o1"))
  expect_equal(colnames(out)[3:4], c("pat", "mat"))
  expect_true(is.na(out["o1", "F_founder"]))
})

test_that("Ne is checked and resolved per population", {
  x <- testset.gl[1:20, 1:50]
  pop(x) <- rep(c("A", "B"), each = 10)
  expect_equal(resolveNe(x, 40, NULL, 0), c(40, 40))
  expect_equal(resolveNe(x, c(40, 60), NULL, 0), c(40, 60))
  expect_error(resolveNe(x, c(1, 2, 3), NULL, 0), "one per population")
  expect_error(resolveNe(x, NULL, NULL, 0), "give Ne, or neest.path")
  expect_error(gl.diagnostics.relatedness(x, Ne = -5, verbose = 0),
               "Ne must be positive")
})

test_that("ancestorPedigree keeps the ids and all their ancestors only", {
  ped <- data.frame(id = c("G1", "G2", "G3", "P1", "P2", "Q1", "C1", "C2"),
                    dad = c(NA, NA, NA, "G1", "G1", "G3", "P1", "P2"),
                    mom = c(NA, NA, NA, "G2", "G2", "G3", NA, NA))
  out <- ancestorPedigree(ped, c("C1", "C2"))
  expect_setequal(out$id, c("G1", "G2", "P1", "P2", "C1", "C2"))
  # cousins C1, C2 linked through unsampled parents and grandparents
  K <- pedigreeKinship(out)
  expect_equal(K["C1", "C2"], 0.0625)
})

test_that("an Ne estimate is unreliable with an infinite upper limit or >10n", {
  expect_equal(neUnreliable(c(36.5, 812, 88.3), c(93, Inf, 150), c(23, 17, 41)),
               c(FALSE, TRUE, FALSE))
  expect_true(neUnreliable(500, 900, 20))
})

test_that("bias tile plot shows each estimator and class present", {
  df <- data.frame(RelDegree = rep(c("full_sibs", "unrelated"), each = 2),
                   rel = c(0.25, 0.25, 0, 0),
                   wang = c(0.2, 0.3, 0.01, -0.03),
                   lynchrd = c(0.25, 0.27, 0, 0.02))
  p <- biasTilePlot(df, c("wang", "lynchrd"))
  expect_s3_class(p, "ggplot")
  expect_equal(nrow(p$data), 4)
  expect_equal(p$data$bias[p$data$estimator == "lynchrd" &
                             p$data$RelDegree == "full_sibs"], 0.01)
})

test_that("inbreeding on high-call-rate loci drops loci with lost calls", {
  x <- gl.filter.monomorphs(gl.filter.callrate(dartR.data::platypus.gl,
                                               threshold = 0.9, verbose = 0),
                            verbose = 0)
  f99 <- inbreedingHighCallrate(x, 0.99, 100)
  f90 <- inbreedingHighCallrate(x, 0.90, 100)
  expect_equal(names(f99), levels(pop(x)))
  expect_true(all(f99 < f90))
  expect_gte(attr(f99, "n.loci"), 100)
  # too few loci at 1: the threshold is lowered until all 150 pass
  f <- inbreedingHighCallrate(x[, 1:150], 1, 150)
  expect_equal(attr(f, "n.loci"), 150)
  expect_lt(attr(f, "threshold"), 1)
})

test_that("a single population simulates without shrinkage or migration", {
  skip_if_not_installed("dartR.coancestry")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- gl.filter.monomorphs(gl.filter.allna(testset.gl[1:24, 1:150],
                                            verbose = 0), verbose = 0)
  pop(x) <- rep("one", nInd(x))
  set.seed(9)
  capture.output(
    res <- gl.diagnostics.relatedness(x, which_tests = "wang", run_sim = TRUE,
                                      Ne = 30, verbose = 0)
  )
  expect_s4_class(res, "finalOutputClass")
  expect_gt(nrow(res@MergedDf[[1]]), 0)
})

test_that("family sizes come from an ind.metrics column, per population", {
  x <- testset.gl[1:12, 1:30]
  pop(x) <- rep(c("A", "B"), each = 6)
  x@other$ind.metrics$fam <- c("f1", "f1", "f1", "f2", "f2", NA,
                               "g1", "g1", "", "u1", "g2", "g2")
  s <- resolveFamilies(x, "fam", NULL, 0)
  expect_equal(s$A, c(3L, 2L))
  expect_equal(s$B, c(2L, 2L))
  expect_error(resolveFamilies(x, "nope", NULL, 0), "name of a column")
  expect_error(resolveFamilies(x, "colony", NULL, 0), "colony.path")
})

test_that("COLONY recovers full-sib families", {
  colony <- path.expand("~/programs")
  skip_if_not(file.exists(file.path(colony, "colony2s.out")), "COLONY not found")
  set.seed(2)
  p <- gl.filter.callrate(dartR.data::platypus.gl[1:4, ], threshold = 1,
                          verbose = 0)
  capture.output(off <- dartR.sim::gl.sim.offspring(p[1, ], p[2, ], 6,
                                                     verbose = 0))
  capture.output(off2 <- dartR.sim::gl.sim.offspring(p[3, ], p[4, ], 5,
                                                      verbose = 0))
  x <- rbind(off, off2)
  indNames(x) <- c(paste0("a", 1:6), paste0("b", 1:5))
  pop(x) <- rep("one", nInd(x))
  x <- gl.compliance.check(x, verbose = 0)
  s <- resolveFamilies(x, "colony", colony, 0)
  expect_equal(s$one, c(6L, 5L))
})

test_that("ExtractParents keeps an individual stored in two generations once", {
  mk <- function(ids, pat, mat) {
    g <- testset.gl[seq_along(ids), 1:10]
    indNames(g) <- ids
    g@other$ind.metrics <- data.frame(sex = "m", phenotype = "c", pat = pat,
                                      mat = mat, row.names = ids)
    g
  }
  sims <- list(list(generation_0 = mk(c("A", "B"), c(NA, NA), c(NA, NA)),
                    generation_1 = mk(c("A", "C"), c(NA, "A"), c(NA, "B"))))
  p <- ExtractParents(sims, 1)
  expect_equal(sort(p$id), c("A", "B", "C"))
})
