# Tests of the reviewed gl.diagnostics.relatedness and its helpers
# (function-review/reports/dartR.captive/gl.diagnostics.relatedness.md).
# Tests that run the relatedness estimators need the non-CRAN package
# 'related'.

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
  skip_if_not_installed("related")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- gl.filter.allna(testset.gl[1:15, ], verbose = 0)
  capture.output(
    res <- gl.diagnostics.relatedness(x, which_tests = "wang", verbose = 0)
  )
  expect_equal(nrow(res@MergedDf[[1]]) / 2, choose(15, 2))
})

test_that("attached pedigree: one row per pair, classes and rel (1, 2, 3)", {
  skip_if_not_installed("related")
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
  skip_if_not_installed("related")
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- gl.filter.allna(testset.gl[1:30, ], verbose = 0)
  set.seed(1)
  capture.output(
    res <- gl.diagnostics.relatedness(x, which_tests = "wang",
                                      run_sim = TRUE, rmseOut = TRUE,
                                      verbose = 0)
  )
  m <- res@MergedDf[[1]]
  expect_equal(nrow(m), choose(nInd(res@SimOutput[[1]]), 2))
  expect_false(anyDuplicated(m[, c("ind1", "ind2")]) > 0)
  expect_true("unrelated" %in% m$RelDegree)
})

test_that("SilicoDArT input errors (6)", {
  skip_if_not_installed("related")
  expect_error(gl.diagnostics.relatedness(testset.gs[1:10, 1:50],
                                          verbose = 0),
               "Only SNP data")
})
