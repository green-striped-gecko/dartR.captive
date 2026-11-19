# Tests of the reviewed gl.filter.parent.offspring
# (function-review/reports/dartR.captive/gl.filter.parent.offspring.md).

removed_from <- function(x, res) setdiff(indNames(x), indNames(res))

test_that("'best' keeps the better-genotyped member of each pair (1, 3)", {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- platypus.gl
  nas <- rowSums(is.na(as.matrix(x)))
  pairs <- gl.report.parent.offspring(x, verbose = 0)
  res <- gl.filter.parent.offspring(x, verbose = 0)
  removed <- removed_from(x, res)
  expect_setequal(removed, c("T3", "T42", "T36", "T13", "SUS36", "T39",
                             "SUS19", "T22"))
  # every pair is broken and no pair loses both members
  in.pair <- cbind(pairs$ind1 %in% removed, pairs$ind2 %in% removed)
  expect_true(all(rowSums(in.pair) == 1))
  # in pairs whose members are in no other pair, the removed one has more
  # missing genotypes
  expect_gt(nas[["T3"]], nas[["T5"]])
  expect_gt(nas[["T42"]], nas[["T28"]])
})

test_that("metadata stay in sync; one history entry (4)", {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  x <- platypus.gl
  res <- gl.filter.parent.offspring(x, rm.monomorphs = TRUE, verbose = 0)
  expect_equal(nrow(res@other$ind.metrics), nInd(res))
  expect_equal(length(pop(res)), nInd(res))
  expect_equal(nrow(res@other$loc.metrics), nLoc(res))
  expect_equal(length(res@other$history), length(x@other$history) + 1)
})

test_that("range is honoured by both methods and 3 does not crash (2)", {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_equal(nInd(gl.filter.parent.offspring(platypus.gl, range = 3,
                                               verbose = 0)), 81)
  set.seed(1)
  expect_equal(nInd(gl.filter.parent.offspring(platypus.gl,
                                               method = "random",
                                               range = 3, verbose = 0)), 81)
  set.seed(1)
  res <- gl.filter.parent.offspring(platypus.gl, method = "random",
                                    verbose = 0)
  expect_equal(nInd(res), 73)
})

test_that("no pairs: object unchanged, quiet at verbose 0 (5)", {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  out <- capture.output(
    res <- gl.filter.parent.offspring(testset.gl[1:10, 1:50], verbose = 0)
  )
  expect_equal(nInd(res), 10)
  expect_length(out, 0)
})

test_that("SilicoDArT and unknown method error (2, 5)", {
  pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_error(gl.filter.parent.offspring(testset.gs[1:20, 1:50],
                                          verbose = 0), "Only SNP data")
  expect_error(gl.filter.parent.offspring(platypus.gl, method = "Best",
                                          verbose = 0), "should be one of")
})
