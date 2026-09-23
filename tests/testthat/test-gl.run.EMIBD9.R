# Tests of the reviewed gl.run.EMIBD9
# (function-review/reports/dartR.captive/gl.run.EMIBD9.md). Needs the EMIBD9 binary; set
# DARTR_EMIBD9_PATH to its folder (default ~/programs/emibd9-v1.0).

emibd9_dir <- function() {
  Sys.getenv("DARTR_EMIBD9_PATH", path.expand("~/programs/emibd9-v1.0"))
}
skip_without_emibd9 <- function() {
  exe <- if (Sys.info()[["sysname"]] == "Windows") "EM_IBD_P.exe" else
    "EM_IBD_P"
  skip_if_not(file.exists(file.path(emibd9_dir(), exe)), "EMIBD9 not found")
}
emibd9_fixture <- function(loc = 1:120) {
  gl.filter.allna(testset.gl[1:10, loc], verbose = 0)
}

test_that("result structure, kinship scale and numeric columns (3, 5, 8)", {
  skip_without_emibd9()
  res <- gl.run.EMIBD9(emibd9_fixture(), emibd9.path = emibd9_dir(),
                       plot.out = FALSE, verbose = 0)
  expect_named(res, c("rel", "raw", "processed", "inbreeding"))
  expect_equal(dim(res$rel), c(10, 10))
  expect_equal(unname(diag(res$rel)), rep(0.5, 10))
  expect_equal(attr(res$rel, "scale"), "kinship")
  expect_equal(nrow(res$processed), 45)
  num.cols <- setdiff(names(res$raw), c("Indiv1", "Indiv2"))
  expect_true(all(vapply(res$raw[num.cols], is.numeric, logical(1))))
  expect_true(is.numeric(res$processed[["r(1,2)"]]))
  # r(1,2) equals the kinship coefficient: D1 + (D3 + D5 + D7) / 2 + D8 / 4
  p <- as.data.frame(res$processed)
  theta <- p$Delta1 + (p$Delta3 + p$Delta5 + p$Delta7) / 2 + p$Delta8 / 4
  expect_equal(theta, p[["r(1,2)"]], tolerance = 1e-3)
})

test_that("change 1/2: OutAlleleFre = TRUE runs and results are fresh", {
  skip_without_emibd9()
  r1 <- gl.run.EMIBD9(emibd9_fixture(1:120), emibd9.path = emibd9_dir(),
                      plot.out = FALSE, verbose = 0)
  r2 <- gl.run.EMIBD9(emibd9_fixture(121:250), emibd9.path = emibd9_dir(),
                      OutAlleleFre = TRUE, plot.out = FALSE, verbose = 0)
  expect_false(isTRUE(all.equal(r1$rel, r2$rel)))
})

test_that("change 1: a failed EMIBD9 run stops with its console output", {
  skip_without_emibd9()
  expect_error(
    gl.run.EMIBD9(emibd9_fixture(), emibd9.path = emibd9_dir(),
                  EM_Method = "bad", plot.out = FALSE, verbose = 0),
    "EMIBD9 did not write its results"
  )
})

test_that("change 7: plot.file works with plot.out = FALSE", {
  skip_without_emibd9()
  d <- withr::local_tempdir()
  res <- gl.run.EMIBD9(emibd9_fixture(), emibd9.path = emibd9_dir(),
                       plot.out = FALSE, plot.file = "p", plot.dir = d,
                       verbose = 0)
  expect_true(file.exists(file.path(d, "p.RDS")))
})

test_that("change 4: SilicoDArT input errors", {
  expect_error(
    gl.run.EMIBD9(testset.gs[1:10, 1:100], emibd9.path = emibd9_dir(),
                  plot.out = FALSE, verbose = 0),
    "Only SNP data"
  )
})

test_that("change 7: missing executable names the file", {
  expect_error(
    gl.run.EMIBD9(emibd9_fixture(), emibd9.path = tempfile(),
                  plot.out = FALSE, verbose = 0),
    "Cannot find EM_IBD_P"
  )
})

test_that("change 7: the output file is copied to outpath", {
  skip_without_emibd9()
  d <- withr::local_tempdir()
  gl.run.EMIBD9(emibd9_fixture(), emibd9.path = emibd9_dir(), outpath = d,
                plot.out = FALSE, verbose = 0)
  expect_true(file.exists(file.path(d, "EMIBD9_Res.ibd9")))
})
