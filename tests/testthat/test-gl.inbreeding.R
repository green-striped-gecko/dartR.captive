# Tests of gl.inbreeding (function-review/reports/dartR.captive/gl.inbreeding.md).
# Needs PLINK 1.9 and GCTA; set DARTR_PLINK_PATH and DARTR_GCTA_PATH to the
# executables or their folders (default ~/programs).

plink_path <- function() {
  Sys.getenv("DARTR_PLINK_PATH", path.expand("~/programs"))
}
gcta_path <- function() {
  Sys.getenv("DARTR_GCTA_PATH", path.expand("~/programs"))
}
# the executable must exist and run on this platform
runs <- function(path, name, banner) {
  exe <- if (file_test("-f", path)) path else file.path(path, name)
  if (!file_test("-f", exe)) {
    return(FALSE)
  }
  out <- tryCatch(suppressWarnings(system2(exe, stdout = TRUE, stderr = TRUE)),
                  error = function(e) "")
  any(grepl(banner, out))
}
skip_without_plink <- function() {
  skip_if_not(runs(plink_path(), "plink", "PLINK"), "PLINK not found")
}
skip_without_gcta <- function() {
  skip_if_not(runs(gcta_path(), "gcta64", "GCTA"), "GCTA not found")
}

# Independent computation of PLINK --het F and GCTA --ibc Fhat1-3 on the
# loci polymorphic in m (individuals x loci, 0/1/2/NA)
ref_inbreeding <- function(m) {
  p <- colMeans(m, na.rm = TRUE) / 2
  k <- !is.na(p) & p > 0 & p < 1
  m <- m[, k, drop = FALSE]
  p <- p[k]
  h <- 2 * p * (1 - p)
  P <- matrix(p, nrow(m), ncol(m), byrow = TRUE)
  H <- matrix(h, nrow(m), ncol(m), byrow = TRUE)
  obs <- !is.na(m)
  O <- rowSums(obs & m != 1)
  E <- as.vector(obs %*% (1 - h))
  N <- rowSums(obs)
  data.frame(N = N,
             F = (O - E) / (N - E),
             Fhat1 = rowMeans((m - 2 * P)^2 / H, na.rm = TRUE) - 1,
             Fhat2 = 1 - rowMeans(m * (2 - m) / H, na.rm = TRUE),
             Fhat3 = rowMeans((m^2 - (1 + 2 * P) * m + 2 * P^2) / H,
                              na.rm = TRUE))
}

test_that("pooled estimates match PLINK and GCTA formulas on the same loci (1, 3)", {
  skip_without_plink()
  skip_without_gcta()
  x <- platypus.gl
  # loci at position 0 were coded chromosome 0 by gl2plink and dropped by
  # GCTA in the original script (F1)
  x$chromosome <- factor(rep("1", nLoc(x)))
  x$position <- rep(c(0L, 100L), length.out = nLoc(x))
  res <- gl.inbreeding(x, plink.path = plink_path(), gcta.path = gcta_path(),
                       verbose = 0)
  ref <- ref_inbreeding(as.matrix(x))
  expect_named(res, c("id", "pop", "O.HOM", "E.HOM", "N.NM", "F",
                      "NOMISS", "Fhat1", "Fhat2", "Fhat3"))
  expect_equal(res$id, indNames(x))
  expect_equal(res$N.NM, ref$N)
  expect_equal(res$NOMISS, ref$N)
  # PLINK prints 4 decimals, GCTA 6 significant digits
  expect_lt(max(abs(res$F - ref$F)), 1e-4)
  expect_lt(max(abs(res$Fhat1 - ref$Fhat1)), 1e-4)
  expect_lt(max(abs(res$Fhat2 - ref$Fhat2)), 1e-4)
  expect_lt(max(abs(res$Fhat3 - ref$Fhat3)), 1e-4)
})

test_that("by.pop takes allele frequencies within each population (2)", {
  skip_without_plink()
  skip_without_gcta()
  x <- platypus.gl
  res <- gl.inbreeding(x, by.pop = TRUE, plink.path = plink_path(),
                       gcta.path = gcta_path(), verbose = 0)
  m <- as.matrix(x)
  for (p in levels(pop(x))) {
    i <- which(pop(x) == p)
    ref <- ref_inbreeding(m[i, , drop = FALSE])
    expect_equal(res$NOMISS[i], ref$N)
    expect_lt(max(abs(res$F[i] - ref$F)), 1e-4)
    expect_lt(max(abs(res$Fhat3[i] - ref$Fhat3)), 1e-4)
  }
  # pooled frequencies give different estimates on multi-population data
  pooled <- gl.inbreeding(x, method = "plink", plink.path = plink_path(),
                          verbose = 0)
  expect_false(isTRUE(all.equal(pooled$F, res$F)))
})

test_that("a population of one individual gets NA with by.pop (2)", {
  skip_without_plink()
  skip_without_gcta()
  x <- platypus.gl
  pop(x) <- factor(c("single", rep("rest", nInd(x) - 1)))
  res <- gl.inbreeding(x, by.pop = TRUE, plink.path = plink_path(),
                       gcta.path = gcta_path(), verbose = 0)
  expect_true(all(is.na(res[1, c("F", "Fhat3")])))
  expect_false(anyNA(res$F[-1]))
  expect_output(gl.inbreeding(x, by.pop = TRUE, method = "plink",
                              plink.path = plink_path(), verbose = 1),
                "one individual")
})

test_that("method = 'plink' runs without GCTA (4)", {
  skip_without_plink()
  res <- gl.inbreeding(platypus.gl, method = "plink",
                       plink.path = plink_path(),
                       gcta.path = "no/such/folder", verbose = 0)
  expect_named(res, c("id", "pop", "O.HOM", "E.HOM", "N.NM", "F"))
  expect_false(anyNA(res$F))
})

test_that("plot is optional and saved only on request (5)", {
  skip_without_plink()
  skip_without_gcta()
  dir <- withr::local_tempdir()
  res <- gl.inbreeding(platypus.gl, plink.path = plink_path(),
                       gcta.path = gcta_path(), plot.display = FALSE,
                       plot.file = "inb", plot.dir = dir, verbose = 0)
  expect_equal(nrow(res), nInd(platypus.gl))
  expect_true(file.exists(file.path(dir, "inb.RDS")))
  expect_s3_class(readRDS(file.path(dir, "inb.RDS")), "patchwork")
})

test_that("input and executables are checked before any work (5)", {
  expect_error(gl.inbreeding(testset.gs, verbose = 0))
  expect_error(gl.inbreeding(platypus.gl, plink.path = "no/such/folder",
                             verbose = 0), "Cannot find PLINK")
  expect_error(gl.inbreeding(platypus.gl, by.pop = NA, verbose = 0),
               "by.pop")
  expect_error(gl.inbreeding(platypus.gl, method = "beagle", verbose = 0))
})
