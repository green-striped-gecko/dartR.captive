# Characterization baseline (2026-09-23) plus tests of the approved review
# changes (function-review/reports/dartR.captive/gl2colony.md).

colony_fixture <- function() {
  x <- testset.gl[1:20, 1:30]
  x@other$ind.metrics$offspring <- rep(c("yes", "no", "no"), c(10, 5, 5))
  x@other$ind.metrics$father <- rep(c("no", "yes", "no"), c(10, 5, 5))
  x@other$ind.metrics$mother <- rep(c("no", "no", "yes"), c(10, 5, 5))
  x
}

test_that("baseline: default export layout", {
  x <- colony_fixture()
  od <- tempdir()
  capture.output(
    f <- gl2colony(x, outfile = "t_base.dat", outpath = od, seed = 1,
                   verbose = 0)
  )
  l <- readLines(f)
  expect_equal(f, file.path(od, "t_base.dat"))
  expect_length(l, 62)
  expect_equal(trimws(l[3]), "10 \t ! No. offspring")
  expect_equal(trimws(l[4]), "30 \t ! No. of loci")
  # offspring row: ID then 2 alleles per locus
  expect_length(strsplit(l[27], " ")[[1]], 1 + 2 * 30)
  expect_equal(strsplit(l[27], " ")[[1]][1:5],
               c("AA010915", "2", "2", "0", "0"))
  expect_true(any(grepl("^5 5 \t ! Number of candidates", l)))
})

test_that("change 5: nothing prints at verbose = 0", {
  x <- colony_fixture()
  out <- capture.output(
    gl2colony(x, outfile = "t_v0.dat", outpath = tempdir(), seed = 1,
              verbose = 0)
  )
  expect_length(out, 0)
})

test_that("change 2: a missing role column keeps the other roles", {
  x <- colony_fixture()
  x@other$ind.metrics$father <- NULL
  capture.output(
    f <- gl2colony(x, outfile = "t_p4.dat", outpath = tempdir(), seed = 1,
                   verbose = 0)
  )
  l <- readLines(f)
  expect_equal(trimws(l[3]), "10 \t ! No. offspring")
  expect_true(any(grepl("^0 5 \t ! Number of candidates", l)))
})

test_that("change 2: role columns and values match ignoring case", {
  x <- colony_fixture()
  im <- x@other$ind.metrics
  names(im)[match(c("offspring", "father", "mother"), names(im))] <-
    c("Offspring", "Father", "Mother")
  im$Mother[16] <- " YES "
  x@other$ind.metrics <- im
  out <- capture.output(
    f <- gl2colony(x, outfile = "t_case.dat", outpath = tempdir(), seed = 1,
                   verbose = 1)
  )
  expect_false(any(grepl("not found", out)))
  expect_true(any(grepl("^5 5 \t ! Number of candidates", readLines(f))))
})

test_that("change 2: no ind.metrics exports everyone as offspring", {
  x <- testset.gl[1:5, 1:5]
  x@other$ind.metrics <- NULL
  capture.output(
    f <- gl2colony(x, outfile = "t_nm.dat", outpath = tempdir(), seed = 1,
                   verbose = 0)
  )
  expect_equal(trimws(readLines(f)[3]), "5 \t ! No. offspring")
})

test_that("change 4: SilicoDArT input errors", {
  x <- colony_fixture()
  s <- testset.gs[1:20, 1:10]
  s@other$ind.metrics <- x@other$ind.metrics[, c("offspring", "father",
                                                  "mother")]
  expect_error(
    gl2colony(s, outfile = "t_gs.dat", outpath = tempdir(), verbose = 0),
    "Only SNP data"
  )
})

test_that("change 1: individual names with whitespace error", {
  x <- colony_fixture()
  indNames(x)[1] <- "ind 1"
  expect_error(
    gl2colony(x, outfile = "t_sp.dat", outpath = tempdir(), verbose = 0),
    "ind 1"
  )
})

test_that("change 3: settings needing extra data blocks error", {
  x <- colony_fixture()
  od <- tempdir()
  expect_error(gl2colony(x, outpath = od, sibship.prior = 2, verbose = 0),
               "sibship.prior")
  expect_error(gl2colony(x, outpath = od, known.allele.freq = 1,
                         verbose = 0), "known.allele.freq")
  expect_error(gl2colony(x, outpath = od,
                         paternity.exclusion.threshold = "1 0", verbose = 0),
               "paternity.exclusion.threshold")
  expect_error(gl2colony(x, outpath = od, excluded.paternity = 1,
                         verbose = 0), "excluded.paternity")
})

test_that("change 3: a single rate without '@' applies to all loci", {
  x <- colony_fixture()
  capture.output(
    f <- gl2colony(x, outfile = "t_at.dat", outpath = tempdir(), seed = 1,
                   allelic.dropout = "0.01", other.typ.err = "0.001",
                   verbose = 0)
  )
  l <- readLines(f)
  expect_equal(trimws(l[25]), "0.01@ \t ! Allelic dropout rate")
  expect_equal(trimws(l[26]), "0.001@ \t ! Other typing error rate")
})
