test_that("gl.relatedness returns the numbers of the dartR_coancestry program", {
  skip_if_no_coancestry()
  skip_if_not_installed("dartR.data")
  gl  <- dartR.data::platypus.gl[1:40, 1:1000]
  est <- c("wang", "lynchrd", "dyadml")
  got <- gl.relatedness(gl, estimators = est, plot.out = FALSE,
                        coancestry.path = coancestry_path(), verbose = 0)
  expect_equal(nrow(got$dyads), choose(40, 2))
  # anchors from the engine validated against Coancestry (trior11); the
  # tolerance covers the rounding difference between processors (< 1e-6)
  rows <- c(1, 2, 400, 780)
  expect_equal(got$dyads$wang[rows],
               c(-0.2113080, -0.2130594, -0.1273998, -0.0872843), tolerance = 1e-6)
  expect_equal(got$dyads$lynchrd[rows],
               c(-0.0518246, -0.0417763, -0.0117164, -0.0192881), tolerance = 1e-6)
  expect_equal(got$dyads$dyadml[rows],
               c(8.562220e-03, 6.765268e-07, 5.350713e-02, 3.051880e-02), tolerance = 1e-6)
})

test_that("a missing program stops with the download instructions", {
  skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:10, 1:100]
  expect_error(gl.relatedness(gl, estimators = "wang", plot.out = FALSE,
                              coancestry.path = withr::local_tempdir(), verbose = 0),
               "gl.download.binary")
})

test_that("a failure of the program is reported with its message", {
  skip_on_os("windows")
  skip_if_not_installed("dartR.data")
  d <- withr::local_tempdir()
  fake <- file.path(d, "dartR_coancestry")
  writeLines(c("#!/bin/sh",
               'if [ "$1" = "--version" ]; then echo "dartR_coancestry 0.2.0"; exit 0; fi',
               'echo "dartR_coancestry: error: boom" >&2', "exit 1"), fake)
  Sys.chmod(fake, "755")
  gl <- dartR.data::platypus.gl[1:10, 1:100]
  wd <- getwd()
  expect_error(gl.relatedness(gl, estimators = "wang", plot.out = FALSE,
                              coancestry.path = d, verbose = 0), "boom")
  expect_identical(getwd(), wd)          # working directory restored after the error
})

test_that("the program's last error line is reported even after a long log", {
  skip_on_os("windows")
  skip_if_not_installed("dartR.data")
  d <- withr::local_tempdir()
  fake <- file.path(d, "dartR_coancestry")
  writeLines(c("#!/bin/sh", "seq 1 20000",
               'echo "dartR_coancestry: error: boom-last" >&2', "exit 1"), fake)
  Sys.chmod(fake, "755")
  gl <- dartR.data::platypus.gl[1:10, 1:100]
  expect_error(gl.relatedness(gl, estimators = "wang", plot.out = FALSE,
                              coancestry.path = d, verbose = 0), "boom-last")
})

test_that("a folder or a file that cannot be run is not taken for the program", {
  skip_on_os("windows")
  d <- withr::local_tempdir()
  dir.create(file.path(d, "dartR_coancestry"))
  expect_error(utils.coancestry.exe(d), "is a folder")
  d2 <- withr::local_tempdir()
  writeLines("#!/bin/sh", file.path(d2, "dartR_coancestry"))
  Sys.chmod(file.path(d2, "dartR_coancestry"), "644")
  expect_error(utils.coancestry.exe(d2), "chmod")
})

test_that("no estimator, or a number beyond the integer range, is refused before the program runs", {
  skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:10, 1:100]
  d <- withr::local_tempdir()           # no program here: the checks must come first
  run <- function(...) gl.relatedness(gl, plot.out = FALSE, coancestry.path = d,
                                      verbose = 0, ...)
  expect_error(run(estimators = character(0)), "at least one estimator")
  expect_error(run(estimators = "wang", rng.seed = 1e10), "rng.seed")
  expect_error(run(estimators = "wang", n.boots = 3e9), "n.boots")
  expect_error(run(estimators = "wang", n.threads = 3e9), "n.threads")
  expect_error(run(estimators = "wang", num.trios = 3e9), "num.trios")
})

test_that("genotypes written in blocks of loci are the same bytes", {
  skip_if_not_installed("dartR.data")
  snp <- as.matrix(dartR.data::platypus.gl[1:12, 1:300])
  storage.mode(snp) <- "integer"
  d <- withr::local_tempdir()
  utils.coancestry.write(snp, file.path(d, "all.bin"))
  for (b in c(1L, 7L, 299L)) {
    utils.coancestry.write(snp, file.path(d, "blk.bin"), block.loci = b)
    expect_identical(readBin(file.path(d, "blk.bin"), "raw", 1e5),
                     readBin(file.path(d, "all.bin"), "raw", 1e5))
  }
  r <- readBin(file.path(d, "all.bin"), "raw", 1e5)
  expect_identical(rawToChar(r[1:8]), "DRCG0001")
  expect_length(r, 16L + 12L * 300L)
  expect_identical(r[17:28], as.raw(ifelse(is.na(snp[, 1]), 3L, snp[, 1])))
})

test_that("results the reader does not fully understand are an error", {
  skip_on_os("windows")
  skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:10, 1:100]
  # stand-in programs that write a manifest the reader must refuse
  fake <- function(manifest.lines, bytes) {
    d <- withr::local_tempdir(.local_envir = parent.frame())
    writeLines(c("#!/bin/sh",
                 sprintf("printf '%s' > res.manifest", paste0(manifest.lines, "\\n", collapse = "")),
                 sprintf("head -c %d /dev/zero > res.dyads.bin", bytes)),
               file.path(d, "dartR_coancestry"))
    Sys.chmod(file.path(d, "dartR_coancestry"), "755")
    d
  }
  run <- function(d) gl.relatedness(gl, estimators = "wang", plot.out = FALSE,
                                    coancestry.path = d, verbose = 0)
  ok.line <- "dyads\\t45\\tind1:i,ind2:i,wang:d"
  expect_error(run(fake(c("another_program 1", ok.line), 45 * 16)), "manifest")
  expect_error(run(fake(c("dartR_coancestry 0.2.0", "dyads\\t45\\tind1:i,ind2:i,wang:x"), 45 * 16)),
               "column type")
  expect_error(run(fake(c("dartR_coancestry 0.2.0", ok.line), 45 * 16 - 1)), "size")
  # a column named twice would silently overwrite the first
  expect_error(run(fake(c("dartR_coancestry 0.2.0", "dyads\\t45\\tind1:i,ind2:i,wang:d,wang:d"),
                        45 * 24)), "malformed")
})

test_that("the program runs from a folder with a space; the working directory is restored", {
  skip_if_no_coancestry()
  skip_if_not_installed("dartR.data")
  exe <- if (.Platform$OS.type == "windows") "dartR_coancestry.exe" else "dartR_coancestry"
  d <- file.path(withr::local_tempdir(), "my programs"); dir.create(d)
  file.copy(file.path(coancestry_path(), exe), d, copy.mode = TRUE)
  gl <- dartR.data::platypus.gl[1:10, 1:100]
  wd <- getwd()
  a <- gl.relatedness(gl, estimators = "wang", plot.out = FALSE,
                      coancestry.path = d, verbose = 0)
  b <- gl.relatedness(gl, estimators = "wang", plot.out = FALSE,
                      coancestry.path = coancestry_path(), verbose = 0)
  expect_identical(a$dyads, b$dyads)
  expect_identical(getwd(), wd)
})

test_that("an interrupt (the Stop button) stops the program and removes its files", {
  skip_on_os("windows")
  skip_on_cran()
  d <- withr::local_tempdir()
  # a stand-in program that records its process and folder, then interrupts R
  # alone, as RStudio's Stop button does, and keeps running
  writeLines(c("#!/bin/sh",
               sprintf('echo $$ > "%s/pid"; pwd > "%s/wd"', d, d),
               "sleep 1",
               'kill -INT "$DARTR_TEST_R_PID"',
               "exec sleep 30"),
             file.path(d, "dartR_coancestry"))
  Sys.chmod(file.path(d, "dartR_coancestry"), "755")
  withr::local_envvar(DARTR_TEST_R_PID = Sys.getpid())
  snp <- matrix(c(0L, 1L, 2L, NA), 4, 3)
  wd <- getwd()
  t0 <- Sys.time()
  got <- tryCatch(utils.coancestry.run(snp, "wang", exe = file.path(d, "dartR_coancestry")),
                  interrupt = function(e) "interrupted")
  expect_identical(got, "interrupted")
  expect_lt(as.numeric(difftime(Sys.time(), t0, units = "secs")), 10)
  Sys.sleep(0.5)
  pid <- readLines(file.path(d, "pid"))
  # ps lists nothing for a process that no longer exists
  expect_length(suppressWarnings(system2("ps", c("-o", "pid=", "-p", pid), stdout = TRUE)), 0)
  expect_false(dir.exists(readLines(file.path(d, "wd"))))
  expect_identical(getwd(), wd)
})

test_that("gl.relatedness rejects an unknown estimator", {
  skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:10, 1:100]
  expect_error(gl.relatedness(gl, estimators = "bogus", verbose = 0),
               regexp = "Unknown estimator")
})

test_that("gl.relatedness returns symmetric per-estimator matrices keyed by name", {
  skip_if_no_coancestry(); skip_if_not_installed("dartR.data")
  gl  <- dartR.data::platypus.gl[1:40, 1:1000]
  got <- gl.relatedness(gl, coancestry.path = coancestry_path(), estimators = "wang", plot.out = FALSE, verbose = 0)
  M <- got$wang
  expect_equal(dim(M), c(40L, 40L))
  expect_identical(rownames(M), indNames(gl))
  expect_identical(colnames(M), indNames(gl))
  expect_true(all(is.na(diag(M))))
  expect_equal(M[lower.tri(M)], t(M)[lower.tri(M)])          # symmetric
  d <- got$dyads
  i <- match(d$ind1[1], indNames(gl)); j <- match(d$ind2[1], indNames(gl))
  expect_equal(M[i, j], d$wang[1])                            # matrix matches long table
  expect_true(all(d$ind1 %in% indNames(gl)))                 # ind cols are names now
})

test_that("gl.relatedness default estimators and optional frames", {
  skip_if_no_coancestry(); skip_if_not_installed("dartR.data")
  gl  <- dartR.data::platypus.gl[1:40, 1:1000]
  got <- gl.relatedness(gl, coancestry.path = coancestry_path(), plot.out = FALSE, verbose = 0)
  expect_true(all(c("wang","lynchli","lynchrd","ritland","quellergt","loiselle")
                  %in% names(got)))
  expect_true("dyads" %in% names(got))
  expect_null(got$delta19); expect_null(got$trio_delta); expect_null(got$inbreeding)
})

test_that("gl.relatedness pre-filters and produces CIs when bootstrapping", {
  skip_if_no_coancestry(); skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:40, 1:1000]   # contains 405 monomorphic/all-NA loci
  # dartR convention: the pre-filter reports through cat(warn(...)) gated on
  # verbose, it never signals an R warning condition.
  expect_no_warning(
    got <- gl.relatedness(gl, coancestry.path = coancestry_path(), estimators = "wang", n.boots = 20,
                          plot.out = FALSE, verbose = 0))
  expect_true(all(c("wang_lo", "wang_hi") %in% names(got$dyads)))
})

test_that("gl.relatedness pre-filter reports only at verbose >= 1", {
  skip_if_no_coancestry(); skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:40, 1:1000]
  expect_output(
    gl.relatedness(gl, coancestry.path = coancestry_path(), estimators = "wang", n.boots = 5,
                   plot.out = FALSE, verbose = 1),
    regexp = "pre-filter")
  quiet <- capture.output(
    gl.relatedness(gl, coancestry.path = coancestry_path(), estimators = "wang", n.boots = 5,
                   plot.out = FALSE, verbose = 0))
  expect_false(any(grepl("pre-filter", quiet)))
})

test_that("gl.relatedness plotting is optional; bad plot.stat errors", {
  skip_if_no_coancestry(); skip_if_not_installed("dartR.data")
  gl  <- dartR.data::platypus.gl[1:40, 1:1000]
  res <- gl.relatedness(gl, coancestry.path = coancestry_path(), estimators = c("wang", "lynchrd"),
                        plot.out = FALSE, verbose = 0)
  expect_true(is.list(res) && "wang" %in% names(res))
  expect_error(
    gl.relatedness(gl, coancestry.path = coancestry_path(), estimators = "wang", plot.stat = "trioml",
                   plot.out = TRUE, verbose = 0),
    regexp = "plot.stat")
})

test_that("gl.relatedness populates and name-maps optional frames", {
  skip_if_no_coancestry(); skip_if_not_installed("dartR.data")
  gl  <- dartR.data::platypus.gl[1:20, 1:300]
  got <- gl.relatedness(gl, coancestry.path = coancestry_path(), estimators = c("dyadml", "inbreeding"),
                        allow.inbreeding = TRUE, plot.out = FALSE, verbose = 0)
  expect_true("dyadml" %in% names(got))
  expect_false(is.null(got$delta19))
  expect_false(is.null(got$inbreeding))
  expect_true(all(got$delta19$ind1 %in% indNames(gl)))   # mapped to names
  expect_true(all(got$inbreeding$ind %in% indNames(gl)))
  expect_true(all(c("LH", "LR", "dyadml_F") %in% names(got$inbreeding)))
})

test_that("gl.relatedness with only 'inbreeding' returns no dyad matrices", {
  skip_if_no_coancestry(); skip_if_not_installed("dartR.data")
  gl  <- dartR.data::platypus.gl[1:20, 1:300]
  got <- gl.relatedness(gl, coancestry.path = coancestry_path(), estimators = "inbreeding", plot.out = FALSE, verbose = 0)
  expect_false(is.null(got$inbreeding))
  expect_true(all(!c("wang","lynchli","lynchrd","ritland","quellergt",
                     "loiselle","dyadml","trioml") %in% names(got)))
})

# Tests added by the function review
# (function-review/reports/dartR.captive/gl.relatedness.md).

test_that("agreement with related::coancestry (Coancestry Fortran)", {
  skip_if_no_coancestry(); skip_if_not_installed("dartR.data")
  skip_if_not_installed("related")
  # related reads and writes files in the working directory
  withr::local_dir(withr::local_tempdir())
  # named through a variable: related is not a declared dependency
  coancestry <- getExportedValue("related", "coancestry")
  gl <- gl.filter.monomorphs(dartR.data::platypus.gl[1:30, ], verbose = 0)
  gl <- gl.filter.callrate(gl, method = "loc", threshold = 1, verbose = 0)
  gl <- gl.filter.monomorphs(gl, verbose = 0)
  est <- c("wang", "lynchli", "lynchrd", "ritland", "quellergt")
  got <- gl.relatedness(gl, coancestry.path = coancestry_path(), estimators = est, plot.out = FALSE, verbose = 0)
  m <- as.matrix(gl)
  gd <- data.frame(id = indNames(gl))
  for (j in seq_len(ncol(m))) {
    gd[[paste0("L", j, "a")]] <- ifelse(m[, j] == 2, 2, 1)
    gd[[paste0("L", j, "b")]] <- ifelse(m[, j] == 0, 1, 2)
  }
  invisible(capture.output(cc <- suppressWarnings(
    coancestry(gd, lynchli = 1, lynchrd = 1, quellergt = 1,
                        ritland = 1, wang = 1))))
  co <- cc$relatedness
  key <- function(a, b) paste(pmin(a, b), pmax(a, b))
  ix <- match(key(got$dyads$ind1, got$dyads$ind2),
              key(co$ind1.id, co$ind2.id))
  # related rounds its output to 4 decimals
  for (e in c("lynchrd", "ritland", "quellergt")) {
    expect_lt(max(abs(got$dyads[[e]] - co[[e]][ix])), 1e-4)
  }
  # wang and lynchli sit above related by an amount that shrinks with n
  # (documented in Details)
  expect_equal(mean(got$dyads$wang - co$wang[ix]), 0.0469, tolerance = 0.01)
  expect_equal(mean(got$dyads$lynchli - co$lynchli[ix]), 0.0183,
               tolerance = 0.01)
})

test_that("output contract: invisible, tagged, NA for pairs without data", {
  skip_if_no_coancestry(); skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:12, 1:300]
  vis <- withVisible(gl.relatedness(gl, coancestry.path = coancestry_path(), estimators = "wang",
                                    plot.out = FALSE, verbose = 0))
  expect_false(vis$visible)                                   # change 7
  expect_equal(attr(vis$value$wang, "scale"), "relatedness")  # change 4
  # an individual with no calls: NA in every estimator (was 0 / NaN; change 2)
  gx <- gl
  gx@gen[[1]] <- new("SNPbin", rep(NA_integer_, nLoc(gx)), ploidy = 2L)
  expect_output(r <- gl.relatedness(gx, coancestry.path = coancestry_path(), estimators = c("wang", "dyadml"),
                                    plot.out = FALSE, verbose = 1),
                "11 pairs share no called locus")
  expect_true(all(is.na(r$wang[1, ])))
  expect_true(all(is.na(r$dyadml[1, ])))
  expect_false(anyNA(r$wang[2:12, 2:12][upper.tri(diag(11))]))
  expect_equal(sum(is.na(r$dyads$wang)), 11)
  # with n.boots > 0 the dropped individual keeps an NA row (was 11 x 11)
  rb <- gl.relatedness(gx, coancestry.path = coancestry_path(), estimators = "wang", n.boots = 3,
                       plot.out = FALSE, verbose = 0)
  expect_equal(dim(rb$wang), c(12L, 12L))
  expect_identical(rownames(rb$wang), indNames(gx))
  expect_true(all(is.na(rb$wang[1, ])))
})

test_that("arguments are validated before the engine runs (change 5)", {
  skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:10, 1:100]
  bad <- list(list(n.boots = NA), list(n.boots = -1), list(num.trios = 0),
              list(n.threads = 0), list(rng.seed = NA),
              list(allow.inbreeding = NA), list(plot.out = "yes"),
              list(plot.stat = "trioml"))
  for (b in bad) {
    expect_error(do.call(gl.relatedness,
                         c(list(gl, estimators = "wang", verbose = 0), b)),
                 paste(names(b), "must be"))
  }
})

test_that("the heatmap is saved without being displayed (change 6)", {
  skip_if_no_coancestry(); skip_if_not_installed("dartR.data")
  gl <- dartR.data::platypus.gl[1:10, 1:200]
  d <- withr::local_tempdir()
  n.dev <- length(grDevices::dev.list())
  out <- capture.output(
    gl.relatedness(gl, coancestry.path = coancestry_path(), estimators = "wang", plot.out = FALSE,
                   plot.file = "h", plot.dir = d, verbose = 0))
  expect_length(list.files(d), 1)
  expect_equal(length(grDevices::dev.list()), n.dev)
  expect_false(any(grepl("gl.colors", out)))
})
