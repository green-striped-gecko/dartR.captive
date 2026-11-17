# Internal helpers that run the standalone program 'dartR_coancestry', a
# closed-source C++ implementation of Jinliang Wang's Coancestry estimators that
# is downloaded separately (see ?gl.relatedness). The program and R exchange
# binary files, all values little-endian:
#   genotypes in: 8 bytes "DRCG0001", int32 individuals, int32 loci, then one
#                 byte per genotype, locus after locus (the order of an R
#                 matrix): 0, 1, 2 = dosage, 3 = missing
#   results out:  <prefix>.manifest, a text file (line 1: program and version;
#                 then one line per table: name TAB rows TAB col:type,... with
#                 type i = int32 and d = float64), written last, and
#                 <prefix>.<table>.bin with the columns one after another

# Full path of the program in the folder coancestry.path; stops with the
# download instructions when it is not there.
utils.coancestry.exe <- function(coancestry.path) {
  exe <- if (.Platform$OS.type == "windows") {
    "dartR_coancestry.exe"
  } else {
    "dartR_coancestry"
  }
  ok <- is.character(coancestry.path) && length(coancestry.path) == 1 &&
    !is.na(coancestry.path)
  if (!ok || !file.exists(file.path(coancestry.path, exe))) {
    stop(error(paste0(
      "  Cannot find the program '", exe, "' in the folder given by ",
      "coancestry.path", if (ok) paste0(": ", coancestry.path), "\n",
      "  gl.relatedness computes the Coancestry estimators with this program,\n",
      "  which is closed-source and downloaded separately. Download it with:\n",
      "    dartRverse::gl.download.binary(software = 'coancestry', ",
      "out.dir = tempdir())\n",
      "  then use coancestry.path = file.path(tempdir(), 'coancestry')\n")))
  }
  normalizePath(file.path(coancestry.path, exe))
}

# Writes the genotype file. The matrix is written in blocks of loci, so the
# extra memory is a few copies of one block and not of the whole matrix.
utils.coancestry.write <- function(snp, file, block.loci = NULL) {
  if (is.null(block.loci)) {
    block.loci <- max(1L, 1e7 %/% nrow(snp))
  }
  con <- file(file, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(charToRaw("DRCG0001"), con)
  writeBin(dim(snp), con, size = 4L, endian = "little")
  for (from in seq(1L, ncol(snp), by = block.loci)) {
    g <- snp[, from:min(ncol(snp), from + block.loci - 1L), drop = FALSE]
    g[is.na(g)] <- 3L
    writeBin(as.raw(g), con)
  }
  invisible(file)
}

# Reads the result tables named in <prefix>.manifest. Anything the reader does
# not fully understand is an error, never a guess: a manifest from another
# program, a column type other than i or d, or a table file whose size does not
# match its rows and columns.
utils.coancestry.read <- function(prefix) {
  man <- readLines(paste0(prefix, ".manifest"))
  if (!length(man) || !startsWith(man[1], "dartR_coancestry ")) {
    stop(error("  The results are not a dartR_coancestry manifest.\n"))
  }
  out <- list()
  for (ln in man[-1]) {
    f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    n <- suppressWarnings(as.integer(f[2]))
    cols <- strsplit(strsplit(f[3], ",", fixed = TRUE)[[1]], ":", fixed = TRUE)
    types <- vapply(cols, `[`, "", 2L)
    if (length(f) != 3 || is.na(n) || anyNA(types) ||
        !all(types %in% c("i", "d"))) {
      stop(error(paste0("  Unknown column type or malformed line in the ",
                        "results of dartR_coancestry: ", ln, "\n")))
    }
    bin <- paste0(prefix, ".", f[1], ".bin")
    want <- n * sum(ifelse(types == "i", 4, 8))
    if (!file.exists(bin) || file.size(bin) != want) {
      stop(error(paste0("  The results of dartR_coancestry are incomplete: ",
                        "table '", f[1], "' has the wrong size.\n")))
    }
    con <- file(bin, "rb")
    df <- list()
    for (cl in cols) {
      df[[cl[1]]] <- if (cl[2] == "i") {
        readBin(con, "integer", n, size = 4L, endian = "little")
      } else {
        readBin(con, "double", n, size = 8L, endian = "little")
      }
    }
    close(con)
    out[[f[1]]] <- as.data.frame(df, stringsAsFactors = FALSE)
  }
  out
}

# Runs the program on an individual x locus dosage matrix (0/1/2, NA missing)
# for all pairs and returns its tables as a list of data frames: $dyads,
# $delta19, $trio_delta and $inbreeding, whichever the estimators produce, with
# individuals as 1-based row numbers of snp. The files live in a temporary
# folder that is removed on exit; the program is run from that folder with
# relative file names, so spaces or accents in the session's paths cannot
# break the command line.
# R waits for the program in steps of 0.2 s rather than in one system2() call,
# because R does not see an interrupt (RStudio's Stop button, Esc, Ctrl-C)
# while system2() waits, and the program would run on to the end. On an
# interrupt, an error or a normal end, the program is ended before its folder
# is removed.
utils.coancestry.run <- function(snp,
                                 estimators,
                                 n_threads = 1L,
                                 n_bootstrap = 0L,
                                 rng_seed = 42L,
                                 allow_inbreeding = FALSE,
                                 num_trios = 100L,
                                 exe) {
  if (any(snp < 0L | snp > 2L, na.rm = TRUE)) {
    stop(error("  Genotypes must be 0, 1, 2 or NA.\n"))
  }
  run.dir <- tempfile("coancestry_")
  dir.create(run.dir)
  old.wd <- setwd(run.dir)
  on.exit({
    setwd(old.wd)
    unlink(run.dir, recursive = TRUE)
  }, add = TRUE)

  utils.coancestry.write(snp, "geno.bin")

  int <- function(v) sprintf("%d", as.integer(v))
  args <- c("--geno", "geno.bin", "--out", "res",
            "--estimators", paste(estimators, collapse = ","),
            "--threads", int(n_threads), "--boots", int(n_bootstrap),
            "--seed", int(rng_seed), "--trios", int(num_trios),
            if (allow_inbreeding) "--inbreeding")
  proc <- NULL
  on.exit(if (!is.null(proc) && proc$is_alive()) proc$kill(),
          add = TRUE, after = FALSE)
  proc <- processx::process$new(exe, args, stdout = "log.txt", stderr = "2>&1")
  while (proc$is_alive()) {
    proc$wait(200)
  }
  status <- proc$get_exit_status()
  if (!identical(status, 0L) || !file.exists("res.manifest")) {
    log <- if (file.exists("log.txt")) readLines("log.txt", warn = FALSE)
    stop(error(paste0("  The program dartR_coancestry failed (status ",
                      status, "):\n    ",
                      paste(log, collapse = "\n    "), "\n")))
  }

  utils.coancestry.read("res")
}
