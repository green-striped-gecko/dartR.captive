#' @name gl.relatedness
#' @title Pairwise relatedness and inbreeding from SNP data (Coancestry estimators)
#' @family relatedness
#'
#' @description
#' Estimates pairwise relatedness (and optionally inbreeding) between individuals
#' in a genlight object using Jinliang Wang's Coancestry estimators, translated to
#' C++/Rcpp. Method-of-moments estimators (Wang, Lynch & Li, Lynch & Ritland,
#' Ritland, Queller & Goodnight, Loiselle) and maximum-likelihood estimators
#' (DyadML, TrioML), plus per-individual inbreeding, are available.
#'
#' The numerical engine is the closed-source, binary-only package
#' \code{dartR.coancestry}, which must be installed separately (see Details).
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param estimators Character vector of estimators; any subset of \code{"wang"},
#'   \code{"lynchli"}, \code{"lynchrd"}, \code{"ritland"}, \code{"quellergt"},
#'   \code{"loiselle"}, \code{"inbreeding"}, \code{"dyadml"}, \code{"trioml"}
#'   [default the six moment estimators].
#' @param allow.inbreeding If TRUE, fit the inbreeding-aware ML model for
#'   \code{dyadml}/\code{trioml} [default FALSE].
#' @param n.boots Bootstrap-over-loci replicates for 95\% CIs; 0 disables CIs
#'   [default 0].
#' @param num.trios Reference trios for \code{trioml} [default 100].
#' @param rng.seed Random seed (matters for \code{dyadml}/\code{trioml} and the
#'   bootstrap) [default 42].
#' @param n.threads Threads for the moment tier (the ML tier is serial) [default 1].
#' @param plot.stat Estimator whose matrix is drawn as a heatmap; NULL = the first
#'   requested dyad estimator [default NULL].
#' @param plot.out If TRUE, display a heatmap of \code{plot.stat}; the heatmap
#'   is saved when \code{plot.file} is given, whether or not it is displayed
#'   [default TRUE].
#' @param plot.colors Divergent palette; NULL uses \code{gl.colors("div")}
#'   [default NULL].
#' @param plot.dir Directory to save the plot RDS [default working dir or tempdir].
#' @param plot.file Filename (no extension) for the saved RDS plot [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'   brief progress messages; 3, progress and results summary; 5, full report
#'   [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' \code{as.matrix(x)} supplies the individual x locus alternate-allele dosage
#' matrix (0/1/2, NA missing) the engine expects; allele frequencies are computed
#' internally. When \code{n.boots > 0}, monomorphic / all-missing loci and
#' all-missing individuals are dropped first (reported at \code{verbose >= 1}):
#' the bootstrap loop, faithful to Coancestry, does not terminate on a pair with
#' no usable polymorphic locus.
#'
#' The engine \code{dartR.coancestry} is closed-source, binary-only, not on CRAN.
#' Download the build matching your operating system, CPU architecture and R
#' version from
#' \url{https://github.com/mijangos81/dartR.coancestry-binaries/releases}, then
#' install it with:
#' \preformatted{install.packages("<path-to-downloaded-file>",
#'   repos = NULL, type = "binary")}
#' Binary R packages are tied to the R minor version they were built against, so
#' a build for R 4.4 will not load under R 4.5.
#'
#' A pair of individuals that shares no called locus has no information on
#' relatedness; every estimate for it is set to NA (the engine returns 0 for the
#' moment estimators and NaN for the likelihood estimators), with a warning at
#' \code{verbose >= 1}. Individuals dropped by the bootstrap pre-filter keep
#' their rows and columns in the matrices, filled with NA.
#'
#' The matrices are on the relatedness scale (twice kinship) and carry
#' \code{attr(, "scale") = "relatedness"}, so the kinship functions of this
#' package (for example \code{gl.report.kin.classes}) halve them when they are
#' passed as \code{kin}.
#'
#' Agreement with Coancestry: on 30 platypus (platypus.gl, 420 loci with no
#' missing data), \code{lynchrd}, \code{ritland} and \code{quellergt} match
#' \code{related::coancestry} (the Coancestry Fortran code) within its
#' 4-decimal rounding. \code{wang} is on average 0.047 higher and
#' \code{lynchli} 0.018 higher. In simulated biallelic SNP data the gap for
#' unrelated pairs shrinks with sample size (0.050 at 20 individuals, 0.023 at
#' 50, 0.013 at 120), which points to a different small-sample correction. In
#' the same simulations neither implementation of \code{wang} or
#' \code{lynchli} recovers r = 0 for unrelated pairs (both about -0.16), while
#' \code{quellergt} does (-0.009 at 120 individuals); for full sibs all give
#' 0.57-0.58.
#'
#' @return Invisibly, a named list: one symmetric individual x individual matrix
#' per requested dyad estimator (e.g. \code{$wang}; rownames/colnames =
#' \code{indNames(x)}, diagonal NA, \code{attr(, "scale") = "relatedness"}); a
#' long per-pair data frame \code{$dyads} (with \code{_lo}/\code{_hi}
#' when \code{n.boots > 0}); and, when requested, \code{$delta19} (DyadML deltas),
#' \code{$trio_delta} (TrioML deltas) and \code{$inbreeding} (per-individual LH/LR
#' and/or ML inbreeding).
#'
#' @author Author(s): Luis Mijangos. Custodian: Luis Mijangos -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \dontrun{
#' gl  <- dartR.data::platypus.gl[1:40, 1:1000]
#' res <- gl.relatedness(gl, estimators = c("wang", "lynchrd"))
#' res$wang[1:5, 1:5]
#' head(res$dyads)
#' }
#'
#' @references
#' Wang, J. (2002). An estimator for pairwise relatedness using molecular markers.
#' Genetics, 160, 1203-1215.
#'
#' @export

gl.relatedness <- function(x,
                           estimators = c("wang", "lynchli", "lynchrd",
                                          "ritland", "quellergt", "loiselle"),
                           allow.inbreeding = FALSE,
                           n.boots = 0,
                           num.trios = 100,
                           rng.seed = 42,
                           n.threads = 1,
                           plot.stat = NULL,
                           plot.out = TRUE,
                           plot.colors = NULL,
                           plot.dir = NULL,
                           plot.file = NULL,
                           verbose = NULL) {

  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # as() adds the fbm slot; class<- only relabels and gives an invalid object
  if (!is(x, "dartR")) x <- methods::as(x, "dartR")

  # FUNCTION-SPECIFIC ERROR CHECKING ---------------------------------
  if (datatype != "SNP")
    stop(error("  gl.relatedness requires SNP data (not SilicoDArT).\n"))
  if (any(ploidy(x) != 2))
    stop(error("  gl.relatedness requires diploid SNP data (ploidy 2).\n"))
  valid.est <- c("wang", "lynchli", "lynchrd", "ritland", "quellergt",
                 "loiselle", "inbreeding", "dyadml", "trioml")
  bad <- setdiff(estimators, valid.est)
  if (length(bad))
    stop(error("  Unknown estimator(s): ", paste(bad, collapse = ", "),
               ". Valid: ", paste(valid.est, collapse = ", "), "\n"))
  dyad.est <- intersect(estimators,
    c("wang", "lynchli", "lynchrd", "ritland", "quellergt", "loiselle",
      "dyadml", "trioml"))
  # Scalar arguments are checked before the engine runs, so a long run is
  # never discarded on a typing error
  is.num1 <- function(v) is.numeric(v) && length(v) == 1 && !is.na(v)
  is.flag <- function(v) is.logical(v) && length(v) == 1 && !is.na(v)
  if (!is.num1(n.boots) || n.boots < 0)
    stop(error("  n.boots must be a single number >= 0.\n"))
  if (!is.num1(num.trios) || num.trios < 1)
    stop(error("  num.trios must be a single number >= 1.\n"))
  if (!is.num1(n.threads) || n.threads < 1)
    stop(error("  n.threads must be a single number >= 1.\n"))
  if (!is.num1(rng.seed))
    stop(error("  rng.seed must be a single number.\n"))
  if (!is.flag(allow.inbreeding))
    stop(error("  allow.inbreeding must be TRUE or FALSE.\n"))
  if (!is.flag(plot.out))
    stop(error("  plot.out must be TRUE or FALSE.\n"))
  if (!is.null(plot.stat) &&
      !(is.character(plot.stat) && length(plot.stat) == 1 &&
        plot.stat %in% dyad.est))
    stop(error("  plot.stat must be one of the requested dyad estimators (",
               paste(dyad.est, collapse = ", "), ").\n"))
  # The engine is not a declared dependency (closed-source, not in a CRAN-like
  # repository); naming it through a variable keeps R CMD check clean
  engine.pkg <- "dartR.coancestry"
  if (!requireNamespace(engine.pkg, quietly = TRUE))
    stop(error(
      "  gl.relatedness needs the 'dartR.coancestry' engine (closed-source, binary-only).\n",
      "  Download the build matching your operating system, CPU architecture and\n",
      "  R version from:\n",
      "    https://github.com/mijangos81/dartR.coancestry-binaries/releases\n",
      "  then install it with:\n",
      "    install.packages('<path-to-downloaded-file>', repos = NULL, type = 'binary')\n"))

  relatedness.cpp <- getExportedValue(engine.pkg, "relatedness_cpp")
  nm.all <- indNames(x)

  # BOOTSTRAP PRE-FILTER: the engine's resample loop (faithful to trior11) does not
  # terminate on a pair with no usable polymorphic locus, so drop monomorphic /
  # all-NA loci and all-NA individuals before any bootstrap.
  if (n.boots > 0) {
    n0.loc <- nLoc(x); n0.ind <- nInd(x)
    x <- gl.filter.monomorphs(x, verbose = 0)
    am <- as.matrix(x)
    keep.ind <- rowSums(!is.na(am)) > 0
    if (!all(keep.ind)) x <- x[keep.ind, ]
    d.loc <- n0.loc - nLoc(x); d.ind <- n0.ind - nInd(x)
    if ((d.loc > 0 || d.ind > 0) && verbose >= 1)
      cat(warn(sprintf(
        "  Bootstrap pre-filter: dropped %d monomorphic/all-NA loci and %d all-NA individuals (prevents the rejection-loop hang).\n",
        d.loc, d.ind)))
  }

  # PREPARE DOSAGE + RUN ENGINE --------------------------------------
  snp <- as.matrix(x)
  storage.mode(snp) <- "integer"
  nd <- choose(nInd(x), 2)
  if (nd > .Machine$integer.max)
    stop(error("  Too many individuals (", nInd(x), "): ", nd,
               " dyads exceeds the integer limit for a single run.\n"))
  max.dyads <- as.integer(nd)
  res <- relatedness.cpp(
    snp, estimators = estimators, max_dyads = max.dyads,
    n_threads = as.integer(n.threads), n_bootstrap = as.integer(n.boots),
    rng_seed = as.integer(rng.seed), allow_inbreeding = allow.inbreeding,
    num_trios = as.integer(num.trios))

  # PAIRS WITHOUT DATA: a pair sharing no called locus carries no
  # information; the engine returns 0 (moment) or NaN (ML) for it
  nm <- indNames(x)
  called <- !is.na(snp)
  shared <- tcrossprod(called * 1)
  no.data <- function(df) shared[cbind(df$ind1, df$ind2)] == 0
  blank <- function(df, rows) {
    num <- setdiff(names(df)[vapply(df, is.numeric, logical(1))],
                   c("ind1", "ind2"))
    df[rows, num] <- NA
    df
  }
  if (!is.null(res$dyads)) {
    nd.rows <- no.data(res$dyads)
    if (any(nd.rows)) {
      res$dyads <- blank(res$dyads, nd.rows)
      if (verbose >= 1)
        cat(warn("  Warning:", sum(nd.rows), "pairs share no called locus;",
                 "their estimates are set to NA\n"))
    }
  }
  if (!is.null(res$delta19)) res$delta19 <- blank(res$delta19, no.data(res$delta19))
  if (!is.null(res$trio_delta)) res$trio_delta <- blank(res$trio_delta, no.data(res$trio_delta))
  if (!is.null(res$inbreeding)) {
    ib.rows <- rowSums(called)[res$inbreeding$ind] == 0
    num <- setdiff(names(res$inbreeding)[vapply(res$inbreeding, is.numeric,
                                                logical(1))], "ind")
    res$inbreeding[ib.rows, num] <- NA
  }

  # POST-PROCESS: map 1-based indices to names; build symmetric matrices over
  # every individual of the input (pre-filter drops become NA rows) ----
  out <- list()
  if (!is.null(res$dyads)) {
    i1 <- nm[res$dyads$ind1]; i2 <- nm[res$dyads$ind2]
    for (e in dyad.est) {
      M <- matrix(NA_real_, length(nm.all), length(nm.all),
                  dimnames = list(nm.all, nm.all))
      v <- res$dyads[[e]]
      M[cbind(i1, i2)] <- v
      M[cbind(i2, i1)] <- v
      attr(M, "scale") <- "relatedness"
      out[[e]] <- M
    }
    dy <- res$dyads; dy$ind1 <- nm[dy$ind1]; dy$ind2 <- nm[dy$ind2]
    out$dyads <- dy
  }
  if (!is.null(res$delta19)) {
    d19 <- res$delta19; d19$ind1 <- nm[d19$ind1]; d19$ind2 <- nm[d19$ind2]
    out$delta19 <- d19
  }
  if (!is.null(res$trio_delta)) {
    td <- res$trio_delta; td$ind1 <- nm[td$ind1]; td$ind2 <- nm[td$ind2]
    out$trio_delta <- td
  }
  if (!is.null(res$inbreeding)) {
    ib <- res$inbreeding; ib$ind <- nm[ib$ind]
    out$inbreeding <- ib
  }

  # PLOT: one heatmap of the chosen estimator (default the first requested);
  # built when it is displayed or saved. The heatmap is base graphics, so a
  # saved-but-not-displayed plot is drawn on a null device ----
  if ((plot.out || !is.null(plot.file)) && length(dyad.est) > 0) {
    stat <- if (is.null(plot.stat)) dyad.est[1] else plot.stat
    pal <- if (is.null(plot.colors)) gl.colors("div", verbose = 0) else plot.colors
    if (!plot.out) {
      grDevices::pdf(NULL)
      on.exit(grDevices::dev.off(), add = TRUE)
    }
    p <- gl.plot.heatmap(out[[stat]], palette.divergent = pal,
                         plot.out = TRUE, verbose = 0)
    if (!is.null(plot.file))
      utils.plot.save(p, dir = plot.dir, file = plot.file, verbose = verbose)
  }

  # FLAG SCRIPT END
  if (verbose >= 1) cat(report("Completed:", funname, "\n"))
  invisible(out)
}
