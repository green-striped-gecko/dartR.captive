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
#' @param plot.out If TRUE, display a heatmap of \code{plot.stat} [default TRUE].
#' @param plot.colors Divergent palette; NULL uses \code{gl.colors("div")}
#'   [default NULL].
#' @param plot.dir Directory to save the plot RDS [default working dir or tempdir].
#' @param plot.file Filename (no extension) for the saved RDS plot [default NULL].
#' @param verbose Verbosity 0..5 [default NULL, uses gl.set.verbosity].
#'
#' @details
#' \code{as.matrix(x)} supplies the individual x locus alternate-allele dosage
#' matrix (0/1/2, NA missing) the engine expects; allele frequencies are computed
#' internally. When \code{n.boots > 0}, monomorphic / all-missing loci and
#' all-missing individuals are dropped first (with a warning): the bootstrap loop,
#' faithful to Coancestry, does not terminate on a pair with no usable polymorphic
#' locus.
#'
#' The engine \code{dartR.coancestry} is closed-source, binary-only, not on CRAN.
#' Install from its GitHub release, e.g.:
#' \preformatted{install.packages(
#'   "https://github.com/mijangos81/dartR.coancestry/releases/download/v0.1.0/dartR.coancestry_0.1.0.tgz",
#'   repos = NULL, type = "binary")}
#' (the repository is private; you need access).
#'
#' @return A named list: one symmetric individual x individual matrix per requested
#' dyad estimator (e.g. \code{$wang}; rownames/colnames = \code{indNames(x)},
#' diagonal NA); a long per-pair data frame \code{$dyads} (with \code{_lo}/\code{_hi}
#' when \code{n.boots > 0}); and, when requested, \code{$delta19} (DyadML deltas),
#' \code{$trio_delta} (TrioML deltas) and \code{$inbreeding} (per-individual LH/LR
#' and/or ML inbreeding).
#'
#' @author Custodian: Luis Mijangos -- Post to
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
  utils.flag.start(func = funname, build = "v.2023.2", verbose = verbose)
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
  if (!requireNamespace("dartR.coancestry", quietly = TRUE))
    stop(error(
      "  gl.relatedness needs the 'dartR.coancestry' engine (closed-source, binary-only).\n",
      "  Install from its GitHub release, e.g.:\n",
      "    install.packages('https://github.com/mijangos81/dartR.coancestry/releases/download/v0.1.0/dartR.coancestry_0.1.0.tgz', repos = NULL, type = 'binary')\n",
      "  (the repository is private; you need access).\n"))

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
    if (d.loc > 0 || d.ind > 0)
      warning(sprintf(
        "Bootstrap pre-filter: dropped %d monomorphic/all-NA loci and %d all-NA individuals (prevents the rejection-loop hang).",
        d.loc, d.ind))
  }

  # PREPARE DOSAGE + RUN ENGINE --------------------------------------
  snp <- as.matrix(x)
  storage.mode(snp) <- "integer"
  nd <- choose(nInd(x), 2)
  if (nd > .Machine$integer.max)
    stop(error("  Too many individuals (", nInd(x), "): ", nd,
               " dyads exceeds the integer limit for a single run.\n"))
  max.dyads <- as.integer(nd)
  res <- dartR.coancestry::relatedness_cpp(
    snp, estimators = estimators, max_dyads = max.dyads,
    n_threads = as.integer(n.threads), n_bootstrap = as.integer(n.boots),
    rng_seed = as.integer(rng.seed), allow_inbreeding = allow.inbreeding,
    num_trios = as.integer(num.trios))

  # POST-PROCESS: map 1-based indices to names; build symmetric matrices ----
  nm <- indNames(x)
  dyad.est <- intersect(estimators,
    c("wang", "lynchli", "lynchrd", "ritland", "quellergt", "loiselle",
      "dyadml", "trioml"))
  out <- list()
  if (!is.null(res$dyads)) {
    i1 <- res$dyads$ind1; i2 <- res$dyads$ind2
    for (e in dyad.est) {
      M <- matrix(NA_real_, nInd(x), nInd(x), dimnames = list(nm, nm))
      v <- res$dyads[[e]]
      M[cbind(i1, i2)] <- v
      M[cbind(i2, i1)] <- v
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

  # PLOT: one heatmap of the chosen estimator (default the first requested) ----
  if (plot.out && length(dyad.est) > 0) {
    stat <- if (is.null(plot.stat)) dyad.est[1] else plot.stat
    if (!stat %in% dyad.est)
      stop(error("  plot.stat '", stat,
                 "' is not among the requested dyad estimators (",
                 paste(dyad.est, collapse = ", "), ").\n"))
    pal <- if (is.null(plot.colors)) gl.colors("div") else plot.colors
    p <- gl.plot.heatmap(out[[stat]], palette.divergent = pal,
                         plot.out = TRUE, verbose = 0)
    if (!is.null(plot.file))
      utils.plot.save(p, dir = plot.dir, file = plot.file, verbose = verbose)
  }

  # FLAG SCRIPT END
  if (verbose >= 1) cat(report("Completed:", funname, "\n"))
  return(out)
}
