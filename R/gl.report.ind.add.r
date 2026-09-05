#' @name gl.report.ind.add
#' @title Reports the gene diversity gained by adding candidate individuals to
#' a target population
#' @family captive management
#'
#' @description
#' Ranks candidate individuals by the change in gene diversity of a target
#' (managed) population that their addition would produce, based on an
#' empirical genomic kinship matrix.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param candidates Character vector of individual names in \code{x} to
#' evaluate as additions; or a single population name in \code{popNames(x)},
#' meaning all its members [required].
#' @param target.pop Name of the managed target population, one of
#' \code{popNames(x)} [required].
#' @param kin Kinship matrix as returned by \code{gl.kin}, with row and column
#' names identical to \code{indNames(x)}; if NULL, computed internally with
#' \code{gl.kin} [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @details
#' This function is the genomic analogue of evaluating the value of potential
#' new founders in the pedigree-management program PMx (Lacy, Ballou & Pollak
#' 2012). PMx must assume that a new founder is unrelated to the living
#' population; here the candidate's kinships to every member of the target
#' population are measured values taken directly from the kinship matrix, so
#' the reported gain is exact and candidates that are cryptic relatives of
#' the existing stock are ranked appropriately lower. It is a natural
#' companion to the \code{gl.assign} suite: having assigned a stray or
#' wild-caught individual to a source population, ask what it would
#' contribute to the managed colony.
#'
#' Gene diversity of the target block is GD = 1 - mean(kin[block, block]),
#' the mean taken over the full sub-matrix including the diagonal. For each
#' candidate c the function computes
#' dgd = GD(target block plus c's row and column from kin) - GD(target block).
#' A positive dgd means adding c raises the target population's gene
#' diversity.
#'
#' If \code{candidates} is a single name that matches a population in
#' \code{popNames(x)}, it is expanded to all members of that population (the
#' population-name interpretation takes precedence over an individual of the
#' same name). Otherwise every element must match an individual name in
#' \code{indNames(x)}. Candidates already belonging to \code{target.pop} are
#' rejected with a fatal error.
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # SNP data
#' kin <- gl.kin(testset.gl)
#' res <- gl.report.ind.add(testset.gl, candidates = "EmmacMaclGeor",
#'                          target.pop = "EmmacCaptBred", kin = kin)
#' head(res)  # wild individuals ranked by the diversity they would add
#' # Tag P/A data (SilicoDArT; kinship computed internally)
#' res.gs <- gl.report.ind.add(testset.gs, candidates = "EmmacMaclGeor",
#'                             target.pop = "EmmacCaptBred")
#' head(res.gs)
#'
#' @seealso \code{\link{gl.kin}}, \code{\link{gl.report.ind.remove}},
#' \code{\link{gl.report.ind.move}}
#'
#' @export
#' @return Invisibly, a data.frame with one row per candidate and columns id,
#' from (the candidate's population), dgd, rank, sorted by dgd descending.

gl.report.ind.add <- function(x,
                              candidates,
                              target.pop,
                              kin = NULL,
                              verbose = NULL) {
  # PRELIMINARIES -- checking ----------------
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname,
                   build = "v.2026.1",
                   verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (missing(candidates) || missing(target.pop)) {
    stop(error("Fatal Error: both candidates and target.pop must be specified\n"))
  }
  if (!is.character(target.pop) || length(target.pop) != 1 ||
      !(target.pop %in% popNames(x))) {
    stop(error("Fatal Error: target.pop must be a single population name in popNames(x)\n"))
  }
  if (!is.character(candidates) || length(candidates) < 1) {
    stop(error("Fatal Error: candidates must be a character vector of individual names, or a single population name\n"))
  }
  if (length(candidates) == 1 && candidates %in% popNames(x)) {
    if (verbose >= 2) {
      cat(report("  Expanding candidates to all members of population",
                 candidates, "\n"))
    }
    candidates <- indNames(x)[pop(x) == candidates]
  } else {
    missing.ids <- setdiff(candidates, indNames(x))
    if (length(missing.ids) > 0) {
      stop(error(
        "Fatal Error: candidate(s) not found in indNames(x): ",
        paste(missing.ids, collapse = ", "), "\n"
      ))
    }
  }
  cand.pop <- as.character(pop(x))[match(candidates, indNames(x))]
  inside <- candidates[cand.pop == target.pop]
  if (length(inside) > 0) {
    stop(error(
      "Fatal Error: candidate(s) already belong to target.pop: ",
      paste(inside, collapse = ", "), "\n"
    ))
  }
  kin <- utils.kin.check(x, kin, verbose = verbose)

  # DO THE JOB ----------------------
  if (verbose >= 2) {
    cat(report("  Evaluating", length(candidates),
               "candidate additions to", target.pop, "\n"))
  }

  t.ids <- indNames(x)[pop(x) == target.pop]
  gd.target <- utils.kin.dgd(kin[t.ids, t.ids, drop = FALSE])

  dgd <- vapply(candidates, function(id) {
    joint <- c(t.ids, id)
    utils.kin.dgd(kin[joint, joint, drop = FALSE]) - gd.target
  }, numeric(1))

  res <- data.frame(
    id = candidates,
    from = cand.pop,
    dgd = dgd,
    stringsAsFactors = FALSE
  )
  res <- res[order(res$dgd, decreasing = TRUE), ]
  res$rank <- seq_len(nrow(res))
  rownames(res) <- NULL

  # Print out the results summary ---------------
  if (verbose >= 3) {
    cat("  Gene diversity of", target.pop, ":", round(gd.target, 4), "\n")
    cat("  Candidates ranked by gene diversity added:\n")
    tmp <- res
    tmp$dgd <- round(tmp$dgd, 4)
    print(tmp, row.names = FALSE)
  }

  # FLAG SCRIPT END ---------------
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN
  return(invisible(res))
}
