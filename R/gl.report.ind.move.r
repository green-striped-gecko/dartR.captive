#' @name gl.report.ind.move
#' @title Reports the effect on gene diversity of moving individuals between
#' populations
#' @family captive management
#'
#' @description
#' Evaluates every possible transfer of an individual from its population to
#' each other population, reporting the change in gene diversity of the source
#' population, of the destination population, and their sum, based on an
#' empirical genomic kinship matrix.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param kin Kinship matrix as returned by \code{gl.kin}, with row and column
#' names identical to \code{indNames(x)}; if NULL, computed internally with
#' \code{gl.kin} [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @details
#' This function is the genomic analogue of the transfer (translocation)
#' analysis of the pedigree-management program PMx (Lacy, Ballou & Pollak
#' 2012), which evaluates moving animals between institutions or management
#' sets and reports the genetic effect on each set. Here the populations of
#' the genlight object play the role of the management sets.
#'
#' Gene diversity of a population block is GD = 1 - mean(kin[block, block]),
#' the mean taken over the full sub-matrix including the diagonal. For each
#' individual i (source population s) and each destination population t
#' (t != s):
#' \itemize{
#' \item dgd.source = GD(s-block without i) - GD(s-block) -- positive when the
#' source population's gene diversity RISES on losing i (i is overrepresented
#' at home);
#' \item dgd.dest = GD(t-block with i's row and column included) - GD(t-block)
#' -- positive when the destination gains diversity by receiving i;
#' \item net = dgd.source + dgd.dest -- the summed effect of the move.
#' }
#' Because the kinship matrix is empirical, the kinships of the moved
#' individual to the destination members are real measured values taken
#' directly from \code{kin} -- no pedigree assumption that a transferred
#' animal is unrelated to its new population is needed.
#'
#' Individuals whose source population contains only themselves receive
#' dgd.source = NA (removing the sole member leaves no population to score);
#' a warning lists the populations concerned. Moves are reported sorted by
#' net change, best first (NA net values last).
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # SNP data: captive colony plus two wild populations
#' ex <- gl.keep.pop(testset.gl,
#'   pop.list = c("EmmacCaptBred", "EmmacMaclGeor", "EmmacBurnBara"),
#'   verbose = 0)
#' res <- gl.report.ind.move(ex)
#' head(res)
#' # Tag P/A data (SilicoDArT; kinship computed internally)
#' ex.gs <- gl.keep.pop(testset.gs,
#'   pop.list = c("EmmacCaptBred", "EmmacMaclGeor", "EmmacBurnBara"),
#'   verbose = 0)
#' res.gs <- gl.report.ind.move(ex.gs)
#' head(res.gs)
#'
#' @seealso \code{\link{gl.kin}}, \code{\link{gl.report.ind.remove}},
#' \code{\link{gl.report.ind.add}}
#'
#' @export
#' @return Invisibly, a data.frame with one row per candidate move and columns
#' id, from, to, dgd.source, dgd.dest, net, sorted by net descending.

gl.report.ind.move <- function(x,
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
  if (nPop(x) < 2) {
    stop(error("Fatal Error: at least two populations are required to assess moves\n"))
  }
  kin <- utils.kin.check(x, kin, verbose = verbose)

  pops <- popNames(x)
  blocks <- split(indNames(x), pop(x))
  blocks <- blocks[pops]

  singletons <- pops[vapply(blocks, length, integer(1)) == 1]
  if (length(singletons) > 0 && verbose >= 1) {
    cat(warn(
      "  Warning: population(s) with a single individual: ",
      paste(singletons, collapse = ", "),
      "; dgd.source is NA for moves out of these populations\n"
    ))
  }

  # DO THE JOB ----------------------
  if (verbose >= 2) {
    cat(report("  Evaluating", nInd(x) * (length(pops) - 1), "candidate moves\n"))
  }

  # Gene diversity of each population block, computed once
  gd.block <- vapply(blocks, function(ids) {
    utils.kin.dgd(kin[ids, ids, drop = FALSE])
  }, numeric(1))

  out <- vector("list", length(pops))
  names(out) <- pops
  for (s in pops) {
    ids.s <- blocks[[s]]
    dests <- setdiff(pops, s)
    rows <- vector("list", length(ids.s))
    for (k in seq_along(ids.s)) {
      id <- ids.s[k]
      # Effect on the source block of removing the individual
      if (length(ids.s) > 1) {
        dgd.source <- utils.kin.dgd(kin[ids.s, ids.s, drop = FALSE],
                                    drop = id) - gd.block[s]
      } else {
        dgd.source <- NA_real_
      }
      # Effect on each destination block of adding the individual's real
      # row/column from the kinship matrix
      dgd.dest <- vapply(dests, function(t) {
        joint <- c(blocks[[t]], id)
        utils.kin.dgd(kin[joint, joint, drop = FALSE]) - gd.block[t]
      }, numeric(1))
      rows[[k]] <- data.frame(
        id = id,
        from = s,
        to = dests,
        dgd.source = unname(dgd.source),
        dgd.dest = unname(dgd.dest),
        stringsAsFactors = FALSE
      )
    }
    out[[s]] <- do.call(rbind, rows)
  }
  res <- do.call(rbind, out)
  res$net <- res$dgd.source + res$dgd.dest
  res <- res[order(res$net, decreasing = TRUE, na.last = TRUE), ]
  rownames(res) <- NULL

  # Print out the results summary ---------------
  if (verbose >= 3) {
    cat("  Top candidate moves by net change in gene diversity:\n")
    tmp <- head(res, 10)
    tmp[, c("dgd.source", "dgd.dest", "net")] <-
      round(tmp[, c("dgd.source", "dgd.dest", "net")], 4)
    print(tmp, row.names = FALSE)
  }

  # FLAG SCRIPT END ---------------
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN
  return(invisible(res))
}
