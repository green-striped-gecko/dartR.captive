#' @name gl.report.kin.sets
#' @title Reports kinship-based management-set statistics by population
#' @family captive management
#'
#' @description
#' Summarises each population (management set) of a genlight object by size,
#' mean kinship, within-set gene diversity and mean inbreeding, and reports
#' pairwise between-set mean kinships and a kinship-based Fst, using an
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
#' This function is the genomic analogue of the Management Sets tab (MetaMK
#' module) of the pedigree-management program PMx (Lacy, Ballou & Pollak 2012),
#' with the populations of the genlight object taking the role of PMx's
#' management sets and the empirical kinship matrix taking the role of
#' pedigree kinship.
#'
#' For each population s the function reports:
#' \itemize{
#' \item n -- number of individuals;
#' \item meanMK -- mean over members of MK_i = rowMeans(kin), the mean kinship
#' of individual i to the whole living population (all populations pooled);
#' \item GD.w -- within-set gene diversity, 1 - mean(kin[s,s]) taken over the
#' full within-set block of the kinship matrix including the diagonal;
#' \item meanF -- mean inbreeding, mean(2*k_ii - 1) over the diagonal entries
#' of the set, since the diagonal of the kinship matrix is 0.5*(1+F).
#' }
#'
#' For each pair of populations s and t the function reports:
#' \itemize{
#' \item MKb[s,t] -- the between-set mean kinship, mean(kin[s-block, t-block]),
#' the mean over the full cross-block (every member of s against every member
#' of t). The diagonal of the MKb matrix holds the within-set mean kinship
#' including the diagonal (equal to 1 - GD.w).
#' \item Fst[s,t] -- the kinship-based analogue of PMx's Fst = 1 - GDt/GDb
#' (Management Sets tab), computed exactly as
#' Fst[s,t] = 1 - GDt.pair/GDb.pair, where GDb.pair = 1 - MKb[s,t] and
#' GDt.pair = 1 - mean(kin) over the joint block of the two populations
#' (both within-blocks and both cross-blocks, diagonal included). The
#' diagonal of the Fst matrix is 0.
#' }
#'
#' Populations of n = 1 are reported but flagged with a warning: their
#' within-set gene diversity is simply 1 minus the individual's
#' self-kinship, and their meanF is that individual's inbreeding coefficient.
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # SNP data
#' kin <- gl.kin(testset.gl)
#' res <- gl.report.kin.sets(testset.gl, kin = kin)
#' res$sets[res$sets$pop == "EmmacCaptBred", ]
#' res$fst[1:4, 1:4]  # captive pop shows elevated within-set kinship
#' # Tag P/A data (SilicoDArT; kinship computed internally)
#' res.gs <- gl.report.kin.sets(testset.gs)
#'
#' @seealso \code{\link{gl.kin}}, \code{\link{gl.report.ind.move}},
#' \code{\link{gl.report.ind.add}}
#'
#' @export
#' @return Invisibly, a list with three components:
#' \itemize{
#' \item sets -- data.frame with columns pop, n, meanMK, GD.w, meanF;
#' \item mkb -- matrix of between-set mean kinships (within-set mean kinship
#' on the diagonal);
#' \item fst -- matrix of pairwise kinship-based Fst (0 on the diagonal).
#' }

gl.report.kin.sets <- function(x,
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
  kin <- utils.kin.check(x, kin, verbose = verbose)

  pops <- popNames(x)
  blocks <- split(indNames(x), pop(x))
  blocks <- blocks[pops]

  singletons <- pops[vapply(blocks, length, integer(1)) == 1]
  if (length(singletons) > 0 && verbose >= 1) {
    cat(warn(
      "  Warning: population(s) with a single individual: ",
      paste(singletons, collapse = ", "),
      "; their within-set gene diversity is just 1 minus self-kinship\n"
    ))
  }
  if (length(pops) < 2 && verbose >= 1) {
    cat(warn("  Warning: only one population; between-set matrices are trivial\n"))
  }

  # DO THE JOB ----------------------
  if (verbose >= 2) {
    cat(report("  Computing per-set statistics for", length(pops), "populations\n"))
  }

  mk <- rowMeans(kin)
  self <- diag(kin)

  sets <- do.call(rbind, lapply(pops, function(s) {
    ids <- blocks[[s]]
    data.frame(
      pop = s,
      n = length(ids),
      meanMK = mean(mk[ids]),
      GD.w = 1 - mean(kin[ids, ids, drop = FALSE]),
      meanF = mean(2 * self[ids] - 1),
      stringsAsFactors = FALSE
    )
  }))
  rownames(sets) <- NULL

  if (verbose >= 2) {
    cat(report("  Computing between-set mean kinship and kinship-based Fst\n"))
  }

  npop <- length(pops)
  mkb <- matrix(NA_real_, nrow = npop, ncol = npop, dimnames = list(pops, pops))
  fst <- matrix(0, nrow = npop, ncol = npop, dimnames = list(pops, pops))

  for (s in seq_len(npop)) {
    ids.s <- blocks[[s]]
    mkb[s, s] <- mean(kin[ids.s, ids.s, drop = FALSE])
    if (s < npop) {
      for (t in (s + 1):npop) {
        ids.t <- blocks[[t]]
        # Between-set mean kinship: full cross-block, every s member x every t member
        mkb.st <- mean(kin[ids.s, ids.t, drop = FALSE])
        mkb[s, t] <- mkb.st
        mkb[t, s] <- mkb.st
        # Kinship-based Fst = 1 - GDt.pair/GDb.pair (PMx: Fst = 1 - GDt/GDb)
        joint <- c(ids.s, ids.t)
        gdt.pair <- 1 - mean(kin[joint, joint, drop = FALSE])
        gdb.pair <- 1 - mkb.st
        fst.st <- 1 - gdt.pair / gdb.pair
        fst[s, t] <- fst.st
        fst[t, s] <- fst.st
      }
    }
  }

  # Print out the results summary ---------------
  if (verbose >= 3) {
    cat("  Per-population management-set statistics:\n")
    tmp <- sets
    tmp[, c("meanMK", "GD.w", "meanF")] <-
      round(tmp[, c("meanMK", "GD.w", "meanF")], 4)
    print(tmp, row.names = FALSE)
    cat("  Between-set mean kinship (mkb) and kinship-based Fst (fst)\n")
    cat("  returned as", npop, "x", npop, "matrices\n")
  }

  # FLAG SCRIPT END ---------------
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN
  return(invisible(list(sets = sets, mkb = mkb, fst = fst)))
}
