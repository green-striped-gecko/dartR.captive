#' @name gl.report.ind.remove
#' @title Reports the effect on gene diversity of removing individuals
#' @family captive management
#'
#' @description
#' Ranks individuals by the change in gene diversity that their removal from
#' the population would produce, based on an empirical genomic kinship matrix,
#' and constructs a greedy set of removals that maximises gene diversity.
#'
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param kin Kinship matrix as returned by \code{gl.kin}, with row and column
#' names identical to \code{indNames(x)}; if NULL, computed internally with
#' \code{gl.kin} [default NULL].
#' @param n.best Maximum number of individuals in the greedy removal set; if
#' NULL, removals continue until no single removal increases gene diversity
#' [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @details
#' Kinship should be estimated on the most inclusive dataset available and the
#' matrix subset to the managed group; computing kinship on a small family
#' group alone inflates estimates (the allele-frequency reference collapses
#' onto the family itself).
#' This function is the genomic analogue of the culling (removal) analysis of
#' the pedigree-management program PMx (Lacy, Ballou & Pollak 2012), in which
#' individuals whose removal would increase gene diversity are flagged as
#' genetically overrepresented -- typically inbred individuals or members of
#' large, well-represented sibships.
#'
#' Gene diversity is GD = 1 - mean(kin), the mean taken over the full kinship
#' matrix including the diagonal. For each individual i the function computes
#' dGD_i = GD(kinship matrix with i's row and column removed) - GD(all), via
#' \code{utils.kin.dgd}. A positive dGD_i means removing i RAISES gene
#' diversity (i is overrepresented); a negative dGD_i means i carries
#' relatively unique genetic material. The ranking table reports each
#' individual's population, mean kinship MK_i = rowMeans(kin), and dGD, sorted
#' with the most expendable individuals first.
#'
#' Because dGD values are not additive (removing one member of a sibship
#' reduces the redundancy of the rest), the function also builds a greedy
#' removal set: at each step the individual whose removal most increases the
#' current gene diversity is removed and all dGD values are recomputed.
#' Removal stops when no candidate yields a gain, or after \code{n.best}
#' removals, whichever comes first (never beyond nInd - 1).
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # SNP data: kinship from the FULL dataset, subset to the captive colony
#' kin <- gl.kin(testset.gl)
#' cb <- gl.keep.pop(testset.gl, pop.list = "EmmacCaptBred", verbose = 0)
#' res <- gl.report.ind.remove(cb, kin = kin[indNames(cb), indNames(cb)])
#' head(res$ranking)  # most-inbred/most-redundant sibs rank first
#' res$removal.set
#' # Tag P/A data (SilicoDArT; kinship computed internally)
#' cb.gs <- gl.keep.pop(testset.gs, pop.list = "EmmacCaptBred", verbose = 0)
#' res.gs <- gl.report.ind.remove(cb.gs, n.best = 3)
#'
#' @seealso \code{\link{gl.kin}}, \code{\link{gl.report.ind.move}},
#' \code{\link{gl.report.ind.add}}
#'
#' @export
#' @return Invisibly, a list with two components:
#' \itemize{
#' \item ranking -- data.frame with columns id, pop, MK, dGD, sorted by dGD
#' descending (most expendable first);
#' \item removal.set -- data.frame with columns step, id, gd.after, the greedy
#' removal sequence and the gene diversity after each removal (zero rows if no
#' removal increases gene diversity).
#' }

gl.report.ind.remove <- function(x,
                                 kin = NULL,
                                 n.best = NULL,
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
  if (nInd(x) < 2) {
    stop(error("Fatal Error: at least two individuals are required to assess removals\n"))
  }
  kin <- utils.kin.check(x, kin, verbose = verbose)
  max.removals <- nInd(x) - 1
  if (!is.null(n.best)) {
    if (!is.numeric(n.best) || length(n.best) != 1 || n.best < 1) {
      if (verbose >= 1) {
        cat(warn("  Warning: n.best must be a single number >= 1; setting to NULL (unlimited)\n"))
      }
    } else if (n.best > nInd(x) - 1) {
      if (verbose >= 1) {
        cat(warn("  Warning: n.best exceeds nInd - 1; clamped to", nInd(x) - 1, "\n"))
      }
    } else {
      max.removals <- as.integer(n.best)
    }
  }

  # DO THE JOB ----------------------
  if (verbose >= 2) {
    cat(report("  Computing per-individual change in gene diversity on removal\n"))
  }

  ids <- indNames(x)
  gd.all <- utils.kin.dgd(kin)
  dgd <- vapply(ids, function(id) {
    utils.kin.dgd(kin, drop = id) - gd.all
  }, numeric(1))

  ranking <- data.frame(
    id = ids,
    pop = as.character(pop(x)),
    MK = rowMeans(kin),
    dGD = dgd,
    stringsAsFactors = FALSE
  )
  ranking <- ranking[order(ranking$dGD, decreasing = TRUE), ]
  rownames(ranking) <- NULL

  # Greedy removal set: recompute gains after each removal
  if (verbose >= 2) {
    cat(report("  Constructing greedy removal set (max", max.removals, "removals)\n"))
  }
  dropped <- character(0)
  remaining <- ids
  gd.current <- gd.all
  set.step <- integer(0)
  set.id <- character(0)
  set.gd <- numeric(0)

  while (length(dropped) < max.removals) {
    gains <- vapply(remaining, function(id) {
      utils.kin.dgd(kin, drop = c(dropped, id))
    }, numeric(1))
    best <- which.max(gains)
    if (gains[best] <= gd.current) {
      break
    }
    dropped <- c(dropped, remaining[best])
    gd.current <- gains[best]
    set.step <- c(set.step, length(dropped))
    set.id <- c(set.id, remaining[best])
    set.gd <- c(set.gd, gd.current)
    remaining <- remaining[-best]
  }

  removal.set <- data.frame(
    step = set.step,
    id = set.id,
    gd.after = set.gd,
    stringsAsFactors = FALSE
  )

  # Print out the results summary ---------------
  if (verbose >= 3) {
    cat("  Gene diversity of the full population:", round(gd.all, 4), "\n")
    cat("  Ranking by dGD on removal (head):\n")
    tmp <- head(ranking)
    tmp[, c("MK", "dGD")] <- round(tmp[, c("MK", "dGD")], 4)
    print(tmp, row.names = FALSE)
    if (nrow(removal.set) > 0) {
      cat("  Greedy removal set:\n")
      tmp <- removal.set
      tmp$gd.after <- round(tmp$gd.after, 4)
      print(tmp, row.names = FALSE)
    } else {
      cat("  No removal increases gene diversity\n")
    }
  }

  # FLAG SCRIPT END ---------------
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN
  return(invisible(list(ranking = ranking, removal.set = removal.set)))
}
