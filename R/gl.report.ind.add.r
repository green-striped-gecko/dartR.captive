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
#' the reported gain is exact given kin, and candidates that are cryptic
#' relatives of the existing stock are ranked appropriately lower. Kinship is
#' relative to the individuals it was estimated on, so estimate kin on the
#' widest dataset available: on testset2.gl, kin estimated on the target and
#' one candidate population only gave gains correlating 0.72 with those from
#' kin estimated on the full dataset. It is a natural
#' companion to the \code{gl.assign} suite: having assigned a stray or
#' wild-caught individual to a source population, ask what it would
#' contribute to the managed colony.
#'
#' Gene diversity of the target block is GD = 1 - mean(kin[block, block]),
#' the mean taken over the full sub-matrix including the diagonal. For each
#' candidate c the function computes
#' dgd = GD(target block plus c's row and column from kin) - GD(target block).
#' A positive dgd means adding c raises the target population's gene
#' diversity. Pairs with missing (NA) kinship are ignored in these means,
#' with a warning.
#'
#' Missing genotypes pull kinship toward 0, because gl.kin fills them with
#' the locus mean, so the ranking depends on call rate: on testset2.gl
#' (captive-bred target, call rates 0.70-0.80) only 5 of the top 10 wild
#' candidates stay in the top 10 after gl.filter.callrate(method = "loc",
#' threshold = 0.95). The function warns when any target or candidate
#' individual has call rate below 0.8; filter on call rate before gl.kin.
#'
#' If \code{candidates} is a single name that matches a population in
#' \code{popNames(x)}, it is expanded to all members of that population (the
#' population-name interpretation takes precedence over an individual of the
#' same name). Otherwise every element must match an individual name in
#' \code{indNames(x)}. Duplicated candidates are evaluated once. Candidates
#' already belonging to \code{target.pop} are rejected with a fatal error.
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # SNP data
#' kin <- gl.kin(testset2.gl)
#' res <- gl.report.ind.add(testset2.gl, candidates = "EmmacMaclGeor",
#'                          target.pop = "EmmacCaptBred", kin = kin)
#' head(res)  # wild individuals ranked by the diversity they would add
#' # Tag P/A data (SilicoDArT; kinship computed internally)
#' res.gs <- gl.report.ind.add(testset2.gs, candidates = "EmmacMaclGeor",
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
  utils.flag.start(func = funname, verbose = verbose)

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
    if (anyDuplicated(candidates)) {
      if (verbose >= 2) {
        cat(report("  Dropping", sum(duplicated(candidates)),
                   "duplicated candidate id(s)\n"))
      }
      candidates <- unique(candidates)
    }
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

  # Mean imputation of missing genotypes (gl.kin) pulls kinship toward 0
  # for individuals with low call rates, which shifts the ranking
  used <- match(c(t.ids, candidates), indNames(x))
  ind.cr <- 1 - vapply(x@gen[used], function(e) length(e@NA.posi),
                       numeric(1)) / nLoc(x)
  if (any(ind.cr < 0.8) && verbose >= 1) {
    cat(warn(paste0("  Warning: ", sum(ind.cr < 0.8),
                    " target or candidate individuals have call rate below ",
                    "0.8 (lowest ", round(min(ind.cr), 3), "); missing ",
                    "genotypes pull their kinship toward 0. Consider ",
                    "filtering on call rate before gl.kin (see Details)\n")))
  }
  sub <- kin[c(t.ids, candidates), c(t.ids, candidates), drop = FALSE]
  n.na <- sum(is.na(sub[upper.tri(sub)]))
  if (n.na > 0 && verbose >= 1) {
    cat(warn("  Warning:", n.na,
             "pairs have missing kinship and are ignored in the gene diversity means\n"))
  }

  gd.target <- utils.kin.dgd(kin[t.ids, t.ids, drop = FALSE], na.rm = TRUE)

  dgd <- vapply(candidates, function(id) {
    joint <- c(t.ids, id)
    utils.kin.dgd(kin[joint, joint, drop = FALSE], na.rm = TRUE) - gd.target
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
    cat(report("  Gene diversity of", target.pop, ":", round(gd.target, 4),
               "\n"))
    cat(report("  Candidates ranked by gene diversity added:\n"))
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
