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
#' @param kin Kinship matrix estimated on a wider set of individuals than x,
#' for example gl.kin() on a dataset that includes the source populations;
#' it may cover more individuals than x and is restricted to indNames(x).
#' Kinship estimated on x alone (including kin = NULL) is an error, because
#' mean kinship over the individuals it was estimated on is 0 by
#' construction [required].
#' @param n.best Maximum number of individuals in the greedy removal set, a
#' single number >= 1 (values above nInd - 1 are clamped); if NULL, removals
#' continue until no single removal increases gene diversity [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#'
#' @details
#' Kinship must come from a wider reference than x (see the kin argument).
#'
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
#' removals, whichever comes first (never beyond nInd - 1). The search keeps
#' running sums of the remaining kinship block, so each step costs O(n^2)
#' rather than recomputing the full mean for every candidate.
#'
#' Pairs with missing (NA) kinship are ignored in the gene diversity means,
#' with a warning. Missing genotypes pull kinship and self-kinship toward 0,
#' because gl.kin fills them with the locus mean, so the removal set depends
#' on call rate: for the testset2.gl captive colony (call rates 0.70-0.80)
#' the greedy set has 7 animals, and 11 (6 of them different) after
#' gl.filter.callrate(method = "loc", threshold = 0.95). The function warns
#' when any individual has call rate below 0.8; filter on call rate before
#' gl.kin.
#'
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # SNP data: kinship from the FULL dataset, subset to the captive colony
#' kin <- gl.kin(testset2.gl)
#' cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
#' res <- gl.report.ind.remove(cb, kin = kin[indNames(cb), indNames(cb)])
#' head(res$ranking)  # most-inbred/most-redundant sibs rank first
#' res$removal.set
#' # Tag P/A data (SilicoDArT): the full-dataset kinship is restricted to cb.gs
#' cb.gs <- gl.keep.pop(testset2.gs, pop.list = "EmmacCaptBred", verbose = 0)
#' res.gs <- gl.report.ind.remove(cb.gs, kin = gl.kin(testset2.gs), n.best = 3)
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
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (nInd(x) < 2) {
    stop(error("Fatal Error: at least two individuals are required to assess removals\n"))
  }
  kin <- utils.kin.check(x, kin, verbose = verbose,
                           need.reference = TRUE)
  max.removals <- nInd(x) - 1
  if (!is.null(n.best)) {
    if (!is.numeric(n.best) || length(n.best) != 1 || is.na(n.best) ||
        n.best < 1) {
      stop(error("Fatal Error: n.best must be NULL or a single number >= 1\n"))
    } else if (n.best > nInd(x) - 1) {
      if (verbose >= 1) {
        cat(warn("  Warning: n.best exceeds nInd - 1; clamped to", nInd(x) - 1, "\n"))
      }
    } else {
      max.removals <- as.integer(n.best)
    }
  }

  # Mean imputation of missing genotypes (gl.kin) pulls kinship and
  # self-kinship toward 0 for individuals with low call rates
  ind.cr <- 1 - vapply(x@gen, function(e) length(e@NA.posi), numeric(1)) /
    nLoc(x)
  if (any(ind.cr < 0.8) && verbose >= 1) {
    cat(warn(paste0("  Warning: ", sum(ind.cr < 0.8),
                    " individuals have call rate below 0.8 (lowest ",
                    round(min(ind.cr), 3), "); missing genotypes pull ",
                    "their kinship toward 0. Consider filtering on call ",
                    "rate before gl.kin (see Details)\n")))
  }
  n.na <- sum(is.na(kin[upper.tri(kin, diag = TRUE)]))
  if (n.na > 0 && verbose >= 1) {
    cat(warn("  Warning:", n.na,
             "missing kinship values (pairs or self-kinships) are ignored in the gene diversity means\n"))
  }

  # DO THE JOB ----------------------
  if (verbose >= 2) {
    cat(report("  Computing per-individual change in gene diversity on removal\n"))
  }

  ids <- indNames(x)
  gd.all <- utils.kin.dgd(kin, na.rm = TRUE)
  dgd <- vapply(ids, function(id) {
    utils.kin.dgd(kin, drop = id, na.rm = TRUE) - gd.all
  }, numeric(1))

  ranking <- data.frame(
    id = ids,
    pop = as.character(pop(x)),
    MK = rowMeans(kin, na.rm = TRUE),
    dGD = dgd,
    stringsAsFactors = FALSE
  )
  ranking <- ranking[order(ranking$dGD, decreasing = TRUE), ]
  rownames(ranking) <- NULL

  # Greedy removal set: recompute gains after each removal
  if (verbose >= 2) {
    cat(report("  Constructing greedy removal set (max", max.removals, "removals)\n"))
  }
  # GD of the remaining set R is 1 - S/C, with S the sum and C the count of
  # its non-missing kinships. Removing candidate i gives
  # S - 2 * rowsum_i(R) + k_ii (and likewise for C), so one matrix-vector
  # product per step scores every candidate
  k0 <- kin
  k0[is.na(k0)] <- 0
  obs <- (!is.na(kin)) * 1
  in.set <- rep(1, length(ids))
  S <- sum(k0)
  C <- sum(obs)
  gd.current <- gd.all
  dropped <- character(0)
  set.step <- integer(0)
  set.id <- character(0)
  set.gd <- numeric(0)

  while (length(dropped) < max.removals) {
    rs <- as.vector(k0 %*% in.set)
    cs <- as.vector(obs %*% in.set)
    S.i <- S - 2 * rs + diag(k0)
    C.i <- C - 2 * cs + diag(obs)
    gains <- 1 - S.i / C.i
    gains[in.set == 0] <- -Inf
    best <- which.max(gains)
    if (gains[best] <= gd.current) {
      break
    }
    S <- S.i[best]
    C <- C.i[best]
    in.set[best] <- 0
    dropped <- c(dropped, ids[best])
    gd.current <- gains[best]
    set.step <- c(set.step, length(dropped))
    set.id <- c(set.id, ids[best])
    set.gd <- c(set.gd, gd.current)
  }

  removal.set <- data.frame(
    step = set.step,
    id = set.id,
    gd.after = set.gd,
    stringsAsFactors = FALSE
  )

  # Print out the results summary ---------------
  if (verbose >= 3) {
    cat(report("  Gene diversity of the full population:", round(gd.all, 4),
               "\n"))
    cat(report("  Ranking by dGD on removal (head):\n"))
    tmp <- head(ranking)
    tmp[, c("MK", "dGD")] <- round(tmp[, c("MK", "dGD")], 4)
    print(tmp, row.names = FALSE)
    if (nrow(removal.set) > 0) {
      cat(report("  Greedy removal set:\n"))
      tmp <- removal.set
      tmp$gd.after <- round(tmp$gd.after, 4)
      print(tmp, row.names = FALSE)
    } else {
      cat(report("  No removal increases gene diversity\n"))
    }
  }

  # FLAG SCRIPT END ---------------
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN
  return(invisible(list(ranking = ranking, removal.set = removal.set)))
}
