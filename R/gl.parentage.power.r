#' @name gl.parentage.power
#' @title Estimate the power of a SNP panel to assign parentage
#' @family captive management
#'
#' @description
#' Quantifies, by simulation, how reliably a SNP panel can assign offspring to
#' their parents given a set of candidate parents. Offspring are simulated by
#' Mendelian inheritance from parent pairs drawn from the supplied genlight,
#' optionally degraded with genotyping error and missing data to mimic a real
#' DArTseq run. Each simulated offspring is then compared against every
#' candidate parent by counting opposite-homozygote mismatches (a candidate
#' that shares no allele at a locus cannot be a parent there), and the true
#' parent's retention, the rate of uniquely correct assignment, and the rate of
#' false compatibility with non-parents are tallied. The whole procedure is
#' repeated across a rarefaction series of locus counts, giving a power-versus-
#' number-of-loci curve.
#'
#' @param x Name of the genlight object holding the candidate parents [required].
#' @param n.offspring Number of offspring to simulate per replicate [default 200].
#' @param engine Assignment engine: "exclusion" for the fast Mendelian
#' opposite-homozygote / trio test used to build the power curve, or "colony"
#' to assign the simulated cohort with COLONY via \code{gl.run.colony} for
#' verification against the exclusion result [default "exclusion"].
#' @param pairs If TRUE, assess parent-pair (trio) power: both parents of each
#' offspring are among the candidates, and the offspring must be consistent
#' with one allele from each of the two. If FALSE, assess single-parent power:
#' only one parent of each offspring was sampled (the other parent is removed
#' from the candidates), and candidates are excluded by opposite-homozygote
#' mismatches [default TRUE].
#' @param error.rate Per-genotype rate at which a call is replaced by a random
#' genotype, mimicking miscalls [default 0].
#' @param missing.rate Per-genotype rate at which a call is set to missing
#' [default 0].
#' @param max.mismatch Maximum number of opposite-homozygote mismatches at which
#' a candidate is still treated as a compatible parent. Must be raised above 0
#' whenever error.rate or missing.rate is non-zero, or the true parent will be
#' wrongly excluded [default 0].
#' @param n.loci.steps Vector of locus counts at which to evaluate power. NULL
#' builds an even series from a floor up to all polymorphic loci [default NULL].
#' @param n.rep Number of independent simulation replicates to average over.
#' With engine = "colony" each replicate simulates a new cohort and is a
#' separate COLONY run [default 1].
#' @param plot.display If TRUE, the power curve is displayed [default TRUE].
#' @param plot.theme Theme for the plot [default theme_dartR()].
#' @param plot.colors List of two color names for the two power lines
#' [default c("#2171B5","#6BAED6")].
#' @param plot.file Name for the RDS binary file to save (base name only)
#' [default NULL].
#' @param plot.dir Directory to save plot RDS files [default tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, unless specified using gl.set.verbosity].
#' @param ... When engine = "colony", additional parameters passed to
#' \code{gl.run.colony} (e.g. colony.path, length.run, num.runs, likelihood).
#'
#' @details
#' The function assumes the candidate panel is a single breeding population,
#' as a hatchery broodstock normally is: parents are drawn from it and every
#' offspring is tested against it, so it should be representative of the real
#' breeding group. Only loci polymorphic in the panel carry information, and
#' monomorphic loci are dropped before simulation.
#'
#' Offspring are generated one true parent pair at a time by transmitting one
#' allele from each parent at each locus (a parent scored as a missing genotype
#' transmits a missing genotype).
#'
#' With \code{engine = "exclusion"} (the default) assignment uses a fast
#' Mendelian test with no external program: in pair mode an offspring must be
#' consistent with one allele from each of the two candidates (the trio test);
#' in single-parent mode a candidate is excluded only by an opposite-homozygote
#' mismatch. In single-parent mode one of the two true parents of each
#' offspring, chosen at random, is the sampled parent, and the other is removed
#' from the candidates, as when the second parent was never genotyped. If both
#' stayed in, they would tie at zero mismatches whenever there is no
#' genotyping error, and no offspring could be assigned to one of them alone.
#' The exclusion engine is run across a rarefaction series of locus counts to
#' give the power curve.
#'
#' Three quantities are reported at each locus count:
#' \itemize{
#'   \item \strong{true.retained} -- the proportion of offspring whose sampled
#'   parent (single mode) or true pair (pair mode) is still within the mismatch
#'   tolerance (sensitivity);
#'   \item \strong{correct.unique} -- the proportion of offspring assigned
#'   uniquely and correctly. In pair mode, the true pair is within the
#'   tolerance and no other pair of candidates is; every pair is checked, not
#'   a sample. In single mode, the sampled parent is within the tolerance and
#'   has fewer mismatches than every other candidate. It never exceeds
#'   true.retained;
#'   \item \strong{false.rate} -- the per-comparison false-positive rate: the
#'   proportion of non-parent candidates (single mode) or of all candidate
#'   pairs other than the true pair (pair mode) wrongly within the mismatch
#'   tolerance. Being a rate, it does not depend on the size of the candidate
#'   pool.
#' }
#'
#' Checking every pair stays fast because a pair can only be compatible when
#' each of its members is: an opposite-homozygote mismatch with one parent is
#' a trio mismatch whatever the other parent's genotype. Only pairs among the
#' candidates that pass the single-parent test are evaluated.
#'
#' With \code{engine = "colony"} the simulated cohort is instead assigned by
#' COLONY through \code{gl.run.colony}, at the full panel, so its accuracy can
#' be checked against the exclusion prediction; this is much slower and needs
#' the COLONY executable. Each replicate is a separate COLONY run in its own
#' folder, \code{plot.dir/colony_run/rep1}, \code{rep2} and so on. COLONY
#' assigns parents rather than testing compatibility, and \code{pairs} is not
#' used, so the three columns mean: true.retained, at least one true parent
#' is among the assigned candidates; correct.unique, the assigned candidates
#' are exactly the two true parents; false.rate, the proportion of offspring
#' assigned at least one non-parent (an assignment error rate, not a
#' per-comparison rate).
#'
#' @author Custodian: Peter J. Unmack -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # A single breeding population from the example set.
#' pop1 <- gl.keep.pop(testset.gl, pop.list = popNames(testset.gl)[1],
#'                     verbose = 0)
#' out <- gl.parentage.power(pop1, n.offspring = 50, verbose = 0)
#'
#' @export
#' @return A list with \code{$power} (the power table), \code{$settings} (the
#' call parameters) and \code{$plot} (the ggplot object), returned invisibly.

gl.parentage.power <- function(x,
                               n.offspring = 200,
                               engine = "exclusion",
                               pairs = TRUE,
                               error.rate = 0,
                               missing.rate = 0,
                               max.mismatch = 0,
                               n.loci.steps = NULL,
                               n.rep = 1,
                               plot.display = TRUE,
                               plot.theme = theme_dartR(),
                               plot.colors = NULL,
                               plot.file = NULL,
                               plot.dir = NULL,
                               verbose = NULL,
                               ...) {

  # PRELIMINARIES -- checking ----------------
  funname <- match.call()[[1]]
  verbose <- gl.check.verbosity(verbose)
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  if (is.null(plot.colors)) {
    plot.colors <- c("#2171B5", "#6BAED6")
  }
  utils.flag.start(func = funname, build = "v.2026.1", verbose = verbose)
  datatype <- utils.check.datatype(x, accept = "genlight", verbose = verbose)
  engine <- match.arg(engine, c("exclusion", "colony"))

  if (nInd(x) < 2) {
    stop(error("Fatal error: at least two candidate parents are required\n"))
  }
  if (error.rate > 0 || missing.rate > 0) {
    if (max.mismatch == 0 && verbose >= 1) {
      cat(warn(
        "  Warning: error.rate or missing.rate is non-zero but max.mismatch",
        "is 0; the true parent may be wrongly excluded. Consider raising",
        "max.mismatch.\n"
      ))
    }
  }

  # DO THE JOB --------------------------------
  # Drop monomorphic loci: they carry no parentage information.
  x <- gl.filter.monomorphs(x, verbose = 0)
  n.loc.total <- nLoc(x)
  if (n.loc.total < 1) {
    stop(error("Fatal error: no polymorphic loci remain after filtering\n"))
  }

  # Candidate dosage matrix: individuals x loci, values in {0, 1, 2}, NA allowed.
  geno <- as.matrix(x)
  n.cand <- nrow(geno)

  # Rarefaction series of locus counts.
  if (is.null(n.loci.steps)) {
    floor.loci <- min(50, n.loc.total)
    n.loci.steps <- unique(round(seq(floor.loci, n.loc.total, length.out = 6)))
    n.loci.steps <- n.loci.steps[n.loci.steps >= 1]
  }
  n.loci.steps <- sort(unique(pmin(n.loci.steps, n.loc.total)))

  # Simulate one offspring's dosage from two parent dosage vectors. Missing
  # parent genotypes are set to a transmission probability of 0 to avoid an
  # rbinom NA warning, then the offspring locus is blanked to NA afterwards.
  sim_offspring <- function(pa, pb) {
    ppa <- pa / 2
    ppb <- pb / 2
    ppa[is.na(ppa)] <- 0
    ppb[is.na(ppb)] <- 0
    off <- stats::rbinom(length(pa), 1, ppa) + stats::rbinom(length(pb), 1, ppb)
    off[is.na(pa) | is.na(pb)] <- NA
    off
  }

  # Degrade a dosage vector with miscalls and missingness.
  degrade <- function(g) {
    if (error.rate > 0) {
      hit <- which(stats::runif(length(g)) < error.rate)
      if (length(hit)) g[hit] <- sample(0:2, length(hit), replace = TRUE)
    }
    if (missing.rate > 0) {
      miss <- which(stats::runif(length(g)) < missing.rate)
      if (length(miss)) g[miss] <- NA
    }
    g
  }

  # COLONY verification engine -----------------
  # Assigns the simulated cohort with COLONY at the full panel, so its accuracy
  # can be compared with the exclusion prediction. Candidates are treated as
  # monoecious (unknown sex), each a potential parent of either sex.
  if (engine == "colony") {
    if (verbose >= 2) {
      cat(report("  Engine = colony: assigning the cohort with COLONY at the",
                 "full panel of", n.loc.total, "loci\n"))
    }
    cand.names <- indNames(x)
    results <- list()
    for (rep.i in seq_len(n.rep)) {
      if (verbose >= 2) {
        cat(report("  Replicate", rep.i, "of", n.rep, "\n"))
      }
      pair.idx <- t(sapply(seq_len(n.offspring), function(i)
        sample.int(n.cand, 2)))
      # Filled row by row so a single locus or a single offspring still gives
      # an offspring x locus matrix.
      off.mat <- matrix(NA_real_, nrow = n.offspring, ncol = n.loc.total,
                        dimnames = list(paste0("off", seq_len(n.offspring)),
                                        colnames(geno)))
      for (i in seq_len(n.offspring)) {
        off.mat[i, ] <- degrade(sim_offspring(geno[pair.idx[i, 1], ],
                                              geno[pair.idx[i, 2], ]))
      }

      comb <- rbind(geno, off.mat)
      g <- new("genlight", comb, ploidy = 2)
      g@other$loc.metrics <- x@other$loc.metrics
      g@loc.all <- x@loc.all
      g@position <- x@position
      g@chromosome <- x@chromosome
      is.off <- indNames(g) %in% rownames(off.mat)
      g@other$ind.metrics <- data.frame(
        id = indNames(g),
        offspring = ifelse(is.off, "yes", "no"),
        father = ifelse(is.off, "no", "yes"),
        mother = ifelse(is.off, "no", "yes"),
        stringsAsFactors = FALSE
      )

      outp <- file.path(plot.dir, "colony_run", paste0("rep", rep.i))
      dir.create(outp, showWarnings = FALSE, recursive = TRUE)
      res <- gl.run.colony(g, outpath = outp, di.mono.ecious = 1,
                           verbose = if (verbose >= 2) 2 else 0, ...)
      bc <- res$best.config

      true.ret <- 0
      correct.uni <- 0
      false.tot <- 0
      for (i in seq_len(n.offspring)) {
        oid <- paste0("off", i)
        rowi <- bc[bc$OffspringID == oid, , drop = FALSE]
        assigned <- character(0)
        if (nrow(rowi) == 1) {
          assigned <- c(rowi$FatherID, rowi$MotherID)
        }
        # Keep only assignments to real candidates (inferred parents are "#n").
        assigned.real <- assigned[assigned %in% cand.names]
        true.names <- cand.names[pair.idx[i, ]]
        if (any(true.names %in% assigned.real)) true.ret <- true.ret + 1
        if (setequal(assigned.real, true.names)) correct.uni <- correct.uni + 1
        if (any(!(assigned.real %in% true.names))) false.tot <- false.tot + 1
      }
      results[[rep.i]] <- data.frame(
        n.loci = n.loc.total,
        true.retained = true.ret / n.offspring,
        correct.unique = correct.uni / n.offspring,
        false.rate = false.tot / n.offspring
      )
    }
    power <- do.call(rbind, results)
    # Average over replicates.
    power <- stats::aggregate(
      cbind(true.retained, correct.unique, false.rate) ~ n.loci,
      data = power, FUN = mean
    )
  } else {

  # Exclusion engine ---------------------------
  results <- list()
  row <- 1
  for (rep.i in seq_len(n.rep)) {
    if (verbose >= 2) {
      cat(report("  Replicate", rep.i, "of", n.rep, "\n"))
    }
    # Draw true parent pairs (two distinct candidates) for each offspring. The
    # order is random, so in single mode the first is the sampled parent.
    pair.idx <- t(sapply(seq_len(n.offspring), function(i)
      sample.int(n.cand, 2)))
    for (nl in n.loci.steps) {
      loc.sub <- sort(sample.int(n.loc.total, nl))
      cand.sub <- geno[, loc.sub, drop = FALSE]
      true.ret <- 0
      correct.uni <- 0
      false.tot <- 0
      for (i in seq_len(n.offspring)) {
        true.pars <- pair.idx[i, ]
        off <- degrade(sim_offspring(cand.sub[true.pars[1], ],
                                     cand.sub[true.pars[2], ]))
        sc <- if (pairs) {
          utils.parentage.score.pair(off, cand.sub, true.pars, max.mismatch)
        } else {
          utils.parentage.score.single(off, cand.sub, target = true.pars[1],
                                       other = true.pars[2], max.mismatch)
        }
        true.ret <- true.ret + sc$retained
        correct.uni <- correct.uni + sc$unique
        false.tot <- false.tot + sc$false.rate
      }
      results[[row]] <- data.frame(
        n.loci = nl,
        true.retained = true.ret / n.offspring,
        correct.unique = correct.uni / n.offspring,
        false.rate = false.tot / n.offspring
      )
      row <- row + 1
    }
  }
  power <- do.call(rbind, results)
  # Average over replicates.
  power <- stats::aggregate(
    cbind(true.retained, correct.unique, false.rate) ~ n.loci,
    data = power, FUN = mean
  )
  } # end exclusion engine

  if (verbose >= 3) {
    cat(report("  Power at the full panel of", n.loc.total, "loci:\n"))
    full <- power[power$n.loci == max(power$n.loci), ]
    false.label <- if (engine == "colony") {
      "offspring assigned a non-parent"
    } else {
      "false-positive rate per comparison"
    }
    cat(report(sprintf(
      "    true parent retained: %.3f; correct unique assignment: %.3f; %s: %.4f\n",
      full$true.retained, full$correct.unique, false.label, full$false.rate
    )))
  }

  # PLOT THE RESULTS --------------------------
  plot.df <- rbind(
    data.frame(n.loci = power$n.loci, value = power$true.retained,
               metric = "true parent retained"),
    data.frame(n.loci = power$n.loci, value = power$correct.unique,
               metric = "correct unique assignment")
  )
  n.loci <- value <- metric <- NULL
  p3 <- ggplot2::ggplot(plot.df,
                        ggplot2::aes(x = n.loci, y = value, colour = metric)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = plot.colors) +
    ggplot2::ylim(0, 1) +
    ggplot2::labs(x = "Number of loci", y = "Power", colour = NULL) +
    plot.theme +
    ggplot2::theme(legend.position = "bottom")

  if (plot.display) {
    suppressWarnings(print(p3))
  }

  # SAVE PLOT (OPTIONAL) ----------------------
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p3, dir = plot.dir, file = plot.file,
                           verbose = verbose)
  }

  # FLAG SCRIPT END ---------------------------
  if (verbose >= 1) {
    cat(report("Completed:", as.character(funname), "\n"))
  }

  # RETURN
  invisible(list(
    power = power,
    settings = list(engine = engine, pairs = pairs,
                    n.offspring = n.offspring,
                    error.rate = error.rate, missing.rate = missing.rate,
                    max.mismatch = max.mismatch, n.rep = n.rep,
                    n.loci.total = n.loc.total),
    plot = p3
  ))
}


###################### Define function utils.parentage.trio.mismatch ###########
## Trio Mendelian incompatibilities between an offspring and one candidate
## pair, as dosage vectors over the same loci. An offspring dosage is possible
## only if it can be formed by transmitting one allele (B-count 0 or 1) from
## each parent. A locus with a missing call counts only when it is
## incompatible whatever the missing genotype is (an offspring 0 with a parent
## 2), which R's three-valued logic gives through na.rm.
utils.parentage.trio.mismatch <- function(off, a, b) {
  a0 <- a <= 1; a1 <- a >= 1          # parent a can transmit a 0 / a 1 allele
  b0 <- b <= 1; b1 <- b >= 1
  ok0 <- a0 & b0                       # offspring 0 achievable
  ok1 <- (a0 & b1) | (a1 & b0)         # offspring 1 achievable
  ok2 <- a1 & b1                       # offspring 2 achievable
  achievable <- (off == 0 & ok0) | (off == 1 & ok1) | (off == 2 & ok2)
  sum(!achievable, na.rm = TRUE)
}
################################################################################


###################### Define function utils.parentage.pair.mismatch ###########
## The same count for every pair of candidates at once: a symmetric matrix
## whose [a, b] cell equals utils.parentage.trio.mismatch(off, cand[a, ],
## cand[b, ]) (the diagonal is meaningless). Per offspring dosage, a pair is
## incompatible at a locus when
##   offspring 0: either candidate is 2 (union of the two candidates' counts);
##   offspring 2: either candidate is 0;
##   offspring 1: both candidates are 0, or both are 2.
## Unions are counted as count(a) + count(b) - count(a and b), with the
## intersections as cross-products.
utils.parentage.pair.mismatch <- function(off, cand) {
  hit <- function(dosage, loci) {
    h <- cand[, loci, drop = FALSE] == dosage
    h[is.na(h)] <- FALSE
    h * 1
  }
  o0 <- which(off == 0)
  o1 <- which(off == 1)
  o2 <- which(off == 2)
  a2 <- hit(2, o0)
  a0 <- hit(0, o2)
  z1 <- hit(0, o1)
  t1 <- hit(2, o1)
  single <- rowSums(a2) + rowSums(a0)
  outer(single, single, "+") - tcrossprod(a2) - tcrossprod(a0) +
    tcrossprod(z1) + tcrossprod(t1)
}
################################################################################


###################### Define function utils.parentage.single.mismatch #########
## Opposite-homozygote mismatches between an offspring and every candidate
## (rows of cand): the single-parent exclusion signal.
utils.parentage.single.mismatch <- function(off, cand) {
  o0 <- which(off == 0)
  o2 <- which(off == 2)
  rowSums(cand[, o0, drop = FALSE] == 2, na.rm = TRUE) +
    rowSums(cand[, o2, drop = FALSE] == 0, na.rm = TRUE)
}
################################################################################


###################### Define function utils.parentage.score.pair ##############
## Pair-mode result for one offspring. Every candidate pair other than the
## true pair is checked. A pair can only be compatible when both members pass
## the single-parent test, because an opposite-homozygote mismatch with one
## parent is a trio mismatch whatever the other parent is, so only pairs among
## those candidates are evaluated.
utils.parentage.score.pair <- function(off, cand, true.pars, max.mismatch) {
  n.alt <- choose(nrow(cand), 2) - 1
  true.mm <- utils.parentage.trio.mismatch(off, cand[true.pars[1], ],
                                           cand[true.pars[2], ])
  retained <- true.mm <= max.mismatch
  ok <- which(utils.parentage.single.mismatch(off, cand) <= max.mismatch)
  n.false <- 0
  if (length(ok) >= 2) {
    mm <- utils.parentage.pair.mismatch(off, cand[ok, , drop = FALSE])
    compat <- mm <= max.mismatch & upper.tri(mm)
    tp <- match(true.pars, ok)
    if (!anyNA(tp)) compat[min(tp), max(tp)] <- FALSE
    n.false <- sum(compat)
  }
  list(retained = retained,
       unique = retained && n.false == 0,
       false.rate = if (n.alt > 0) n.false / n.alt else 0)
}
################################################################################


###################### Define function utils.parentage.score.single ############
## Single-mode result for one offspring. Only the target parent was sampled:
## the other true parent is not a candidate. The target is assigned when it is
## within the tolerance and has fewer mismatches than every other candidate.
utils.parentage.score.single <- function(off, cand, target, other,
                                         max.mismatch) {
  mm <- utils.parentage.single.mismatch(off, cand)
  pool <- setdiff(seq_len(nrow(cand)), other)
  best <- pool[mm[pool] == min(mm[pool])]
  nonpar <- setdiff(pool, target)
  retained <- mm[target] <= max.mismatch
  list(retained = retained,
       unique = retained && length(best) == 1 && best == target,
       false.rate = if (length(nonpar)) mean(mm[nonpar] <= max.mismatch) else 0)
}
################################################################################
