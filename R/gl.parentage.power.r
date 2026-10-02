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
#' @param pairs If TRUE, assess parent-pair (trio) power -- whether the
#' offspring is consistent with one allele from each of the two candidates; if
#' FALSE, assess single-parent power by opposite-homozygote exclusion
#' [default TRUE].
#' @param n.false.pairs In pair mode, the number of non-true candidate pairs
#' sampled per offspring to estimate the false-pair rate [default 200].
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
#' @param n.rep Number of independent simulation replicates to average over
#' [default 1].
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
#' mismatch. The exclusion engine is run across a rarefaction series of locus
#' counts to give the power curve. With \code{engine = "colony"} the simulated
#' cohort is instead assigned by COLONY through \code{gl.run.colony}, at the
#' full panel, so its accuracy can be checked against the exclusion prediction;
#' this is much slower and needs the COLONY executable.
#'
#' Three quantities are reported at each locus count:
#' \itemize{
#'   \item \strong{true.retained} -- the proportion of offspring whose true
#'   parent (single mode) or true pair (pair mode) is still within the mismatch
#'   tolerance (sensitivity);
#'   \item \strong{correct.unique} -- the proportion of offspring assigned
#'   uniquely and correctly to the true parent (single mode) or true pair
#'   (pair mode);
#'   \item \strong{false.rate} -- the per-comparison false-positive rate: the
#'   proportion of non-parent candidates (single mode) or sampled non-true
#'   pairs (pair mode) wrongly within the mismatch tolerance. Being a rate, it
#'   does not depend on the size of the candidate pool.
#' }
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
                               n.false.pairs = 200,
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
    pair.idx <- t(sapply(seq_len(n.offspring), function(i)
      sample.int(n.cand, 2)))
    off.mat <- t(sapply(seq_len(n.offspring), function(i)
      degrade(sim_offspring(geno[pair.idx[i, 1], ], geno[pair.idx[i, 2], ]))))
    rownames(off.mat) <- paste0("off", seq_len(n.offspring))
    colnames(off.mat) <- colnames(geno)

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

    outp <- file.path(plot.dir, "colony_run")
    dir.create(outp, showWarnings = FALSE, recursive = TRUE)
    res <- gl.run.colony(g, outpath = outp, di.mono.ecious = 1,
                         verbose = if (verbose >= 2) 2 else 0, ...)
    bc <- res$best.config
    cand.names <- indNames(x)

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
    power <- data.frame(
      n.loci = n.loc.total,
      true.retained = true.ret / n.offspring,
      correct.unique = correct.uni / n.offspring,
      false.rate = false.tot / n.offspring
    )
  } else {

  # Exclusion engine ---------------------------
  # Opposite-homozygote mismatches between an offspring and one candidate,
  # over the locus subset. Single-parent exclusion signal.
  single_mismatch <- function(off, cand) {
    sum((off == 0 & cand == 2) | (off == 2 & cand == 0), na.rm = TRUE)
  }

  # Trio Mendelian incompatibilities between an offspring and a candidate pair.
  # An offspring dosage o is possible only if it can be formed by transmitting
  # one allele (B-count 0 or 1) from each parent; a locus with any NA is skipped.
  trio_mismatch <- function(off, a, b) {
    a0 <- a <= 1; a1 <- a >= 1          # parent a can transmit a 0 / a 1 allele
    b0 <- b <= 1; b1 <- b >= 1
    ok0 <- a0 & b0                       # offspring 0 achievable
    ok1 <- (a0 & b1) | (a1 & b0)         # offspring 1 achievable
    ok2 <- a1 & b1                       # offspring 2 achievable
    achievable <- (off == 0 & ok0) | (off == 1 & ok1) | (off == 2 & ok2)
    sum(!achievable, na.rm = TRUE)
  }

  n.pairs.total <- choose(n.cand, 2)

  results <- list()
  row <- 1
  for (rep.i in seq_len(n.rep)) {
    if (verbose >= 2) {
      cat(report("  Replicate", rep.i, "of", n.rep, "\n"))
    }
    # Draw true parent pairs (two distinct candidates) for each offspring.
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
        pa <- geno[true.pars[1], loc.sub]
        pb <- geno[true.pars[2], loc.sub]
        off <- degrade(sim_offspring(pa, pb))

        if (pairs) {
          # Parent-pair (trio) assessment. The true pair is tested for
          # retention; a Monte-Carlo sample of other pairs estimates the
          # per-pair false-positive rate.
          true.mm <- trio_mismatch(off, pa, pb)
          if (true.mm <= max.mismatch) true.ret <- true.ret + 1
          k <- min(n.false.pairs, max(n.pairs.total - 1, 0))
          false.compat <- 0
          if (k > 0) {
            for (s in seq_len(k)) {
              repeat {
                pr <- sample.int(n.cand, 2)
                if (!setequal(pr, true.pars)) break
              }
              mm.s <- trio_mismatch(off, geno[pr[1], loc.sub],
                                    geno[pr[2], loc.sub])
              if (mm.s <= max.mismatch) false.compat <- false.compat + 1
            }
          }
          false.tot <- false.tot + (if (k > 0) false.compat / k else 0)
          # correct unique: true pair compatible and no sampled pair compatible.
          if (true.mm <= max.mismatch && false.compat == 0) {
            correct.uni <- correct.uni + 1
          }
        } else {
          # Single-parent (opposite-homozygote) assessment against every
          # candidate; argmin assignment. The false rate is the proportion of
          # non-parent candidates within the mismatch tolerance.
          mm <- vapply(seq_len(n.cand),
                       function(j) single_mismatch(off, cand.sub[j, ]),
                       numeric(1))
          if (any(mm[true.pars] <= max.mismatch)) true.ret <- true.ret + 1
          best <- which(mm == min(mm))
          if (length(best) == 1 && best %in% true.pars) {
            correct.uni <- correct.uni + 1
          }
          within <- which(mm <= max.mismatch)
          n.nonpar <- n.cand - length(unique(true.pars))
          false.tot <- false.tot +
            (if (n.nonpar > 0) sum(!(within %in% true.pars)) / n.nonpar else 0)
        }
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
    cat(report(sprintf(
      "    true parent retained: %.3f; correct unique assignment: %.3f; false-positive rate per comparison: %.4f\n",
      full$true.retained, full$correct.unique, full$false.rate
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
                    n.offspring = n.offspring, n.false.pairs = n.false.pairs,
                    error.rate = error.rate, missing.rate = missing.rate,
                    max.mismatch = max.mismatch, n.rep = n.rep,
                    n.loci.total = n.loc.total),
    plot = p3
  ))
}
