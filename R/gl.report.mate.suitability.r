#' @name gl.report.mate.suitability
#' @title Reports a mate suitability matrix for all male x female pairings
#' @family captive management
#' @description
#' Scores every potential male x female pairing in a genlight object against
#' four criteria derived from a pairwise kinship matrix and from genotype
#' completeness, and combines them into a Mate Suitability Index (MSI) ranging
#' from 1 (very beneficial to the population) to 6 (very detrimental), with
#' 'NoWay' flagging pairings that should not be bred at all. This is a
#' SNP-based analogue of the MateRx screen in the PMx studbook-management
#' software, implementing the full published MateRx rule set (the 'Tulsa'
#' rules) rather than a simplification.
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param kin Pairwise kinship matrix with row and column names identical to
#' indNames(x), as produced by gl.kin(); if NULL, computed internally with
#' gl.kin() [default NULL].
#' @param f.noway No Way point on the prospective offspring inbreeding
#' coefficient (the pair's kinship); pairings at or above this value are scored
#' 'NoWay' [default 0.125, the inbreeding of offspring of a half-sib mating].
#' Automatically raised to 0.25, 0.5 or 1.0 when the mean pairwise kinship of
#' the candidate pairs exceeds 0.125, 0.25 or 0.5 respectively, per the PMx
#' rule for inbred populations.
#' @param unknown.breaks Either NULL to disable the unknown-genome penalty, or
#' four strictly increasing thresholds on the unknown (uncalled) fraction of
#' the prospective offspring genome at which the MSI is raised to at least 4,
#' at least 5, exactly 6, and 'NoWay' respectively. The PMx pedigree defaults
#' are c(0.0625, 0.125, 0.25, 0.5); they are calibrated for pedigree
#' unknownness and on routine SNP missingness they dominate the score, so the
#' default here is NULL [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#' @details
#' Kinship should be estimated on the most inclusive dataset available and the
#' matrix subset to the managed group; computing kinship on a small family
#' group alone inflates estimates (the allele-frequency reference collapses
#' onto the family itself).
#'
#' This function implements the MateRx capability of PMx (Lacy, Ballou &
#' Pollak 2012, Methods in Ecology and Evolution 3:433-437; PMx Users Manual
#' v1.0 pp. 84-86 and 107-112; original rules in Ballou et al. 1999), which
#' assigns every potential pairing a Mate Suitability Index of 1 (very
#' beneficial genetically), 2 (moderately beneficial), 3 (slightly
#' beneficial), 4 (slightly detrimental), 5 (detrimental, only if
#' demographically necessary) or 6 (very detrimental), plus a do-not-breed
#' category (PMx displays '-'; here labelled 'NoWay'). The pedigree-derived
#' inputs of PMx are replaced by their SNP-derived analogues, all taken from
#' the pairwise kinship matrix.
#'
#' Four component matrices (males in rows, females in columns) are computed
#' over the sexed individuals:
#' \itemize{
#' \item F.off -- the inbreeding coefficient of the prospective offspring,
#' equal to the parents' pairwise kinship kin[m, f].
#' \item deltaGD -- the change in gene diversity of the population were the
#' pair to produce one offspring, computed analytically via utils.kin.dgd()
#' (offspring kinships are parental means), relative to the current gene
#' diversity.
#' \item MKdiff -- the absolute difference in the parents' mean kinships,
#' abs(MK_m - MK_f). Pairing a genetically over-represented individual with an
#' under-represented one squanders the rarer genome, so smaller is better.
#' \item completeness -- 1 minus the unknown fraction of the prospective
#' offspring genome, taken as the mean of the parents' per-individual call
#' rates (each parent contributes half the offspring genome). This is the
#' SNP-data stand-in for PMx's percent-known-ancestry component.
#' }
#'
#' Scoring follows the published PMx defaults. Pairings with F.off at or above
#' the No Way point (or, when unknown.breaks is supplied, with more than
#' unknown.breaks[4] of the offspring genome uncalled) are scored 'NoWay' and
#' excluded from rating (Manual pp. 86, 108-109). Each remaining pair is ranked 1-6 on each genetic
#' component (Manual pp. 108-110): deltaGD around a break point of zero, with
#' bins 1-3 each spanning a third of the observed range above zero and bins
#' 4-6 a third of the range below; MKdiff around a break point at the mean
#' MKdiff of the rated pairs, bins 1-3 spanning thirds of the range from the
#' smallest observed value to the break point and bins 4-6 thirds of the range
#' from the break point to the largest; F with bin 1 pinned to F <= 0, bins 2
#' and 3 each spanning half the range from 0 to the break point (the mean
#' pairwise kinship of all candidate pairs), and bins 4-6 each a third of the
#' range from the break point to the No Way point. The three ranks are then
#' combined by the Tulsa logic rules (Manual p. 111): if all ranks are below
#' 4, equal triples map to themselves, a rank sum below 5 gives MSI 1, a rank
#' sum below 6 or a deltaGD rank of 1 gives MSI 2, and anything else gives
#' MSI 3; if any rank is 4 or more, an F rank of 6 or an MKdiff rank of 5 or 6
#' gives MSI 6, otherwise a deltaGD rank of 1-4 gives MSI 4, a deltaGD rank of
#' 6 gives MSI 5, and a deltaGD rank of 5 gives MSI 5 when the F rank is 5 and
#' MSI 4 otherwise. When unknown.breaks is supplied, the unknown-genome
#' penalty is then applied by the PMx default 'MSI Minimum' method (Manual
#' p. 111): offspring unknown fraction above unknown.breaks[1] raises the MSI
#' to at least 4, above unknown.breaks[2] to at least 5, above
#' unknown.breaks[3] to 6 (above unknown.breaks[4] the pair is already
#' 'NoWay'); by default (unknown.breaks = NULL) no penalty is applied.
#'
#' Genomic departures from the pedigree formulation: bin 1 of F is F <= 0
#' rather than F = 0 exactly (genomic kinship estimates are continuous and can
#' be negative); and when the mean pairwise kinship is not strictly between 0
#' and the No Way point, the F break point falls back to half the No Way
#' point. The published proration thresholds for unknown ancestry are
#' calibrated for pedigree unknownness, a quantity not comparable to SNP
#' no-call rate: routine missingness (10-30%) would dominate the score while
#' leaving the kinship estimates essentially unaffected. The penalty is
#' therefore off by default and opt-in via unknown.breaks.
#'
#' Sexes are read from x@@other$ind.metrics$sex (values 'Male', 'Female',
#' 'Unknown'); individuals of unknown sex are excluded from the pools with a
#' warning (the PMx default treatment, Manual p. 112).
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # Estimate kinship on the FULL dataset, then subset to the captive colony
#' kin <- gl.kin(testset.gl)
#' cb <- gl.keep.pop(testset.gl, pop.list = "EmmacCaptBred", verbose = 0)
#' kin.cb <- kin[indNames(cb), indNames(cb)]
#' res <- gl.report.mate.suitability(cb, kin = kin.cb)
#' table(res$msi)   # spans 1-6; kin pairings (F.off >= 0.125) score NoWay
#' res$f.off["CB_AB_01", "CB_AB_02"]  # full sibs: ~0.23 -> 'NoWay'
#' # A low offspring F does not guarantee a good MSI: gene-diversity change
#' # and mean-kinship difference also count
#' res$f.off["CB_AB_01", "CB_CD_02"]  # across families: F.off ~0.01
#' res$msi["CB_AB_01", "CB_CD_02"]    # yet detrimental: both well represented
#' @seealso \code{\link{gl.kin}}, \code{\link{gl.select.pairs}},
#' \code{\link{gl.report.repro.targets}}
#' @export
#' @return Invisibly, a list of ten components, the matrices all males x
#' females: $msi (character, '1'..'6' or 'NoWay'), $f.off, $dgd, $mkdiff and
#' $completeness (numeric), $rank.dgd, $rank.mkdiff and $rank.f (integer bin
#' ranks 1-6, NA for 'NoWay' pairs), $msi.base (numeric Tulsa score before the
#' unknown-genome penalty, NA for 'NoWay' pairs), and $settings (list: the
#' No Way point and break points used, the component range limits, and
#' unknown.breaks).
#'
# ----------------------
# Function
gl.report.mate.suitability <- function(x,
                                       kin = NULL,
                                       f.noway = 0.125,
                                       unknown.breaks = NULL,
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

# FUNCTION SPECIFIC ERROR CHECKING ----------------
    # Kinship matrix -- auto-computed via gl.kin when NULL, validated against
    # indNames(x), by the shared series validator
    kin <- utils.kin.check(x, kin, verbose = verbose)

    # Sexes (PMx default 'For Unknown Sexes: Exclude Unknowns')
    if (is.null(x@other$ind.metrics$sex)) {
        stop(error(
            "Fatal Error: mate suitability requires a 'sex' column (values 'Male', 'Female', 'Unknown') in x@other$ind.metrics; none found\n"
        ))
    }
    ids <- indNames(x)
    sex <- as.character(x@other$ind.metrics$sex)
    males <- ids[sex == "Male"]
    females <- ids[sex == "Female"]
    n.unknown <- sum(!(sex %in% c("Male", "Female")))
    if (n.unknown > 0 && verbose >= 1) {
        cat(warn("  Warning:", n.unknown,
                 "individuals of unknown sex excluded from the male/female pools\n"))
    }
    if (length(males) == 0 || length(females) == 0) {
        stop(error(
            "Fatal Error: at least one Male and one Female are required to score pairings\n"
        ))
    }

    # f.noway
    if (!is.numeric(f.noway) || length(f.noway) != 1 ||
        f.noway <= 0 || f.noway > 0.5) {
        if (verbose >= 1) {
            cat(warn("  Warning: f.noway must be a single value in (0, 0.5]; resetting to 0.125\n"))
        }
        f.noway <- 0.125
    }

    # unknown.breaks: NULL disables the unknown-genome penalty
    if (!is.null(unknown.breaks)) {
        if (!is.numeric(unknown.breaks) || length(unknown.breaks) != 4 ||
            is.unsorted(unknown.breaks, strictly = TRUE) ||
            any(unknown.breaks <= 0) || any(unknown.breaks > 1)) {
            stop(error(
                "Fatal Error: unknown.breaks must be NULL or four strictly increasing values in (0, 1]\n"
            ))
        }
    }

# DO THE JOB ----------------------
    n.m <- length(males)
    n.f <- length(females)

    # Component matrices ------------------------------------

    # (a) Prospective offspring inbreeding
    f.off <- kin[males, females, drop = FALSE]

    # (b) Delta gene diversity per pairing (baseline computed once)
    if (verbose >= 2) {
        cat(report("  Computing delta gene diversity for", n.m * n.f,
                   "candidate pairings\n"))
    }
    gd.base <- utils.kin.dgd(kin)
    dgd <- matrix(NA_real_, nrow = n.m, ncol = n.f,
                  dimnames = list(males, females))
    for (m in males) {
        for (f in females) {
            dgd[m, f] <- utils.kin.dgd(kin,
                                       add.pairs = matrix(c(m, f), nrow = 1)) - gd.base
        }
    }

    # (c) Mean kinship difference (PMx defaults: Mean Kinship, Absolute Diffs)
    mk <- rowMeans(kin)
    mkdiff <- abs(outer(mk[males], mk[females], "-"))
    dimnames(mkdiff) <- list(males, females)

    # (d) Completeness: the prospective offspring draws half its genome from
    # each parent, so its unknown (uncalled) fraction is the mean of the
    # parents' missing-call fractions; completeness = 1 - that
    cr <- rowMeans(!is.na(as.matrix(x)))
    names(cr) <- ids
    completeness <- outer(cr[males], cr[females], "+") / 2
    dimnames(completeness) <- list(males, females)
    unknown <- 1 - completeness

    # No Way point (raised automatically for inbred populations) -----------
    f.mean <- mean(f.off)
    noway <- f.noway
    if (f.mean > 0.5) {
        noway <- max(noway, 1.0)
    } else if (f.mean > 0.25) {
        noway <- max(noway, 0.5)
    } else if (f.mean > 0.125) {
        noway <- max(noway, 0.25)
    }
    if (noway > f.noway && verbose >= 1) {
        cat(warn("  Warning: mean pairwise kinship", round(f.mean, 4),
                 "exceeds 0.125; No Way point raised from", f.noway, "to",
                 noway, "per the PMx rule for inbred populations\n"))
    }

    noway.mask <- f.off >= noway
    if (!is.null(unknown.breaks)) {
        noway.mask <- noway.mask | (unknown > unknown.breaks[4])
    }
    rated <- !noway.mask
    if (!any(rated) && verbose >= 1) {
        cat(warn("  Warning: every candidate pairing is 'NoWay'; no pairings rated\n"))
    }

    # Bin ranks 1-6 per component (PMx default bin definitions) ------------
    rank.dgd <- matrix(NA_integer_, nrow = n.m, ncol = n.f,
                       dimnames = list(males, females))
    rank.mkdiff <- rank.f <- rank.dgd
    msi.base <- matrix(NA_real_, nrow = n.m, ncol = n.f,
                       dimnames = list(males, females))
    dgd.lo <- dgd.hi <- mkd.lo <- mkd.hi <- mkd.bp <- f.bp <- NA_real_

    if (any(rated)) {
        # Range limits and break points from the rated pairings; the F break
        # point is the mean kinship over all candidate pairs (PMx default
        # 'F Break Point: Use Average F')
        dgd.lo <- min(dgd[rated])
        dgd.hi <- max(dgd[rated])
        mkd.lo <- min(mkdiff[rated])
        mkd.hi <- max(mkdiff[rated])
        mkd.bp <- mean(mkdiff[rated])
        f.bp <- f.mean
        if (f.bp <= 0 || f.bp >= noway) {
            f.bp <- noway / 2
            if (verbose >= 2) {
                cat(report("  F break point (mean pairwise kinship) outside (0, No Way); using No Way/2 =",
                           round(f.bp, 4), "\n"))
            }
        }

        # deltaGD: higher is better; break point 0; bins 1-3 thirds of the
        # range above 0, bins 4-6 thirds of the range below 0
        bin.dgd <- function(v) {
            if (v >= 0) {
                if (dgd.hi <= 0) return(3L)
                as.integer(min(floor(3 * (dgd.hi - v) / dgd.hi) + 1, 3))
            } else {
                if (dgd.lo >= 0) return(4L)
                as.integer(min(floor(3 * v / dgd.lo) + 4, 6))
            }
        }

        # MKdiff: lower is better; break point = mean of rated pairs; bins
        # 1-3 thirds of [best, break], bins 4-6 thirds of (break, worst]
        bin.mkd <- function(v) {
            if (v <= mkd.bp) {
                w <- mkd.bp - mkd.lo
                if (w <= 0) return(3L)
                as.integer(min(floor(3 * (v - mkd.lo) / w) + 1, 3))
            } else {
                w <- mkd.hi - mkd.bp
                if (w <= 0) return(4L)
                as.integer(min(floor(3 * (v - mkd.bp) / w) + 4, 6))
            }
        }

        # F: bin 1 pinned to F <= 0; bins 2-3 halves of (0, break]; bins 4-6
        # thirds of (break, No Way)
        bin.f <- function(v) {
            if (v <= 0) return(1L)
            if (v <= f.bp) {
                if (v <= f.bp / 2) 2L else 3L
            } else {
                as.integer(min(floor(3 * (v - f.bp) / (noway - f.bp)) + 4, 6))
            }
        }

        idx <- which(rated)
        rank.dgd[idx] <- vapply(dgd[idx], bin.dgd, integer(1))
        rank.mkdiff[idx] <- vapply(mkdiff[idx], bin.mkd, integer(1))
        rank.f[idx] <- vapply(f.off[idx], bin.f, integer(1))

        # Tulsa logic rules ---------------------------------
        tulsa <- function(r.gd, r.mk, r.f) {
            if (r.gd < 4 && r.mk < 4 && r.f < 4) {
                if (r.gd == r.mk && r.mk == r.f) return(as.numeric(r.gd))
                s <- r.gd + r.mk + r.f
                if (s < 5) return(1)
                if (s < 6 || r.gd == 1) return(2)
                return(3)
            }
            if (r.f == 6 || r.mk >= 5) return(6)
            if (r.gd <= 4) return(4)
            if (r.gd == 6) return(5)
            if (r.f == 5) return(5)
            4
        }
        msi.base[idx] <- mapply(tulsa,
                                rank.dgd[idx], rank.mkdiff[idx], rank.f[idx])
    }

    # Unknown-genome proration ('MSI Minimum' method); off when
    # unknown.breaks is NULL (the default)
    msi.num <- msi.base
    if (!is.null(unknown.breaks)) {
        pen <- rated & (unknown > unknown.breaks[1])
        msi.num[pen] <- pmax(msi.num[pen], 4)
        pen <- rated & (unknown > unknown.breaks[2])
        msi.num[pen] <- pmax(msi.num[pen], 5)
        pen <- rated & (unknown > unknown.breaks[3])
        msi.num[pen] <- 6
    }

    msi <- matrix(as.character(msi.num), nrow = n.m, ncol = n.f,
                  dimnames = list(males, females))
    msi[noway.mask] <- "NoWay"

# Printing outputs -----------
    if (verbose >= 3) {
        cat(report("  Mate Suitability Index summary (", n.m, "males x",
                   n.f, "females; No Way for F =", noway, "):\n"))
        print(table(factor(msi, levels = c(as.character(1:6), "NoWay"))))
    }

# FLAG SCRIPT END ---------------
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
# ----------------------

    # RETURN
    invisible(list(msi = msi,
                   f.off = f.off,
                   dgd = dgd,
                   mkdiff = mkdiff,
                   completeness = completeness,
                   rank.dgd = rank.dgd,
                   rank.mkdiff = rank.mkdiff,
                   rank.f = rank.f,
                   msi.base = msi.base,
                   settings = list(msi.method = "Tulsa",
                                   f.noway = noway,
                                   f.breakpoint = f.bp,
                                   dgd.range = c(lo = dgd.lo, hi = dgd.hi),
                                   mkdiff.range = c(lo = mkd.lo, hi = mkd.hi),
                                   mkdiff.breakpoint = mkd.bp,
                                   unknown.breaks = unknown.breaks)))
}
