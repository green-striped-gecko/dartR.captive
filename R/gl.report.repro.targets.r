#' @name gl.report.repro.targets
#' @title Reports per-individual reproductive targets to equalise founder representation
#' @family captive management
#' @description
#' Calculates a lifetime offspring target for each sexed individual in a
#' genlight object that lowers the expected mean kinship of the next
#' generation: individuals with low mean kinship (carrying under-represented
#' genomes) receive high targets, individuals with high mean kinship receive
#' low targets, and each sex's targets sum to the total number of offspring
#' sought. This is a SNP-based analogue of the Repro Goals screen
#' in the PMx studbook-management software.
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param kin Kinship matrix estimated on a wider set of individuals than x,
#' for example gl.kin() on a dataset that includes the source populations;
#' it may cover more individuals than x and is restricted to indNames(x).
#' Kinship estimated on x alone (including kin = NULL) is an error, because
#' mean kinship over the individuals it was estimated on is 0 by
#' construction [required].
#' @param n.target Total number of offspring sought in the next generation;
#' each sex's individual targets sum to this value [default NULL, resolving to
#' nInd(x), i.e. one-for-one replacement of the current population].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' brief progress messages; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' Kinship should be estimated on the most inclusive dataset available and the
#' matrix subset to the managed group; computing kinship on a small family
#' group alone inflates estimates (the allele-frequency reference collapses
#' onto the family itself).
#' This function mirrors the Repro Goals capability of PMx (Lacy, Ballou &
#' Pollak 2012, Methods in Ecology and Evolution 3:433-437; PMx Users Manual
#' v1.0 pp. 90-91), which recommends how many offspring each individual should
#' contribute so that genome representation is equalised across founders in the
#' next generation. PMx derives its goals by an iterative simulation of pairing
#' and culling; here the same intent is implemented as a direct
#' mean-kinship-proportional allocation, which is closed-form and
#' deterministic.
#'
#' Mean kinship MK_i is the row mean of the kinship matrix (self included). An
#' individual's target is proportional to (max(MK) - MK_i + eps), where max(MK)
#' is taken within the individual's sex pool and eps = 0.001 keeps every weight
#' positive: the higher an individual's mean kinship (the more over-represented
#' its genome), the fewer offspring it should contribute. Because every
#' offspring requires exactly one sire and one dam, the weights are normalised
#' within each sex so that the male targets and the female targets each sum to
#' n.target. Fractional targets are integerised by the largest-remainder
#' method, so each sex's targets sum exactly to n.target; remainder ties are
#' broken in favour of lower mean kinship, then by individual id.
#'
#' The allocation uses each individual's mean kinship only, not the kinship
#' between the individuals it favours, so it does not reach the optimum. On
#' the testset2.gl captive colony (unfiltered, n.target = nInd(x)), the
#' expected mean kinship of the next generation is 0.0582, against 0.0667
#' for equal contributions and 0.0560 for optimal (non-negative)
#' contributions.
#'
#' Missing genotypes change the targets. gl.kin fills missing genotypes with
#' the locus mean, which pulls kinship toward 0 for individuals with low call
#' rates. On the same colony (call rates 0.70-0.80), filtering loci with
#' gl.filter.callrate(method = "loc", threshold = 0.95) changes the target
#' of 9 of 24 individuals (one goes from 2 offspring to 0). The function
#' warns when any individual has call rate below 0.8; filter on call rate
#' before gl.kin. Missing (NA) kinship values are ignored in the mean
#' kinship, with a warning.
#'
#' Sexes are read from x@@other$ind.metrics$sex (values 'Male', 'Female',
#' 'Unknown'); individuals of any other or missing sex are excluded (with a
#' warning) and
#' receive no target, although the default n.target = nInd(x) still counts them
#' as part of the population being replaced.
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # Filter on call rate, estimate kinship on the FULL dataset, and pass it
#' # for the captive colony (it is restricted to the colony's individuals)
#' tf <- gl.filter.callrate(testset2.gl, method = "loc", threshold = 0.95,
#'                          verbose = 0)
#' kin <- gl.kin(tf)
#' cb <- gl.keep.pop(tf, pop.list = "EmmacCaptBred", verbose = 0)
#' tg <- gl.report.repro.targets(cb, kin = kin, n.target = 12)
#' head(tg[order(-tg$target), ])
#' @seealso \code{\link{gl.kin}}, \code{\link{gl.select.pairs}},
#' \code{\link{gl.report.mate.suitability}}
#' @export
#' @return Invisibly, a data frame with one row per sexed individual: $id,
#' $pop, $sex, $MK (mean kinship, self included) and $target (integer offspring
#' target; targets sum to n.target within each sex).
#'
# ----------------------
# Function
gl.report.repro.targets <- function(x,
                                    kin = NULL,
                                    n.target = NULL,
                                    verbose = NULL) {
# PRELIMINARIES -- checking ----------------
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)

    # FLAG SCRIPT START
    funname <- match.call()[[1]]
    utils.flag.start(func = funname, verbose = verbose)

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)

# FUNCTION SPECIFIC ERROR CHECKING ----------------
    # Kinship matrix
    kin <- utils.kin.check(x, kin, verbose = verbose,
                           need.reference = TRUE)

    # Sexes
    if (is.null(x@other$ind.metrics$sex)) {
        stop(error(
            "Fatal Error: reproductive targets require a 'sex' column (values 'Male', 'Female', 'Unknown') in x@other$ind.metrics; none found\n"
        ))
    }
    ids <- indNames(x)
    sex <- as.character(x@other$ind.metrics$sex)
    # %in% rather than ==, so that an NA sex counts as unknown
    males <- ids[sex %in% "Male"]
    females <- ids[sex %in% "Female"]
    n.unknown <- sum(!(sex %in% c("Male", "Female")))
    if (n.unknown > 0 && verbose >= 1) {
        cat(warn("  Warning:", n.unknown,
                 "individuals of unknown sex excluded; they receive no target\n"))
    }
    if (length(males) == 0 || length(females) == 0) {
        stop(error(
            "Fatal Error: at least one Male and one Female are required to allocate reproductive targets\n"
        ))
    }

    # n.target
    if (is.null(n.target)) {
        n.target <- nInd(x)
    } else if (!is.numeric(n.target) || length(n.target) != 1 ||
               is.na(n.target) || n.target < 1) {
        stop(error(
            "Fatal Error: n.target must be a single number >= 1, or NULL for nInd(x)\n"
        ))
    }
    if (n.target != floor(n.target)) {
        if (verbose >= 1) {
            cat(warn("  Warning: n.target rounded down to", floor(n.target), "\n"))
        }
        n.target <- floor(n.target)
    }

# DO THE JOB ----------------------
    if (verbose >= 2) {
        cat(report("  Allocating", n.target, "offspring across",
                   length(males), "males and", length(females),
                   "females by mean-kinship-proportional weights\n"))
    }
    # Individual call rates; gl.kin fills missing genotypes with the locus
    # mean, which pulls kinship toward 0 for individuals with low call rates
    cr <- 1 - vapply(x@gen, function(e) length(e@NA.posi), numeric(1)) /
        nLoc(x)
    if (any(cr < 0.8) && verbose >= 1) {
        cat(warn(paste0("  Warning: ", sum(cr < 0.8),
                        " individuals have call rate below 0.8 (lowest ",
                        round(min(cr), 3), "); missing genotypes pull ",
                        "their kinship toward 0 and change the targets. ",
                        "Consider filtering on call rate before gl.kin ",
                        "(see Details)\n")))
    }
    n.na <- sum(is.na(kin[upper.tri(kin, diag = TRUE)]))
    if (n.na > 0 && verbose >= 1) {
        cat(warn("  Warning:", n.na, "missing (NA) kinship values ignored",
                 "in the mean kinship\n"))
    }
    mk <- rowMeans(kin, na.rm = TRUE)

    # Mean-kinship-proportional allocation with largest-remainder rounding;
    # ties broken toward lower MK, then id, for determinism
    allocate <- function(mk.pool, n.target) {
        eps <- 0.001
        w <- max(mk.pool) - mk.pool + eps
        raw <- w / sum(w) * n.target
        tgt <- floor(raw)
        rem <- n.target - sum(tgt)
        if (rem > 0) {
            frac <- raw - tgt
            ord <- order(-frac, mk.pool, names(mk.pool))
            tgt[ord[seq_len(rem)]] <- tgt[ord[seq_len(rem)]] + 1
        }
        tgt
    }
    tgt.m <- allocate(mk[males], n.target)
    tgt.f <- allocate(mk[females], n.target)
    target <- c(tgt.m, tgt.f)[c(males, females)]

    keep <- ids %in% c(males, females)
    out <- data.frame(
        id = ids[keep],
        pop = as.character(pop(x))[keep],
        sex = sex[keep],
        MK = mk[ids[keep]],
        target = target[ids[keep]],
        stringsAsFactors = FALSE
    )
    rownames(out) <- NULL

# Printing outputs -----------
    if (verbose >= 3) {
        cat(report("  Reproductive targets (sorted by target, descending):\n"))
        print(out[order(-out$target, out$MK, out$id), ], row.names = FALSE)
    }

# FLAG SCRIPT END ---------------
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
# ----------------------

    # RETURN
    invisible(out)
}
