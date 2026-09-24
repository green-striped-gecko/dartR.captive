#' @name gl.report.kin.groups
#' @title Reports group-level kinship, mean kinship and gene diversity
#' @family captive management

#' @description
#' Aggregates a pairwise genomic kinship matrix to management groups, reporting
#' the group-by-group kinship matrix, each group's size-weighted mean kinship,
#' mean inbreeding of its members, and the group-level gene diversity of the
#' whole managed population.

#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param group.col Name of a column in x@@other$ind.metrics defining the
#' management group of each individual; if NULL, groups are taken from pop(x)
#' [default NULL].
#' @param kin Kinship matrix estimated on a wider set of individuals than x,
#' for example gl.kin() on a dataset that includes the source populations;
#' it may cover more individuals than x and is restricted to indNames(x).
#' Kinship estimated on x alone (including kin = NULL) is an error, because
#' group mean kinship over the individuals it was estimated on is 0 and gene
#' diversity 1 by construction [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @details
#' Kinship must be estimated on a wider reference than x, for example the most
#' inclusive dataset available, and is restricted to the managed individuals.
#' gl.kin centres kinship on the individuals it is estimated on, so over
#' those individuals every row of the matrix averages 0: every group MK would
#' be 0 and GD 1, whatever the data. Such a matrix (including kin = NULL) is
#' rejected with an error.
#'
#' Missing genotypes pull kinship and self-kinship toward 0, because gl.kin
#' fills them with the locus mean. Groups of individuals with low call rates
#' therefore show negative meanF, lower MK and higher GD: on the testset2.gl
#' captive cohorts (call rates 0.70-0.80) meanF is -0.25 to -0.38, and -0.15
#' to +0.01 after gl.filter.callrate(method = "loc", threshold = 0.95). The
#' function warns when any individual's call rate is below 0.8; filter on
#' call rate before gl.kin. Pairs with missing (NA) kinship are ignored in
#' the block means, with a warning.
#'
#' Many species are managed as groups (herds, flocks, tanks, enclosures) within
#' which parentage is unobserved, so individual-level pedigree management as
#' implemented in PMx (Lacy, Ballou & Pollak 2012) is not possible.
#' Jimenez-Mena et al. (2016) extended mean-kinship management to populations
#' managed as groups, treating groups as the units of management (their
#' MERGE/SPLIT/EXTRACT formalism) and propagating group kinships through a
#' group pedigree. This function is the genomic analogue: the group-pedigree
#' recursion is replaced by direct genomic estimates, so no group pedigree is
#' required.
#'
#' Groups are defined by the nominated ind.metrics column (fatal error if the
#' column is absent), or by pop(x) when group.col is NULL. Individuals with a
#' missing group value are dropped with a warning.
#'
#' The group kinship matrix is the block average of the individual kinship
#' matrix: f[g,h] = mean(kin[members of g, members of h]) for g != h, and the
#' group self-kinship f[g,g] is the mean of the group's block INCLUDING the
#' diagonal (so it reflects both within-group relatedness and member
#' inbreeding). The mean kinship of group g is the size-weighted average over
#' all groups h (including g itself), MK_g = sum_h n_h * f[g,h] / sum_h n_h,
#' and the group-level gene diversity is GD = 1 - sum_g n_g * MK_g / sum_g n_g.
#' Because the weights telescope, this GD equals 1 - mean(kin) over the full
#' individual kinship matrix -- consistent with the gene diversity reported by
#' the individual-level functions in this family.
#'
#' The mean inbreeding coefficient reported per group is derived from the
#' kinship diagonal, F_i = 2 * kin[i,i] - 1, averaged over group members.
#' Under the dominant (SilicoDArT) kinship estimator the diagonal is fixed at
#' 0.5, so meanF is 0 by construction.

#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' # Captive cohorts as management groups; kinship estimated on the FULL
#' # dataset and subset to the colony (see Details on the reference population)
#' kin <- gl.kin(testset2.gl)
#' cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
#' res <- gl.report.kin.groups(cb, group.col = "cohort",
#'                             kin = kin[indNames(cb), indNames(cb)])
#' res$groups
#' res$kin.groups  # F1_AB vs F1_AE elevated (shared sire)

#' @seealso \code{\link{gl.kin}}, \code{\link{gl.report.gd.projection}}

#' @export
#' @return Invisibly, a list with three components: groups, a data frame with
#' columns group, n, MK (size-weighted group mean kinship) and meanF (mean
#' inbreeding of members; 0 by construction for SilicoDArT data); kin.groups, the group-by-group kinship matrix; and
#' gd, the group-level gene diversity of the managed population.
#'
# ----------------------
# Function
gl.report.kin.groups <- function(x,
                                 group.col = NULL,
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

    # KINSHIP MATRIX -- MK and GD are means over all individuals, so a
    # self-referenced matrix (rows averaging 0) would give MK 0 and GD 1
    kin <- utils.kin.check(x, kin, verbose = verbose,
                           need.reference = TRUE)

    # Mean imputation of missing genotypes (gl.kin) pulls kinship and
    # self-kinship toward 0 for individuals with low call rates
    ind.cr <- 1 - vapply(x@gen, function(e) length(e@NA.posi), numeric(1)) /
        nLoc(x)
    if (any(ind.cr < 0.8) && verbose >= 1) {
        cat(warn(paste0("  Warning: ", sum(ind.cr < 0.8),
                        " individuals have call rate below 0.8 (lowest ",
                        round(min(ind.cr), 3),
                        "); missing genotypes pull their kinship and ",
                        "inbreeding toward 0. Consider filtering on call ",
                        "rate before gl.kin (see Details)\n")))
    }

    # FUNCTION SPECIFIC ERROR CHECKING -- group membership
    if (is.null(group.col)) {
        if (verbose >= 2) {
            cat(report("  No group column nominated; using pop(x) as groups\n"))
        }
        grp <- as.character(pop(x))
    } else {
        im <- x@other$ind.metrics
        if (is.null(im) || !(group.col %in% names(im))) {
            stop(error(paste("Fatal Error: Column", group.col,
                             "not found in x@other$ind.metrics\n")))
        }
        grp <- as.character(im[[group.col]])
    }

    keep <- !is.na(grp) & grp != ""
    if (sum(keep) == 0) {
        stop(error("Fatal Error: No individuals with a non-missing group value\n"))
    }
    if (any(!keep)) {
        if (verbose >= 1) {
            cat(warn("  Warning:", sum(!keep),
                     "individual(s) with missing group value dropped\n"))
        }
        kin <- kin[keep, keep, drop = FALSE]
        grp <- grp[keep]
    }

    grp <- factor(grp)
    levs <- levels(grp)
    G <- length(levs)
    n.g <- as.integer(table(grp)[levs])
    N <- sum(n.g)

# DO THE JOB ----------------------
    if (verbose >= 2) {
        cat(report("  Aggregating kinship over", G, "groups,", N, "individuals\n"))
    }

    n.na <- sum(is.na(kin[upper.tri(kin)]))
    if (n.na > 0 && verbose >= 1) {
        cat(warn("  Warning:", n.na,
                 "pairs have missing kinship and are ignored in the block means\n"))
    }

    # Group kinship matrix: block means of the individual kinship matrix
    # (self blocks include the diagonal)
    f <- matrix(NA_real_, nrow = G, ncol = G, dimnames = list(levs, levs))
    idx <- split(seq_along(grp), grp)
    for (g in seq_len(G)) {
        for (h in seq_len(G)) {
            f[g, h] <- mean(kin[idx[[levs[g]]], idx[[levs[h]]], drop = FALSE],
                            na.rm = TRUE)
        }
    }

    # Group mean kinship (size-weighted over all groups incl. self) and GD
    MK <- as.vector(f %*% (n.g / N))
    gd <- 1 - sum(n.g * MK) / N

    # Mean inbreeding per group from the kinship diagonal
    Fi <- 2 * diag(kin) - 1
    meanF <- as.numeric(tapply(Fi, grp, mean, na.rm = TRUE)[levs])

    groups <- data.frame(group = levs,
                         n = n.g,
                         MK = MK,
                         meanF = meanF,
                         stringsAsFactors = FALSE)

    # Printing outputs -----------
    if (verbose >= 3) {
        cat(report("  Group-level kinship summary\n"))
        cat(report("    Groups            :", G, "\n"))
        cat(report("    Individuals       :", N, "\n"))
        cat(report("    Gene diversity    :", round(gd, 4), "\n"))
        cat(report("  Per-group mean kinship and inbreeding\n"))
        print(data.frame(group = groups$group,
                         n = groups$n,
                         MK = round(groups$MK, 4),
                         meanF = round(groups$meanF, 4)))
    }

# FLAG SCRIPT END ---------------
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
# ----------------------

    # RETURN
    invisible(list(groups = groups, kin.groups = f, gd = gd))
}
