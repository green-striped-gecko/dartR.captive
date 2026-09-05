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
#' @param kin Kinship matrix with both dimnames identical to indNames(x), as
#' produced by gl.kin; if NULL, computed internally with gl.kin [default NULL].
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

#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' # Captive cohorts as management groups; kinship estimated on the FULL
#' # dataset and subset to the colony (see Details on the reference population)
#' kin <- gl.kin(testset.gl)
#' cb <- gl.keep.pop(testset.gl, pop.list = "EmmacCaptBred", verbose = 0)
#' res <- gl.report.kin.groups(cb, group.col = "cohort",
#'                             kin = kin[indNames(cb), indNames(cb)])
#' res$groups
#' res$kin.groups  # F1_AB vs F1_AE elevated (shared sire)

#' @seealso \code{\link{gl.kin}}, \code{\link{gl.report.gd.projection}}

#' @export
#' @return Invisibly, a list with three components: groups, a data frame with
#' columns group, n, MK (size-weighted group mean kinship) and meanF (mean
#' inbreeding of members); kin.groups, the group-by-group kinship matrix; and
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
    utils.flag.start(func = funname,
                     build = "v.2026.1",
                     verbose = verbose)

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)

    # KINSHIP MATRIX
    kin <- utils.kin.check(x, kin, verbose = verbose)

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

    # Group kinship matrix: block means of the individual kinship matrix
    # (self blocks include the diagonal)
    f <- matrix(NA_real_, nrow = G, ncol = G, dimnames = list(levs, levs))
    idx <- split(seq_along(grp), grp)
    for (g in seq_len(G)) {
        for (h in seq_len(G)) {
            f[g, h] <- mean(kin[idx[[levs[g]]], idx[[levs[h]]], drop = FALSE])
        }
    }

    # Group mean kinship (size-weighted over all groups incl. self) and GD
    MK <- as.vector(f %*% (n.g / N))
    gd <- 1 - sum(n.g * MK) / N

    # Mean inbreeding per group from the kinship diagonal
    Fi <- 2 * diag(kin) - 1
    meanF <- as.numeric(tapply(Fi, grp, mean)[levs])

    groups <- data.frame(group = levs,
                         n = n.g,
                         MK = MK,
                         meanF = meanF,
                         stringsAsFactors = FALSE)

    # Printing outputs -----------
    if (verbose >= 3) {
        cat("  Group-level kinship summary\n")
        cat(paste("    Groups            :", G), "\n")
        cat(paste("    Individuals       :", N), "\n")
        cat(paste("    Gene diversity    :", round(gd, 4)), "\n")
        cat("  Per-group mean kinship and inbreeding\n")
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
