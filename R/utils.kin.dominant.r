#' @name utils.kin.dominant
#' @title Estimates a kinship matrix from dominant (presence/absence) data
#' @family captive management
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' Estimates pairwise kinship for presence/absence (SilicoDArT) data from the
#' standardized covariance of band states, serving as the dominant-marker
#' engine behind gl.kin(method='dominant'). The returned matrix follows the
#' kinship contract of the captive management series: dimnames = indNames(x),
#' off-diagonal = pairwise kinship estimate, diagonal fixed at 0.5.
#' @param x Name of the genlight object containing the presence/absence
#' (SilicoDArT) data [required].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#' @details
#' For each locus l, the band frequency mu_l is computed over all individuals
#' (missing values removed). Only loci with 0 < mu_l < 1 are informative and
#' retained. For a pair of individuals i and j, the standardized band
#' covariance is
#'
#' r_ij = mean over pairwise-complete loci of
#'        (x_il - mu_l)(x_jl - mu_l) / (mu_l (1 - mu_l))
#'
#' and the kinship estimate is r_ij / 2. The computation is done by matrix
#' algebra: the band matrix is centred and scaled by column, missing entries
#' are zero-filled so cross-products accumulate only over shared scored loci,
#' and the denominators are the cross-products of the not-missing indicator
#' (pairwise-complete locus counts). Pairs of individuals sharing no scored
#' loci receive NA with a warning.
#'
#' The diagonal is fixed at 0.5, the self-kinship of a non-inbred individual,
#' because individual inbreeding cannot be estimated from single dominant
#' bands (a band-present phenotype does not distinguish homozygous from
#' heterozygous carriers).
#'
#' This estimator is a rank-valid approximation. Dominant markers
#' underestimate absolute kinship by a factor that depends on the band
#' frequency spectrum of the loci (cf. Lynch & Milligan 1994, Molecular
#' Ecology 3:91-99): rankings of individuals or pairs and relative comparisons
#' (for example, choosing the lower-kinship pairing) are valid, but absolute
#' values should not be compared against SNP-based or pedigree-based kinship
#' thresholds. As with all methods in this series, estimates are relative to
#' the band frequencies of the full dataset supplied; datasets composed only
#' of close relatives will yield inflated relative values.
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # Examples for testing
#' # kin <- utils.kin.dominant(testset.gs)
#' # kin[1:5, 1:5]
#' # 1 - mean(kin)  # gene diversity over the full matrix
#' @seealso \code{\link{gl.kin}}, \code{\link{gl.report.kinship}}
#' @keywords internal
# @export
#' @return A square numeric kinship matrix with dimnames = indNames(x),
#' off-diagonal r_ij/2 and diagonal 0.5, returned invisibly.
#'
# ----------------------
# Function
utils.kin.dominant <- function(x,
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
    datatype <- utils.check.datatype(x,
                                     accept = "SilicoDArT",
                                     verbose = verbose)

    # DO THE JOB ----------------------

    mat <- as.matrix(x)

    # Per-locus band frequency over all individuals
    mu <- colMeans(mat, na.rm = TRUE)
    polym <- !is.na(mu) & mu > 0 & mu < 1
    if (sum(polym) == 0) {
        stop(error(paste0("Fatal Error: no loci with band frequency strictly ",
                          "between 0 and 1; kinship cannot be estimated\n")))
    }
    if (verbose >= 2) {
        cat(report("  Using", sum(polym),
                   "loci with 0 < band frequency < 1 (",
                   sum(!polym), "monomorphic or all-missing loci excluded)\n"))
    }
    mat <- mat[, polym, drop = FALSE]
    mu <- mu[polym]

    # Centre and scale by column so that the cross-product of two individuals
    # at locus l is (x_il - mu_l)(x_jl - mu_l) / (mu_l (1 - mu_l))
    Z <- sweep(mat, 2, mu, "-")
    Z <- sweep(Z, 2, sqrt(mu * (1 - mu)), "/")

    # Pairwise-complete means via cross-products on a zero-filled matrix;
    # denominators from the cross-product of the not-missing indicator
    obs <- !is.na(Z)
    Z0 <- Z
    Z0[!obs] <- 0
    num <- tcrossprod(Z0)        # sums of scaled cross-products over shared loci
    den <- tcrossprod(obs * 1L)  # pairwise-complete locus counts
    r <- num / den               # standardized band covariance r_ij

    if (any(den == 0)) {
        r[den == 0] <- NA
        if (verbose >= 1) {
            cat(warn(paste0("  Warning: some pairs of individuals share no ",
                            "scored loci; their kinship is set to NA\n")))
        }
    }

    kin <- r / 2
    # Individual inbreeding is not estimable from single dominant bands
    diag(kin) <- 0.5
    dimnames(kin) <- list(indNames(x), indNames(x))

    # Results summary -----------
    if (verbose >= 3) {
        off <- kin[upper.tri(kin)]
        cat(report("  Dominant-marker kinship summary\n"))
        cat(report("    Individuals:", nInd(x), "| Informative loci:",
                   sum(polym), "\n"))
        cat(report("    Mean off-diagonal kinship:",
                   round(mean(off, na.rm = TRUE), 4), "\n"))
        cat(report("    Off-diagonal range:",
                   round(min(off, na.rm = TRUE), 4), "to",
                   round(max(off, na.rm = TRUE), 4), "\n"))
    }

    # FLAG SCRIPT END ---------------
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    # ----------------------

    # RETURN
    invisible(kin)
}
