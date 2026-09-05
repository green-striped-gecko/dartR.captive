#' @name gl.report.kin.confidence
#' @title Reports bootstrap confidence of genomic kinship estimates
#' @family captive management

#' @description
#' Estimates the sampling uncertainty of pairwise genomic kinship estimates by
#' bootstrapping over loci, reporting a per-pair standard error, a per-pair
#' confidence-interval width, and summary statistics of the standard errors.

#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param kin Kinship matrix with both dimnames identical to indNames(x), as
#' produced by gl.kin; used only for the verbose diagnostic of how many pairs
#' are resolvable given their standard errors; if NULL, computed internally
#' with gl.kin [default NULL].
#' @param nboots Number of bootstrap resamples of the loci [default 100].
#' @param conf Confidence level for the reported interval widths [default 0.95].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @details
#' The PMx pedigree-management software cautions that molecular estimates of
#' kinship should be weighted by how much they can be trusted, and notes that
#' quantifying that trust is unresolved (PMx Users Manual pp. 83, 106; Lacy,
#' Ballou & Pollak 2012). This function supplies that quantity for genomic
#' kinships: the dispersion of each pairwise estimate under resampling of loci,
#' which is the dominant source of sampling error when individuals are typed at
#' a finite number of markers.
#'
#' For SNP data, loci (columns of the dosage matrix) are resampled with
#' replacement nboots times and the kinship matrix is recomputed on each
#' resample using a light internal VanRaden (2008) genomic relationship matrix:
#' allele frequencies p are taken from column means of dosage/2 on the ORIGINAL
#' data, the dosage matrix is centred as Z = dosage - 2p (missing values mean
#' imputed to zero after centring), G = ZZ'/(2*sum(p(1-p))) over the resampled
#' columns, and kinship via the identical Goudet-style centring gl.kin applies
#' for its grm method: MS = mean(diag(G) - 1), off-diagonals G/2 - MS,
#' diagonal G/2 -- so the bootstrap distribution sits on the scale of the grm
#' method. This is a close approximation of the rrBLUP A.mat engine used by
#' gl.kin (which differs mainly in its treatment of missing data), and is
#' adequate for estimating the dispersion of the estimates while being fast
#' enough to run inside a bootstrap loop.
#'
#' For presence/absence (SilicoDArT) data, band columns are resampled with
#' replacement and the dominant-marker estimator recomputed on each resample:
#' kinship = correlation of band profiles between individuals divided by 2.
#'
#' Confidence-interval widths are computed from the bootstrap percentiles at
#' (1-conf)/2 and 1-(1-conf)/2. This requires holding all bootstrap kinship
#' matrices in memory (nInd x nInd x nboots doubles, about 60 MB for 274
#' individuals and 100 bootstraps); reduce nboots or subset individuals if
#' memory is limiting.
#'
#' The kin matrix itself (supplied, or computed via gl.kin when NULL) plays no
#' part in the standard errors; it is used only in the verbose >= 3 summary to
#' report the proportion of pairs whose kinship differs from the dataset
#' baseline (median off-diagonal kinship) by more than two standard errors,
#' i.e. the pairs whose relatedness the data can actually resolve.

#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' # SNP data -- bootstrap on the FULL dataset (kinship references collapse on
#' # small family groups), reduced bootstraps for speed
#' res <- gl.report.kin.confidence(testset.gl, nboots = 20)
#' res$summary
#' # SE for a known full-sib pair in the captive colony
#' res$se["CB_AB_01", "CB_AB_02"]
#' # Tag P/A data
#' res.gs <- gl.report.kin.confidence(testset.gs, nboots = 20)
#' res.gs$summary

#' @seealso \code{\link{gl.kin}}, \code{\link{gl.report.kin.classes}}

#' @export
#' @return Invisibly, a list with three components: se, a matrix of per-pair
#' bootstrap standard deviations of the kinship estimates; ci.width, a matrix
#' of bootstrap-percentile confidence-interval widths at level conf; and
#' summary, a one-row data frame with columns mean.se, median.se and q90.se
#' computed over the off-diagonal pairs.
#'
# ----------------------
# Function
gl.report.kin.confidence <- function(x,
                                     kin = NULL,
                                     nboots = 100,
                                     conf = 0.95,
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

    # FUNCTION SPECIFIC ERROR CHECKING
    if (!is.numeric(nboots) || length(nboots) != 1 || nboots < 2) {
        stop(error("Fatal Error: nboots must be a single number of 2 or more\n"))
    }
    nboots <- as.integer(round(nboots))
    if (!is.numeric(conf) || length(conf) != 1 || conf <= 0 || conf >= 1) {
        stop(error("Fatal Error: conf must be a single value in (0, 1)\n"))
    }

    nI <- nInd(x)
    nL <- nLoc(x)
    ids <- indNames(x)

# DO THE JOB ----------------------
    if (verbose >= 2) {
        cat(report("  Bootstrapping kinship over", nL, "loci,",
                   nboots, "resamples\n"))
    }

    boots <- array(NA_real_, dim = c(nI, nI, nboots))
    mat <- as.matrix(x)

    if (datatype == "SNP") {
        # DO THE JOB -- SNP data: light internal VanRaden GRM per resample
        p <- colMeans(mat, na.rm = TRUE) / 2
        # Loci with no called genotypes have undefined allele frequency and
        # would propagate NA into the resampling denominator; drop them
        called <- !is.na(p)
        if (!all(called)) {
            if (verbose >= 1) {
                cat(warn("  Warning:", sum(!called),
                    "loci with no called genotypes dropped from the bootstrap
"))
            }
            mat <- mat[, called, drop = FALSE]
            p <- p[called]
            nL <- ncol(mat)
        }
        Z <- sweep(mat, 2, 2 * p)
        Z[is.na(Z)] <- 0  # mean imputation after centring

        for (b in seq_len(nboots)) {
            cols <- sample.int(nL, nL, replace = TRUE)
            denom <- 2 * sum(p[cols] * (1 - p[cols]))
            if (denom <= 0) {
                stop(error(
                    "Fatal Error: Resampled loci are all monomorphic; kinship undefined\n"
                ))
            }
            Gb <- tcrossprod(Z[, cols, drop = FALSE]) / denom
            # Centring copied exactly from gl.kin's grm method (Goudet et al.
            # 2018 as in gl.grm.network): MS = mean(diag(G) - 1) is the
            # reference, off-diagonals G/2 - MS, diagonal G/2 -- so bootstrap
            # SEs are on the same scale as the estimator users see
            MS.b <- mean(diag(Gb) - 1)
            kb <- Gb / 2 - MS.b
            diag(kb) <- diag(Gb) / 2
            boots[, , b] <- kb
            if (verbose >= 2 && b %% 25 == 0) {
                cat(report("  Completed", b, "of", nboots, "bootstraps\n"))
            }
        }
    } else {
        # DO THE JOB -- Tag P/A data: dominant estimator per resample
        for (b in seq_len(nboots)) {
            cols <- sample.int(nL, nL, replace = TRUE)
            kb <- suppressWarnings(
                stats::cor(t(mat[, cols, drop = FALSE]),
                           use = "pairwise.complete.obs")
            ) / 2
            diag(kb) <- 0.5
            boots[, , b] <- kb
            if (verbose >= 2 && b %% 25 == 0) {
                cat(report("  Completed", b, "of", nboots, "bootstraps\n"))
            }
        }
    }

    # Per-pair dispersion and interval widths
    if (verbose >= 2) {
        cat(report("  Summarising per-pair standard errors and interval widths\n"))
    }
    se <- apply(boots, c(1, 2), function(v) stats::sd(v, na.rm = TRUE))
    alpha <- (1 - conf) / 2
    qs <- apply(boots, c(1, 2), function(v) {
        v <- v[!is.na(v)]
        if (length(v) < 2) return(c(NA_real_, NA_real_))
        stats::quantile(v, probs = c(alpha, 1 - alpha), names = FALSE)
    })
    ci.width <- qs[2, , ] - qs[1, , ]
    dimnames(se) <- list(ids, ids)
    dimnames(ci.width) <- list(ids, ids)

    off <- row(se) != col(se)
    sumry <- data.frame(
        mean.se = mean(se[off], na.rm = TRUE),
        median.se = stats::median(se[off], na.rm = TRUE),
        q90.se = stats::quantile(se[off], probs = 0.90, na.rm = TRUE, names = FALSE)
    )

    # Printing outputs -----------
    if (verbose >= 3) {
        cat("  Bootstrap confidence of kinship estimates\n")
        cat(paste("    Individuals, loci, bootstraps:", nI, ",", nL, ",", nboots), "\n")
        cat(paste("    Mean SE (off-diagonal)  :", round(sumry$mean.se, 5)), "\n")
        cat(paste("    Median SE               :", round(sumry$median.se, 5)), "\n")
        cat(paste("    90th percentile SE      :", round(sumry$q90.se, 5)), "\n")
        cat(paste("    Mean", conf, "CI width       :",
                  round(mean(ci.width[off], na.rm = TRUE), 5)), "\n")
        baseline <- stats::median(kin[row(kin) != col(kin)])
        resolvable <- mean(abs(kin[off] - baseline) > 2 * se[off], na.rm = TRUE)
        cat(paste("    Pairs resolvable from baseline at 2 SE:",
                  round(100 * resolvable, 1), "%"), "\n")
    }

# FLAG SCRIPT END ---------------
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
# ----------------------

    # RETURN
    invisible(list(se = se, ci.width = ci.width, summary = sumry))
}
