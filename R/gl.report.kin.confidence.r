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
#' produced by gl.kin; used only for the verbose >= 3 diagnostic of how many
#' pairs are resolvable given their standard errors; if NULL and verbose >= 3,
#' computed internally with gl.kin. The standard errors are always those of
#' the gl.kin 'grm' (SNP) or 'dominant' (SilicoDArT) estimator, whatever
#' method produced kin; a kin from another method (e.g. 'emibd9') gives a
#' warning [default NULL].
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
#' columns, and kinship = G/2, as gl.kin returns for its 'grm' method. On the
#' full set of loci this reproduces gl.kin (via gl.grm and rrBLUP A.mat)
#' closely -- on testset2.gl off-diagonal values differ by at most 8e-5
#' against a spread (SD) of 0.037, the diagonal by at most 0.007 -- and it
#' is fast enough to run inside a bootstrap loop.
#'
#' For presence/absence (SilicoDArT) data, band columns are resampled with
#' replacement and gl.kin's 'dominant' estimator (utils.kin.dominant) is
#' recomputed on each resample: band covariance standardised by the band
#' frequencies of the ORIGINAL data, over the loci scored in both individuals,
#' divided by 2, with monomorphic loci excluded and the diagonal fixed at 0.5.
#'
#' Missing data: both estimators fill missing calls with the locus mean, which
#' pulls an individual's kinships toward 0 and also narrows their spread
#' across resamples. An individual with a low call rate therefore gets
#' SMALLER standard errors, although its kinships are the least reliable; the
#' bootstrap measures sampling error over loci, not this bias. The function
#' warns when any individual's call rate is below 0.8; filter individuals on
#' call rate (gl.filter.callrate(method = "ind")) before interpreting the
#' standard errors.
#'
#' Confidence-interval widths are computed from the bootstrap percentiles at
#' (1-conf)/2 and 1-(1-conf)/2. This requires holding all bootstrap kinship
#' matrices in memory (nInd x nInd x nboots doubles, about 60 MB for 274
#' individuals and 100 bootstraps); reduce nboots or subset individuals if
#' memory is limiting.
#'
#' The kin matrix itself (supplied, or computed via gl.kin when NULL at
#' verbose >= 3) plays no part in the standard errors; it is used only in the
#' verbose >= 3 summary to
#' report the proportion of pairs whose kinship differs from the dataset
#' baseline (median off-diagonal kinship) by more than two standard errors,
#' i.e. the pairs whose relatedness the data can actually resolve.

#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' # SNP data -- bootstrap on the FULL dataset (kinship references collapse on
#' # small family groups), reduced bootstraps for speed
#' res <- gl.report.kin.confidence(testset2.gl, nboots = 20)
#' res$summary
#' # SE for a known full-sib pair in the captive colony; the captive-bred
#' # individuals have call rates of 0.70-0.80, so the function warns that
#' # their SEs are underestimated (see Details)
#' res$se["CB_AB_01", "CB_AB_02"]
#' # Tag P/A data
#' res.gs <- gl.report.kin.confidence(testset2.gs, nboots = 20)
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
    utils.flag.start(func = funname, verbose = verbose)

    # CHECK DATATYPE
    datatype <- utils.check.datatype(x, verbose = verbose)

    # KINSHIP MATRIX -- used only by the verbose >= 3 summary, so it is
    # computed only then; a supplied matrix is always validated
    if (!is.null(kin)) {
        kin.method <- attr(kin, "method")
        boot.method <- if (datatype == "SNP") "grm" else "dominant"
        if (!is.null(kin.method) && kin.method != boot.method &&
            verbose >= 1) {
            cat(warn("  Warning: kin was estimated with method '", kin.method,
                     "', but the standard errors are those of the '",
                     boot.method, "' estimator\n", sep = ""))
        }
    }
    if (!is.null(kin) || verbose >= 3) {
        kin <- utils.kin.check(x, kin, verbose = verbose)
    }

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

    # Low call rates shrink both the kinships and their standard errors
    # (mean imputation), so the least reliable individuals look the most
    # precise
    ind.cr <- 1 - rowMeans(is.na(mat))
    if (any(ind.cr < 0.8) && verbose >= 1) {
        cat(warn(paste0("  Warning: ", sum(ind.cr < 0.8),
                        " individuals have call rate below 0.8 (lowest ",
                        round(min(ind.cr), 3),
                        "); their standard errors are underestimated. ",
                        "Consider gl.filter.callrate(method = 'ind') first\n")))
    }

    if (datatype == "SNP") {
        # DO THE JOB -- SNP data: light internal VanRaden GRM per resample
        p <- colMeans(mat, na.rm = TRUE) / 2
        # Loci with no called genotypes have undefined allele frequency and
        # would propagate NA into the resampling denominator; drop them
        called <- !is.na(p)
        if (!all(called)) {
            if (verbose >= 1) {
                cat(warn("  Warning:", sum(!called),
                         "loci with no called genotypes dropped from the bootstrap\n"))
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
            # kinship = G/2, exactly as gl.kin's grm method returns it
            boots[, , b] <- Gb / 2
            if (verbose >= 2 && b %% 25 == 0) {
                cat(report("  Completed", b, "of", nboots, "bootstraps\n"))
            }
        }
    } else {
        # DO THE JOB -- Tag P/A data: gl.kin's dominant estimator
        # (utils.kin.dominant) per resample, band frequencies from the
        # original data
        mu <- colMeans(mat, na.rm = TRUE)
        polym <- !is.na(mu) & mu > 0 & mu < 1
        if (sum(polym) == 0) {
            stop(error("Fatal Error: no loci with band frequency strictly between 0 and 1; kinship cannot be estimated\n"))
        }
        mat <- mat[, polym, drop = FALSE]
        mu <- mu[polym]
        nL <- ncol(mat)
        Z <- sweep(mat, 2, mu, "-")
        Z <- sweep(Z, 2, sqrt(mu * (1 - mu)), "/")
        obs <- !is.na(Z) * 1
        Z[is.na(Z)] <- 0

        for (b in seq_len(nboots)) {
            cols <- sample.int(nL, nL, replace = TRUE)
            # pairwise-complete mean of standardised cross-products
            den <- tcrossprod(obs[, cols, drop = FALSE])
            kb <- tcrossprod(Z[, cols, drop = FALSE]) / den / 2
            kb[den == 0] <- NA
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
        cat(report("  Bootstrap confidence of kinship estimates\n"))
        cat(report("    Individuals, loci, bootstraps:", nI, ",", nL, ",",
                   nboots, "\n"))
        cat(report("    Mean SE (off-diagonal)  :", round(sumry$mean.se, 5),
                   "\n"))
        cat(report("    Median SE               :", round(sumry$median.se, 5),
                   "\n"))
        cat(report("    90th percentile SE      :", round(sumry$q90.se, 5),
                   "\n"))
        cat(report("    Mean", conf, "CI width       :",
                   round(mean(ci.width[off], na.rm = TRUE), 5), "\n"))
        baseline <- stats::median(kin[row(kin) != col(kin)], na.rm = TRUE)
        resolvable <- mean(abs(kin[off] - baseline) > 2 * se[off], na.rm = TRUE)
        cat(report("    Pairs resolvable from baseline at 2 SE:",
                   round(100 * resolvable, 1), "%\n"))
    }

# FLAG SCRIPT END ---------------
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
# ----------------------

    # RETURN
    invisible(list(se = se, ci.width = ci.width, summary = sumry))
}
