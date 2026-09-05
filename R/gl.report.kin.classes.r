#' @name gl.report.kin.classes
#' @title Classifies pairwise relationships from genomic kinship
#' @family captive management

#' @description
#' Assigns each pair of individuals to a relationship class (parent-offspring,
#' full-sib, second-degree, third-degree, unrelated) from baseline-adjusted
#' genomic kinship, and cross-checks the genomic classes against any recorded
#' sire/dam links in the individual metadata.

#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param kin Kinship matrix with both dimnames identical to indNames(x), as
#' produced by gl.kin; if NULL, computed internally with gl.kin [default NULL].
#' @param oh.thresh Maximum opposite-homozygote rate for a first-degree pair to
#' be classed as parent-offspring rather than full-sib; allows for genotyping
#' error [default 0.005].
#' @param all.pairs If TRUE, the returned pairs table includes pairs classified
#' as unrelated; if FALSE, only related pairs are returned [default FALSE].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @details
#' The PMx pedigree-management software (Lacy, Ballou & Pollak 2012) reasons
#' about individuals in terms of categorical pedigree relationships; this
#' function recovers those categories directly from genomic kinship, which also
#' provides the pedigree-verification capability PMx delegates to external
#' molecular analyses (PMx Users Manual pp. 83, 106).
#'
#' Kinship estimates from frequency-based estimators are RELATIVE to the sample
#' used as the allele-frequency reference: the average pair in the sample sits
#' near zero (or a common baseline) regardless of its true kinship. Before
#' classification, pairwise kinships are therefore adjusted by subtracting the
#' dataset median off-diagonal kinship, which anchors the "unrelated" class at
#' zero provided the majority of pairs in the dataset are indeed unrelated. If
#' most pairs are close relatives (e.g. a single family), the baseline is
#' biased and classes will be shifted downwards -- include unrelated
#' individuals in the dataset where possible.
#'
#' Adjusted kinships are assigned to the nearest of the theoretical class
#' values 0.25 (first-degree), 0.125 (second-degree), 0.0625 (third-degree)
#' and 0 (unrelated), using the midpoints between adjacent values as
#' thresholds (0.1875, 0.09375, 0.03125).
#'
#' For SNP data, first-degree pairs are further split into parent-offspring
#' and full-sib using opposite homozygotes: a true parent and offspring cannot
#' be opposite homozygotes (one scored 0, the other 2) at any locus, barring
#' genotyping error, whereas full sibs are expected to be opposite homozygotes
#' at a proportion of loci. A first-degree pair is classed parent-offspring if
#' its opposite-homozygote rate (count of opposite-homozygote loci over loci
#' called in both individuals) is below oh.thresh, otherwise full-sib. For
#' presence/absence (SilicoDArT) data homozygotes cannot be distinguished, so
#' this split is not possible and first-degree pairs retain the class
#' "first-degree".
#'
#' If the individual metadata (x@@other$ind.metrics) contains sire and dam
#' columns, each recorded parent-offspring link for which both individuals are
#' present in the data is cross-checked against the genomic class. Links whose
#' genomic class is not parent-offspring (SNP) or first-degree (SilicoDArT)
#' are returned in the conflicts table with their kinship and genomic class,
#' flagging possible pedigree errors, sample mix-ups, or genotyping problems.

#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' # SNP data -- run on the FULL dataset so the median baseline is set by
#' # unrelated pairs (a family-only matrix biases it; see Details)
#' kin <- gl.kin(testset.gl)
#' res <- gl.report.kin.classes(testset.gl, kin = kin)
#' table(res$pairs$class)
#' # The known captive full sibs classify as full-sib
#' res$pairs[res$pairs$id1 == "CB_AB_01" & res$pairs$id2 == "CB_AB_02", ]
#' # Recorded parent-offspring links whose genomic class disagrees; a few
#' # near the class boundary are expected at this number of loci
#' head(res$conflicts)
#' # Tag P/A data (no PO/FS split possible)
#' res.gs <- gl.report.kin.classes(testset.gs, kin = gl.kin(testset.gs))

#' @seealso \code{\link{gl.kin}}, \code{\link{gl.report.kin.confidence}},
#' \code{\link{gl.report.parent.offspring}}

#' @export
#' @return Invisibly, a list with two components: pairs, a data frame with
#' columns id1, id2, kinship, kinship.adj and class (restricted to pairs with
#' class != "unrelated" unless all.pairs = TRUE); and conflicts, a data frame
#' of recorded parent-offspring links whose genomic class disagrees (columns
#' offspring, parent, kinship, kinship.adj, class), or NULL if there are no
#' conflicts or no sire/dam columns in ind.metrics.
#'
# ----------------------
# Function
gl.report.kin.classes <- function(x,
                                  kin = NULL,
                                  oh.thresh = 0.005,
                                  all.pairs = FALSE,
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
    if (!is.numeric(oh.thresh) || length(oh.thresh) != 1 ||
        oh.thresh < 0 || oh.thresh >= 1) {
        stop(error("Fatal Error: oh.thresh must be a single value in [0, 1)\n"))
    }
    if (nInd(x) < 2) {
        stop(error("Fatal Error: At least two individuals are required\n"))
    }

    ids <- indNames(x)

# DO THE JOB ----------------------
    # Baseline adjustment: frequency-based estimators are relative to the
    # sample reference, so anchor "unrelated" at the median off-diagonal value
    baseline <- stats::median(kin[row(kin) != col(kin)])
    if (verbose >= 2) {
        cat(report("  Baseline (median off-diagonal kinship):",
                   round(baseline, 4), "\n"))
    }

    ut <- which(upper.tri(kin), arr.ind = TRUE)
    pairs.all <- data.frame(
        id1 = ids[ut[, 1]],
        id2 = ids[ut[, 2]],
        kinship = kin[ut],
        kinship.adj = kin[ut] - baseline,
        stringsAsFactors = FALSE
    )

    # Classify to nearest of {0.25, 0.125, 0.0625, 0} via midpoint thresholds
    pairs.all$class <- as.character(cut(
        pairs.all$kinship.adj,
        breaks = c(-Inf, 0.03125, 0.09375, 0.1875, Inf),
        labels = c("unrelated", "third-degree", "second-degree", "first-degree"),
        right = FALSE
    ))

    # DO THE JOB -- SNP data: split first-degree into PO vs FS
    if (datatype == "SNP") {
        fd <- which(pairs.all$class == "first-degree")
        if (length(fd) > 0) {
            if (verbose >= 2) {
                cat(report("  Splitting", length(fd),
                           "first-degree pairs into PO vs FS by opposite homozygotes\n"))
            }
            mat <- as.matrix(x)
            for (k in fd) {
                gi <- mat[pairs.all$id1[k], ]
                gj <- mat[pairs.all$id2[k], ]
                ok <- !is.na(gi) & !is.na(gj)
                n.ok <- sum(ok)
                if (n.ok == 0) {
                    if (verbose >= 1) {
                        cat(warn("  Warning: No shared called loci for pair",
                                 pairs.all$id1[k], "x", pairs.all$id2[k],
                                 "; class left as first-degree\n"))
                    }
                    next
                }
                oh <- sum((gi[ok] == 0 & gj[ok] == 2) |
                          (gi[ok] == 2 & gj[ok] == 0))
                pairs.all$class[k] <-
                    if (oh / n.ok < oh.thresh) "parent-offspring" else "full-sib"
            }
        }
    } else {
        # DO THE JOB -- Tag P/A data: PO/FS split not possible
        if (verbose >= 2) {
            cat(report(
                "  Tag P/A data (SilicoDArT): PO vs FS split not possible, first-degree retained\n"
            ))
        }
    }

    # Cross-check against recorded sire/dam links -----------
    conflicts <- NULL
    im <- x@other$ind.metrics
    if (!is.null(im) && all(c("sire", "dam") %in% names(im))) {
        recorded <- data.frame(
            offspring = rep(ids, 2),
            parent = c(as.character(im$sire), as.character(im$dam)),
            stringsAsFactors = FALSE
        )
        recorded <- recorded[!is.na(recorded$parent) &
                             recorded$parent != "" &
                             recorded$parent %in% ids, ]
        if (nrow(recorded) > 0) {
            if (verbose >= 2) {
                cat(report("  Cross-checking", nrow(recorded),
                           "recorded parent-offspring links\n"))
            }
            key.all <- paste(pmin(pairs.all$id1, pairs.all$id2),
                             pmax(pairs.all$id1, pairs.all$id2))
            key.rec <- paste(pmin(recorded$offspring, recorded$parent),
                             pmax(recorded$offspring, recorded$parent))
            hit <- match(key.rec, key.all)
            recorded$kinship <- pairs.all$kinship[hit]
            recorded$kinship.adj <- pairs.all$kinship.adj[hit]
            recorded$class <- pairs.all$class[hit]
            expected <- if (datatype == "SNP") "parent-offspring" else "first-degree"
            bad <- which(recorded$class != expected)
            if (length(bad) > 0) {
                conflicts <- recorded[bad, , drop = FALSE]
                rownames(conflicts) <- NULL
            }
        }
    } else {
        if (verbose >= 2) {
            cat(report("  No sire/dam columns in ind.metrics; pedigree cross-check skipped\n"))
        }
    }

    # Printing outputs -----------
    if (verbose >= 3) {
        cat("  Relationship class counts (all pairs)\n")
        tbl <- table(pairs.all$class)
        for (cl in names(tbl)) {
            cat(paste0("    ", format(cl, width = 16), ": ", tbl[[cl]]), "\n")
        }
        if (is.null(conflicts)) {
            cat("  Pedigree cross-check: no conflicts detected\n")
        } else {
            cat(paste("  Pedigree cross-check:", nrow(conflicts),
                      "conflicting recorded parent-offspring link(s)\n"))
            print(conflicts)
        }
    }

    pairs.out <- if (all.pairs) {
        pairs.all
    } else {
        pairs.all[pairs.all$class != "unrelated", , drop = FALSE]
    }
    rownames(pairs.out) <- NULL

# FLAG SCRIPT END ---------------
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
# ----------------------

    # RETURN
    invisible(list(pairs = pairs.out, conflicts = conflicts))
}
