#' @name gl.kin
#' @title Calculates a standardized kinship matrix from SNP or presence/absence data
#' @family captive management
#' @description
#' Produces the single standardized kinship matrix that underpins the captive
#' management (PMx-style) series of functions in dartR.captive. Whatever the
#' estimation method, the returned matrix has the same currency: a base numeric
#' matrix with row and column names set to the individual labels, diagonal
#' elements equal to the self-kinship 0.5*(1+F) where F is the individual
#' inbreeding coefficient, and off-diagonal elements equal to the pairwise
#' kinship coefficient between individuals.
#'
#' The method is selected automatically from the datatype if not specified:
#' SNP data default to method 'grm' (genomic relationship matrix via rrBLUP);
#' presence/absence (SilicoDArT) data are forced to method 'dominant'
#' (band-sharing kinship via utils.kin.dominant). Method 'emibd9' (SNP only)
#' uses the external program EMIBD9 via gl.run.EMIBD9.
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param method Method for estimating kinship, one of 'grm' (SNP only),
#' 'emibd9' (SNP only, requires the external EMIBD9 binary) or 'dominant'
#' (SilicoDArT only) [default NULL, resolved to 'grm' for SNP data and
#' 'dominant' for SilicoDArT data].
#' @param emibd9.path Path to the folder containing the EMIBD9 executable,
#' passed to gl.run.EMIBD9; used only with method 'emibd9'
#' [default getwd()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#' @details
#' The kinship contract: the returned matrix kin has dimnames equal to
#' indNames(x); kin[i,i] = 0.5*(1+F_i), the self-kinship of individual i (0.5
#' for a non-inbred individual); kin[i,j] = f_ij, the kinship between
#' individuals i and j. Gene diversity for any set of individuals is
#' GD = 1 - mean(kin) taken over the full matrix including the diagonal, and
#' mean kinship for individual i is MK_i = rowMeans(kin)[i], including self.
#' These are the conventions PMx uses for pedigree-based management, with one
#' difference: PMx measures kinship from the founders, whereas here the base
#' is the supplied dataset. Because genomic kinship is centred on the sample,
#' GD of the full dataset is about 1 by construction (exactly 1 for method
#' 'grm'); GD values are meaningful as comparisons between subsets of
#' individuals, or before and after removals and additions, not as an
#' absolute measure of diversity.
#'
#' Method 'grm' computes the genomic relationship matrix G with
#' gl.grm (rrBLUP::A.mat, diagonal ~ 1+F, off-diagonal ~ twice the kinship)
#' and returns kin = G/2, the same scale gl.grm.network uses by default: the
#' diagonal is 0.5*(1+F) and the off-diagonal elements are pairwise kinship,
#' both relative to the allele frequencies of the supplied dataset.
#' gl.grm imputes missing genotypes with the locus mean, which pulls kinship
#' toward 0: on testset2.gl (15% missing) parent-offspring pairs average
#' 0.19, and 0.25 after gl.filter.callrate(threshold = 0.95). Filter on call
#' rate before gl.kin when kinship values are compared against thresholds or
#' relationship classes.
#'
#' Method 'emibd9' runs the external EMIBD9 program (Wang 2022, Methods in
#' Ecology and Evolution 13:2443-2462) via gl.run.EMIBD9 and symmetrises the
#' populated triangle of the returned $rel matrix. $rel is EMIBD9's r(1,2),
#' the kinship coefficient, so it is used without rescaling; where EMIBD9
#' leaves the diagonal unpopulated it defaults to the non-inbred self-kinship
#' 0.5. gl.run.EMIBD9 runs with Inbreed = FALSE, so the diagonal is 0.5.
#' EMIBD9 must be installed separately; give its folder with emibd9.path (see
#' gl.run.EMIBD9 for download details).
#'
#' Method 'dominant' calls utils.kin.dominant, a standardized band-sharing
#' covariance for presence/absence data. Dominant-marker kinship is a
#' rank-valid approximation: rankings and relative comparisons are sound but
#' absolute values are approximate (see utils.kin.dominant for details), and
#' the diagonal is fixed at 0.5 because individual inbreeding is not estimable
#' from single dominant bands.
#'
#' Kinship between two different individuals cannot exceed 0.5 except for
#' duplicate samples or clones; if any pair exceeds 0.5 a warning reports how
#' many (see gl.report.replicates in dartR.base to identify duplicates).
#'
#' Important caveat: with all methods, kinship values are relative to the
#' allele (or band) frequencies of the full dataset supplied, which acts as
#' the reference population. Estimates from a dataset containing only close
#' relatives (for example, a single captive family) will be inflated. Where
#' possible, estimate kinship on a dataset that includes the wider source
#' population, then subset the resulting matrix.
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # SNP data
#' kin <- gl.kin(testset2.gl)
#' kin[1:5, 1:5]
#' # Tag P/A data
#' kins <- gl.kin(testset2.gs)
#' \dontrun{
#' # Requires the external EMIBD9 binary (see gl.run.EMIBD9)
#' kin.ibd <- gl.kin(testset2.gl, method = "emibd9",
#'                   emibd9.path = "path/to/EMIBD9")
#' }
#' @seealso \code{\link{gl.report.kinship}}, \code{\link{gl.grm}},
#' \code{\link{gl.run.EMIBD9}}
#' @export
#' @return A square numeric kinship matrix with dimnames = indNames(x),
#' diagonal 0.5*(1+F), off-diagonal pairwise kinship, and attributes 'method',
#' 'datatype', 'nLoc' and 'scale' (= "kinship", which gl.grm.network reads to
#' plot the matrix without rescaling), returned invisibly.
#'
# ----------------------
# Function
gl.kin <- function(x,
                   method = NULL,
                   emibd9.path = getwd(),
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
                                     accept = c("genlight", "SNP", "SilicoDArT"),
                                     verbose = verbose)

    # FUNCTION SPECIFIC ERROR CHECKING
    if (is.null(method)) {
        method <- ifelse(datatype == "SilicoDArT", "dominant", "grm")
        if (verbose >= 2) {
            cat(report("  No method specified, defaulting to method '",
                       method, "' for ", datatype, " data\n", sep = ""))
        }
    }
    method <- tolower(as.character(method[1]))
    if (!method %in% c("grm", "emibd9", "dominant")) {
        stop(error(paste0("Fatal Error: method must be one of 'grm', 'emibd9' ",
                          "or 'dominant', not '", method, "'\n")))
    }
    if (method %in% c("grm", "emibd9") && datatype == "SilicoDArT") {
        stop(error(paste0("Fatal Error: method '", method, "' requires SNP ",
                          "genotypes, but the supplied object holds ",
                          "presence/absence (SilicoDArT) data. Use ",
                          "method='dominant' (the SilicoDArT default)\n")))
    }
    if (method == "dominant" && datatype == "SNP") {
        stop(error(paste0("Fatal Error: method 'dominant' is for ",
                          "presence/absence (SilicoDArT) data, but the ",
                          "supplied object holds SNP genotypes. Use ",
                          "method='grm' (the SNP default) or method='emibd9'\n")))
    }

    # DO THE JOB ----------------------

    if (method == "grm") {
        if (verbose >= 2) {
            cat(report("  Computing kinship from the genomic relationship matrix (rrBLUP)\n"))
        }
        G <- gl.grm(x, plotheatmap = FALSE, verbose = 0)
        if (!is.matrix(G)) {
            stop(error(paste0("Fatal Error: gl.grm did not return a matrix; ",
                              "check that package rrBLUP is installed\n")))
        }
        # G is on the relatedness scale (diagonal 1 + F, off-diagonal ~2x
        # kinship). Subtracting the mean inbreeding from pairs would shift
        # every kinship by the Wahlund and missing-data effects in diag(G)
        kin <- G / 2
    }

    if (method == "emibd9") {
        if (verbose >= 2) {
            cat(report("  Computing kinship with the external program EMIBD9 (gl.run.EMIBD9)\n"))
        }
        res <- gl.run.EMIBD9(x, emibd9.path = emibd9.path,
                             plot.out = FALSE, verbose = 0)
        m <- res$rel
        if (is.null(m) || !is.matrix(m)) {
            stop(error("Fatal Error: gl.run.EMIBD9 did not return a $rel matrix\n"))
        }
        # Merge the populated triangle(s) into a symmetric matrix
        kin <- pmax(m, t(m), na.rm = TRUE)
        # $rel is EMIBD9's r(1,2), already the kinship coefficient (see
        # gl.run.EMIBD9), so no rescaling. Diagonal 0.5*(1+F) where EMIBD9
        # populated it, else the non-inbred self-kinship
        d <- diag(kin)
        d[is.na(d)] <- 0.5
        diag(kin) <- d
    }

    if (method == "dominant") {
        if (verbose >= 2) {
            cat(report("  Computing dominant-marker (band-sharing) kinship\n"))
        }
        kin <- utils.kin.dominant(x, verbose = verbose)
    }

    # Standardize the currency of the matrix
    dimnames(kin) <- list(indNames(x), indNames(x))
    attr(kin, "method") <- method
    attr(kin, "datatype") <- datatype
    attr(kin, "nLoc") <- nLoc(x)
    attr(kin, "scale") <- "kinship"

    # Pairs above 0.5 are beyond any relationship other than identity
    n.high <- sum(kin[upper.tri(kin)] > 0.5, na.rm = TRUE)
    if (n.high > 0 && verbose >= 1) {
        cat(warn("  Warning:", n.high, "pair(s) of individuals have kinship",
                 "above 0.5; they may be duplicate samples or clones",
                 "(see gl.report.replicates)\n"))
    }

    # Results summary -----------
    if (verbose >= 3) {
        cat(report("  Kinship matrix summary\n"))
        cat(report("    Method:", method, "| Datatype:", datatype,
                   "| Loci:", nLoc(x), "| Individuals:", nInd(x), "\n"))
        # Not GD = 1 - mean(kin): kinship is centred on this dataset, so
        # that is about 1 whatever the diversity (see Details)
        cat(report("    Mean pairwise kinship (off-diagonal):",
                   round(mean(kin[upper.tri(kin)], na.rm = TRUE), 4), "\n"))
        cat(report("    Self-kinship (diagonal) range:",
                   round(min(diag(kin), na.rm = TRUE), 4), "to",
                   round(max(diag(kin), na.rm = TRUE), 4), "\n"))
    }

    # FLAG SCRIPT END ---------------
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    # ----------------------

    # RETURN
    invisible(kin)
}
