#' @name utils.kin.check
#' @title Checks (and if necessary computes) a kinship matrix for a genlight object
#' @family captive management
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' The shared kin-validation gate of the captive management series. If no
#' kinship matrix is supplied (kin = NULL), one is computed with gl.kin using
#' the default method for the datatype, unless the caller needs a reference
#' (see need.reference). A supplied matrix tagged
#' attr(kin, "scale") = "relatedness" (such as gl.grm output) is halved to
#' kinship; an untagged matrix is assumed to be kinship. A matrix covering
#' more individuals than x (for example, kinship estimated on a dataset that
#' includes the source populations) is restricted to indNames(x). The result
#' is validated against the series contract: a base numeric matrix with row
#' and column names both identical to indNames(x). Consumer functions call
#' this once, immediately after their standard preamble, in place of a local
#' validation stanza.
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param kin A kinship matrix following the series contract (numeric matrix,
#' dimnames = indNames(x), or a superset of them), as produced by gl.kin
#' [default NULL, computed internally with gl.kin].
#' @param need.reference TRUE for callers whose results depend on mean
#' kinship or on the gene diversity of the whole of x; the matrix must then
#' have been estimated on a wider set of individuals than x [default FALSE].
#' @param verbose Verbosity, already resolved by the caller's
#' gl.check.verbosity(); the only gated output is the verbose >= 2 progress
#' line announcing the internal gl.kin call [default 0].
#' @details
#' This is a lightweight internal validator intended to be called from inside
#' functions that have already run the standard preamble
#' (gl.check.verbosity(), utils.flag.start(), utils.check.datatype()), so it
#' runs none of that machinery itself; the caller passes its resolved verbose
#' value through. Validation failure is a fatal error: kin must be a base
#' numeric matrix whose rownames and colnames are each identical() to
#' indNames(x) — the same order, not merely the same set — so that all
#' downstream row/column indexing by individual id is sound.
#'
#' Genomic kinship is centred on the individuals it was estimated on: each
#' row averages about 0 over that set (exactly 0 for gl.kin method 'grm'),
#' so mean kinship and whole-set gene diversity computed on the same set
#' carry no information. With need.reference = TRUE the matrix is rejected
#' as self-referenced when kin is NULL, when its 'ref.ids' attribute (set by
#' gl.kin) equals indNames(x), or when its row means over x are all within
#' 1e-10 of 0.
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # Examples for testing
#' # kin <- utils.kin.check(testset2.gl)              # computes via gl.kin
#' # kin2 <- utils.kin.check(testset2.gl, kin = kin)  # validates and passes through
#' @seealso \code{\link{gl.kin}}, \code{\link{gl.report.kinship}}
#' @keywords internal
# @export
#' @return The validated kinship matrix (the supplied kin, or the matrix
#' computed with gl.kin when kin is NULL), returned visibly for immediate
#' reassignment by the caller.
#'
# ----------------------
# Function
utils.kin.check <- function(x,
                            kin = NULL,
                            verbose = 0,
                            need.reference = FALSE) {
    ref.msg <- paste0(
        "Fatal Error: this analysis needs kinship estimated on a wider ",
        "reference than the individuals analysed. Kinship is centred on the ",
        "individuals it is estimated on, so their mean kinship is 0 and their ",
        "gene diversity 1 by construction. Estimate it on a dataset that ",
        "includes the source populations and pass it as kin, e.g. ",
        "kin = gl.kin(full.dataset); it is restricted to indNames(x)\n")

    # Auto-compute when no kinship matrix is supplied ----------
    if (is.null(kin)) {
        if (need.reference) {
            stop(error(ref.msg))
        }
        if (verbose >= 2) {
            cat(report("  No kinship matrix supplied; computing with gl.kin\n"))
        }
        kin <- gl.kin(x, verbose = 0)
    }

    # Convert by scale tag (relatedness -> kinship; unknown tag is fatal)
    kin <- utils.kin.as.kinship(kin, verbose = verbose)

    if (!is.matrix(kin) || !is.numeric(kin)) {
        stop(error("Fatal Error: kin must be a numeric matrix\n"))
    }
    ref.ids <- attr(kin, "ref.ids")

    # Restrict a larger (reference) matrix to the individuals of x ----------
    ids <- indNames(x)
    if (!identical(rownames(kin), ids) || !identical(colnames(kin), ids)) {
        if (!all(ids %in% rownames(kin)) || !all(ids %in% colnames(kin))) {
            stop(error("Fatal Error: kin must be a numeric matrix whose row and column names include every indNames(x)\n"))
        }
        if (verbose >= 2) {
            cat(report("  Kinship matrix of", nrow(kin), "individuals",
                       "restricted to the", length(ids), "in x\n"))
        }
        # an untagged larger matrix was estimated on at least its own rows
        if (is.null(ref.ids)) {
            ref.ids <- rownames(kin)
        }
        keep <- attributes(kin)[c("method", "datatype", "nLoc", "scale")]
        keep <- keep[!vapply(keep, is.null, logical(1))]
        kin <- kin[ids, ids, drop = FALSE]
        attributes(kin)[names(keep)] <- keep
    }
    attr(kin, "ref.ids") <- ref.ids

    # Reject self-referenced kinship when the caller needs a reference ------
    if (need.reference) {
        self <- (!is.null(ref.ids) && setequal(ref.ids, ids)) ||
            all(abs(rowMeans(kin, na.rm = TRUE)) < 1e-10)
        if (self) {
            stop(error(ref.msg))
        }
    }

    # RETURN
    return(kin)
}
