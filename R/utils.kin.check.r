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
#' the default method for the datatype. The matrix (supplied or computed) is
#' then validated against the series contract: a base numeric matrix with row
#' and column names both identical to indNames(x). Consumer functions call
#' this once, immediately after their standard preamble, in place of a local
#' validation stanza.
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param kin A kinship matrix following the series contract (numeric matrix,
#' dimnames = indNames(x)), as produced by gl.kin [default NULL, computed
#' internally with gl.kin].
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
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # Examples for testing
#' # kin <- utils.kin.check(testset.gl)              # computes via gl.kin
#' # kin2 <- utils.kin.check(testset.gl, kin = kin)  # validates and passes through
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
                            verbose = 0) {
    # Auto-compute when no kinship matrix is supplied ----------
    if (is.null(kin)) {
        if (verbose >= 2) {
            cat(report("  No kinship matrix supplied; computing with gl.kin\n"))
        }
        kin <- gl.kin(x, verbose = 0)
    }

    # Validate against the series contract ----------
    if (!is.matrix(kin) || !is.numeric(kin) ||
        !identical(rownames(kin), indNames(x)) ||
        !identical(colnames(kin), indNames(x))) {
        stop(error("Fatal Error: kin must be a numeric matrix with row and column names identical to indNames(x)\n"))
    }

    # RETURN
    return(kin)
}
