#' @name utils.kin.as.kinship
#' @title Converts a relationship matrix to kinship using its scale tag
#' @family captive management
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' The single place where the captive management series reads the scale of a
#' relationship matrix. Matrices produced in dartR.captive carry
#' attr(m, "scale"): "relatedness" (gl.grm; diagonal 1 + F, off-diagonal
#' about twice the kinship) or "kinship" (gl.run.EMIBD9, gl.kin; diagonal
#' 0.5 * (1 + F)).
#' @param m A square numeric relationship matrix [required].
#' @param verbose Verbosity, already resolved by the caller; the only gated
#' output is the verbose >= 2 line announcing a conversion [default 0].
#' @details
#' A matrix tagged "relatedness" is halved; one tagged "kinship" is returned
#' unchanged. An untagged matrix is assumed to be kinship, the scale the
#' series works in, and returned unchanged. Any other tag is a fatal error.
#' The returned matrix is tagged "kinship"; its other attributes are kept.
#' @author Author(s): Luis Mijangos. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # Examples for testing
#' # G <- gl.grm(testset2.gl, plotheatmap = FALSE)  # tagged "relatedness"
#' # kin <- utils.kin.as.kinship(G)                  # G / 2, tagged "kinship"
#' @seealso \code{\link{gl.kin}}, \code{\link{utils.kin.check}}
#' @keywords internal
# @export
#' @return The matrix on the kinship scale, with attr(, "scale") = "kinship".
#'
# ----------------------
# Function
utils.kin.as.kinship <- function(m,
                                 verbose = 0) {
    scale <- attr(m, "scale")
    if (is.null(scale) || identical(scale, "kinship")) {
        attr(m, "scale") <- "kinship"
        return(m)
    }
    if (identical(scale, "relatedness")) {
        if (verbose >= 2) {
            cat(report("  Matrix tagged as relatedness; halved to kinship\n"))
        }
        m <- m / 2
        attr(m, "scale") <- "kinship"
        return(m)
    }
    stop(error(paste0("Fatal Error: unknown matrix scale '",
                      paste(scale, collapse = " "), "'; expected ",
                      "'kinship' or 'relatedness'\n")))
}
