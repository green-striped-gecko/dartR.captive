#' @name utils.kin.dgd
#' @title Gene diversity of a kinship matrix after removals and virtual offspring
#' @family captive management
#' @description
#' WARNING: UTILITY SCRIPTS ARE FOR INTERNAL USE ONLY AND SHOULD NOT BE USED
#' BY END USERS AS THEIR USE OUT OF CONTEXT COULD LEAD TO UNPREDICTABLE
#' OUTCOMES.
#'
#' The shared delta-GD engine of the captive management series. Given a
#' kinship matrix following the series contract (square, both dimnames set to
#' individual ids, diagonal = self-kinship 0.5*(1+F), off-diagonal = pairwise
#' kinship), it returns the gene diversity GD = 1 - mean(kin) over the full
#' modified matrix (diagonal included) after optionally (a) removing named
#' individuals and (b) appending virtual offspring of nominated parent pairs.
#' Management functions call this repeatedly to score candidate pairings,
#' removals or releases by their effect on GD.
#' @param kin A square numeric kinship matrix with identical row and column
#' names (individual ids), as produced by gl.kin [required].
#' @param drop Character vector of individual ids to remove from the matrix
#' before computing GD; unknown ids are a fatal error [default NULL].
#' @param add.pairs A 2-column character matrix of parent ids, one virtual
#' offspring appended per row, processed in row order; parents may be
#' previously added virtual offspring from earlier rows. Row names, if
#' supplied, become the ids of the virtual offspring, otherwise ids
#' 'offspring_1', 'offspring_2', ... are generated [default NULL].
#' @details
#' Removals are applied first, then virtual offspring are appended one row of
#' add.pairs at a time. For a virtual offspring v of parents a and b, kinship
#' with every individual j already in the matrix (including previously added
#' virtual offspring) is
#'
#' kin[v, j] = (kin[a, j] + kin[b, j]) / 2
#'
#' and its self-kinship is
#'
#' kin[v, v] = 0.5 * (1 + kin[a, b])
#'
#' so kinship between two virtual offspring arises from the same recursion as
#' rows accumulate. Because rows are processed in order, a parent id must be
#' either a retained individual or a virtual offspring defined in an earlier
#' row; anything else is a fatal error. The returned value is the single
#' numeric GD = 1 - mean(kin) over the full modified matrix including the
#' diagonal, per the series convention. NA entries in kin (for example,
#' individual pairs sharing no scored loci under the dominant estimator)
#' propagate to an NA gene diversity; callers should resolve missing kinship
#' values before optimisation.
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # Examples for testing
#' # kin <- gl.kin(testset.gl, verbose = 0)
#' # utils.kin.dgd(kin)                                # current gene diversity
#' # utils.kin.dgd(kin, drop = indNames(testset.gl)[1:5])
#' # off <- matrix(c("CB_AB_01", "CB_X_01"), ncol = 2) # one virtual offspring
#' # utils.kin.dgd(kin, add.pairs = off)
#' @seealso \code{\link{gl.kin}}, \code{\link{gl.report.kinship}}
#' @keywords internal
# @export
#' @return A single numeric value: the gene diversity GD = 1 - mean of the
#' modified kinship matrix (full matrix including the diagonal), returned
#' invisibly.
#'
# ----------------------
# Function
utils.kin.dgd <- function(kin,
                          drop = NULL,
                          add.pairs = NULL) {
    # ARGUMENT VALIDATION (pure engine -- no verbosity machinery) ----------
    if (!is.matrix(kin) || !is.numeric(kin)) {
        stop(error("Fatal Error: kin must be a numeric matrix\n"))
    }
    if (nrow(kin) != ncol(kin)) {
        stop(error("Fatal Error: kin must be a square matrix\n"))
    }
    if (is.null(rownames(kin)) || is.null(colnames(kin)) ||
        !identical(rownames(kin), colnames(kin))) {
        stop(error(paste0("Fatal Error: kin must have identical row and ",
                          "column names (individual ids)\n")))
    }

    # (a) Remove nominated individuals ----------
    if (!is.null(drop)) {
        drop <- as.character(drop)
        unknown <- setdiff(drop, rownames(kin))
        if (length(unknown) > 0) {
            stop(error(paste0("Fatal Error: id(s) in drop not present in the ",
                              "kinship matrix: ",
                              paste(unknown, collapse = ", "), "\n")))
        }
        keep <- setdiff(rownames(kin), drop)
        if (length(keep) == 0) {
            stop(error(paste0("Fatal Error: dropping all individuals leaves ",
                              "an empty kinship matrix\n")))
        }
        kin <- kin[keep, keep, drop = FALSE]
    }

    # (b) Append virtual offspring, one per row of add.pairs, in row order ----
    if (!is.null(add.pairs)) {
        if (is.data.frame(add.pairs)) {
            add.pairs <- as.matrix(add.pairs)
        }
        if (!is.matrix(add.pairs) || ncol(add.pairs) != 2) {
            stop(error(paste0("Fatal Error: add.pairs must be a 2-column ",
                              "character matrix of parent ids\n")))
        }
        mode(add.pairs) <- "character"
        n0 <- nrow(kin)
        n.add <- nrow(add.pairs)

        # Ids for the virtual offspring: rownames(add.pairs) if supplied,
        # else generated
        off.ids <- rownames(add.pairs)
        if (is.null(off.ids)) {
            off.ids <- paste0("offspring_", seq_len(n.add))
        }
        ids <- c(rownames(kin), off.ids)
        if (any(duplicated(ids))) {
            stop(error(paste0("Fatal Error: virtual offspring ids duplicate ",
                              "existing ids in the kinship matrix or each ",
                              "other\n")))
        }

        big <- matrix(NA_real_, nrow = n0 + n.add, ncol = n0 + n.add,
                      dimnames = list(ids, ids))
        big[seq_len(n0), seq_len(n0)] <- kin

        for (r in seq_len(n.add)) {
            k <- n0 + r
            # Parents must already exist: retained individuals or virtual
            # offspring from earlier rows only
            a <- match(add.pairs[r, 1], ids[seq_len(k - 1)])
            b <- match(add.pairs[r, 2], ids[seq_len(k - 1)])
            if (is.na(a) || is.na(b)) {
                bad <- setdiff(add.pairs[r, ], ids[seq_len(k - 1)])
                stop(error(paste0("Fatal Error: parent id(s) in add.pairs ",
                                  "row ", r, " not found among retained ",
                                  "individuals or previously added ",
                                  "offspring: ",
                                  paste(bad, collapse = ", "), "\n")))
            }
            v <- (big[a, seq_len(k - 1)] + big[b, seq_len(k - 1)]) / 2
            big[k, seq_len(k - 1)] <- v
            big[seq_len(k - 1), k] <- v
            big[k, k] <- 0.5 * (1 + big[a, b])
        }
        kin <- big
    }

    # GENE DIVERSITY over the full modified matrix, diagonal included ------
    gd <- 1 - mean(kin)

    # RETURN
    invisible(gd)
}
