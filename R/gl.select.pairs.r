#' @name gl.select.pairs
#' @title Selects a set of breeding pairs to maximise retained gene diversity
#' @family captive management
#' @description
#' Chooses a set of male x female breeding pairs from a genlight object using a
#' pairwise kinship matrix, under constraints on the number of pairings per
#' sire and per dam, a ceiling on prospective offspring inbreeding, and a
#' target number of pairs. Three selection schemes are offered, including a
#' dynamic scheme that re-evaluates projected gene diversity as each virtual
#' offspring accrues. This is a SNP-based analogue of the Auto Pair capability
#' of the PMx studbook-management software.
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param kin Pairwise kinship matrix with row and column names identical to
#' indNames(x), as produced by gl.kin(); if NULL, computed internally with
#' gl.kin() [default NULL].
#' @param n.pairs Number of pairs to select [default NULL, resolving to
#' min(#males x max.per.sire, #females x max.per.dam), capped at the number of
#' distinct male x female pairings].
#' @param scheme Pair selection scheme, one of 'dynamic', 'static' or 'ranked';
#' see Details [default 'dynamic'].
#' @param max.per.sire Maximum number of pairings per male; may be Inf
#' [default 1].
#' @param max.per.dam Maximum number of pairings per female; may be Inf
#' [default 1].
#' @param f.max Ceiling on the prospective offspring inbreeding coefficient
#' (the pair's kinship); pairs at or above this value are never selected
#' [default 0.125].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#' @details
#' Kinship should be estimated on the most inclusive dataset available and the
#' matrix subset to the managed group; computing kinship on a small family
#' group alone inflates estimates (the allele-frequency reference collapses
#' onto the family itself).
#' This function mirrors the Auto Pair capability of PMx (Lacy, Ballou & Pollak
#' 2012, Methods in Ecology and Evolution 3:433-437; PMx Users Manual v1.0
#' pp. 87-92), including its static, ranked and dynamic pairing schemes and its
#' dynamic recalculation of population genetic values as virtual offspring
#' accrue -- here made analytic through the kinship matrix (an offspring's
#' kinships are the means of its parents' kinships; see utils.kin.dgd), so no
#' simulation is needed.
#'
#' A pair (m, f) is feasible when both parents have remaining capacity
#' (max.per.sire / max.per.dam) and kin[m, f] < f.max. Schemes:
#' \itemize{
#' \item 'static' -- all male x female pairings are ranked once, ascending, by
#' mean parental mean kinship (MK_m + MK_f)/2, and accepted greedily in that
#' order while feasible, until n.pairs is reached.
#' \item 'ranked' -- repeatedly take the female with remaining capacity and
#' lowest mean kinship who has at least one feasible mate, and pair her with
#' the feasible male of lowest mean kinship.
#' \item 'dynamic' -- at each step, every feasible remaining pair's projected
#' gene diversity is evaluated via utils.kin.dgd() given the pairs already
#' selected (each selected pair contributing one virtual offspring), and the
#' pair maximising projected gene diversity is accepted; the evaluation is
#' repeated from scratch each round.
#' }
#' Mean kinships are row means of the kinship matrix (self included). Ties are
#' broken by sire then dam name, so all schemes are deterministic. If no
#' feasible pair remains before n.pairs is reached, the function warns and
#' stops early with the pairs selected so far.
#'
#' Sexes are read from x@@other$ind.metrics$sex (values 'Male', 'Female',
#' 'Unknown'); individuals of unknown sex are excluded from the pools with a
#' warning.
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # Estimate kinship on the FULL dataset, then subset to the captive colony
#' kin <- gl.kin(testset.gl)
#' cb <- gl.keep.pop(testset.gl, pop.list = "EmmacCaptBred", verbose = 0)
#' kin.cb <- kin[indNames(cb), indNames(cb)]
#' pairs <- gl.select.pairs(cb, kin = kin.cb, n.pairs = 5)
#' pairs  # avoids sib x sib pairings (f.off < 0.125)
#' @seealso \code{\link{gl.kin}}, \code{\link{gl.report.mate.suitability}},
#' \code{\link{gl.report.repro.targets}}
#' @export
#' @return Invisibly, a data frame with one row per selected pair, in order of
#' selection: $sire, $dam, $f.off (prospective offspring inbreeding,
#' kin[sire, dam]), $mk.sire, $mk.dam, $dgd.cum (projected gene diversity of
#' the population after this and all previously selected pairs each contribute
#' one offspring). Attributes 'gd.start' (gene diversity before pairing) and
#' 'gd.projected' (after all selected pairs) are attached.
#'
# ----------------------
# Function
gl.select.pairs <- function(x,
                            kin = NULL,
                            n.pairs = NULL,
                            scheme = "dynamic",
                            max.per.sire = 1,
                            max.per.dam = 1,
                            f.max = 0.125,
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

# FUNCTION SPECIFIC ERROR CHECKING ----------------
    # Kinship matrix
    kin <- utils.kin.check(x, kin, verbose = verbose)

    # Sexes
    if (is.null(x@other$ind.metrics$sex)) {
        stop(error(
            "Fatal Error: pair selection requires a 'sex' column (values 'Male', 'Female', 'Unknown') in x@other$ind.metrics; none found\n"
        ))
    }
    ids <- indNames(x)
    sex <- as.character(x@other$ind.metrics$sex)
    males <- ids[sex == "Male"]
    females <- ids[sex == "Female"]
    n.unknown <- sum(!(sex %in% c("Male", "Female")))
    if (n.unknown > 0 && verbose >= 1) {
        cat(warn("  Warning:", n.unknown,
                 "individuals of unknown sex excluded from the male/female pools\n"))
    }
    if (length(males) == 0 || length(females) == 0) {
        stop(error(
            "Fatal Error: at least one Male and one Female are required to select pairs\n"
        ))
    }
    n.m <- length(males)
    n.f <- length(females)

    # Scheme
    if (!is.character(scheme) || length(scheme) != 1 ||
        !(scheme %in% c("dynamic", "static", "ranked"))) {
        if (verbose >= 1) {
            cat(warn("  Warning: scheme must be one of 'dynamic', 'static' or 'ranked'; resetting to 'dynamic'\n"))
        }
        scheme <- "dynamic"
    }

    # Capacities
    if (!is.numeric(max.per.sire) || length(max.per.sire) != 1 || max.per.sire < 1) {
        if (verbose >= 1) {
            cat(warn("  Warning: max.per.sire must be a single value >= 1; resetting to 1\n"))
        }
        max.per.sire <- 1
    }
    if (!is.numeric(max.per.dam) || length(max.per.dam) != 1 || max.per.dam < 1) {
        if (verbose >= 1) {
            cat(warn("  Warning: max.per.dam must be a single value >= 1; resetting to 1\n"))
        }
        max.per.dam <- 1
    }

    # f.max
    if (!is.numeric(f.max) || length(f.max) != 1 || f.max <= 0 || f.max > 1) {
        if (verbose >= 1) {
            cat(warn("  Warning: f.max must be a single value in (0, 1]; resetting to 0.125\n"))
        }
        f.max <- 0.125
    }

    # n.pairs
    n.pairs.default <- min(n.m * max.per.sire, n.f * max.per.dam, n.m * n.f)
    if (is.null(n.pairs)) {
        n.pairs <- n.pairs.default
    } else if (!is.numeric(n.pairs) || length(n.pairs) != 1 || n.pairs < 1) {
        if (verbose >= 1) {
            cat(warn("  Warning: n.pairs must be a single value >= 1; resetting to the default",
                     n.pairs.default, "\n"))
        }
        n.pairs <- n.pairs.default
    }
    n.pairs <- floor(n.pairs)

# DO THE JOB ----------------------
    mk <- rowMeans(kin)
    gd.start <- utils.kin.dgd(kin)
    if (verbose >= 2) {
        cat(report("  Selecting up to", n.pairs, "pairs by the '", scheme,
                   "' scheme; starting gene diversity", round(gd.start, 4), "\n"))
    }

    cap.m <- stats::setNames(rep(max.per.sire, n.m), males)
    cap.f <- stats::setNames(rep(max.per.dam, n.f), females)
    sel <- matrix(character(0), ncol = 2)
    rows <- list()

    accept.pair <- function(m, f) {
        sel <<- rbind(sel, c(m, f))
        cap.m[m] <<- cap.m[m] - 1
        cap.f[f] <<- cap.f[f] - 1
        rows[[length(rows) + 1]] <<- data.frame(
            sire = m,
            dam = f,
            f.off = kin[m, f],
            mk.sire = mk[m],
            mk.dam = mk[f],
            dgd.cum = utils.kin.dgd(kin, add.pairs = sel),
            stringsAsFactors = FALSE
        )
    }
    warn.early <- function() {
        if (verbose >= 1) {
            cat(warn("  Warning: no feasible pair remains after",
                     nrow(sel), "of", n.pairs,
                     "requested pairs; stopping early\n"))
        }
    }

    if (scheme == "static") {
        # Rank all pairings once by mean parental mean kinship, ascending
        grid <- expand.grid(sire = males, dam = females,
                            stringsAsFactors = FALSE)
        grid$mk.mean <- (mk[grid$sire] + mk[grid$dam]) / 2
        grid <- grid[order(grid$mk.mean, grid$sire, grid$dam), ]
        for (i in seq_len(nrow(grid))) {
            if (nrow(sel) >= n.pairs) break
            m <- grid$sire[i]
            f <- grid$dam[i]
            if (cap.m[m] > 0 && cap.f[f] > 0 && kin[m, f] < f.max) {
                accept.pair(m, f)
            }
        }
        if (nrow(sel) < n.pairs) warn.early()

    } else if (scheme == "ranked") {
        # Alternate: lowest-MK available female, best feasible male
        repeat {
            if (nrow(sel) >= n.pairs) break
            fem.ord <- females[cap.f[females] > 0]
            fem.ord <- fem.ord[order(mk[fem.ord], fem.ord)]
            accepted <- FALSE
            for (f in fem.ord) {
                m.ok <- males[cap.m[males] > 0 & kin[males, f] < f.max]
                if (length(m.ok) > 0) {
                    m <- m.ok[order(mk[m.ok], m.ok)][1]
                    accept.pair(m, f)
                    accepted <- TRUE
                    break
                }
            }
            if (!accepted) {
                warn.early()
                break
            }
        }

    } else {
        # 'dynamic': each round, accept the feasible pair maximising
        # projected gene diversity given the pairs already selected
        repeat {
            if (nrow(sel) >= n.pairs) break
            cand <- expand.grid(sire = males[cap.m[males] > 0],
                                dam = females[cap.f[females] > 0],
                                stringsAsFactors = FALSE)
            if (nrow(cand) > 0) {
                cand <- cand[kin[cbind(cand$sire, cand$dam)] < f.max, ,
                             drop = FALSE]
            }
            if (nrow(cand) == 0) {
                warn.early()
                break
            }
            cand <- cand[order(cand$sire, cand$dam), ]
            proj <- vapply(seq_len(nrow(cand)), function(i) {
                utils.kin.dgd(kin,
                              add.pairs = rbind(sel,
                                                c(cand$sire[i], cand$dam[i])))
            }, numeric(1))
            best <- which.max(proj)
            accept.pair(cand$sire[best], cand$dam[best])
            if (verbose >= 2) {
                cat(report("    Pair", nrow(sel), ":", cand$sire[best], "x",
                           cand$dam[best], "; projected GD",
                           round(proj[best], 4), "\n"))
            }
        }
    }

    # Assemble the pair table
    if (length(rows) > 0) {
        out <- do.call(rbind, rows)
    } else {
        out <- data.frame(sire = character(0),
                          dam = character(0),
                          f.off = numeric(0),
                          mk.sire = numeric(0),
                          mk.dam = numeric(0),
                          dgd.cum = numeric(0),
                          stringsAsFactors = FALSE)
    }
    rownames(out) <- NULL
    attr(out, "gd.start") <- gd.start
    attr(out, "gd.projected") <- if (nrow(out) > 0) {
        out$dgd.cum[nrow(out)]
    } else {
        gd.start
    }

# Printing outputs -----------
    if (verbose >= 3) {
        cat(report("  Selected", nrow(out), "pairs; gene diversity",
                   round(gd.start, 4), "->",
                   round(attr(out, "gd.projected"), 4), "\n"))
        print(out, row.names = FALSE)
    }

# FLAG SCRIPT END ---------------
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
# ----------------------

    # RETURN
    invisible(out)
}
