#' @name gl.report.kinship
#' @title Reports a PMx-style kinship overview for a genlight object
#' @family captive management
#' @description
#' Produces the standard genetics overview used in captive management, in the
#' style of the PMx software, from a molecular kinship matrix: a
#' per-individual table of mean kinship (MK) within the individual's
#' population, MK rank within population and sex, and inbreeding coefficient
#' F, and a population summary (overall and per population) of sample size,
#' mean kinship, gene diversity (GD), founder genome equivalents (FGE) and
#' mean F. If no kinship matrix is supplied, one
#' is computed with gl.kin using the default method for the datatype (SNP:
#' 'grm'; SilicoDArT: 'dominant').
#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data [required].
#' @param kin A kinship matrix following the series contract (dimnames =
#' indNames(x), diagonal = 0.5*(1+F), off-diagonal = pairwise kinship), as
#' produced by gl.kin [default NULL, computed internally with gl.kin].
#' @param plot.display If TRUE, resultant plots are displayed in the plot
#' window [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors List of two color names for the borders and fill of the
#' plots [default c("#2171B5","#6BAED6")].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].
#' @details
#' Per-individual statistics follow the PMx conventions, where each
#' population is treated as a managed population: MK_i is the mean kinship
#' of individual i with the members of its own population, including itself
#' (the row mean of that population's block of the kinship matrix); MKrank
#' ranks individuals within their population and sex (ascending, so rank 1
#' is the lowest MK and hence the most genetically valuable individual of
#' that sex in that population; ties broken by order of appearance;
#' individuals of unknown sex -- NA, empty or 'unknown' in any case -- are
#' ranked in one group); and F = 2*diag(kin) - 1 recovers the inbreeding
#' coefficient from the self-kinship. The kinship matrix itself is best
#' estimated on the full dataset (see gl.kin), which is why MK is taken from
#' blocks of it rather than re-estimated per population. A mean kinship over
#' the whole dataset is not reported: genomic kinship is centred on the
#' supplied individuals, so each row of the matrix averages about 0
#' (exactly 0 for method 'grm'), whatever the individual. Sex is taken
#' from the ind.metrics slot (column 'sex', matched case-insensitively) and
#' is NA if absent. Note that under the
#' dominant (SilicoDArT) kinship estimator the diagonal is fixed at 0.5, so F
#' is 0 by construction.
#'
#' The population summary reports, for the overall dataset and for each
#' population: n, meanMK = mean of that population's kinship block (full
#' block including the diagonal, equal to the mean of the MK column for that
#' population), GD = 1 - meanMK, FGE = 1/(2*meanMK), and meanF. GD and FGE
#' are the PMx identities. In the 'overall' row GD and FGE are NA: genomic
#' kinship is centred on the supplied individuals, so the mean of the full
#' matrix is about 0 by construction (exactly 0 for method 'grm'), and GD =
#' 1 or an unbounded FGE would carry no information. Where a population
#' block's mean kinship is zero or negative (within 1e-8), FGE is undefined
#' and reported as NA.
#'
#' Missing kinship values (for example, pairs sharing no scored loci under
#' the dominant estimator) are ignored in the means, with a warning giving
#' the number of such pairs.
#'
#' The plot is a histogram of the off-diagonal pairwise kinship values with a
#' dashed line at their mean. A color vector can be obtained with
#' gl.select.colors() and passed with the plot.colors parameter. Themes can
#' be obtained from \url{https://ggplot2.tidyverse.org/reference/ggtheme.html}.
#' If a plot.file is given, the ggplot is saved as an RDS binary file using
#' saveRDS(); it can be reloaded with readRDS(). If a plot directory
#' (plot.dir) is specified, the ggplot binary is saved to that directory;
#' otherwise to the tempdir().
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' # SNP data
#' kin <- gl.kin(testset2.gl)
#' res <- gl.report.kinship(testset2.gl, kin = kin)
#' head(res$ind)
#' res$pop
#' # Tag P/A data -- dominant kinship computed internally
#' res.gs <- gl.report.kinship(testset2.gs)
#' # Captive-bred individuals, most valuable first within each sex
#' cb <- res$ind[res$ind$pop == "EmmacCaptBred", ]
#' cb[order(cb$sex, cb$MKrank), ]
#' @seealso \code{\link{gl.kin}}, \code{\link{gl.grm}},
#' \code{\link{gl.run.EMIBD9}}
#' @export
#' @return A list with two data frames, returned invisibly: $ind, the
#' per-individual table (id, pop, sex, MK, MKrank, F); and $pop, the
#' population summary (pop, n, meanMK, GD, FGE, meanF), with 'overall' as the
#' first row (GD and FGE NA, see Details).
#'
# ----------------------
# Function
gl.report.kinship <- function(x,
                              kin = NULL,
                              plot.display = TRUE,
                              plot.theme = theme_dartR(),
                              plot.colors = NULL,
                              plot.file = NULL,
                              plot.dir = NULL,
                              verbose = NULL) {
    # PRELIMINARIES -- checking ----------------
    # SET VERBOSITY
    verbose <- gl.check.verbosity(verbose)
    if (verbose == 0) {
        plot.display <- FALSE
    }

    # SET WORKING DIRECTORY
    plot.dir <- gl.check.wd(plot.dir, verbose = 0)

    # SET COLOURS
    if (is.null(plot.colors)) {
        plot.colors <- c("#2171B5", "#6BAED6")
    } else {
        if (length(plot.colors) > 2) {
            if (verbose >= 2) {
                cat(warn("  More than 2 colors specified, only the first 2 are used\n"))
            }
            plot.colors <- plot.colors[1:2]
        }
    }

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
    kin <- utils.kin.check(x, kin, verbose = verbose)

    # Assign a population if none is specified
    if (is.null(pop(x)) | is.na(length(pop(x))) | length(pop(x)) <= 0) {
        if (verbose >= 2) {
            cat(important("  Population assignments not detected, individuals",
                          "assigned to a single population labelled 'pop1'\n"))
        }
        pop(x) <- factor(rep("pop1", nInd(x)))
    }
    # With one population, MK is taken over the set kinship was estimated on
    # and is ~0 for everyone unless kin comes from a wider reference
    if (nPop(x) == 1 && verbose >= 1) {
        ref.ids <- attr(kin, "ref.ids")
        self <- (!is.null(ref.ids) && setequal(ref.ids, indNames(x))) ||
            all(abs(rowMeans(kin, na.rm = TRUE)) < 1e-10)
        if (self) {
            cat(warn("  Warning: one population and kinship estimated on these",
                     "individuals only, so MK is ~0 for everyone and MKrank is",
                     "not informative; supply kin estimated on a wider dataset",
                     "(e.g. kin = gl.kin(full.dataset))\n"))
        }
    }

    # DO THE JOB ----------------------

    # Per-individual table
    if (verbose >= 2) {
        cat(report("  Computing per-individual mean kinship, rank within sex and F\n"))
    }
    n.na <- sum(is.na(kin[upper.tri(kin)]))
    if (n.na > 0 && verbose >= 1) {
        cat(warn("  Warning:", n.na, "pair(s) of individuals have missing",
                 "kinship; they are ignored in the means\n"))
    }
    pops <- as.character(pop(x))
    # PMx MK is relative to the managed population: use each population's
    # block of the (full-dataset) kinship matrix
    MK <- rep(NA_real_, nInd(x))
    for (p in unique(pops)) {
        idx <- which(pops == p)
        MK[idx] <- rowMeans(kin[idx, idx, drop = FALSE], na.rm = TRUE)
    }
    Fi <- 2 * diag(kin) - 1

    sex <- rep(NA_character_, nInd(x))
    if (!is.null(x@other$ind.metrics)) {
        sex.col <- which(tolower(names(x@other$ind.metrics)) == "sex")
        if (length(sex.col) >= 1) {
            sex <- as.character(x@other$ind.metrics[[sex.col[1]]])
        }
    }
    sex.group <- ifelse(is.na(sex) | sex == "" | tolower(sex) == "unknown",
                        "unknown", sex)
    MKrank <- stats::ave(MK, pops, sex.group,
                         FUN = function(v) rank(v, ties.method = "first"))

    ind.df <- data.frame(
        id = indNames(x),
        pop = pops,
        sex = sex,
        MK = unname(MK),
        MKrank = as.integer(MKrank),
        F = unname(Fi),
        stringsAsFactors = FALSE,
        row.names = NULL
    )

    # Population summary (overall + per population)
    if (verbose >= 2) {
        cat(report("  Summarising kinship by population (n, meanMK, GD, FGE, meanF)\n"))
    }
    pop.block <- function(idx, label, overall = FALSE) {
        sub <- kin[idx, idx, drop = FALSE]
        mk <- mean(sub, na.rm = TRUE)
        # the full matrix averages ~0 by construction (centred kinship), so
        # overall GD and FGE carry no information; a mean within rounding
        # of 0 would give an astronomically large FGE
        ok <- !overall && is.finite(mk) && mk > 1e-8
        data.frame(
            pop = label,
            n = length(idx),
            meanMK = mk,
            GD = if (overall) NA_real_ else 1 - mk,
            FGE = if (ok) 1 / (2 * mk) else NA_real_,
            meanF = mean(2 * diag(sub) - 1, na.rm = TRUE),
            stringsAsFactors = FALSE
        )
    }
    blocks <- lapply(levels(pop(x)), function(p) {
        pop.block(which(pops == p), p)
    })
    pop.df <- rbind(pop.block(seq_len(nInd(x)), "overall", overall = TRUE),
                    do.call(rbind, blocks))
    row.names(pop.df) <- NULL

    if (any(is.na(pop.df$FGE[-1])) && verbose >= 1) {
        cat(warn(paste0("  Warning: mean kinship <= 0 for one or more ",
                        "populations (possible with centred GRM estimates); ",
                        "FGE reported as NA for those rows\n")))
    }

    # Printing outputs -----------
    if (verbose >= 3) {
        cat(report("  Per-individual summary (first 10 rows):\n"))
        print(head(ind.df, 10))
        cat(report("  Population summary:\n"))
        print(pop.df)
    }

    # Plot the results -----------------
    vals <- kin[upper.tri(kin)]
    df.plot <- data.frame(kinship = vals)
    df.plot <- df.plot[!is.na(df.plot$kinship), , drop = FALSE]
    kinship <- NULL  # silence R CMD check note for ggplot NSE

    method.lab <- attr(kin, "method")
    if (is.null(method.lab)) {
        method.lab <- "supplied"
    }
    p1 <- ggplot(data = df.plot, aes(x = kinship)) +
        geom_histogram(bins = 50,
                       color = plot.colors[1],
                       fill = plot.colors[2]) +
        geom_vline(xintercept = mean(df.plot$kinship),
                   color = plot.colors[1],
                   linetype = "dashed") +
        xlab("Pairwise kinship (off-diagonal)") +
        ylab("Count") +
        ggtitle(paste0("Pairwise Kinship [method: ", method.lab, "]")) +
        plot.theme
    if (plot.display) {
        print(p1)
    }

    # Optionally save the plot ---------------------
    if (!is.null(plot.file) && exists("p1")) {
        tmp <- utils.plot.save(p1,
                               dir = plot.dir,
                               file = plot.file,
                               verbose = verbose)
    }

    # FLAG SCRIPT END ---------------
    if (verbose >= 1) {
        cat(report("Completed:", funname, "\n"))
    }
    # ----------------------

    # RETURN
    invisible(list(ind = ind.df, pop = pop.df))
}
