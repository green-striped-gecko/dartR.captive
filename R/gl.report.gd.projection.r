#' @name gl.report.gd.projection
#' @title Projects retention of gene diversity under drift, PMx Goals-screen style
#' @family captive management

#' @description
#' Projects the proportion of current gene diversity (expected heterozygosity)
#' retained over a management horizon given an effective population size, and
#' solves for the effective population size required to meet a retention target
#' and the years until the target retention is crossed.

#' @param x Name of the genlight object containing the SNP or presence/absence
#' (SilicoDArT) data; alternative to gd.now [default NULL].
#' @param gd.now Current gene diversity, supplied directly as a value in (0, 1];
#' alternative to x, and takes precedence if both are given [default NULL].
#' @param ne Effective population size assumed for the projection [required].
#' @param n Census population size, used only to report the Ne/N ratio
#' [default NULL].
#' @param gen.length Generation length in years [default 1].
#' @param years Projection horizon in years [default 100].
#' @param gd.target Target proportion of current gene diversity to retain
#' [default 0.9].
#' @param kin Kinship matrix with both dimnames identical to indNames(x), as
#' produced by gl.kin; if NULL and x is supplied, computed internally with
#' gl.kin [default NULL].
#' @param plot.display If TRUE, resultant plots are displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot.colors List of two color names, the first for the projection
#' line, the second for the target line [default c("#2171B5","#6BAED6")].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default NULL, adopting the global verbosity set by gl.set.verbosity(),
#' or 2 if no global is set].

#' @details
#' This function is the genomic analogue of the Goals screen in the PMx
#' pedigree-management software (Lacy, Ballou & Pollak 2012; PMx Users Manual
#' pp. 116-118), which projects gene-diversity retention from the pedigree and
#' solves for the population parameters needed to meet a retention goal
#' (conventionally 90 percent over 100 years). Here the starting gene diversity
#' is taken from genomic kinships rather than the pedigree: when a genlight
#' object is supplied, gd.now = 1 - mean(kin), the mean taken over the full
#' kinship matrix including the diagonal, with the kinship matrix computed by
#' gl.kin if not supplied.
#'
#' The projection applies the standard drift expectation
#' GD_t = gd.now * (1 - 1/(2*ne))^(t/gen.length) for t = 0, 1, ..., years.
#'
#' Two solve-for outputs are reported: ne.required, the effective size at which
#' the proportion retained equals gd.target exactly at the end of the horizon,
#' ne.required = 1/(2*(1 - gd.target^(gen.length/years))); and years.to.target,
#' the (continuous) time at which the proportion retained falls to gd.target at
#' the given ne, years.to.target = gen.length * log(gd.target) /
#' log(1 - 1/(2*ne)). The first whole year in the projection table at which the
#' proportion retained drops below gd.target is reported at verbose >= 3, and
#' is NA if the target is never crossed within the horizon.
#'
#' If both x and gd.now are supplied, gd.now takes precedence and the genlight
#' object is ignored, with a warning. If gd.now is supplied directly, no
#' genlight machinery (datatype check, kinship computation) is invoked.
#'
#' If a plot.file is given, the ggplot arising from this function is saved as
#' an "RDS" binary file using saveRDS(); can be reloaded with readRDS(). If a
#' plot directory (plot.dir) is specified, the ggplot binary is saved to that
#' directory; otherwise to the tempdir().

#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#' \url{https://groups.google.com/d/forum/dartr}

#' @examples
#' # From a genlight object (kinship computed internally)
#' proj <- gl.report.gd.projection(testset.gl, ne = 50, years = 100)
#' proj$summary
#' # Direct gene-diversity input, no genlight required
#' gl.report.gd.projection(gd.now = 0.3, ne = 25, years = 50)

#' @seealso \code{\link{gl.kin}}, \code{\link{gl.report.kin.groups}}

#' @export
#' @return Invisibly, a list with two components: projection, a data frame with
#' columns year, gd and prop.retained; and summary, a named list with elements
#' gd.now, ne, years.to.target and ne.required.
#'
# ----------------------
# Function
gl.report.gd.projection <- function(x = NULL,
                                    gd.now = NULL,
                                    ne,
                                    n = NULL,
                                    gen.length = 1,
                                    years = 100,
                                    gd.target = 0.9,
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
    if (verbose == 0) { plot.display <- FALSE }

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

    # RESOLVE THE SOURCE OF CURRENT GENE DIVERSITY
    # The datatype check is deferred deliberately: it runs only when a genlight
    # object is actually consumed. A direct gd.now bypasses all genlight
    # machinery.
    if (is.null(x) && is.null(gd.now)) {
        stop(error(
            "Fatal Error: Supply either a genlight object (x) or a current gene diversity (gd.now)\n"
        ))
    }
    if (!is.null(x) && !is.null(gd.now)) {
        if (verbose >= 1) {
            cat(warn(
                "  Warning: Both x and gd.now supplied; using gd.now and ignoring the genlight object\n"
            ))
        }
        x <- NULL
    }
    if (is.null(gd.now)) {
        # CHECK DATATYPE
        datatype <- utils.check.datatype(x, verbose = verbose)

        # KINSHIP MATRIX
        if (is.null(kin)) {
            if (verbose >= 2) {
                cat(report("  No kinship matrix supplied; computing with gl.kin\n"))
            }
            kin <- gl.kin(x, verbose = 0)
        }
        if (!is.matrix(kin) || nrow(kin) != ncol(kin) ||
            !identical(rownames(kin), indNames(x)) ||
            !identical(colnames(kin), indNames(x))) {
            stop(error(
                "Fatal Error: kin must be a square matrix with both dimnames identical to indNames(x)\n"
            ))
        }
        gd.now <- 1 - mean(kin)
    }

    # FUNCTION SPECIFIC ERROR CHECKING
    if (missing(ne) || !is.numeric(ne) || length(ne) != 1 || ne <= 0.5) {
        stop(error("Fatal Error: ne must be a single number greater than 0.5\n"))
    }
    if (!is.numeric(gd.now) || length(gd.now) != 1 || gd.now <= 0 || gd.now > 1) {
        stop(error("Fatal Error: gd.now must be a single value in (0, 1]\n"))
    }
    if (!is.numeric(gen.length) || length(gen.length) != 1 || gen.length <= 0) {
        stop(error("Fatal Error: gen.length must be a single positive number\n"))
    }
    if (!is.numeric(years) || length(years) != 1 || years < 1) {
        stop(error("Fatal Error: years must be a single number of 1 or more\n"))
    }
    if (years != round(years)) {
        if (verbose >= 1) {
            cat(warn("  Warning: years is not a whole number, truncating to", floor(years), "\n"))
        }
        years <- floor(years)
    }
    if (!is.numeric(gd.target) || length(gd.target) != 1 ||
        gd.target <= 0 || gd.target >= 1) {
        stop(error("Fatal Error: gd.target must be a single value in (0, 1)\n"))
    }
    if (!is.null(n)) {
        if (!is.numeric(n) || length(n) != 1 || n <= 0) {
            stop(error("Fatal Error: n, if supplied, must be a single positive number\n"))
        }
        if (ne > n && verbose >= 1) {
            cat(warn("  Warning: ne exceeds the census size n; check inputs\n"))
        }
    }

# DO THE JOB ----------------------
    if (verbose >= 2) {
        cat(report("  Projecting gene diversity over", years, "years at Ne =", ne, "\n"))
    }

    lambda <- 1 - 1 / (2 * ne)
    t <- 0:years
    gd <- gd.now * lambda^(t / gen.length)
    projection <- data.frame(year = t,
                             gd = gd,
                             prop.retained = gd / gd.now)

    # First whole year at which retention drops below target (NA if never)
    cross.year <- projection$year[projection$prop.retained < gd.target][1]

    # Solve-for outputs
    years.to.target <- gen.length * log(gd.target) / log(lambda)
    ne.required <- 1 / (2 * (1 - gd.target^(gen.length / years)))

    sumry <- list(gd.now = gd.now,
                  ne = ne,
                  years.to.target = years.to.target,
                  ne.required = ne.required)

    # Printing outputs -----------
    if (verbose >= 3) {
        cat("  Gene diversity projection\n")
        cat(paste("    Current gene diversity     :", round(gd.now, 4)), "\n")
        cat(paste("    Effective size (Ne)        :", ne), "\n")
        if (!is.null(n)) {
            cat(paste("    Census size (N), Ne/N      :", n, ",",
                      round(ne / n, 3)), "\n")
        }
        cat(paste("    Generation length (years)  :", gen.length), "\n")
        cat(paste("    Target retention           :", gd.target), "\n")
        cat(paste("    Retained at year", years, "        :",
                  round(projection$prop.retained[years + 1], 4)), "\n")
        if (is.na(cross.year)) {
            cat(paste("    Year target crossed        : not within horizon (exact:",
                      round(years.to.target, 1), "years )"), "\n")
        } else {
            cat(paste("    Year target crossed        :", cross.year,
                      "( exact:", round(years.to.target, 1), "years )"), "\n")
        }
        cat(paste("    Ne required for target     :",
                  round(ne.required, 1)), "\n")
        if (ne >= ne.required) {
            cat(paste("    The nominated Ne meets the retention target over the horizon\n"))
        } else {
            cat(paste("    The nominated Ne does NOT meet the retention target over the horizon\n"))
        }
    }

    # Plot the results -----------
    p1 <- ggplot(projection, aes(x = year, y = prop.retained)) +
        geom_line(color = plot.colors[1], linewidth = 1) +
        geom_hline(yintercept = gd.target,
                   linetype = "dashed",
                   color = plot.colors[2]) +
        xlab("Years") +
        ylab("Proportion of gene diversity retained") +
        ggtitle(paste0("Gene diversity retention (GD now = ",
                       round(gd.now, 3), ", Ne = ", ne, ")")) +
        plot.theme

    if (!is.na(cross.year)) {
        p1 <- p1 + geom_vline(xintercept = years.to.target,
                              linetype = "dotted",
                              color = plot.colors[2])
    }

    if (plot.display) { print(p1) }

    # Optionally save the plot ---------------------
    if (!is.null(plot.file)) {
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
    invisible(list(projection = projection, summary = sumry))
}
