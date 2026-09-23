#' @name gl.filter.parent.offspring
#' @title Filters putative parent offspring within a population
#' @description
#' This script removes individuals suspected of being related as
#' parent-offspring,using the output of the function
#' \code{\link{gl.report.parent.offspring}}, which examines the frequency of
#' pedigree inconsistent loci, that is, those loci that are homozygotes in the
#' parent for the reference allele, and homozygous in the offspring for the
#' alternate allele. This condition is not consistent with any pedigree,
#' regardless of the (unknown) genotype of the other parent.
#' The proportion of pedigree inconsistent loci is used as an indication of
#' whether or not it is reasonable to propose the two individuals are in a
#' parent-offspring relationship.
#' @param x Name of the genlight object containing the SNP genotypes [required].
#' @param min.rdepth Minimum read depth to include in analysis [default 12].
#' @param min.reproducibility Minimum reproducibility to include in analysis
#' [default 1].
#' @param range Specifies the range to extend beyond the interquartile range for
#'  delimiting outliers [default 1.5 interquartile ranges].
#' @param method Method of selecting the individual to remove from each pair
#' in a parent offspring relationship: 'best' removes the individual with more
#' missing genotypes (lower call rate) and keeps the other; 'random' removes
#' one at random [default 'best'].
#' @param rm.monomorphs If TRUE, remove monomorphic loci after filtering
#' individuals [default FALSE].
#' @param plot_theme Theme for the plot. See Details for options
#' [default theme_dartR()].
#' @param plot_colors List of two color names for the borders and fill of the
#'  plots [default NULL, which uses gl.colors(2)].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()]
#' @param plot.file Name for the RDS binary file to save (base name only, exclude extension) [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#'  progress log ; 3, progress and results summary; 5, full report
#'  [default 2, unless specified using gl.set.verbosity].
#' @details
#' If two individuals are in a parent offspring relationship, the true number of
#' pedigree inconsistent loci should be zero, but SNP calling is not infallible.
#' Some loci will be miscalled. The problem thus becomes one of determining if
#' the two focal individuals have a proportion of pedigree inconsistent loci
#' lower than would be expected of typical unrelated individuals. There are some quite
#' sophisticated software packages available to formally apply likelihoods to
#' the decision, but we use a simple outlier comparison.
#' 
#' To reduce the frequency of miscalls, and so emphasize the difference
#' between true parent-offspring pairs and unrelated pairs, the data can be
#' filtered on read depth. Typically minimum read depth is set to 5x, but you
#' can examine the distribution of read depths with the function
#' \code{\link[dartR.base]{gl.report.rdepth}} and push this up with an acceptable loss of
#' loci. 12x might be a good minimum for this particular analysis. It is
#' sensible also to push the minimum reproducibility up to 1, if that does not
#' result in an unacceptable loss of loci. Reproducibility is stored in the slot
#'  \code{@other$loc.metrics$RepAvg} and is defined as the proportion of
#'  technical replicate assay pairs for which the marker score is consistent.
#' You can examine the distribution of reproducibility with the function
#' \code{\link[dartR.base]{gl.report.reproducibility}}.
#' 
#' Note that the null expectation is not well defined, and the power reduced, if
#' the population from which the putative parent-offspring pairs are drawn
#' contains many sibs. Note also that if an individual has been genotyped twice
#' in the dataset, the replicate pair will be assessed by this script as being
#' in a parent-offspring relationship.
#' 
#' You should run \code{\link{gl.report.parent.offspring}} before filtering. Use
#' this report to decide min.rdepth and min.reproducibility and assess impact on
#' your dataset.
#' 
#' Note that if your dataset does not contain RepAvg or rdepth among the locus
#' metrics, the filters for reproducibility and read depth are not used.
#'
#' The pairs are identified with \code{\link{gl.report.parent.offspring}}
#' using the same arguments. An individual that belongs to several pairs is
#' removed first, because that resolves all of its pairs at once; the
#' remaining pairs are resolved in order of evidence (lowest proportion of
#' inconsistent loci first), and a pair drops out as soon as either member
#' has been removed. This keeps the number of removed individuals low.
#' 
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#' @return The genlight object without one individual from each putative
#' parent-offspring pair. If no pairs are found, the object is returned
#' unchanged.
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges -- Post to
#'  \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' out <- gl.filter.parent.offspring(testset.gl[1:10, 1:50])
#' @seealso  \code{\link[dartR.base]{gl.report.rdepth}} , \code{\link[dartR.base]{gl.report.reproducibility}},
#' \code{\link{gl.report.parent.offspring}}
#' @family filter functions
#' @importFrom stats median IQR setNames
#' @importFrom utils combn
#' @import patchwork
#' @export

gl.filter.parent.offspring <- function(x,
                                       min.rdepth = 12,
                                       min.reproducibility = 1,
                                       range = 1.5,
                                       method = "best",
                                       rm.monomorphs = FALSE,
                                       plot_theme = theme_dartR(),
                                       plot_colors = NULL,
                                       plot.file = NULL,
                                       plot.dir = NULL,
                                       verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)

  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  utils.flag.start(func = funname, verbose = verbose)

  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  method <- match.arg(method, c("best", "random"))

  # DO THE JOB
  hold <- x
  # Pairs come from gl.report.parent.offspring, so the filter always removes
  # individuals from the pairs the report lists (it also draws the plots)
  pairs <- gl.report.parent.offspring(
    x,
    min.rdepth = min.rdepth,
    min.reproducibility = min.reproducibility,
    range = range,
    plot_theme = plot_theme,
    plot_colors = plot_colors,
    plot.dir = plot.dir,
    plot.file = plot.file,
    verbose = 0
  )

  ind_to_remove <- character(0)
  if (nrow(pairs) > 0) {
    if (verbose >= 2) {
      cat(report(
        "  Selecting one individual from each pair",
        if (method == "best") "based on call rate\n" else "at random\n"
      ))
    }
    fbm <- .fbm_or_null(hold)
    missing.count <- function(ind) {
      i <- which(indNames(hold) == ind)[1]
      if (is.null(fbm)) {
        sum(is.na(as.matrix(hold[i, ])))
      } else {
        sum(is.na(hold@fbm[i, ]))
      }
    }
    # An individual in several pairs is removed first, because removing it
    # resolves all of them; otherwise pairs are resolved in order of evidence
    # (pairs is ordered by the lowest proportion of inconsistent loci first).
    # A pair drops out once either member has been removed.
    unresolved <- pairs[, c("ind1", "ind2")]
    while (nrow(unresolved) > 0) {
      n.pairs <- table(c(unresolved$ind1, unresolved$ind2))
      if (max(n.pairs) > 1) {
        candidates <- names(n.pairs)[n.pairs == max(n.pairs)]
      } else {
        candidates <- c(unresolved$ind1[1], unresolved$ind2[1])
      }
      if (method == "best") {
        # remove the candidate with more missing genotypes (lower call
        # rate); on a tie, the first candidate
        drop <- candidates[which.max(vapply(candidates, missing.count,
                                            numeric(1)))]
      } else {
        drop <- candidates[sample(length(candidates), 1)]
      }
      ind_to_remove <- c(ind_to_remove, drop)
      unresolved <- unresolved[unresolved$ind1 != drop &
                                 unresolved$ind2 != drop, , drop = FALSE]
    }

    # record only this call in the history, not the internal calls
    history <- hold@other$history
    hold <- gl.drop.ind(hold, ind.list = ind_to_remove, verbose = 0)
    if (rm.monomorphs == TRUE) {
      hold <- gl.filter.monomorphs(hold, verbose = 0)
    }
    hold@other$history <- history

    # REPORT THE RESULTS
    if (verbose >= 2) {
      cat(report("  Initial number of individuals:", nInd(x), "\n"))
      cat(report("  Individuals removed:", length(ind_to_remove), "\n"))
      cat("   ", ind_to_remove, sep = "\n    ")
      cat("\n")
    }
    if (verbose >= 3) {
      cat(report("  Pairs of individuals in a parent offspring relationship:\n"))
      print(pairs)
    }
  } else {
    if (verbose >= 1) {
      cat(important(
        "  No individuals were found to be in a parent offspring",
        "relationship, therefore the genlight object is returned unchanged\n"
      ))
    }
  }

  # ADD ACTION TO HISTORY

  nh <- length(hold@other$history)
  hold@other$history[[nh + 1]] <- match.call()

  # FLAG SCRIPT END

  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN
  invisible(hold)
}
