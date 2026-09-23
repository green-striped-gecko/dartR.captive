#' @name gl.report.parent.offspring
#' @title Identifies putative parent offspring within a population
#' @description
#' This script examines the proportion of pedigree inconsistent loci, that is,
#' those loci that are homozygotes in the parent for the reference allele, and
#' homozygous in the offspring for the alternate allele. This condition is not
#' consistent with any pedigree, regardless of the (unknown) genotype of the
#' other parent. The pedigree inconsistent loci are counted, as a proportion of
#' the loci genotyped in both individuals, as an indication of
#' whether or not it is reasonable to propose the two individuals are in a
#' parent-offspring relationship.
#' @param x Name of the genlight object containing the SNP genotypes [required].
#' @param min.rdepth Minimum read depth to include in analysis [default 12].
#' @param min.reproducibility Minimum reproducibility to include in analysis
#' [default 1].
#' @param range Specifies the range to extend beyond the interquartile range for
#' delimiting outliers [default 1.5 interquartile ranges].
#' @param plot.filters Whether to show the plots of filters within the function 
#' [default FALSE].
#' @param plot_theme Theme for the plot. See Details for options
#'  [default theme_dartR()].
#' @param plot_colors List of two color names for the borders and fill of the
#'  plots [default NULL, which uses gl.colors(2)].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()]
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' If two individuals are in a parent offspring relationship, the true number of
#' pedigree inconsistent loci should be zero, but SNP calling is not infallible.
#' Some loci will be miscalled. The problem thus becomes one of determining
#' if the two focal individuals have a proportion of pedigree inconsistent
#' loci lower than would be expected of typical unrelated individuals. The
#' proportion is taken over the loci genotyped in both individuals, so
#' individuals with much missing data are not flagged merely because fewer
#' loci could be compared. Pairs are flagged when their proportion lies
#' strictly below the first quartile minus range times the interquartile
#' range of all pairs. There are some quite
#' sophisticated software packages available to formally apply likelihoods to
#' the decision, but we use a simple outlier comparison.
#' 
#' To reduce the frequency of miss-calls, and so emphasize the difference
#' between true parent-offspring pairs and unrelated pairs, the data can be
#' filtered on read depth.
#' 
#' Typically minimum read depth is set to 5x, but you can examine the
#' distribution of read depths with the function \code{\link[dartR.base]{gl.report.rdepth}}
#' and push this up with an acceptable loss of loci. 12x might be a good minimum
#' for this particular analysis. It is sensible also to push the minimum
#' reproducibility up to 1, if that does not result in an unacceptable loss of
#' loci. Reproducibility is stored in the slot \code{@other$loc.metrics$RepAvg}
#' and is defined as the proportion of technical replicate assay pairs for which
#' the marker score is consistent. You can examine the distribution of
#'  reproducibility with the function \code{\link[dartR.base]{gl.report.reproducibility}}.
#'  
#' Note that the null expectation is not well defined, and the power reduced, if
#' the population from which the putative parent-offspring pairs are drawn
#' contains many sibs. Note also that if an individual has been genotyped twice
#' in the dataset, the replicate pair will be assessed by this script as being
#' in a parent-offspring relationship.
#' 
#' The function \code{\link{gl.filter.parent.offspring}} will filter out those
#' individuals in a parent offspring relationship.
#' 
#' Note that if your dataset does not contain RepAvg or rdepth among the locus
#' metrics, the filters for reproducibility and read depth are not used.
#'  Examples of other themes that can be used can be consulted in \itemize{
#'  \item \url{https://ggplot2.tidyverse.org/reference/ggtheme.html} and \item
#'  \url{https://yutannihilation.github.io/allYourFigureAreBelongToUs/ggthemes/}
#'  }
#' @return A data frame of putative parent-offspring pairs, ordered by prop,
#' with columns: Outlier, the number of pedigree inconsistent loci; ind1 and
#' ind2, the two individuals; n.loci, the number of loci genotyped in both;
#' prop, Outlier / n.loci; zscore, the standardised prop; and p, the
#' one-sided (lower tail) normal probability of zscore. An empty data frame
#' with the same columns if no parent-offspring relationships were found.
#' @author Author(s): Arthur Georges. Custodian: Arthur Georges (Post to
#' \url{https://groups.google.com/d/forum/dartr})
#' @examples
#' if (isTRUE(getOption("dartR_fbm"))) testset.gl <- gl.gen2fbm(testset.gl)
#' out <- gl.report.parent.offspring(testset.gl[1:10, 1:100])
#' @seealso \code{\link[dartR.base]{gl.report.rdepth}} ,\code{\link[dartR.base]{gl.report.reproducibility}},
#'  \code{\link{gl.filter.parent.offspring}}
#' @family report functions
#' @importFrom stats median IQR
#' @import patchwork
#' @export

gl.report.parent.offspring <- function(x,
                                       min.rdepth = 12,
                                       min.reproducibility = 1,
                                       range = 1.5,
                                       plot.filters = FALSE,
                                       plot_theme = theme_dartR(),
                                       plot_colors = NULL,
                                       plot.dir = NULL,
                                       plot.file = NULL,
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
  # presence/absence data cannot show opposing homozygotes (0 vs 2)
  if (datatype == "SilicoDArT") {
    stop(error(
      "  Only SNP data are supported; x contains SilicoDArT data\n"
    ))
  }
  if (is.null(plot_colors)) {
    plot_colors <- gl.colors(2, verbose = 0)
  }

  # DO THE JOB

  if (verbose >= 2) {
    cat(report(
      "  Generating null expectation for the distribution of pedigree",
      "incompatibility\n"
    ))
  }
  # Filter stringently on reproducibility to minimize miscalls
  if (is.null(x@other$loc.metrics$RepAvg)) {
    if (verbose >= 1) {
      cat(warn(
        "  Dataset does not include RepAvg among the locus metrics,",
        "therefore the reproducibility filter was not used\n"
      ))
    }
  } else {
    x <-
      gl.filter.reproducibility(x,
        threshold = min.reproducibility,
        verbose = 0,
        plot.display = plot.filters
      )
  }
  # Filter stringently on read depth, to further minimize miscalls
  if (is.null(x@other$loc.metrics$rdepth)) {
    if (verbose >= 1) {
      cat(warn(
        "  Dataset does not include rdepth among the locus metrics,",
        "therefore the read depth filter was not used\n"
      ))
    }
  } else {
    x <- gl.filter.rdepth(x, lower = min.rdepth, verbose = 0,
                          plot.display = plot.filters)
  }

  # Pedigree-inconsistent loci: one individual homozygous for the reference
  # allele (0) and the other for the alternate allele (2). Counted for all
  # pairs at once with matrix products, and divided by the number of loci
  # typed in both individuals, so pairs with missing data are not flagged
  # merely because fewer loci could be compared.
  genmat <- as.matrix(x)
  typed <- !is.na(genmat)
  hom.ref <- (genmat == 0) & typed
  hom.alt <- (genmat == 2) & typed
  hom.ref[is.na(hom.ref)] <- FALSE
  hom.alt[is.na(hom.alt)] <- FALSE
  storage.mode(hom.ref) <- storage.mode(hom.alt) <- "double"
  storage.mode(typed) <- "double"
  count.mat <- tcrossprod(hom.ref, hom.alt) + tcrossprod(hom.alt, hom.ref)
  n.mat <- tcrossprod(typed)
  dimnames(count.mat) <- dimnames(n.mat) <- list(indNames(x), indNames(x))

  pairs.idx <- which(lower.tri(count.mat), arr.ind = TRUE)
  pairs <- data.frame(
    Outlier = count.mat[pairs.idx],
    ind1 = rownames(count.mat)[pairs.idx[, "row"]],
    ind2 = colnames(count.mat)[pairs.idx[, "col"]],
    n.loci = n.mat[pairs.idx],
    stringsAsFactors = FALSE
  )
  pairs$prop <- ifelse(pairs$n.loci > 0, pairs$Outlier / pairs$n.loci,
                       NA_real_)

  if (verbose >= 2) {
    cat(report(
      "  Identifying pairs with lower than expected proportions of",
      "pedigree inconsistent loci\n"
    ))
  }
  # lower outliers: strictly below the first quartile minus range x IQR
  prop <- pairs$prop
  cutoff <- stats::quantile(prop, 0.25, na.rm = TRUE) -
    range * stats::IQR(prop, na.rm = TRUE)
  pairs$zscore <- (prop - mean(prop, na.rm = TRUE)) / sd(prop, na.rm = TRUE)
  # zscore is standardised, so the p-value comes from the standard normal
  pairs$p <- stats::pnorm(pairs$zscore)
  df <- pairs[!is.na(prop) & prop < cutoff,
              c("Outlier", "ind1", "ind2", "n.loci", "prop", "zscore", "p")]
  df <- df[order(df$prop, df$Outlier), ]
  rownames(df) <- NULL

  title <- "SNP data (DArTSeq)\nProportion of pedigree incompatible loci per pair"
  counts_plot <- data.frame(prop = prop[!is.na(prop)])

  # Boxplot
  p1 <-
    ggplot(counts_plot, aes(y = prop)) +
    geom_boxplot(color = plot_colors[1], fill = plot_colors[2],
                 coef = range) +
    coord_flip() +
    plot_theme +
    xlim(range = c(-1, 1)) +
    ylab(" ") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    ggtitle(title)

  # Histogram
  p2 <-
    ggplot(counts_plot, aes(x = prop)) +
    geom_histogram(bins = 50,
                   color = plot_colors[1],
                   fill = plot_colors[2]) +
    geom_vline(xintercept = cutoff,
               color = "red",
               linewidth = 1) +
    xlab("Proportion pedigree incompatible") +
    ylab("Count") +
    plot_theme

  if (nrow(df) == 0) {
    if (verbose >= 1) cat(important("  No outliers detected\n"))
  } else {
    if (verbose >= 3) {
      print(df)
    }
  }

  # PRINTING OUTPUTS
    # using package patchwork
    p3 <- (p1 / p2) + plot_layout(heights = c(1, 4))
    print(p3)

    # Optionally save the plot ---------------------
    if (!is.null(plot.file)) {
  tmp <- utils.plot.save(p3,
    dir = plot.dir,
    file = plot.file,
    verbose = verbose
  )
  }

  # FLAG SCRIPT END

  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN
  return(df)
}
