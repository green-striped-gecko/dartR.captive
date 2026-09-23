#' @name gl.sim.relatedness
#' @title Simulate kinship estimates for known relationships
#' @description
#' Simulates pairs of individuals with a known relationship (full siblings,
#' half siblings or first cousins) from the individuals in a genlight object,
#' estimates their kinship with EMIBD9 (Wang, J. 2022) and summarises the
#' distribution of the estimates over replicates. The result shows which
#' kinship values that relationship produces in this dataset, to guide the
#' choice of kinship thresholds.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param rel The relationship to simulate. One of 'full.sib', 'half.sib',
#' 'first.cousin' [default 'full.sib'].
#' @param nboots The number of simulated pairs (replicates); use at least 100
#' for a stable interval [default 10].
#' @param emibd9.path The location of all necessary files to run EMIBD9
#' (read more at gl.run.EMIBD9) [required].
#' @param conf The proportion of simulated kinship values contained in the
#' reported interval, and the confidence level of the interval for the mean
#' [default 0.95].
#' @param OutAlleleFre Whether to write (1) or not (0) the EMIBD9 allele
#' frequency file [default 0].
#' @param EM_Method EMIBD9 expectation maximization method: 1, standard; 2,
#' quasi-Newton acceleration; 3, SQUAREM acceleration (see gl.run.EMIBD9)
#' [default 1].
#' @param Inbreed Whether EMIBD9 allows inbreeding when estimating IBD
#' coefficients [default FALSE].
#' @param ISeed Seed for the EMIBD9 random number generator. The simulation of
#' offspring uses R's random number generator; set it with set.seed() for
#' reproducible results [default 42].
#' @param parallel Use the parallel (MPI) version of EMIBD9. Only works for Mac
#' and Linux at the moment [default FALSE].
#' @param ncores How many cores should be used [default 1].
#' @param plot.out A boolean that indicates whether to plot the results
#' [default TRUE].
#' @param plot.dir Directory to save the plot RDS files [default as specified
#' by the global working directory or tempdir()]
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL]
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#' @details
#' In each replicate, parents are drawn at random, without replacement, from
#' the individuals in x, and the relatives are simulated with
#' \code{dartR.sim::gl.sim.offspring}:
#' \itemize{
#' \item full.sib: two offspring of the same two parents.
#' \item half.sib: one offspring of each of two mothers with the same father.
#' \item first.cousin: two full siblings, each mated with a further
#' (different) individual; the two resulting offspring are first cousins.
#' }
#' The simulated individuals are added to x and kinship is estimated with
#' \code{\link{gl.run.EMIBD9}}, so allele frequencies come from the sample.
#'
#' The interval reported is the range that contains the central conf
#' proportion of the simulated kinship values (e.g. the 2.5\% and 97.5\%
#' quantiles for conf = 0.95). It describes the spread of individual pairs and
#' is the one to compare with thresholds. The confidence interval of the mean
#' is reported separately; it narrows as nboots grows.
#'
#' Below is a table modified from Speed & Balding (2015) showing kinship
#' values, and their confidence intervals (CI), for different relationships.
#'
#' \tabular{lll}{
#'   \strong{Relationship} \tab \strong{Kinship} \tab \strong{95\% CI} \cr
#'   Identical twins / clones / same individual \tab 0.5   \tab -              \cr
#'   Sibling / Parent-Offspring                \tab 0.25  \tab (0.204, 0.296)\cr
#'   Half-sibling                              \tab 0.125 \tab (0.092, 0.158)\cr
#'   First cousin                              \tab 0.062 \tab (0.038, 0.089)\cr
#'   Half-cousin                               \tab 0.031 \tab (0.012, 0.055)\cr
#'   Second cousin                             \tab 0.016 \tab (0.004, 0.031)\cr
#'   Half-second cousin                        \tab 0.008 \tab (0.001, 0.020)\cr
#'   Third cousin                              \tab 0.004 \tab (0.000, 0.012)\cr
#'   Unrelated                                 \tab 0     \tab -              \cr
#' }
#'
#' @return Invisibly, a named list:
#' \itemize{
#' \item values -- the simulated kinship values, one per replicate.
#' \item mean -- their mean.
#' \item interval -- the lower and upper limits containing the central conf
#' proportion of the values.
#' \item ci.mean -- the confidence interval (level conf) of the mean.
#' \item plot -- a histogram of the values with the mean and the interval.
#' }
#' @author Author(s): Sam Amini. Custodian: Sam Amini -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#' @examples
#' \dontrun{
#' # To run this function needs EMIBD9 installed in your computer
#' t1 <- gl.filter.allna(platypus.gl)
#' set.seed(1)
#' res <- gl.sim.relatedness(t1, rel = "half.sib", nboots = 100,
#'                           emibd9.path = "path/to/emibd9")
#' res$interval
#' }
#'
#' @references
#' \itemize{
#' \item Wang, J. (2022). A joint likelihood estimator of relatedness and allele
#'  frequencies from a small sample of individuals. Methods in Ecology and
#'  Evolution, 13(11), 2443-2462.
#' \item Speed, D., Balding, D. (2015). Relatedness in the post-genomic era: is
#'  it still useful? Nature Reviews Genetics 16, 33-44.
#' }
#' @family captive management
#' @importFrom stringr str_split
#' @importFrom dartR.sim gl.sim.offspring
#' @importFrom stats confint lm na.omit
#' @export

gl.sim.relatedness <- function(x,
                               rel = "full.sib",
                               nboots = 10,
                               emibd9.path = getwd(),
                               conf = 0.95,
                               OutAlleleFre = 0,
                               EM_Method = 1,
                               Inbreed = FALSE,
                               ISeed = 42,
                               parallel = FALSE,
                               ncores = 1,
                               plot.out = TRUE,
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
  rel <- match.arg(rel, c("full.sib", "half.sib", "first.cousin"))
  if (datatype == "SilicoDArT") {
    stop(error(
      "  Only SNP data are supported; x contains SilicoDArT data\n"
    ))
  }
  # first cousins need four distinct parents
  if (nInd(x) < 4) {
    stop(error("  x must contain at least 4 individuals\n"))
  }

  # DO THE JOB

  # simulated individuals get their own names, so they can be found in the
  # EMIBD9 output and never clash with the names in x
  offspring <- function(father, mother, n, prefix) {
    # gl.sim.offspring prints a warning about missing data at every call
    utils::capture.output(
      off <- dartR.sim::gl.sim.offspring(father, mother, noffpermother = n,
                                         sexratio = 0.5, verbose = 0)
    )
    indNames(off) <- paste0(prefix, seq_len(nInd(off)))
    off
  }

  # kinship between two simulated relatives, estimated on the sample plus the
  # simulated individuals so that allele frequencies come from the sample
  kinship.pair <- function(ppoff, id1, id2) {
    res <- gl.run.EMIBD9(ppoff,
                         emibd9.path = emibd9.path,
                         OutAlleleFre = OutAlleleFre,
                         EM_Method = EM_Method,
                         Inbreed = Inbreed,
                         ISeed = ISeed,
                         parallel = parallel,
                         ncores = ncores,
                         plot.out = FALSE,
                         verbose = 0)
    res$rel[id1, id2]
  }

  simulate.pair <- function() {
    if (rel == "full.sib") {
      parents <- sample(nInd(x), 2)
      sibs <- offspring(x[parents[1], ], x[parents[2], ], 2, "simrel_fs")
      kinship.pair(rbind(x, sibs), "simrel_fs1", "simrel_fs2")
    } else if (rel == "half.sib") {
      parents <- sample(nInd(x), 3)
      off1 <- offspring(x[parents[3], ], x[parents[1], ], 1, "simrel_hsA")
      off2 <- offspring(x[parents[3], ], x[parents[2], ], 1, "simrel_hsB")
      kinship.pair(rbind(x, off1, off2), "simrel_hsA1", "simrel_hsB1")
    } else {
      # the two cousins' other parents are distinct from the grandparents
      parents <- sample(nInd(x), 4)
      sibs <- offspring(x[parents[1], ], x[parents[2], ], 2, "simrel_sib")
      cousin1 <- offspring(sibs[1, ], x[parents[3], ], 1, "simrel_fcA")
      cousin2 <- offspring(x[parents[4], ], sibs[2, ], 1, "simrel_fcB")
      kinship.pair(rbind(x, sibs, cousin1, cousin2),
                   "simrel_fcA1", "simrel_fcB1")
    }
  }

  values <- numeric(nboots)
  for (b in seq_len(nboots)) {
    if (verbose >= 2) {
      cat(report("  Simulating", rel, "pair", b, "of", nboots, "\n"))
    }
    values[b] <- simulate.pair()
  }

  mean_kin <- mean(values)
  probs <- c((1 - conf) / 2, 1 - (1 - conf) / 2)
  interval <- stats::quantile(values, probs, names = FALSE)
  names(interval) <- c("lower", "upper")
  if (nboots > 1) {
    ci.mean <- confint(lm(values ~ 1), level = conf)[1, ]
  } else {
    ci.mean <- c(NA_real_, NA_real_)
  }
  names(ci.mean) <- c("lower", "upper")

  # PRINTING OUTPUTS
  rel.label <- c(full.sib = "full sibling", half.sib = "half sibling",
                 first.cousin = "first cousin")[[rel]]
  kinship <- NULL
  p1 <- ggplot(data.frame(kinship = values), aes(x = kinship)) +
    geom_histogram(binwidth = 0.005, colour = "black") +
    geom_vline(xintercept = mean_kin, color = "red", linewidth = 1) +
    geom_vline(xintercept = interval, color = "green", linewidth = 1,
               linetype = 2) +
    labs(y = "Count", x = "Kinship") +
    ggtitle(paste0("Simulated ", rel.label, " kinship\n(mean in red, ",
                   conf * 100, "% of values between the green lines)")) +
    theme_dartR() +
    theme(plot.title = element_text(hjust = 0.5))
  if (plot.out) {
    print(p1)
  }

  # Optionally save the plot
  if (!is.null(plot.file)) {
    tmp <- utils.plot.save(p1,
                           dir = plot.dir,
                           file = plot.file,
                           verbose = verbose)
  }

  if (verbose >= 3) {
    cat(report("  Simulated", rel.label, "kinship,", nboots, "replicates:\n"))
    cat("    Mean:", round(mean_kin, 4), "\n")
    cat("   ", paste0(conf * 100, "% of values between:"),
        round(interval, 4), "\n")
    cat("   ", paste0(conf * 100, "% CI of the mean:"), round(ci.mean, 4),
        "\n")
  }

  # FLAG SCRIPT END
  if (verbose >= 1) {
    cat(report("Completed:", funname, "\n"))
  }

  # RETURN
  invisible(list(values = values,
                 mean = mean_kin,
                 interval = interval,
                 ci.mean = ci.mean,
                 plot = p1))
}
