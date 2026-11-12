#' @name gl.report.errorrate
#' @title Estimates the genotyping error rate from replicated individuals
#' @family reporting functions
#' @description
#' Estimates the per-allele genotyping error rate (epsilon) of SNP data from
#' pairs of replicate samples of the same individual (technical replicates,
#' re-genotyped samples, clones). The estimate is made for every replicate
#' pair and for all pairs pooled, and can be used directly as the epsilon
#' of \code{\link{gl.test.offspring}}.
#'
#' @param x Name of the genlight object containing the SNP data [required].
#' @param pairs A data frame with two columns giving the replicate pairs,
#' as individual names or row indices of x. Use either pairs or replicates
#' [default NULL].
#' @param replicates A vector with one entry per individual in x giving its
#' replicate group (e.g. the sample ID shared by all replicates of one
#' individual). All pairs within each group are compared; individuals with
#' NA, and groups with one member, are ignored [default NULL].
#' @param min.loci Minimum number of SNPs genotyped in both samples of a pair
#' for the pair to be used [default 100].
#' @param alpha The confidence intervals have coverage 1 - alpha
#' [default 0.05].
#' @param plot.display If TRUE, the plot is displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plot [default theme_dartR()].
#' @param plot.colors Vector of two colours: per-pair estimates and pooled
#' estimate [default c('#2171B5', '#D7301F')].
#' @param plot.file Name for the RDS binary file to save (base name only,
#' exclude extension) [default NULL].
#' @param plot.dir Directory to save the plot RDS files [default as
#' specified by the global working directory or tempdir()].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2, unless specified using gl.set.verbosity].
#'
#' @details
#' \strong{Model.} Each allele copy is miscalled independently with
#' probability epsilon, so a true homozygote is observed as a heterozygote
#' with probability 2 epsilon (1 - epsilon) and a true heterozygote as either
#' homozygote with probability epsilon (1 - epsilon). This is the error model
#' of gl.test.offspring. Two replicates share the same true genotype, so the
#' expected frequency of each observed genotype combination is
#' P(g1, g2) = sum over t of f(t) P(g1 | t) P(g2 | t), where f(t) is the
#' frequency of true genotype t across SNPs. Epsilon and f are estimated
#' jointly by maximum likelihood from the 3 x 3 table of observed genotype
#' combinations. f is estimated from the replicates themselves, not from the
#' whole data set, because observed genotype frequencies include errors and
#' the replicated individuals may come from a different population.
#'
#' \strong{Estimates.} For every pair, epsilon is estimated from that pair's
#' table, with a profile-likelihood confidence interval. The pooled estimate
#' uses the summed table of all pairs (pairs with more SNPs get more weight).
#' The mean and median of the per-pair estimates are also reported; a mean
#' clearly above the pooled estimate, or pairs far outside the pooled
#' interval, point to a few bad samples. When a group has three or more
#' replicates, its pairs share samples and are not independent, so the
#' pooled interval is too narrow.
#'
#' \strong{Model fit.} The fit of the pooled table is tested with a G-test
#' (2 df). A poor fit has two usual causes: errors that are not symmetric
#' per allele, typically allelic dropout, which turns heterozygotes into
#' homozygotes much more often than the reverse (compare the observed and
#' expected tables); or pairs with clearly different error rates (check the
#' per-pair estimates and the plot). With many SNPs the test also picks up
#' small departures that do not matter in practice.
#'
#' \strong{Missing data.} Loci called in one replicate but missing in the
#' other are not errors in this model but are reported per pair as
#' missing_discordance (proportion of loci called in at least one sample).
#'
#' @return Returns invisibly a list with the components:
#' \itemize{
#' \item per.pair: data frame with one row per pair: ind1, ind2, group,
#' loci (SNPs called in both), discordant (SNPs with different genotypes),
#' discordance_rate, missing_discordance, epsilon, lower, upper.
#' \item pooled: data frame with the pooled epsilon and its confidence
#' interval, the mean and median of the per-pair estimates, and the number
#' of pairs and SNP comparisons.
#' \item genotype.freq: estimated true genotype frequencies (0, 1, 2).
#' \item observed: pooled 3 x 3 table of observed genotype combinations
#' (made symmetric, so each pair counts in both orders, halved).
#' \item expected: expected table under the pooled estimate.
#' \item fit: G statistic, df and p-value of the model fit.
#' \item plot: the ggplot.
#' }
#'
#' @author Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' \donttest{
#' # treat the first two individuals as replicates
#' res <- gl.report.errorrate(bandicoot.gl, pairs = data.frame(1, 2))
#' # replicate groups from a sample ID in the individual names, e.g.
#' # "AB123_rep1", "AB123_rep2"
#' # reps <- sub("_rep[0-9]+$", "", indNames(gl))
#' # res <- gl.report.errorrate(gl, replicates = reps)
#' }
#'
#' @seealso \code{\link{gl.test.offspring}}
#' @importFrom stats optim optimize uniroot qchisq pchisq median
#' @importFrom utils combn
#' @export

gl.report.errorrate <- function(x,
                                pairs = NULL,
                                replicates = NULL,
                                min.loci = 100,
                                alpha = 0.05,
                                plot.display = TRUE,
                                plot.theme = theme_dartR(),
                                plot.colors = c('#2171B5', '#D7301F'),
                                plot.file = NULL,
                                plot.dir = NULL,
                                verbose = NULL) {
  # SET VERBOSITY
  verbose <- gl.check.verbosity(verbose)
  # SET WORKING DIRECTORY
  plot.dir <- gl.check.wd(plot.dir, verbose = 0)
  # FLAG SCRIPT START
  funname <- match.call()[[1]]
  funname <- utils.flag.start(func = funname, build = "v.2023.3",
                              verbose = verbose)
  # CHECK DATATYPE
  datatype <- utils.check.datatype(x, accept = "SNP", verbose = verbose)

  # FUNCTION SPECIFIC ERROR CHECKING
  if (length(plot.colors) < 2)
    stop(error("  plot.colors must contain two colours\n"))
  if (length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop(error("  alpha must be between 0 and 1\n"))
  if (length(min.loci) != 1 || !is.finite(min.loci) || min.loci < 1)
    stop(error("  min.loci must be a positive number\n"))
  if (is.null(pairs) == is.null(replicates))
    stop(error("  Provide either pairs or replicates\n"))
  ind <- adegenet::indNames(x)
  n.ind <- adegenet::nInd(x)

  # Build the pair table
  if (!is.null(replicates)) {
    if (length(replicates) != n.ind)
      stop(error(paste0("  replicates must have one entry per individual (",
                        n.ind, ")\n")))
    grp <- as.character(replicates)
    groups <- split(seq_len(n.ind), grp)        # NA groups are dropped
    groups <- groups[lengths(groups) >= 2]
    if (!length(groups))
      stop(error("  No replicate group has two or more members\n"))
    pairs <- do.call(rbind, lapply(names(groups), function(g) {
      cc <- utils::combn(groups[[g]], 2)  # vectors of length >= 2 only
      data.frame(i1 = cc[1, ], i2 = cc[2, ], group = g,
                 stringsAsFactors = FALSE)
    }))
  } else {
    pairs <- as.data.frame(pairs, stringsAsFactors = FALSE)
    if (ncol(pairs) < 2 || nrow(pairs) < 1)
      stop(error("  pairs must be a data frame with two columns\n"))
    to.index <- function(v) {
      if (is.factor(v)) v <- as.character(v)
      if (is.character(v)) {
        idx <- match(v, ind)
        if (anyNA(idx))
          stop(error(paste0("  Individuals not found in x: ",
                            paste(utils::head(v[is.na(idx)], 5),
                                  collapse = ", "), "\n")))
      } else {
        idx <- as.integer(v)
        if (anyNA(idx) || any(idx < 1 | idx > n.ind))
          stop(error(paste0("  Individual indices must be between 1 and ",
                            n.ind, "\n")))
      }
      idx
    }
    pairs <- data.frame(i1 = to.index(pairs[[1]]), i2 = to.index(pairs[[2]]),
                        group = NA_character_, stringsAsFactors = FALSE)
  }
  if (any(pairs$i1 == pairs$i2))
    stop(error("  A pair contains the same individual twice\n"))
  if (verbose >= 2)
    cat(report("  Comparing", nrow(pairs), "replicate pairs\n"))

  # ---------------------------------------------------------------- tables
  gmat <- as.matrix(x)
  tabs <- matrix(0, nrow(pairs), 9)   # counts of (g1, g2), g1 * 3 + g2 + 1
  n.loci <- n.miss <- integer(nrow(pairs))
  for (r in seq_len(nrow(pairs))) {
    g1 <- gmat[pairs$i1[r], ]
    g2 <- gmat[pairs$i2[r], ]
    both <- !is.na(g1) & !is.na(g2)
    tabs[r, ] <- tabulate(g1[both] * 3 + g2[both] + 1, nbins = 9)
    n.loci[r] <- sum(both)
    n.miss[r] <- sum(xor(is.na(g1), is.na(g2)))
  }
  # Order within a pair is arbitrary, so make each table symmetric
  sym <- function(v) { m <- matrix(v, 3, 3, byrow = TRUE); (m + t(m)) / 2 }

  # ---------------------------------------------------------------- model
  conf <- function(eps)
    matrix(c((1 - eps)^2, 2 * eps * (1 - eps), eps^2,
             eps * (1 - eps), (1 - eps)^2 + eps^2, eps * (1 - eps),
             eps^2, 2 * eps * (1 - eps), (1 - eps)^2), 3, byrow = TRUE)
  expected.p <- function(eps, f) { E <- conf(eps); t(E) %*% (f * E) }
  softmax <- function(b) { z <- exp(c(0, b)); z / sum(z) }
  # negative log-likelihood for fixed eps, profiled over f
  nll.f <- function(b, eps, tab) -sum(tab * log(pmax(expected.p(eps, softmax(b)), 1e-300)))
  prof <- function(eps, tab, b0) {
    o <- stats::optim(b0, nll.f, eps = eps, tab = tab, method = "BFGS")
    list(value = o$value, b = o$par)
  }
  eps.max <- 0.4999
  fit.table <- function(tab) {
    marg <- (rowSums(tab) + 0.5) / (sum(tab) + 1.5)
    b0 <- log(marg[2:3] / marg[1])
    f.eps <- function(eps) prof(eps, tab, b0)$value
    o <- stats::optimize(f.eps, c(0, eps.max), tol = 1e-9)
    # optimize() never evaluates the end points; check eps = 0 (no errors)
    if (f.eps(0) <= o$objective) o <- list(minimum = 0, objective = f.eps(0))
    best <- prof(o$minimum, tab, b0)
    crit <- stats::qchisq(1 - alpha, 1) / 2
    g <- function(eps) f.eps(eps) - best$value - crit
    lower <- if (o$minimum == 0 || g(0) <= 0) 0 else
      stats::uniroot(g, c(0, o$minimum), tol = 1e-10)$root
    upper <- if (g(eps.max) <= 0) eps.max else
      stats::uniroot(g, c(o$minimum, eps.max), tol = 1e-10)$root
    list(eps = o$minimum, lower = lower, upper = upper,
         f = softmax(best$b), nll = best$value)
  }

  # ---------------------------------------------------------------- per pair
  use <- n.loci >= min.loci
  if (!any(use))
    stop(error("  No pair has", min.loci, "or more SNPs called in both samples\n"))
  if (any(!use) && verbose >= 1)
    cat(warn("  Warning:", sum(!use), "pairs have fewer than", min.loci,
             "SNPs called in both samples and are not used\n"))
  per <- data.frame(
    ind1 = ind[pairs$i1], ind2 = ind[pairs$i2], group = pairs$group,
    loci = n.loci,
    discordant = n.loci - as.integer(round(tabs[, c(1, 5, 9)] %*% c(1, 1, 1))),
    stringsAsFactors = FALSE)
  per$discordance_rate <- per$discordant / per$loci
  per$missing_discordance <- n.miss / (n.loci + n.miss)
  per$epsilon <- per$lower <- per$upper <- NA_real_
  for (r in which(use)) {
    ft <- fit.table(sym(tabs[r, ]))
    per$epsilon[r] <- ft$eps
    per$lower[r] <- ft$lower
    per$upper[r] <- ft$upper
  }
  per <- per[, c("ind1", "ind2", "group", "loci", "discordant",
                 "discordance_rate", "missing_discordance", "epsilon",
                 "lower", "upper")]

  # ---------------------------------------------------------------- pooled
  pool.tab <- sym(colSums(tabs[use, , drop = FALSE]))
  pf <- fit.table(pool.tab)
  exp.tab <- expected.p(pf$eps, pf$f) * sum(pool.tab)
  dimnames(pool.tab) <- dimnames(exp.tab) <-
    list(sample1 = c("0", "1", "2"), sample2 = c("0", "1", "2"))
  # G-test on the 6 distinct cells of the symmetric table (5 df - 3 params)
  ut <- upper.tri(pool.tab, diag = TRUE)
  o6 <- (pool.tab + t(pool.tab) * upper.tri(pool.tab))[ut]
  e6 <- (exp.tab + t(exp.tab) * upper.tri(exp.tab))[ut]
  G <- 2 * sum(ifelse(o6 > 0, o6 * log(o6 / e6), 0))
  fit <- data.frame(G = G, df = 2, p_value = stats::pchisq(G, 2, lower.tail = FALSE))
  pooled <- data.frame(
    epsilon = pf$eps, lower = pf$lower, upper = pf$upper,
    mean_pairs = mean(per$epsilon[use]),
    median_pairs = stats::median(per$epsilon[use]),
    n_pairs = sum(use), n_comparisons = sum(per$loci[use]))
  gfreq <- stats::setNames(pf$f, c("0", "1", "2"))

  # ---------------------------------------------------------------- plot
  # Columns used in aes(); defined to avoid R CMD check notes
  label <- epsilon <- lower <- upper <- NULL
  pd <- per[use, , drop = FALSE]
  pd$label <- paste(pd$ind1, "-", pd$ind2)
  pd$label <- factor(make.unique(pd$label),
                     levels = make.unique(pd$label)[order(pd$epsilon)])
  ci.txt <- paste0(round(100 * (1 - alpha)), "% CI")
  p <- ggplot2::ggplot(pd, ggplot2::aes(y = label)) +
    ggplot2::annotate("rect", xmin = pf$lower, xmax = pf$upper,
                      ymin = -Inf, ymax = Inf,
                      fill = plot.colors[2], alpha = 0.15) +
    ggplot2::geom_vline(xintercept = pf$eps, colour = plot.colors[2],
                        linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = pooled$mean_pairs,
                        colour = plot.colors[2], linetype = "dashed") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = lower, xmax = upper),
                           width = 0, colour = plot.colors[1], alpha = 0.6) +
    ggplot2::geom_point(ggplot2::aes(x = epsilon), colour = plot.colors[1],
                        size = 2.3) +
    ggplot2::labs(
      title = "Genotyping error rate from replicate pairs",
      subtitle = paste0(
        "Pooled epsilon = ", signif(pf$eps, 3), " (", ci.txt, " ",
        signif(pf$lower, 3), "-", signif(pf$upper, 3), "; solid line and band)",
        "\nMean of pairs = ", signif(pooled$mean_pairs, 3), " (dashed)"),
      x = paste0("Per-allele error rate (epsilon), points with ", ci.txt),
      y = "Replicate pair") +
    plot.theme
  if (nrow(pd) > 60)
    p <- p + ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                            axis.ticks.y = ggplot2::element_blank())

  # ---------------------------------------------------------------- report
  if (verbose >= 3) {
    cat("  Pooled per-allele error rate (epsilon):", signif(pf$eps, 4),
        paste0("(", ci.txt, " ", signif(pf$lower, 4), " - ",
               signif(pf$upper, 4), ")\n"))
    cat("  Mean / median over", sum(use), "pairs:",
        signif(pooled$mean_pairs, 4), "/", signif(pooled$median_pairs, 4), "\n")
    cat("  SNP comparisons:", sum(per$loci[use]), "\n")
    cat("  Model fit: G =", round(G, 2), "on 2 df, p =",
        signif(fit$p_value, 3), "\n")
    if (fit$p_value < 0.05)
      cat(warn("  The symmetric per-allele error model fits poorly: check for",
               "allelic dropout ($observed vs $expected) and for pairs with",
               "different error rates (plot)\n"))
  }
  if (plot.display) print(p)
  if (!is.null(plot.file))
    tmp <- utils.plot.save(p, dir = plot.dir, file = plot.file,
                           verbose = verbose)

  # FLAG SCRIPT END
  if (verbose >= 1) cat(report("Completed:", funname, "\n"))

  invisible(list(per.pair = per, pooled = pooled, genotype.freq = gfreq,
                 observed = pool.tab, expected = exp.tab, fit = fit,
                 plot = p))
}
