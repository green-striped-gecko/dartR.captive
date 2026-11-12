#' @name gl.test.offspring
#' @title Tests candidate parent pairs for offspring using Mendelian
#' mismatches and LOD scores
#' @family relatedness functions
#' @description
#' Tests whether two candidate parents can be the parents of an offspring
#' from biallelic SNP data. For each trio (parent 1, parent 2, offspring)
#' the number of SNPs whose genotypes are Mendelian-incompatible is compared
#' with the number expected from genotyping error alone, and a LOD score
#' ranks the candidate parent pairs of each offspring.
#'
#' @param x Genlight object for exploratory mode: every individual is tested
#' as offspring of every unordered pair of the other individuals. Use either
#' x or p1, p2 and off [default NULL].
#' @param p1 Genlight object of candidate first parents (e.g. fathers)
#' [default NULL].
#' @param p2 Genlight object of candidate second parents (e.g. mothers)
#' [default NULL].
#' @param off Genlight object of offspring to test [default NULL].
#' @param epsilon Per-allele genotyping error probability, applied equally
#' to parents and offspring; must be in [0, 0.5) [default 0.01].
#' @param allele.freq Allele frequencies (of the allele counted by the
#' genotype 0/1/2 coding): a single value for all SNPs, a vector with one
#' value per SNP in locus order, or a vector named by locus. If NULL, they
#' are estimated from all unique individuals in the inputs [default NULL].
#' @param alpha Significance level of the mismatch test and coverage of the
#' prediction interval (1 - alpha) [default 0.05].
#' @param min.loci Minimum number of SNPs genotyped in all three individuals
#' for a trio to be tested; trios with fewer SNPs get NA [default 100].
#' @param return.loci If TRUE, the names of the mismatching SNPs of each
#' trio are returned [default TRUE].
#' @param max.tests Maximum number of trios to test; the function stops if
#' more would be needed [default 100000].
#' @param plot.display If TRUE, the plots are displayed in the plot window
#' [default TRUE].
#' @param plot.theme Theme for the plots [default theme_dartR()].
#' @param plot.colors Vector of two colours. Checked but currently not used
#' by the plots [default c('#2171B5', '#D7301F')].
#' @param plot.file Name of the file to save the plots to with ggsave(),
#' e.g. "parentage.png" (".png" is added if no extension is given). The
#' heatmap is saved with the suffix "_heatmap" [default NULL].
#' @param plot.dir Directory to save the plot files to [default tempdir()].
#' @param plot.top Number of highest-LOD trios per offspring shown in the
#' candidate plot [default 3].
#' @param plot.max Maximum number of trios shown in the candidate plot
#' [default 30].
#' @param verbose Verbosity: 0, silent or fatal errors; 1, begin and end; 2,
#' progress log; 3, progress and results summary; 5, full report
#' [default 2].
#'
#' @details
#' \strong{Modes.} In candidate mode (p1, p2 and off given) every
#' combination of a parent from p1, a parent from p2 and an offspring from
#' off is tested, excluding combinations in which the same individual
#' (by name) appears twice. If an individual is in both p1 and p2 (e.g.
#' when sexes are unknown and p1 = p2), each unordered parent pair is tested
#' once per offspring, because swapping the two parents gives identical
#' results; such individuals must have identical genotypes in p1 and p2.
#' In exploratory mode (x given) each individual
#' is tested as offspring of each unordered pair of the remaining
#' individuals, i.e. n(n - 1)(n - 2)/2 trios, which grows quickly with n.
#' The three objects must contain the same SNPs; if all have locus names,
#' loci are matched by name, otherwise by position.
#'
#' \strong{Mismatch test.} A SNP is Mendelian-incompatible if the observed
#' offspring genotype cannot arise from the observed parental genotypes
#' (e.g. both parents 0 and offspring 1). True parents can still show
#' mismatches because of genotyping error. Errors are modelled per allele
#' with probability epsilon, independently in parents and offspring. Given
#' the observed parental genotypes, Hardy-Weinberg priors and the error
#' model, the probability that each SNP shows a mismatch under the null
#' hypothesis "these are the parents" is computed. The number of mismatches
#' then follows a Poisson-binomial distribution, which gives the expected
#' number of mismatches, an exact (1 - alpha) prediction interval and the
#' p-value P(mismatches >= observed). Trios with p-value < alpha are
#' excluded (reject_H0 = TRUE).
#'
#' The Poisson-binomial is evaluated only over the range of counts that
#' carries probability (the tail beyond it is < 1e-15), with one FFT per
#' trio, and trios are processed in vectorised batches.
#'
#' \strong{LOD.} The LOD score is the log10 ratio of the likelihood of the
#' offspring genotypes given both observed parents to their likelihood as
#' an unrelated individual drawn from the population, both under the same
#' error model. Positive values favour parentage. delta_LOD is the
#' difference to the highest LOD of the same offspring (0 for the best
#' trio). LOD ranks candidates of an offspring; it is not a universal
#' assignment threshold.
#'
#' \strong{Allele frequencies.} If allele.freq is NULL, frequencies are
#' estimated from the unique individuals in the inputs, with Beta(1,1)
#' shrinkage and a correction for genotyping error. With few individuals,
#' or when the inputs are a small family group, these estimates are
#' imprecise and affect both the expected mismatches and LOD; supplying
#' frequencies estimated from a larger population sample is recommended.
#'
#' \strong{Plots.} Two ggplots are produced. The candidate plot shows, for
#' the top plot.top LOD trios of each offspring (plus all non-excluded
#' trios, up to plot.max), the observed mismatches (coloured by LOD, shape
#' by test result), the expected mismatches (white point) and the
#' prediction interval (grey bar). The heatmap shows observed mismatches for
#' all evaluable trios, offspring by parental pair, with non-excluded trios
#' circled and the highest-LOD trio of each offspring starred.
#'
#' @return A list with the components:
#' \itemize{
#' \item results: data frame with one row per trio: test_id, offspring,
#' parent1, parent2, loci (SNPs genotyped in all three),
#' observed_mismatches, expected_mismatches, prediction_lower,
#' prediction_upper, expected_mismatch_rate, p_value, LOD, reject_H0 and
#' delta_LOD.
#' \item mismatch.loci: list with the mismatching SNP names of each trio
#' (NULL if return.loci = FALSE).
#' \item plot: the candidate ggplot (NULL if no trio could be tested).
#' \item heatmap: the heatmap ggplot (NULL if no trio has a finite LOD).
#' \item heatmap.info: list with the numbers of trios in total, with a
#' finite LOD, and plotted in the heatmap.
#' \item allele.freq: named vector of the allele frequencies used.
#' \item settings: list with epsilon, alpha, the source of the allele
#' frequencies, the mode, min.loci and the number of trios tested.
#' \item plot.info: list with the numbers of trios tested and shown in the
#' candidate plot.
#' }
#'
#' @author Author(s): Bernd Gruber. Custodian: Bernd Gruber -- Post to
#' \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' # three candidate fathers and three candidate mothers
#' x <- platypus.gl[1:6, 1:500]
#' dads <- x[1:3, ]
#' mums <- x[4:6, ]
#' # offspring of the first father with the first two mothers
#' kids <- gl.sim.offspring(fathers = dads[1, ], mothers = mums[1:2, ],
#'                          noffpermother = 2, verbose = 0)
#' res <- gl.test.offspring(p1 = dads, p2 = mums, off = kids, verbose = 0)
#' # best candidate pairs first
#' head(res$results[order(-res$results$LOD), ])
#'
#' \donttest{
#' # exploratory mode: all individuals against all pairs of the others
#' res2 <- gl.test.offspring(x = rbind(x[c(1, 4, 5)], kids[1:2]), verbose = 0)
#' }
#'
#' @seealso \code{\link{gl.report.parent.offspring}},
#' \code{\link{gl.filter.parent.offspring}}
#' @importFrom stats setNames ave mvfft nextn
#' @importFrom utils combn head
#' @export

gl.test.offspring <- function(
    x = NULL, p1 = NULL, p2 = NULL, off = NULL,
    epsilon = 0.01, allele.freq = NULL, alpha = 0.05,
    min.loci = 100L, return.loci = TRUE, max.tests = 100000L,
    plot.display = TRUE, plot.theme = NULL,
    plot.colors = c('#2171B5', '#D7301F'),
    plot.file = NULL, plot.dir = tempdir(),
    plot.top = 3L, plot.max = 30L, verbose = 2L) {
  
  if (length(epsilon) != 1L || !is.finite(epsilon) ||
      epsilon < 0 || epsilon >= 0.5)
    stop('epsilon must be a per-allele error probability in [0, 0.5)')
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop('alpha must be in (0, 1)')
  if (length(min.loci) != 1L || !is.finite(min.loci) ||
      min.loci < 1 || min.loci != as.integer(min.loci))
    stop('min.loci must be a positive integer')
  if (length(max.tests) != 1L || !is.finite(max.tests) ||
      max.tests < 1 || max.tests != as.integer(max.tests))
    stop('max.tests must be a positive integer')
  
  exploratory <- !is.null(x)
  if (exploratory) {
    if (!is.null(p1) || !is.null(p2) || !is.null(off))
      stop('Use x OR p1/p2/off, not both')
    p1 <- p2 <- off <- x
  } else if (is.null(p1) || is.null(p2) || is.null(off)) {
    stop('Provide x, or all three of p1, p2, off')
  }
  objs <- list(p1, p2, off)
  if (!all(vapply(objs, inherits, logical(1), 'genlight')))
    stop('All inputs must be genlight objects')
  if (any(vapply(objs, function(z) any(adegenet::ploidy(z) != 2L), logical(1))))
    stop('Only diploid SNPs supported')
  # Matrix dimensions, not locNames(), determine the SNP count.
  G <- lapply(objs, function(z) {
    m <- as.matrix(z)
    if (is.null(dim(m))) m <- matrix(m, nrow = 1L)
    storage.mode(m) <- 'numeric'
    m
  })
  if (any(vapply(G, nrow, integer(1)) == 0L))
    stop('Empty individual set')
  if (length(unique(vapply(G, ncol, integer(1)))) != 1L)
    stop('All objects must have the same number of SNPs')
  L <- ncol(G[[1]])
  if (L < 1L) stop('No SNPs in input')
  
  # Missing locus names are common for simulated genlight objects.
  raw.loci <- lapply(objs, adegenet::locNames)
  valid.loci <- vapply(raw.loci, function(z)
    length(z) == L && !anyNA(z) && all(nzchar(z)) &&
      !anyDuplicated(z), logical(1))
  if (all(valid.loci)) {
    loci <- raw.loci[[1]]
    if (!all(vapply(raw.loci, function(z) setequal(z, loci), logical(1))))
      stop('Locus names differ between input objects')
    for (j in 1:3)
      G[[j]] <- G[[j]][, match(loci, raw.loci[[j]]), drop = FALSE]
  } else if (!any(valid.loci)) {
    loci <- paste0('loc_', seq_len(L))
  } else {
    stop('Some objects have locus names and others do not; supply consistent names')
  }
  nm <- lapply(seq_along(objs), function(j) {
    z <- adegenet::indNames(objs[[j]])
    if (length(z) != nrow(G[[j]]) || anyNA(z) || any(!nzchar(z)))
      z <- paste0('ind_', seq_len(nrow(G[[j]])))
    if (anyDuplicated(z)) stop('Individual names must be unique within objects')
    z
  })
  if (any(vapply(G, function(z) any(!is.na(z) & !(z %in% 0:2)), logical(1))))
    stop('Genotypes must be 0, 1, 2 or NA')
  
  # Estimate allele frequencies from unique individuals across inputs.
  # Beta(1,1) shrinkage avoids extreme 0/1 estimates in small samples.
  if (is.null(allele.freq)) {
    allG <- do.call(rbind, G)
    allNames <- unlist(nm, use.names = FALSE)
    allG <- allG[!duplicated(allNames), , drop = FALSE]
    ncall <- colSums(!is.na(allG))
    qobs <- (colSums(allG, na.rm = TRUE) + 1) / (2 * ncall + 2)
    q <- (qobs - epsilon) / (1 - 2 * epsilon)
    freq.source <- 'estimated (pooled unique individuals; Beta(1,1) shrinkage)'
  } else {
    if (!is.numeric(allele.freq) || !length(allele.freq) ||
        anyNA(allele.freq) || any(!is.finite(allele.freq)))
      stop('allele.freq must contain finite numeric frequencies')
    if (length(allele.freq) == 1L) {
      q <- rep(unname(allele.freq), L)
    } else if (!is.null(names(allele.freq))) {
      if (anyDuplicated(names(allele.freq)) ||
          !all(loci %in% names(allele.freq)))
        stop('Named allele.freq must include every SNP exactly once')
      q <- unname(allele.freq[loci])
    } else if (length(allele.freq) == L) {
      q <- unname(allele.freq)
    } else {
      stop('allele.freq has ', length(allele.freq),
           ' entries; expected 1 or ', L)
    }
    freq.source <- 'user-supplied'
  }
  if (any(!is.finite(q) | q < 0 | q > 1))
    stop('Allele frequencies must be finite and in [0,1]')
  # Avoid a zero prior probability for any genotype.
  q <- pmin(pmax(q, 1e-6), 1 - 1e-6)
  
  if (exploratory) {
    n <- nrow(G[[1]])
    if (n < 3L) stop('Exploratory mode requires >= 3 individuals')
    ntests <- as.double(n) * (n - 1) * (n - 2) / 2
    if (ntests > max.tests) stop('Too many trios: ', ntests)
    cmb <- do.call(rbind, lapply(seq_len(n), function(k) {
      z <- utils::combn(setdiff(seq_len(n), k), 2)
      data.frame(i = z[1, ], j = z[2, ], k = k)
    }))
  } else {
    ntests <- prod(vapply(G, nrow, integer(1)) * 1.0)
    if (ntests > max.tests) stop('Too many candidate trios: ', ntests)
    cmb <- expand.grid(i = seq_len(nrow(G[[1]])),
                       j = seq_len(nrow(G[[2]])),
                       k = seq_len(nrow(G[[3]])))
    distinct <- nm[[1]][cmb$i] != nm[[2]][cmb$j] &
      nm[[1]][cmb$i] != nm[[3]][cmb$k] &
      nm[[2]][cmb$j] != nm[[3]][cmb$k]
    cmb <- cmb[distinct, , drop = FALSE]
    if (!nrow(cmb)) stop('No distinct-individual trios')
    # Individuals in both p1 and p2: the trios (A, B, off) and (B, A, off)
    # give identical results (the model is symmetric in the two parents),
    # so each unordered parent pair is tested only once.
    shared <- intersect(nm[[1]], nm[[2]])
    if (length(shared)) {
      same.geno <- vapply(shared, function(s)
        identical(G[[1]][match(s, nm[[1]]), ], G[[2]][match(s, nm[[2]]), ]),
        logical(1))
      if (!all(same.geno))
        stop('Individuals with the same name in p1 and p2 have different ',
             'genotypes: ', paste(head(shared[!same.geno], 5), collapse = ', '))
      n1 <- nm[[1]][cmb$i]; n2 <- nm[[2]][cmb$j]
      key <- paste(pmin(n1, n2), pmax(n1, n2), cmb$k, sep = '\r')
      dup <- duplicated(key)
      if (any(dup)) {
        cmb <- cmb[!dup, , drop = FALSE]
        if (verbose >= 2L)
          message(sum(dup), ' trios with swapped parents (individuals in both ',
                  'p1 and p2) removed; each parent pair is tested once')
      }
    }
  }
  if (verbose > 0) message('Testing ', nrow(cmb), ' trios at ', length(loci), ' SNPs')
  
  # Rows = true genotype (0,1,2); columns = observed genotype.
  e <- epsilon
  E <- matrix(c((1-e)^2, 2*e*(1-e), e^2,
                e*(1-e), (1-e)^2 + e^2, e*(1-e),
                e^2, 2*e*(1-e), (1-e)^2), nrow = 3, byrow = TRUE)
  # Hardy-Weinberg prior probabilities for the true genotype.
  prior <- rbind((1-q)^2, 2*q*(1-q), q^2)

  # Missing genotypes: work on 0-filled copies plus a 0/1 "typed" mask, so
  # all trio calculations below are plain matrix arithmetic.
  typed <- lapply(G, function(z) (!is.na(z)) * 1)
  G0 <- lapply(G, function(z) { z[is.na(z)] <- 0; z })

  # Posterior(true genotype | observed genotype, allele frequency) for every
  # individual and SNP: a list of three matrices (true genotype 0, 1, 2),
  # each individuals x SNPs.
  posterior <- function(g0) {
    n.ind <- nrow(g0)
    w <- lapply(1:3, function(t)
      matrix(E[t, g0 + 1], n.ind) * rep(prior[t, ], each = n.ind))
    tot <- w[[1]] + w[[2]] + w[[3]]
    lapply(w, function(z) z / tot)
  }
  post <- if (exploratory) rep(list(posterior(G0[[1]])), 3) else
    lapply(G0, posterior)

  # log10 P(observed offspring genotype | unrelated, HWE + error): LOD
  # denominator, individuals x SNPs (offspring object only).
  null.obs <- t(E) %*% prior
  log.null <- matrix(log10(null.obs[cbind(as.vector(G0[[3]]) + 1,
                                          rep(seq_len(L), each = nrow(G0[[3]])))]),
                     nrow(G0[[3]]))

  # Mendelian true-offspring probabilities for each true-parent pair.
  T <- array(0, c(3, 3, 3))
  incompatible <- array(FALSE, c(3, 3, 3))
  for (a in 0:2) for (b in 0:2) {
    pa <- a / 2; pb <- b / 2
    T[a+1, b+1, ] <- c((1-pa)*(1-pb),
                       pa*(1-pb)+(1-pa)*pb, pa*pb)
    incompatible[a+1, b+1, ] <- T[a+1, b+1, ] == 0
  }
  # P(observed offspring genotype | true parental genotypes).
  offspring.obs <- array(0, c(3, 3, 3))
  for (a in 1:3) for (b in 1:3)
    offspring.obs[a, b, ] <- as.vector(T[a, b, ] %*% E)
  # Incompatibility as a 9 x 3 table: row = 3 * parent1 + parent2 + 1.
  inc9 <- matrix(0, 9, 3)
  for (a in 0:2) for (b in 0:2) inc9[3 * a + b + 1, ] <- incompatible[a+1, b+1, ]

  # Trios are processed in batches: each batch is a set of
  # batch x SNP matrices, so the work is vectorised over trios and SNPs.
  ntrio <- nrow(cmb)
  batch <- max(1L, floor(1e6 / L))   # ~8 MB per batch matrix
  n.loci <- observed <- pi.low <- pi.high <- rep(NA_integer_, ntrio)
  exp.n <- pv <- lod <- rep(NA_real_, ntrio)
  bad.loci <- if (return.loci) vector('list', ntrio) else NULL

  for (b0 in seq(1L, ntrio, by = batch)) {
    rr <- b0:min(ntrio, b0 + batch - 1L)
    i <- cmb$i[rr]; j <- cmb$j[rr]; k <- cmb$k[rr]
    keep <- typed[[1]][i, , drop = FALSE] * typed[[2]][j, , drop = FALSE] *
      typed[[3]][k, , drop = FALSE]
    nk <- as.integer(rowSums(keep))
    n.loci[rr] <- nk
    ok <- nk >= min.loci
    if (!any(ok)) next
    rr <- rr[ok]; i <- i[ok]; j <- j[ok]; k <- k[ok]
    keep <- keep[ok, , drop = FALSE]

    ga <- G0[[1]][i, , drop = FALSE]
    gb <- G0[[2]][j, , drop = FALSE]
    go <- G0[[3]][k, , drop = FALSE]
    ab <- 3 * ga + gb + 1                 # observed parental combination
    nb <- nrow(ga)

    # Observed Mendelian mismatches.
    bad <- (matrix(inc9[cbind(as.vector(ab), as.vector(go) + 1)], nb) == 1) &
      keep == 1
    observed[rr] <- as.integer(rowSums(bad))
    if (return.loci)
      bad.loci[rr] <- lapply(seq_len(nb), function(z) loci[bad[z, ]])

    # P(observed offspring genotype g | observed parents, H0), g = 0, 1, 2.
    A <- lapply(post[[1]], function(z) z[i, , drop = FALSE])
    B <- lapply(post[[2]], function(z) z[j, , drop = FALSE])
    probs <- list(0, 0, 0)
    for (u in 1:3) for (v in 1:3) {
      AB <- A[[u]] * B[[v]]
      for (g in 1:3) probs[[g]] <- probs[[g]] + AB * offspring.obs[u, v, g]
    }
    # Per-SNP probability of a mismatch (relative to the OBSERVED parental
    # genotypes); 0 at SNPs not typed in all three.
    pp <- 0
    for (g in 1:3) pp <- pp + probs[[g]] * matrix(inc9[as.vector(ab), g], nb)
    pp <- pmin(pmax(pp * keep, 0), 1)
    exp.n[rr] <- rowSums(pp)

    # LOD: trio likelihood / unrelated-offspring likelihood.
    num <- (go == 0) * probs[[1]] + (go == 1) * probs[[2]] + (go == 2) * probs[[3]]
    zero <- rowSums(keep == 1 & num <= 0) > 0
    num[keep == 0] <- 1
    lb <- rowSums(log10(num)) -
      rowSums(log.null[k, , drop = FALSE] * keep)
    lb[zero] <- NA_real_
    lod[rr] <- lb

    # Poisson-binomial distribution of the number of mismatches.
    pmf <- .pb_pmf(pp)
    cdf <- apply(pmf, 2, cumsum)
    if (is.null(dim(cdf))) cdf <- matrix(cdf, nrow = 1)
    M <- nrow(cdf)
    qlow <- alpha / 2; qhigh <- 1 - alpha / 2
    # Smallest count with cdf >= quantile (same convention as poibin::qpoibin).
    pi.low[rr]  <- as.integer(colSums(cdf < qlow))
    pi.high[rr] <- as.integer(colSums(cdf < qhigh))
    obs <- observed[rr]
    pvb <- rep(0, nb)                     # beyond the computed range: < 1e-15
    in.range <- obs <= M
    pvb[obs == 0] <- 1
    idx <- which(in.range & obs > 0)
    if (length(idx)) pvb[idx] <- 1 - cdf[cbind(obs[idx], idx)]
    pv[rr] <- pmin(1, pmax(0, pvb))
  }

  out <- data.frame(
    test_id = seq_len(ntrio), offspring = nm[[3]][cmb$k],
    parent1 = nm[[1]][cmb$i], parent2 = nm[[2]][cmb$j], loci = n.loci,
    observed_mismatches = observed, expected_mismatches = exp.n,
    prediction_lower = pi.low, prediction_upper = pi.high,
    expected_mismatch_rate = exp.n / ifelse(n.loci >= min.loci, n.loci, NA),
    p_value = pv, LOD = lod,
    reject_H0 = pv < alpha, stringsAsFactors = FALSE)
  rownames(out) <- NULL
  # Best trio for each offspring has delta_LOD = 0.
  # Lower-ranked trios have negative delta_LOD.
  fin <- is.finite(out$LOD)
  out$delta_LOD <- NA_real_
  if (any(fin))
    out$delta_LOD[fin] <- out$LOD[fin] -
    ave(out$LOD[fin], out$offspring[fin], FUN = max)
  if (return.loci) {
    bad.loci <- lapply(bad.loci, function(z) if (is.null(z)) character(0) else z)
    names(bad.loci) <- paste0('test_', out$test_id)
  }

  # Display only the most promising trios per offspring by default.
  # The full, unfiltered table is always returned in $results.
  valid <- out[!is.na(out$p_value), , drop = FALSE]
  p <- NULL
  # Columns used in ggplot aes(); defined to avoid R CMD check notes.
  label <- prediction_lower <- prediction_upper <- expected_mismatches <-
    observed_mismatches <- LOD <- status <- pair <- offspring <- NULL
  plot.info <- list(total = nrow(valid), n.plotted = 0L,
                    subsetted = FALSE)
  if (nrow(valid)) {
    if (length(plot.top) != 1L || is.na(plot.top) ||
        plot.top < 1 || plot.top != floor(plot.top))
      stop("plot.top must be a positive integer")
    if (length(plot.max) != 1L || is.na(plot.max) ||
        plot.max < 1 || plot.max != floor(plot.max))
      stop("plot.max must be a positive integer")
    if (length(plot.colors) < 2L)
      stop("plot.colors must contain two colours")
    
    # Show top LOD candidates per offspring; also retain any
    # non-excluded trio (up to plot.max).
    groups <- split(seq_len(nrow(valid)), valid$offspring)
    chosen <- unique(unlist(lapply(groups, function(ii) {
      ii <- ii[order(valid$LOD[ii], decreasing = TRUE, na.last = TRUE)]
      head(ii, as.integer(plot.top))
    }), use.names = FALSE))
    chosen <- unique(c(chosen, which(!valid$reject_H0)))
    if (length(chosen) > plot.max) {
      chosen <- head(chosen[order(valid$LOD[chosen], decreasing = TRUE,
                                  na.last = TRUE)], as.integer(plot.max))
    }
    shown <- valid[chosen, , drop = FALSE]
    shown <- shown[order(shown$LOD, decreasing = FALSE), , drop = FALSE]
    shown$label <- paste0(shown$offspring, "  |  ",
                          shown$parent1, " x ", shown$parent2)
    shown$label <- make.unique(shown$label)
    shown$label <- factor(shown$label, levels = shown$label)
    shown$status <- factor(
      ifelse(shown$reject_H0, "Excluded", "Not excluded"),
      levels = c("Not excluded", "Excluded"))
    plot.info <- list(total = nrow(valid), n.plotted = nrow(shown),
                      subsetted = nrow(shown) < nrow(valid),
                      plot.top = plot.top, plot.max = plot.max)
    
    if (is.null(plot.theme)) {
      plot.theme <- if (exists("theme_dartR", mode = "function"))
        get("theme_dartR", mode = "function")() else ggplot2::theme_bw()
    }
    
    # LOD is mapped to a continuous colour scale, with the
    # exclusion decision indicated separately by point shape.
    p <- ggplot2::ggplot(shown, ggplot2::aes(y = label)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = prediction_lower, xend = prediction_upper,
                     yend = label), linewidth = 3.3,
        colour = "#AAB4BF", lineend = "round") +
      ggplot2::geom_point(
        ggplot2::aes(x = expected_mismatches),
        shape = 21, size = 2.7, fill = "white",
        colour = "#343A40", stroke = 0.8) +
      ggplot2::geom_point(
        ggplot2::aes(x = observed_mismatches, colour = LOD,
                     shape = status), size = 3.6, stroke = 1.1) +
      ggplot2::scale_colour_gradient2(
        low = "#C0392B", mid = "#A0A0A0", high = "#2166AC",
        midpoint = 0, name = "Trio LOD") +
      ggplot2::scale_shape_manual(
        values = c("Not excluded" = 16, "Excluded" = 4),
        name = "Mismatch test") +
      ggplot2::labs(
        title = "Parentage candidates: mismatch test and LOD",
        subtitle = paste0(
          "Coloured = observed; white = expected; grey = ",
          round(100 * (1 - alpha), 1), "% exact prediction interval"),
        caption = paste0(
          "Showing ", nrow(shown), " of ", nrow(valid),
          " trios; top ", plot.top,
          " LOD candidates per offspring plus non-excluded trios. ",
          "Shared per-allele error = ", epsilon,
          ". Full results in $results."),
        x = "Number of Mendelian-incompatible SNPs",
        y = "Offspring | parent 1 x parent 2") +
      ggplot2::scale_x_continuous(expand =
                                    ggplot2::expansion(mult = c(0.02, 0.06))) +
      plot.theme +
      ggplot2::theme(
        legend.position = "bottom",
        panel.grid.major.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 9),
        plot.caption = ggplot2::element_text(hjust = 0, size = 9))
    if (plot.display) print(p)
    if (!is.null(plot.file)) {
      dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
      fname <- if (grepl("\\.[[:alnum:]]+$", plot.file)) plot.file else
        paste0(plot.file, ".png")
      ggplot2::ggsave(file.path(plot.dir, fname), plot = p,
                      width = 11,
                      height = max(5, 0.32 * nrow(shown) + 2.3),
                      dpi = 300, limitsize = FALSE)
    }
  }
  # Second plot: one matrix containing ALL evaluable trios.
  # Rows are offspring, columns are unique unordered parental pairs.
  # Fill = observed Mendelian mismatches, ring = not excluded,
  # star = highest-LOD trio for each offspring.
  heatmap <- NULL
  heatmap.info <- list(total = nrow(out),
                       valid = sum(is.finite(out$LOD)),
                       plotted = 0L)
  heat.df <- out[is.finite(out$LOD), , drop = FALSE]
  if (nrow(heat.df)) {
    # Canonicalise parental pair labels, so swapped parents share a column.
    # If both orientations exist, display the best-LOD orientation.
    par1 <- as.character(heat.df$parent1)
    par2 <- as.character(heat.df$parent2)
    heat.df$pair <- paste(pmin(par1, par2), pmax(par1, par2), sep = " x ")
    heat.df <- heat.df[order(-heat.df$LOD), , drop = FALSE]
    heat.df <- heat.df[!duplicated(
      paste(heat.df$offspring, heat.df$pair, sep = "\r")), , drop = FALSE]
    heat.df$not.excluded <- !heat.df$reject_H0
    # Rows are sorted by LOD, so the first row of each offspring is its best.
    heat.df$best <- !duplicated(heat.df$offspring)
    # Most promising pairs first: number of non-excluded offspring,
    # then maximum LOD.
    pair.names <- unique(heat.df$pair)
    pair.score <- tapply(heat.df$not.excluded, heat.df$pair, sum)[pair.names]
    pair.lod <- tapply(heat.df$LOD, heat.df$pair, max)[pair.names]
    pair.names <- pair.names[order(-pair.score, -pair.lod)]
    heat.df$pair <- factor(heat.df$pair, levels = pair.names)
    heat.df$offspring <- factor(heat.df$offspring,
                                levels = rev(unique(out$offspring)))
    
    if (is.null(plot.theme)) {
      plot.theme <- if (exists("theme_dartR", mode = "function"))
        get("theme_dartR", mode = "function")() else ggplot2::theme_bw()
    }
    heatmap <- ggplot2::ggplot(
      heat.df, ggplot2::aes(x = pair, y = offspring)) +
      ggplot2::geom_tile(
        ggplot2::aes(fill = observed_mismatches),
        colour = "white", linewidth = 0.4) +
      ggplot2::geom_point(
        data = heat.df[heat.df$not.excluded, , drop = FALSE],
        shape = 21, size = 2.7, stroke = 0.8,
        fill = NA, colour = "black") +
      ggplot2::geom_point(
        data = heat.df[heat.df$best, , drop = FALSE],
        shape = 8, size = 3.2, colour = "black") +
      ggplot2::scale_fill_gradient(
        low = "#E5F4E3", high = "#D73027",
        name = "Observed mismatches") +
      ggplot2::labs(
        title = "Parentage overview: all candidate trios",
        subtitle = paste0(
          "Colour = observed Mendelian mismatches; ",
          "circle = not excluded; star = highest LOD per offspring"),
        caption = paste0(
          "LOD selects the best trio per offspring; it is not a universal ",
          "assignment threshold. All ", nrow(heat.df),
          " distinct evaluable offspring-parent-pair combinations shown."),
        x = "Candidate parental pair", y = "Offspring") +
      plot.theme +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(
          angle = 60, hjust = 1, vjust = 1, size = 8),
        axis.text.y = ggplot2::element_text(size = 9),
        panel.grid = ggplot2::element_blank(),
        plot.caption = ggplot2::element_text(hjust = 0))
    heatmap.info$plotted <- nrow(heat.df)
    if (plot.display) print(heatmap)
    if (!is.null(plot.file)) {
      dir.create(plot.dir, recursive = TRUE, showWarnings = FALSE)
      stem <- sub("\\.[[:alnum:]]+$", "", plot.file)
      ext <- if (grepl("\\.[[:alnum:]]+$", plot.file))
        sub("^.*\\.", "", plot.file) else "png"
      ggplot2::ggsave(
        file.path(plot.dir, paste0(stem, "_heatmap.", ext)),
        plot = heatmap,
        width = max(10, min(24, 0.35 * length(pair.names) + 4)),
        height = max(5, 0.38 * length(unique(heat.df$offspring)) + 3),
        dpi = 300, limitsize = FALSE)
    }
  }
  if (verbose >= 2L) message('Valid: ', nrow(valid),
                             '; not rejected: ', sum(!valid$reject_H0),
                             '; rejected: ', sum(valid$reject_H0))
  list(results = out, mismatch.loci = bad.loci, plot = p,
       heatmap = heatmap, heatmap.info = heatmap.info,
       allele.freq = setNames(q, loci),
       settings = list(epsilon = epsilon, alpha = alpha,
                       allele.freq.source = freq.source,
                       mode = if (exploratory) 'exploratory' else 'candidate',
                       min.loci = min.loci, n.tests = nrow(cmb)),
       plot.info = plot.info)
}

# Internal: Poisson-binomial probability mass functions for many trios.
# pp: trios x SNPs matrix of per-SNP success probabilities (0 = not used).
# Returns an M x trios matrix: P(X = 0), ..., P(X = M - 1). M is chosen so
# the probability of more than M - 1 successes is below 1e-15 for every
# trio, so only the useful part of the distribution is computed, instead of
# all n + 1 values. Exact up to floating-point error (checked against
# poibin::dpoibin).
# Method: characteristic function at M points on the unit circle, then one
# FFT. For SNPs with p <= 0.1, log(1 + p(z - 1)) is summed through its
# power series, so only the power sums S_k = sum(p^k) are needed (k <= 24,
# truncation error < 1e-17 per SNP). The rare SNPs with p > 0.1 are added
# directly.
.pb_pmf <- function(pp, tol = 1e-15, small = 0.1, nterm = 24L) {
  if (is.null(dim(pp))) pp <- matrix(pp, nrow = 1)
  nb <- nrow(pp)
  mu <- rowSums(pp)
  ntot <- rowSums(pp > 0)
  # Chernoff bound for sums of Bernoullis: P(X >= k) <= exp(-mu) (e mu / k)^k
  kcut <- pmax(ceiling(mu) + 1, 5)
  repeat {
    lb <- -mu + kcut * (1 + log(pmax(mu, 1e-300)) - log(kcut))
    more <- lb > log(tol) & kcut < ntot
    if (!any(more)) break
    kcut[more] <- ceiling(kcut[more] * 1.15) + 1
  }
  kcut <- pmin(kcut, ntot)
  M <- stats::nextn(max(kcut) + 1L)
  zm1 <- exp(2i * pi * (0:(M - 1)) / M) - 1

  big <- pp > small
  ps <- pp
  ps[big] <- 0
  # Power sums S_k for the small-p SNPs, trios x nterm
  S <- matrix(0, nb, nterm)
  pk <- ps
  for (k in seq_len(nterm)) {
    S[, k] <- rowSums(pk)
    if (k < nterm) pk <- pk * ps
  }
  coef <- outer(seq_len(nterm), zm1,
                function(k, z) (-1)^(k + 1) * z^k / k)
  logphi <- S %*% coef                       # trios x M
  # Large-p SNPs, added directly
  if (any(big)) {
    for (r in which(rowSums(big) > 0)) {
      p <- pp[r, big[r, ]]
      logphi[r, ] <- logphi[r, ] + colSums(log(1 + outer(p, zm1)))
    }
  }
  pmf <- Re(stats::mvfft(t(exp(logphi)))) / M  # M x trios
  pmf[pmf < 0] <- 0
  pmf
}
