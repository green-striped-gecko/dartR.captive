# dartR.captive NEWS

## dartR.captive 1.2.7

### New functions

* New `gl.test.offspring()` tests whether two candidate parents can be the
  parents of an offspring. For each trio the number of Mendelian-incompatible
  SNPs is compared with the number expected from genotyping error alone
  (per-allele error `epsilon`, Poisson-binomial distribution, exact
  prediction interval and p-value), and a LOD score ranks the candidate
  parent pairs of each offspring. Candidate mode tests given fathers,
  mothers and offspring (an individual in both parent sets is tested once
  per unordered pair); exploratory mode tests every individual against all
  pairs of the others. Returns the full table, the mismatching SNPs of each
  trio, a candidate plot and a heatmap of all trios.

* New `gl.report.errorrate()` estimates the per-allele genotyping error
  rate (the `epsilon` of `gl.test.offspring()`) from replicate samples of
  the same individual, given as pairs or as replicate groups. Error rate
  and true genotype frequencies are estimated jointly by maximum
  likelihood, per pair and pooled, with profile-likelihood confidence
  intervals; a G-test of the model fit flags allelic dropout and pairs with
  different error rates. Missing-data discordance is reported separately.

* New `gl.inbreeding()` (contributed by Luis Mijangos) estimates the
  inbreeding coefficient of each individual with PLINK 1.9 (`--het`, F) and
  GCTA (`--ibc`, Fhat1, Fhat2 and Fhat3), and plots PLINK F and Fhat3 by
  population. Both programs read the same loci: monomorphic loci are
  removed first, and all loci are treated as autosomal. Allele frequencies
  come from all individuals pooled, or from each population with
  `by.pop = TRUE`. `method = "plink"` runs without GCTA. `ind.labels = TRUE`
  shows the names of the individuals on the vertical axis.

* New `gl.parentage.power()` estimates, by simulation, how reliably a SNP
  panel assigns offspring to candidate parents across a range of locus
  counts, by Mendelian exclusion or with COLONY. Pair mode checks every
  candidate pair, so a correct unique assignment means that no other pair
  is compatible. Single-parent mode (`pairs = FALSE`) models offspring whose
  other parent was not sampled. With `engine = "colony"`, `n.rep` runs
  COLONY once per replicate.

* New `gl.relatedness()` estimates pairwise relatedness, and optionally
  inbreeding, with Jinliang Wang's Coancestry estimators: the moment
  estimators of Wang, Lynch & Li, Lynch & Ritland, Ritland, Queller &
  Goodnight and Loiselle, and the likelihood estimators DyadML and TrioML,
  with optional bootstrap confidence intervals. The numerical engine is the
  closed-source package `dartR.coancestry`, installed separately. A pair of
  individuals that shares no called locus gets `NA` for every estimator,
  with a warning. The matrices are on the relatedness scale (twice kinship)
  and are tagged `attr(, "scale") = "relatedness"`. `wang` and `lynchli`
  sit above `related::coancestry()` by a small-sample term (+0.047 and
  +0.018 on 30 platypus), while `lynchrd`, `ritland` and `quellergt` match
  it.

### Changes and bug fixes

* History: functions called through `do.call(fun, list(x = gl, ...))` no
  longer store the whole input object in `@other$history` (the calls are
  compacted with `dartR.base::utils.history.call()`), and no longer stop
  with "cannot coerce type 'closure'" at `verbose >= 1`. Needs dartR.base
  1.2.7. dartR.data (>= 1.2.5) is now required, for `testset2.gl` in the
  examples.

* `gl.diagnostics.relatedness()`:
  * Simulations now mirror the input when no variable files are given: the
    loci of `x` with their allele frequencies per population (same number
    of loci, no extra neutral loci), the populations of `x` with its sample
    sizes, and parents that mate with several partners so that half sibs
    occur. The previous defaults simulated 100 loci at frequency 0.5 in one
    population of 50, whatever `x` was. `run_sim = TRUE` works without
    variable files.
  * Simulates the effective population size of `x`: new arguments `Ne`
    (one value or one per population) and `neest.path` (estimates Ne with
    `dartR.popgen::gl.LDNe()` when `Ne` is not given; the function stops
    when Ne is neither given nor estimable). Each population gets a census
    size above 2 Ne and each generation stores a sample with the sample
    sizes of `x`. The pedigree kinship is built from the full simulated
    pedigree, so relatives linked through unsampled ancestors are
    recognised. A warning flags `gl.LDNe` estimates of Ne that are
    unreliable (infinite jackknife upper limit or more than 10 times the
    sample size), and a failure of `gl.LDNe` is reported as a request for
    `Ne`.
  * Can recreate samples made of families, a source of bias when allele
    frequencies come from a few families: new `families` (an ind.metrics
    column of full-sib family identifiers, or "colony" to reconstruct them
    with `gl.run.colony()`, with `colony.path`) and `familyParents` (`x`
    also holds the parents). Every stored generation from generation 1 on
    is then sampled with the family sizes of `x`, and generation 0 is left
    out of the analysis.
  * Estimates relatedness within each stored generation by default (new
    `analysisUnit = "generation"`), so the estimators take their allele
    frequencies from a sample the size of `x`, as they do on the real data;
    `analysisUnit = "pooled"` keeps the previous behaviour and is needed
    for pairs across generations. `@MergedDf` gains a `generation` column.
    A single population no longer stops the simulation.
  * Computes the `which_tests` estimators with `gl.relatedness()` instead
    of `related::coancestry()`, so it needs the `dartR.coancestry` engine
    rather than the `related` package. Every estimator except the slow
    `trioml` runs by default. `lynchrd`, `quellergt` and rrBLUP results are
    unchanged; `wang` differs (on 15 individuals of testset.gl, 0.045 higher
    on average on the kinship scale). `loiselle` is now available, and an
    unknown estimator in `which_tests` is an error.
  * Accuracy statistics are computed against the exact pedigree kinship of
    each pair (new column `rel`; inbreeding included), with one
    relationship class per pair (`RelDegree`). New classes
    grandparent_grandchild, avuncular, great_grandparent_grandchild,
    half_avuncular, other_relatives and unrelated (pedigree kinship 0);
    these pairs were labelled unrelated or discarded. Previously pairs
    below 0.05 were discarded before any summary, founders in an attached
    pedigree were labelled half siblings, some pairs were counted under two
    relationships, and "RMSE" was the mean absolute error. RMSE and
    variance values change. New `biasOut` reports the bias of each
    estimator by relationship class, and a tile plot of bias and RMSE is
    added to `@plotList`.
  * The simulated genotypes take the missing-data pattern of `x` (new
    `simMissing`, default TRUE), and `@SimOutput` keeps the parents and
    generation of each individual. Founder inbreeding is measured on loci
    called in at least 99% of individuals and carried into the simulation;
    with a dartR.sim version that stores the founders, they are kept as
    generation 0 and their inbreeding enters the pedigree kinship.
    SilicoDArT input is an error.

* `gl.grm()`:
  * Passes `min.MAF = 1/(2n) - 1e-10` to `rrBLUP::A.mat()` unless `min.MAF`
    is given. With the default `1/(2n)`, a locus with a single minor-allele
    copy sat exactly on the cut-off, and rounding decided whether it was
    kept: Apple Silicon Macs dropped it, and Linux and Windows kept it
    unless the minor allele was the one counted as 2. Such loci are now
    kept on every platform. Results change mostly on Apple Silicon (on
    `platypus.gl[1:12, 1:200]` the T27 diagonal goes from 0.943 to 0.986).
    Same fix as `dartR.spatial::gl.grm2()`.
  * Output is tagged `attr(, "scale") = "relatedness"`. Documentation
    corrected: the diagonal can be below 1, and the off-diagonal mean is
    `-mean(diag)/(n - 1)` with `n` individuals, not loci.
  * Errors on SilicoDArT (presence/absence) input instead of returning a
    meaningless matrix.
  * Fixed a crash (`object 'p3' not found`) when `plotheatmap = FALSE` was
    combined with a `plot.file`; the matrix is returned with a warning that
    nothing was saved. Corrected the documented default of `legendy` (0.5).

* `gl.grm()` and `gl.grm.network()`: colours given to `palette_discrete`
  with population names now go to the populations they name; they were
  assigned by position.

* `gl.grm.network()`:
  * Kinship is computed as G / 2 when `standardise = FALSE` (default);
    relatedness values were treated as kinship, so thresholds and
    categories were off by a factor of two. A `G` tagged
    `attr(G, "scale") = "kinship"` (as returned by `gl.run.EMIBD9()`) is
    used without halving.
  * `standardise = TRUE` follows Goudet et al. (2018): kinship relative to
    the average pair, `(theta - mean) / (1 - mean)`. It subtracted the mean
    inbreeding coefficient from every kinship, so after call-rate filtering
    parent-offspring pairs came out near 0.1 instead of 0.25.
  * Fixed crashes when exactly one pair cleared `kinship.threshold`, and
    when `categorise = TRUE` was combined with a `kinship.threshold` below
    0.1 (the undocumented "First Cousins" category is removed). An
    individual related to two or more others is plotted once instead of
    as overlapping duplicates. `categorise` colours follow the documented
    order.
  * `G` must match `indNames(x)`; an invalid `method` is an error. Returns a
    named list (`plot`, `kinship`). Corrected the documented defaults of
    `node.label.size` and `title`.

* `gl.plot.network()`: new argument `type = c("similarity", "distance")`;
  with `type = "distance"` the closest pairs (lowest values) are drawn.
  Previously the largest values were always drawn, i.e. the most distant
  pairs of a distance matrix; the default keeps the previous behaviour for
  similarity matrices such as `gl.grm()`. Link widths reflect each drawn
  link's own value (they were taken from other pairs). `x = NULL` works as
  documented, D is checked against `indNames(x)`, and an invalid `method`
  is an error. Returns the drawn igraph network invisibly (was NULL). Much
  faster for many individuals.

* `gl.sim.relatedness()`: `rel = "full.sib"` returned the average
  parent-offspring kinship; it now simulates two full siblings. The
  reported interval was the confidence interval of the mean; it is now the
  range holding the central `conf` proportion of simulated kinship values,
  with the CI of the mean returned separately. Returns a list (`values`,
  `mean`, `interval`, `ci.mean`, `plot`) instead of printing, no longer
  requires the `related` package, and runs EMIBD9 quietly.

* `gl.report.parent.offspring()`: pairs are assessed on the proportion of
  pedigree-inconsistent loci among loci genotyped in both individuals (new
  columns `n.loci` and `prop`), so individuals with much missing data are
  no longer flagged as relatives of many others. p-values were computed on
  the wrong scale and are now `pnorm(zscore)` (larger than before). `range`
  below 1.5 flags more pairs, as documented. SilicoDArT input is an error.
  Faster on large data sets; nothing printed at `verbose = 0`.

* `gl.filter.parent.offspring()`: `method = "best"` removed the member of
  each pair with fewer missing genotypes (counts compared as text); it now
  keeps the better-genotyped individual, as documented. Pairs come from
  `gl.report.parent.offspring()`, `range` is honoured by both methods (and
  `range = 3` no longer errors), and SilicoDArT input is an error.
  Individuals in several pairs are removed first, so fewer individuals are
  removed. One history entry per call; an unknown `method` is an error.

* `gl.run.EMIBD9()`: a run in which EMIBD9 fails stops with EMIBD9's
  message; previously the results of the previous call in the same session
  could be returned. `OutAlleleFre = TRUE` no longer crashes EMIBD9. `raw`
  and `processed` columns are numeric (were character). SilicoDArT input is
  an error. `rel` is documented as the kinship coefficient (it always was)
  and carries `attr(rel, "scale") = "kinship"`. Building the input file is
  faster for many loci, and EMIBD9's console output prints only at
  `verbose >= 2`.

* `gl.run.colony()`: COLONY's output files are written to `outpath` instead
  of the R working directory, and the function returns
  `list(files, best.config)`, where `best.config` is the COLONY best
  configuration as a data frame (it returned `outpath`). The default
  `outpath = NULL` works, paths with spaces are supported, and a run in
  which COLONY rejects its input stops with COLONY's error message.
  COLONY's console output prints only at `verbose >= 2`. Individual names
  longer than COLONY's 20-character limit or containing spaces are replaced
  by short IDs for the run and restored in `best.config`.

* `gl2colony()`: stops with an error, instead of writing a file COLONY
  rejects or misreads, when individual names contain whitespace, the data
  are SilicoDArT, or `sibship.prior`, `known.allele.freq`, or any
  known/excluded count is non-zero. A missing offspring/mother/father
  column no longer resets the other two; column names and values are
  matched ignoring case. A single `allelic.dropout`/`other.typ.err` value
  without `@` applies to all loci. Progress messages follow `verbose`.

* `gl.assign.grm()`: errors clearly when `pop(x)` is `NULL` and on duplicate
  individual names (two individuals were folded into population
  'unknown'). The documentation now describes the actual output, a named
  numeric vector of mean pairwise relatedness scores, not assignment
  probabilities. Prints a summary at `verbose >= 3` and a completion
  message at `verbose >= 1`.

### Moved and removed functions

* `gl.assign.pca()`, `gl.assign.pa()` and `gl.assign.on.genotype()` moved
  to dartR.popgen, which has the reviewed versions.

* `gl.assign.mahal()` and `gl.assign.mahalanobis()` removed from
  dartR.captive: two versions of the same function. The reviewed
  `gl.assign.mahalanobis()` is in dartR.popgen.

* `utils.assignment()`, `utils.assignment_2()`, `utils.assignment_3()` and
  `utils.assignment_4()` removed: four revisions of the same helper, called
  by no function.
