# dartR.captive NEWS

## Unreleased

* `gl.relatedness()` (new on `dev`; Coancestry estimators through the
  closed-source `dartR.coancestry` engine): the engine is fetched with
  `getExportedValue()`, so R CMD check no longer warns about an
  undeclared dependency. A pair of individuals that shares no called
  locus now gets `NA` for every estimator (the engine returned 0, i.e.
  "unrelated", for the moment estimators and `NaN` for the likelihood
  ones), with a warning; with `n.boots > 0`, individuals dropped by the
  pre-filter keep `NA` rows, so the matrices always cover `indNames(x)`.
  The matrices are tagged `attr(, "scale") = "relatedness"`, so the
  kinship functions halve them when they are passed as `kin`. Arguments
  are checked before the engine runs (an invalid `plot.stat` used to stop
  only after the run). The heatmap is saved whenever `plot.file` is
  given; the result is returned invisibly. Details document that `wang`
  and `lynchli` sit above `related::coancestry` by a small-sample term
  (+0.047 and +0.018 on 30 platypus), while `lynchrd`, `ritland` and
  `quellergt` match it.
* `gl.select.pairs()`: a warning when any individual has call rate below
  0.8, and Details on why: missing genotypes shrink kinship while `f.max`
  is fixed, so on the `testset2.gl` captive colony one of the 10 pairs
  selected unfiltered has kinship 0.185 after filtering loci at call rate
  0.95; the example now filters first. Missing kinship values are ignored
  in the means with a warning, and a pair with `NA` kinship is never
  selected (one `NA` crashed the `dynamic` scheme and gave `NA` gene
  diversity and a wrong order in `static` and `ranked`); an `NA` sex
  counts as unknown (it crashed). Invalid `scheme`, `max.per.sire`,
  `max.per.dam`, `f.max` or `n.pairs` are now errors (they were reset to
  defaults). The `dynamic` scheme scores candidates in closed form, with
  identical selections: 0.7 s instead of 164 s for 200 individuals.
  Details now say that `dgd.cum` can fall as pairs are added.
* `gl.report.repro.targets()`: a warning when any individual has call
  rate below 0.8 (filtering loci at call rate 0.95 changes 9 of 24
  targets on the `testset2.gl` captive colony; the example now filters
  first). Missing kinship values are ignored in the mean kinship with a
  warning, and an `NA` sex counts as unknown (each gave "missing value
  where TRUE/FALSE needed"). An invalid `n.target` (0, negative, `NA`,
  non-numeric) is now an error (it was reset to `nInd(x)`). The
  description now says the targets lower the expected mean kinship of
  the next generation (0.0582 on the colony, against 0.0667 for equal
  contributions and 0.0560 for optimal contributions) rather than
  equalising founder representation.
* `gl.report.mate.suitability()`: a warning when any individual has call
  rate below 0.8, and Details on why: missing genotypes shrink kinship
  while the No Way point is fixed, so on the `testset2.gl` captive colony
  17 of the 48 pairings that are 'NoWay' after filtering loci at call rate
  0.95 are rated 4-6 unfiltered; the example now filters first. Missing
  kinship values are ignored in the means with a warning, and a pairing
  with `NA` kinship gets MSI `NA` (one `NA` gave "missing value where
  TRUE/FALSE needed"); an `NA` sex counts as unknown (it gave "subscript
  out of bounds"); an out-of-range `f.noway` is now an error (it was
  reset to 0.125). Delta gene diversity uses a closed form with identical
  results: 0.04 s instead of 1.7 s for 200 individuals.
* `gl.report.kin.sets()`: a single-population `x` now needs kinship
  estimated on a wider reference (`kin = NULL` or kinship estimated on
  `x` gave `GD.w = 1` and `meanMK = 0` by construction; multi-population
  calls are unchanged). Missing kinship values are ignored in the block
  means with a warning (one `NA` pair voided a population's statistics
  and its 60 Fst values). A warning when individuals have call rate
  below 0.8. Details now describe the upward self-kinship term in the
  kinship Fst (random halves of one population give 0.0098) and the
  SilicoDArT `meanF` of 0; `@param kin` follows the #113 wording.
* `gl.report.gd.projection()`: `summary` gains `years.to.target.source`
  and `ne.required.source`, the PMx goal read as a proportion of SOURCE
  gene diversity (after Soule et al. 1986). The existing
  `years.to.target` and `ne.required` are relative to current gene
  diversity and do not depend on the data. For the `testset2.gl` captive
  colony at Ne = 50 the source-relative outputs are 3.7 years and
  Ne 1360 (current-relative: 10.5 years, Ne 475); they are `NA`, with a
  warning, when gene diversity is already at or below the target. The
  plot adds the source-relative target line. Missing kinship values are
  ignored in `gd.now` with a warning (one `NA` gave an opaque error), a
  warning when individuals in `x` have call rate below 0.8, and scalar
  arguments are checked before any kinship work.
* `gl.report.ind.remove()`: missing kinship values are ignored in the
  gene diversity means with a warning (one `NA` crashed the function
  with "missing value where TRUE/FALSE needed"); an invalid `n.best`
  (0, negative, non-numeric) is now an error (it removed without limit:
  `n.best = 0` returned 7 removals on the captive colony); a warning when
  any individual has call rate below 0.8 (the greedy set goes from 7 to
  11 animals after filtering loci at call rate 0.95). The greedy search
  now keeps running sums, with identical results: 0.08 s instead of 1.6 s
  for 200 individuals.
* `gl.report.ind.move()`: missing kinship values are ignored in the gene
  diversity means with a warning (one `NA` pair inside a population made
  every move into or out of it `NA`: 70 of 92 moves in the example); a
  warning when any individual has call rate below 0.8. Details now say
  what `net` measures: an unweighted sum of two populations' GD changes
  that favours small destinations and genetically divergent or admixed
  movers, and says nothing about genetic integrity.

* `gl.report.ind.add()`: missing kinship values are ignored in the gene
  diversity means with a warning (one `NA` in the target block made
  every `dgd` `NA`); duplicated candidates are evaluated once; a warning
  when any target or candidate individual has call rate below 0.8
  (on `testset2.gl`, 5 of the top 10 candidates change after filtering
  loci at call rate 0.95). `utils.kin.dgd()` gains `na.rm` (default
  `FALSE`, so its other callers are unchanged).

* `gl.report.kin.groups()` now needs kinship estimated on a wider
  reference than `x` (as `gl.report.ind.remove()` and the other
  functions in #113): group MK averages over all individuals, so with
  kinship estimated on `x` itself (including `kin = NULL`) every group
  had MK = 0 and GD = 1. The #113 audit listed this function as
  unaffected; that was wrong. Missing kinship values are ignored in the
  block means with a warning (one NA made MK and GD NA), and individuals
  with call rate below 0.8 give a warning: missing genotypes pull kinship
  and inbreeding toward 0 (captive cohorts in `testset2.gl`: meanF -0.25
  to -0.38, and -0.15 to +0.01 after filtering loci at call rate 0.95).
* `gl.grm()` passes `min.MAF = 1/(2n) - 1e-10` to `rrBLUP::A.mat()` unless
  `min.MAF` is given. With the default `1/(2n)`, a locus with a single
  minor-allele copy sat exactly on the cut-off, and rounding decided whether
  it was kept: Apple Silicon Macs dropped it, and Linux and Windows kept it
  unless the minor allele was the one counted as 2. Such loci are now kept
  on every platform. Results change mostly on Apple Silicon (on
  `platypus.gl[1:12, 1:200]` the T27 diagonal goes from 0.943 to 0.986),
  and on Linux and Windows only for single-copy loci whose minor allele is
  counted as 2 (2 of the 13 single-copy loci in `testset2.gl`).
  Same fix as `dartR.spatial::gl.grm2()`.

* `gl.report.kin.confidence()`: the standard errors now describe the
  estimator `gl.kin()` returns. For SNP data each resample is `G/2`; it
  still subtracted `mean(diag(G) - 1)`, the centring #108 removed from
  `gl.kin()`, which made standard errors about 17% too large (median
  0.0258 against 0.0220 on `testset2.gl`). For SilicoDArT data each
  resample uses the `dominant` estimator of `gl.kin()`; it used the
  correlation of band profiles / 2, a different estimator (correlation
  0.75 with `gl.kin()`). New warnings: individuals with call rate below
  0.8, whose standard errors mean imputation shrinks (on `testset2.gl`,
  all 24 captive-bred individuals), and a `kin` from a method other than
  the one bootstrapped. `kin` is now computed only when `verbose >= 3`
  uses it.

* `gl.report.kin.classes()`: a missing (`NA`) kinship no longer turns
  every pair's class into `NA`; the median baseline ignores `NA`, only
  the affected pairs are classed `NA` (returned with `all.pairs = TRUE`),
  and a warning gives the count. `conflicts` is now a zero-row data frame
  rather than `NULL` when there are no conflicts or no sire/dam columns.
  At least three individuals are required (with two, the baseline is the
  pair's own kinship, so it was always "unrelated"), and fewer than ten
  give a warning. With more than one population in `x` the function now
  warns that population structure inflates within-population classes and
  deflates between-population ones. It also warns when any individual's
  call rate is below 0.8: missing genotypes pull kinship toward 0, which
  is the main reason 19 of 48 recorded parent-offspring links in
  `testset2.gl` are called second-degree (5 after filtering loci at call
  rate 0.95). An earlier version of this entry attributed those 19 to
  population structure.
* Kinship series: mean kinship (MK) and the gene diversity of a whole
  group need kinship estimated on a wider reference than the group, because
  `gl.kin()` centres kinship on the individuals it is given (each row
  averages 0). `gl.report.ind.remove()`, `gl.report.repro.targets()`,
  `gl.select.pairs()`, `gl.report.mate.suitability()` and
  `gl.report.gd.projection()` (when given `x`) now stop with an explanation
  when `kin` is `NULL` or was estimated on `x` alone; previously their MK
  was ~1e-17 for everyone and they ranked, allocated and paired on rounding
  noise (or, for removals, on inbreeding alone). Every series function now
  accepts a `kin` covering more individuals than `x`, e.g.
  `gl.select.pairs(captive, kin = gl.kin(full.dataset))`. `gl.kin()`
  output records the individuals it was estimated on (`attr(, "ref.ids")`).
  `gl.report.kin.sets()` reports `meanMK` as the within-set mean kinship
  (it averaged the ~0 dataset-wide values). `gl.report.kinship()` warns
  when a single population is analysed with self-referenced kinship.

* `gl.report.kinship()`: `MK` and `MKrank` are now computed within each
  individual's population (PMx's managed population), ranked within
  population and sex; they were taken over the whole dataset, so a
  captive animal's rank depended on its kinship to wild individuals, and
  with `gl.kin()` method `"grm"` every dataset-wide MK was 0 up to
  rounding (each row of the centred matrix sums to 0), so the rank sorted
  rounding noise. In the `overall`
  row `GD` and `FGE` are NA: the full-matrix mean kinship is about 0 by
  construction, which gave `GD = 1` and `FGE` around 1e17. NA, empty and
  "Unknown" sexes form one rank group, and missing kinship values are
  ignored in the means with a warning.

* Relationship matrices now carry their scale: `gl.grm()` output is tagged
  `attr(, "scale") = "relatedness"` (`gl.run.EMIBD9()` and `gl.kin()`
  already tag `"kinship"`). The kinship series (`gl.report.kin*`,
  `gl.report.ind.*`, `gl.select.pairs`, ...) halves a `kin` tagged
  `"relatedness"`, so passing `gl.grm()` output gives the same results as
  `gl.kin()`; previously it was used as kinship, twice too large. A
  matrix with any other tag is an error; untagged matrices are still
  taken as kinship. `gl.grm()` documentation corrected: the diagonal can
  be below 1, and the off-diagonal mean is `-mean(diag)/(n - 1)` with `n`
  individuals, not loci.

* `gl.grm.network(standardise = TRUE)` now follows Goudet et al. (2018):
  kinship relative to the average pair, `(theta - mean) / (1 - mean)`.
  It subtracted the mean inbreeding coefficient from every kinship, so
  after call-rate filtering parent-offspring pairs came out near 0.1
  instead of 0.25. Standardised kinship values, categories and the links
  drawn above `kinship.threshold` change.
* `gl.kin()`: `method = "grm"` (the SNP default) now returns `G / 2`.
  It returned `G / 2 - mean(diag(G) - 1)` off the diagonal, which shifted
  every pairwise kinship by the mean inbreeding of the sample; after
  call-rate filtering, parent-offspring pairs came out near 0.1 instead of
  0.25. `method = "emibd9"` no longer halves EMIBD9 kinship (it returned
  half the kinship and a diagonal of 0.25) and has a new `emibd9.path`
  argument. All kinship-series results computed from `gl.kin()` change
  accordingly. The output carries `attr(, "scale") = "kinship"`, a warning
  reports pairs above 0.5 (possible duplicates), and the verbose summary
  prints mean pairwise kinship instead of a gene diversity that is about 1
  by construction. dartR.data (>= 1.2.5) is now required, for
  `testset2.gl` in the examples.
* `gl.diagnostics.relatedness()`: accuracy statistics are now computed
  against the exact pedigree kinship of each pair (new column `rel`;
  inbreeding included), with one relationship class per pair
  (`RelDegree`, including a new "unrelated" class). Previously pairs
  below 0.05 were discarded before any summary, founders in an attached
  pedigree were labelled half siblings, some pairs were counted under
  two relationships, and "RMSE" was the mean absolute error. RMSE and
  variance values change. `run_sim = TRUE` now works without supplying
  variable files (dartR.sim's defaults are used). SilicoDArT input now
  errors.

* `gl.plot.network()`: new argument `type = c("similarity", "distance")`;
  with `type = "distance"` the closest pairs (lowest values) are drawn.
  Previously the largest values were always drawn, i.e. the most distant
  pairs of a distance matrix; the default keeps the previous behaviour
  for similarity matrices such as `gl.grm()`. Link widths now reflect
  each drawn link's own value (they were taken from other pairs).
  `x = NULL` works as documented. D is checked against `indNames(x)`; an
  invalid `method` errors. Returns the drawn igraph network invisibly
  (was NULL). Much faster for many individuals.
* `gl.sim.relatedness()`: `rel = "full.sib"` returned the average
  parent-offspring kinship; it now simulates two full siblings. The
  reported interval was the confidence interval of the mean; it is now
  the range holding the central `conf` proportion of simulated kinship
  values, with the CI of the mean returned separately. The function now
  returns a list (`values`, `mean`, `interval`, `ci.mean`, `plot`)
  instead of printing, no longer requires the `related` package, and
  runs EMIBD9 quietly without per-replicate heatmaps.

* `gl.filter.parent.offspring()`: `method = "best"` removed the member of
  each pair with fewer missing genotypes (counts compared as text); it
  now keeps the better-genotyped individual, as documented. Pairs now
  come from `gl.report.parent.offspring()`, so the filter uses the same
  proportion-based pairs, `range` is honoured by both methods (and
  `range = 3` no longer errors), and SilicoDArT input errors.
  Individuals in several pairs are removed first, so fewer individuals
  are removed. One history entry per call; an unknown `method` errors.

* `gl.report.parent.offspring()`: pairs are now assessed on the
  proportion of pedigree-inconsistent loci among loci genotyped in both
  individuals (new columns `n.loci` and `prop`), so individuals with much
  missing data are no longer flagged as relatives of many others.
  p-values were computed on the wrong scale and are now
  `pnorm(zscore)` (larger than before). `range` below 1.5 now flags more
  pairs, as documented. SilicoDArT input now errors. Faster on large
  datasets; nothing printed at `verbose = 0`.
* `gl.run.EMIBD9()`: a run in which EMIBD9 fails now stops with EMIBD9's
  message; previously the results of the previous call in the same
  session could be returned. `OutAlleleFre = TRUE` no longer crashes
  EMIBD9. `raw` and `processed` columns are numeric (were character).
  SilicoDArT input now errors. `rel` is documented as the kinship
  coefficient (it always was) and carries `attr(rel, "scale") =
  "kinship"`. Building the input file is faster for large numbers of
  loci. EMIBD9's console output prints only at `verbose >= 2`.

* `gl.grm.network()`: a `G` tagged `attr(G, "scale") = "kinship"` (as
  returned by `gl.run.EMIBD9()`) is used as kinship without halving.
  Since the previous change EMIBD9 kinship was halved, because EMIBD9's
  `r(1,2)` is kinship, not relatedness.
* `gl.run.colony()`: COLONY's output files are now written to `outpath`
  instead of the R working directory, and the function returns
  `list(files, best.config)`, where `best.config` is the COLONY best
  configuration as a data frame (previously it returned `outpath`). The
  default `outpath = NULL` now works, paths with spaces are supported,
  and a run in which COLONY rejects its input now stops with COLONY's
  error message instead of returning normally. COLONY's console output
  prints only at `verbose >= 2`. Individual names longer than COLONY's
  20-character limit or containing spaces are replaced by short IDs for
  the run and restored in `best.config`.

* `gl.grm.network()`: kinship is now computed as G / 2 when
  `standardise = FALSE` (default); previously relatedness values were
  treated as kinship, so thresholds and categories were off by a factor
  of two. `categorise` colours now follow the documented order. `G` must
  match `indNames(x)`; an invalid `method` now errors. The return value
  is a named list (`plot`, `kinship`).

* `gl2colony()`: stops with an error, instead of writing a file COLONY
  rejects or misreads, when individual names contain whitespace, the data
  are SilicoDArT, or `sibship.prior`, `known.allele.freq`, or any
  known/excluded count is non-zero. A missing offspring/mother/father
  column no longer resets the other two; column names and values are
  matched ignoring case. A single `allelic.dropout`/`other.typ.err` value
  without `@` now applies to all loci. Progress messages follow
  `verbose`.

* `gl.assign.pca()` moved out of dartR.captive: the successor in
  dartR.popgen is the newer revision (April 2026 bug fixes and CL
  ellipse documentation caveats); the copy here had diverged and its
  only post-divergence change was a commented-out example line.

* `gl.assign.pa()` moved out of dartR.captive: the reviewed successor
  lives in dartR.popgen (arriving with the assignment-suite
  migration).

* `gl.assign.on.genotype()` moved out of dartR.captive: the reviewed
  successor lives in dartR.popgen (arriving with the assignment-suite
  migration).
* `utils.assignment()`, `utils.assignment_2()`, `utils.assignment_3()`
  and `utils.assignment_4()` REMOVED: four side-by-side revisions of
  the same helper ("Population assignment probabilities", identical
  signatures), all exported yet called by nothing anywhere in the
  family - development history versioned by filename rather than git.

* `gl.assign.mahal()` and `gl.assign.mahalanobis()` REMOVED (family
  consolidation): two side-by-side versions of the same function
  (identical signatures and titles). The reviewed and corrected
  `gl.assign.mahalanobis()` becomes the single implementation across
  the verse, arriving with the assignment-suite migration to
  dartR.popgen; dartR.base's `gl.mahal.assign()` is removed in the
  companion PR. No callers existed in the family.

* `gl.grm()`: now errors on SilicoDArT (presence/absence) input instead of
  silently returning a numerically meaningless matrix — the additive
  relationship algorithm and the documented diagonal range (1 to 2) only
  hold for SNP dosage data. Callers passing SilicoDArT data will need to
  filter to SNP data before calling `gl.grm()`.
* `gl.grm()`: fixed a crash (`object 'p3' not found`) when
  `plotheatmap = FALSE` was combined with a non-`NULL` `plot.file`; that
  combination now computes and returns the matrix and warns that nothing
  was saved, since no plot is generated.
* `gl.grm()`: corrected the documented default for `legendy` (was `1`,
  actual default is `0.5`).
* `gl.grm.network()`: fixed a crash (`` `breaks` and `labels` have
  different lengths``) whenever exactly one pair of individuals cleared
  `kinship.threshold` — previously the common case of confirming a single
  suspected relationship never returned a plot.
* `gl.grm.network()`: fixed a crash (`Insufficient values in manual
  scale`) when `categorise = TRUE` was combined with a
  `kinship.threshold` below `0.1`. The undocumented "First Cousins"
  bucket this produced has been removed; `categorise = TRUE` now always
  shows the 3 documented kinship categories, and pairs below `0.1` are
  left uncategorised.
* `gl.grm.network()`: an individual related (above `kinship.threshold`)
  to 2 or more other individuals was plotted as a duplicate, overlapping
  node; each individual is now plotted once.
* `gl.grm.network()`: corrected the documented defaults for
  `node.label.size` (was `3`, actual default is `2`) and `title` (was
  `'Network of similarity matrix'`, actual default is `'Network of a
  similarity matrix'`).
* `gl.assign.grm()`: now errors clearly when `pop(x)` is `NULL`, instead of
  failing deep inside `order()` with an opaque `argument 1 is not a
  vector` error.
* `gl.assign.grm()`: now errors on genlight objects with duplicate
  individual names, instead of silently folding two different individuals
  into population `'unknown'` and returning a numerically wrong result.
* `gl.assign.grm()`: corrected `@description`/`@return` to describe the
  function's actual output (a named numeric vector of mean pairwise
  relatedness scores) — previously documented as returning assignment
  probabilities in a `data.frame`, which the code never did.
* `gl.assign.grm()`: now prints a results summary (best-matching
  population) at `verbose >= 3` and a completion message at
  `verbose >= 1`, matching the documented verbosity levels.
