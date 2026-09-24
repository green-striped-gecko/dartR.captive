# dartR.captive NEWS

## Unreleased

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
