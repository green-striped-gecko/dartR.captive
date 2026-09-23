# dartR.captive NEWS

## Unreleased

* `gl.plot.network()`: new argument `type = c("similarity", "distance")`;
  with `type = "distance"` the closest pairs (lowest values) are drawn.
  Previously the largest values were always drawn, i.e. the most distant
  pairs of a distance matrix; the default keeps the previous behaviour
  for similarity matrices such as `gl.grm()`. Link widths now reflect
  each drawn link's own value (they were taken from other pairs).
  `x = NULL` works as documented. D is checked against `indNames(x)`; an
  invalid `method` errors. Returns the drawn igraph network invisibly
  (was NULL). Much faster for many individuals.

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
