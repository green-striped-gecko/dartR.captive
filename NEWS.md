# dartR.captive NEWS

## Unreleased

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
