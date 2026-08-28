# Review: gl.grm.network (dartR.captive)
- Family mode: analysis
- Date: 2026-08-28
- Reviewer: Claude (claude-sonnet-5), dartr-function-review v1.0.0
- Package commit: 7a84d3f
- Datasets: testset.gl (constructed small `G` matrices), possums.gl (documented example)
- Baseline: tests/testthat/test-gl.grm.network.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — the canonical anatomy and verbosity/dependency
guards are in place, but three code paths crash or silently duplicate plot
data under ordinary parameter choices, and two documented defaults don't
match the code.
**Spec: Rework** — `categorise` and the default kinship-gradient legend both
fail on inputs the function's own documentation tells users to try (a low
`kinship.threshold` to catch distant relatives; a dataset with one detected
pair), not on edge cases outside its stated purpose.

## Findings

**F1 [HIGH, confidence: high] — the plot crashes whenever exactly one
kinship value clears the threshold**
`R/gl.grm.network.r:346-348,401-407` — `breaks <- pretty(round(seq(min(edges$size),
max(edges$size), 0.1), 1), 4)` is fed into
`scale_colour_gradientn(breaks = breaks, labels = as.character(breaks))`.
When `edges` has a single row (`min == max`), `pretty()` returns break
points outside that single-point range; ggplot drops the out-of-range
break internally when rendering the legend but does not drop the
correspondingly-indexed `labels` entry, and `scale_colour_gradientn()`
errors.
Failure scenario:
```r
G <- diag(1, 3); dimnames(G) <- list(nm, nm)
G[nm[1], nm[2]] <- G[nm[2], nm[1]] <- 0.26
gl.grm.network(G, sub, verbose = 0)
```
errors `` `breaks` and `labels` have different lengths`` inside the
function's own `print(p1)` call — i.e. the call never returns. A dataset
with exactly one relative pair above `kinship.threshold` is not a corner
case; it's the headline use of this function (confirm a single suspected
relationship).
Proposed change: drop the manual `breaks`/`labels` computation and let
`scale_colour_gradientn()` derive its own breaks (verified to render
correctly for both the single-edge case and the documented
`possums.gl` example).

**F2 [HIGH, confidence: high] — `categorise = TRUE` crashes once a 4th
kinship bucket appears (DOC5)**
`R/gl.grm.network.r:360-367` bins `edges$size` into four labelled buckets —
"Same Individual", "Full Siblings/Parent-Offspring", "Half Siblings", and
an undocumented "First Cousins" (0.038-0.1) — then
`scale_color_manual(values = color.categories)` (line 378-381) is handed
only the 3 colors the `@param color.categories` doc and `@details` table
describe. Whenever a call produces all 4 buckets (any `kinship.threshold`
below 0.1 — exactly what the function's own kinship-value reference table
recommends for catching cousins), `ggplot_build()` errors `Insufficient
values in manual scale. 4 needed but only 3 provided.`
Failure scenario: five individuals with pairwise kinships spanning
0.06-0.35 and `categorise = TRUE, kinship.threshold = 0.05` — a
straightforward "show me everything down to first-cousin" call — crashes
on `print(p1)` inside the function.
Proposed change: make the bucket count and `color.categories` length
agree — either drop the undocumented "First Cousins" bucket so
categorisation always yields the 3 documented categories, or extend the
`color.categories` default to 4 colors and document the new bucket in
`@param color.categories` and the `@details` table. Needs the custodian's
choice of which side is the intended behavior.

**F3 [MEDIUM, confidence: high] — an individual in 2+ above-threshold
pairs is plotted as a duplicate node**
`R/gl.grm.network.r:243-249,305` — `links_plot` is built by splitting each
qualifying pair into a `from` row and a `to` row, with no deduplication by
individual. `merge(plotcord, links_plot, by = "label.node", all.x = TRUE)`
then produces one `plotcord` row per matching `links_plot` row, so any
individual related (above `kinship.threshold`) to 2+ others gets
duplicated in `plotcord`.
Failure scenario: individual 1 related to both individual 2 and
individual 3 above threshold — `nInd(sub)` is 4, but the rendered point
layer has 5 rows; individual 1's node and label are drawn twice,
overlapping at the same layout coordinates. Any hub individual (a parent
with multiple offspring, a member of a sibling group) triggers this.
Proposed change: deduplicate `links_plot` to one row per `label.node`
before the merge, e.g. `aggregate(kinship ~ label.node, links_plot, max)`
(verified to keep the correct — highest — kinship for the node's alpha
shading while producing exactly one row per individual).

**F4 [LOW, confidence: high] — fragile character/numeric coercion in the
`standardise` branches (STY1)**
`R/gl.grm.network.r:218-233` — the self-loop rows appended at
`rbind(links_tmp, cbind(from = indNames(x), to = indNames(x), kinship = 0))`
coerce the whole `kinship` column to character (`cbind()` on character
vectors forces the numeric side to character too). Every later consumer
(`links_matrix` via `apply(..., as.numeric)`, and `plotcord$kinship` via
`as.numeric(plotcord$kinship)`) re-parses the string back to a number, so
this currently self-heals and produces correct output for realistic
kinship magnitudes — confirmed by brute-force comparison of
`as.character(x) >= as.character(threshold)` against the numeric
comparison across 400k sampled values in `[-1, 1]`: 0 false negatives (a
value that should pass the threshold silently failing it), only harmless
false positives that are corrected by the later reparse. It's still a
maintenance hazard: it depends on kinship values never being large enough
for R to switch to scientific notation, and the `isFALSE(standardise)`
branch (lines 227-233) sets `links_tmp$kinship <- as.character(...)` then
immediately overwrites `links_tmp` two lines later, discarding that
conversion — dead code that misleads a future reader into thinking the
branches are typed differently than they are.
Proposed change: build the self-loop rows with
`data.frame(from = indNames(x), to = indNames(x), kinship = 0)` instead of
`cbind()`, which keeps `kinship` numeric through the `rbind()`; delete the
dead `as.character()` line in the `isFALSE(standardise)` branch. No
behavior change intended — the round-trip already produces the same
numbers; this only removes the fragility and the dead code.

**F5 [LOW, confidence: high] — two documented defaults don't match the
code (DOC1)**
`R/gl.grm.network.r:19,142` — `@param node.label.size` documents
`[default 3]`; the signature default is `node.label.size = 2`.
`R/gl.grm.network.r:27,147` — `@param title` documents `[default 'Network
of similarity matrix']`; the signature default is `"Network of a
similarity matrix"` (note "a").
Failure scenario: a user reading `?gl.grm.network` sets no `title` and
expects the documented string; the plot shows different text than
documented. Harmless to results, but a direct doc/code contradiction.
Proposed change: reconcile each pair (change the roxygen default or the
code default, whichever the custodian intends), then run
`devtools::document()`.

**F6 [LOW, confidence: high] — inconsistent style within the same
function (STY1)**
`R/gl.grm.network.r:423` uses `if (node.label == T)` — `T` is a rebindable
base-R variable, not a keyword; `isTRUE(node.label)` is the safe idiom
already used elsewhere in this file (`isTRUE(standardise)`, line 218).
`R/gl.grm.network.r:250,340` access `x$ind.names` and `x$pop` directly,
while the rest of the function uses the `indNames(x)` and `pop(x)`
accessors (lines 175-186, 235-236, 300).
Proposed change: replace `== T` with `isTRUE()`; use `indNames(x)`/`pop(x)`
in the two direct-slot-access spots for consistency with the rest of the
file.

**F7 [INFO, confidence: high] — `@author` has no separate `Author(s):`
line (proposed rule DOC7)**
`R/gl.grm.network.r:103-104` — `@author Custodian: Arthur Georges -- Post
to ...` names only a custodian. The same gap was raised and approved as
F7 in the `gl.grm` review (PR #85, same `@family`, same custodian).
Proposed change: prefix with `Author(s): Arthur Georges. ` before the
existing `Custodian:` text, matching the fix already applied to `gl.grm`.

**F8 [INFO, confidence: medium] — duplicated ~50-line plot-building block
(STY3)**
`R/gl.grm.network.r:369-421` — the `categorise` and non-`categorise`
branches build near-identical `ggplot()` chains (segment geom, point geom,
theme, title), differing only in the color scale/mapping
(`scale_color_manual` + `cat` vs. `scale_colour_gradientn` + `size`).
Proposed change: build the shared base plot (segments without color,
points, theme, title) once, then add the `categorise`-dependent
color-scale layer conditionally. No behavior change intended; verify via
the F1-fixed snapshot that plot output is unchanged for both branches.

## Proposed changes

1. Remove the manual `breaks`/`labels` computation feeding
   `scale_colour_gradientn()`; let it derive breaks automatically (F1).
2. Reconcile the `categorise` bucket count with `color.categories` length
   — drop the undocumented "First Cousins" bucket, or extend
   `color.categories` to 4 colors and document the 4th bucket (F2).
   **Needs the custodian's choice of which side is correct** — this
   changes what the categorised plot's legend shows for
   `kinship.threshold < 0.1`, once it stops crashing.
3. Deduplicate `links_plot` to one row per individual before the
   `plotcord` merge (F3).
4. Build the self-loop rows as a typed `data.frame()` instead of
   `cbind()`, and delete the dead `as.character()` line in the
   `isFALSE(standardise)` branch (F4). No behavior change intended.
5. Reconcile the `node.label.size` and `title` documented defaults with
   the code, then run `devtools::document()` (F5).
6. Replace `node.label == T` with `isTRUE(node.label)`; use
   `indNames(x)`/`pop(x)` instead of `x$ind.names`/`x$pop` (F6).
7. Add the `Author(s):` label to `@author`, matching `gl.grm` (F7,
   proposed rule).
8. Refactor the two near-duplicate `ggplot()` blocks into one shared base
   plot plus a conditional color-scale layer (F8).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behavior vs roxygen/`@details`, tested empirically via
  `devtools::load_all()` with constructed `G` matrices on `testset.gl`
  subsets — run, 3 defects found (F1, F2, F3)
- Analysis-family numerical check: the returned `links_matrix` was
  verified unaffected by F4 (round-trips correctly through
  `as.numeric()`); the char/numeric coercion was further checked by
  brute-force comparison against numeric comparison across 400k sampled
  values — no false negatives found, confirming F4 is a robustness issue
  rather than a live output-correctness bug
- DEP1 (`igraph` guard returns `-1` instead of `stop()`): matches the
  codebase-wide idiom (`gl.grm.r`, `gl.assign.pca.r`, `gl.plot.network.r`,
  `gl.sim.relatedness.R`, and others); not flagged as a per-function
  finding, consistent with the sibling `gl.grm` review (PR #85), which
  left the same pattern in place
  (`function-review/reports/dartR.captive/gl.grm.md`)
- `@examples` block run end-to-end on `possums.gl` — succeeds; the
  documented example doesn't hit F1/F2/F3 because it has multiple
  qualifying pairs spanning a wide value range and `categorise = FALSE`
- FBM path (DAT6): not applicable — the function operates on a
  precomputed similarity matrix `G`, not the genlight's genotype matrix

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Bernd Gruber | |
| 2 | approved | Bernd Gruber | option 2a: drop the undocumented "First Cousins" bucket (recommended option; user said "apply all" without objecting to the recommendation) |
| 3 | approved | Bernd Gruber | |
| 4 | approved | Bernd Gruber | |
| 5 | approved | Bernd Gruber | reconcile docs to match code's current defaults (no code default change) |
| 6 | approved | Bernd Gruber | |
| 7 | approved | Bernd Gruber | |
| 8 | approved | Bernd Gruber | |

## Outcome

All eight approved changes applied to `R/gl.grm.network.r` on 2026-08-28.

- Change 1 (F1): removed the manual `breaks`/`labels` computation feeding
  `scale_colour_gradientn()`; the scale now derives its own breaks.
  Verified: a 3-individual case with exactly one pair above threshold now
  returns a plot and `ggplot2::ggplot_build()` on it runs silently (was:
  `` `breaks` and `labels` have different lengths``).
- Change 2 (F2, option 2a): removed the undocumented "First Cousins"
  (0.038-0.1) bucket, so `categorise = TRUE` always yields at most the 3
  documented buckets, matching the 3-color `color.categories` default.
  Verified: 5 individuals with pairwise kinships spanning 0.06-0.35,
  `categorise = TRUE, kinship.threshold = 0.05`, now builds without error
  (was: `Insufficient values in manual scale. 4 needed but only 3
  provided.`); the formerly-crashing 0.06 pair is now uncategorised (NA)
  rather than crashing.
- Change 3 (F3): `links_plot` is now deduplicated to one row per
  individual (`aggregate(kinship ~ label.node, ., max)`) before the
  `plotcord` merge. Verified: an individual in 2 above-threshold pairs
  (4 individuals total) now produces exactly 4 point-layer rows in
  `ggplot_build()$data` (was 5).
- Change 4 (F4): the self-loop rows are now built with
  `data.frame(from =, to =, kinship = 0)` instead of `cbind()`, keeping
  `kinship` numeric through `rbind()`; the dead `as.character()` line in
  the `isFALSE(standardise)` branch is removed.
- Change 5 (F5): `@param node.label.size` and `@param title` roxygen
  defaults corrected to match the code (`2` and `'Network of a similarity
  matrix'`); `devtools::document()` run, `man/gl.grm.network.Rd` updated.
- Change 6 (F6): `node.label == T` replaced with `isTRUE(node.label)`;
  `x$ind.names`/`x$pop` replaced with `indNames(x)`/`pop(x)`.
- Change 7 (F7): `@author` now states `Author(s): Arthur Georges.
  Custodian: Arthur Georges -- ...`, matching the fix applied to `gl.grm`.
- Change 8 (F8): the `categorise` TRUE/FALSE branches now share one base
  `ggplot()` call, branching only on a `color_layer` (aes mapping + scale)
  built once per branch. Along the way, `aes(color = edges$cat)` /
  `aes(color = edges$size)` were written as `aes(color = cat)` /
  `aes(color = size)` (bare column names) to resolve a
  "Use of `edges$cat` is discouraged" ggplot2 lint warning surfaced while
  restructuring — not a separately-approved finding, but a direct
  consequence of applying F8 with no behavior change.

Snapshot verification: `tests/testthat/test-gl.grm.network.R` re-run
against the patched function — all 4 tests pass. The 3 tests that
characterized F1-F3's pre-fix crashes/duplication were updated to assert
the new, approved behavior; all 3 diffs map to approved changes 1-3. The
baseline structural test (dimensions, row names, self-loop diagonal) is
unchanged. `tests/testthat/test-gl.grm.R` (sibling function) re-run
unaffected (12/12 pass, this PR touches no shared code). The documented
`@examples` block re-run end to end at `verbose = 3` on `possums.gl` —
completes without error or warning.
PR: not yet opened.

```json
{
  "function": "gl.grm.network",
  "package": "dartR.captive",
  "family": "analysis",
  "skill_version": "1.0.0",
  "commit": "7a84d3f",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "principle:plot-legend-crash", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "principle:node-dedup", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "STY1", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "STY1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "INFO", "confidence": "high", "rule": "DOC7 (proposed)", "status": "approved", "change": 7},
    {"id": "F8", "severity": "INFO", "confidence": "medium", "rule": "STY3", "status": "approved", "change": 8}
  ],
  "coverage_skipped": [],
  "status": "pr-open",
  "pr": null
}
```
