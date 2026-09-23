# Review: gl.grm.network (dartR.captive) — second pass
- Family mode: analysis
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 1b6a37e (dev_luis, includes PR #87 fixes from the first review)
- Datasets: testset.gl (constructed small `G` matrices), possums.gl (`gl.grm()` output for scale check)
- Baseline: tests/testthat/test-gl.grm.network.R (first-review tests plus four "baseline r2" tests captured pre-review)
- First review: `function-review/reports/dartR.captive/gl.grm.network.md` (PR #87)

## Verdict

**Standards: Needs work** — the first-review fixes hold, but one of them
(the node deduplication in PR #87) reintroduced the empty-network crash
that commit eea754d had fixed, and there is no input check that `G`
belongs to `x`.
**Spec: Rework** — with the default `standardise = FALSE`, the function
treats `gl.grm()` values (relationship, about twice the kinship) as
kinship, so the default threshold and the `categorise` labels are off by
a factor of two; the categorise colours also do not follow the documented
order.

## Findings

**F1 [HIGH, confidence: high] — crash when no pair clears the threshold
(regression)**
`R/gl.grm.network.r:250` — `aggregate(kinship ~ label.node, ...)` (added
in PR #87 to fix duplicate nodes) errors on a zero-row data frame. The
empty-edge handling from eea754d (2026-06-04, lines 347-353) is never
reached.
Failure scenario: `G <- diag(1, 5)` with names of `testset.gl[1:5, ]`, or
any real dataset without close relatives at the default threshold:
`Error in aggregate.data.frame(...): no rows to aggregate`. No test
covered the empty case, so the regression passed the first review.
Proposed change: run `aggregate()` only when `links_plot` has rows;
otherwise use an empty data frame with the same columns. Add a test for
the empty case.

**F2 [HIGH, confidence: high] — `categorise` colours do not follow the
documented order (DOC5)**
`R/gl.grm.network.r:365` — `scale_color_manual(values = color.categories)`
gets an unnamed vector, so ggplot assigns colours to categories in
alphabetical order of the levels present, not in the documented order
(Same Individual, Full Siblings, Half Siblings).
Failure scenario: all three categories present — "Same Individual" is
drawn cyan (`#3ED2E6`, documented for Half Siblings) and "Half Siblings"
yellow. With only two categories present the colours shift again, so the
same relationship gets different colours in two plots.
Proposed change: name the vector:
`values = setNames(color.categories, c("Same Individual",
"Full Siblings\nParent-Offspring", "Half Siblings"))`.
**Consequence: the colours of categorised plots change for every user of
`categorise = TRUE`.**

**F3 [HIGH, confidence: medium] — default scale is relationship, not
kinship (DOC5, API1)**
`R/gl.grm.network.r:227-230` — with `standardise = FALSE` (default), the
off-diagonal of `G` is used directly as "kinship". `gl.grm()` returns
`rrBLUP::A.mat`, whose off-diagonal is the additive relationship, about
2 x kinship (the function's own baseline test comments `0.6 # kinship
0.3`). On `possums.gl`, `gl.grm()` off-diagonals reach 1.32 and diagonals
average 1.27.
Failure scenario: a half-sib pair (kinship 0.125, relationship 0.25) with
`categorise = TRUE` is labelled "Full Siblings/Parent-Offspring"; the
default `kinship.threshold = 0.125` shows links down to kinship 0.0625
(first cousins), not half-sibs as the `@param` and table suggest. The
legend title says "Kinship" throughout.
A related question on `standardise = TRUE` (lines 218-225): it computes
`G/2 - mean(diag(G) - 1)`, subtracting mean inbreeding from kinship. The
Goudet et al. (2018) estimator rescales against the mean kinship among
pairs, `(phi_ij - phi_bar) / (1 - phi_bar)`. The `@details` claim that the
method "aligns with Goudet et al." needs the custodian's check.
Proposed change (custodian's choice): (a) divide `G` by 2 when
`standardise = FALSE` so the plotted value is kinship, or (b) keep values
as they are, relabel the legend "Relationship", and document that
thresholds and categories must be doubled for `gl.grm()` output.
**Consequence of (a): links, node shading, categories, and the returned
matrix change for every call with the default `standardise = FALSE`.**

**F4 [MEDIUM, confidence: high] — no check that `G` matches `x` (FS5,
DAT5)**
`R/gl.grm.network.r:171-208` — nothing compares `rownames(G)` with
`indNames(x)`.
Failure scenarios: `G` with no dimnames or a name not in `x` errors
`length of 'dimnames' [1] not equal to array extent`; `G` with fewer
individuals than `x` runs without warning, adds the missing individuals
as isolated nodes, and returns a 5 x 5 matrix from a 4 x 4 `G`.
Proposed change: after the datatype check, `stop(error(...))` unless `G`
is square, has dimnames, and its names equal `indNames(x)` as a set;
reorder `x` or `G` to match if only the order differs.

**F5 [LOW, confidence: high] — return value undocumented (DOC5)**
`R/gl.grm.network.r:102,433` — `@return` says "A network plot"; the
function returns an unnamed `list(p1, links_matrix)`. `links_matrix` is
lower-triangular (upper triangle `NA`), its diagonal is set to 0
regardless of `G`, and rows follow the row order of `G`.
Failure scenario: a user calls `res[[2]]` expecting a symmetric kinship
matrix and gets `NA` for half of the pairs and 0 for all inbreeding
values.
Proposed change: document both elements and name them
(`list(plot = p1, kinship = links_matrix)`); positional access
(`res[[1]]`, `res[[2]]`) keeps working.

**F6 [LOW, confidence: high] — housekeeping (FS3, VRB3, STY1)**
`R/gl.grm.network.r:166` passes the outdated `build = "Jody"` to
`utils.flag.start()`. Lines 200-208 compare `method` with `||`, so
`method = c("fr", "kk")` errors `'length = 2' in coercion to
'logical(1)'`; an invalid method is silently reset to `"fr"` instead of
failing. Line 180 splits a message across source lines, printing a run of
spaces into the console.
Proposed change: drop `build =`; use `method <- match.arg(method,
c("fr", "kk", "gh", "mds"))` (an invalid value now errors instead of
falling back to "fr"); put the message on one line.

**F7 [INFO, confidence: high] — non-ASCII characters in roxygen
(proposed rule DOC6)**
`R/gl.grm.network.r:90-98` — en dashes and non-breaking hyphens in the
Speed & Balding table (`–`, `‑`).
Failure scenario: PDF manual build with pdflatex may fail on these glyphs.
Proposed change: replace with ASCII `-`, then `devtools::document()`.

## Proposed changes

1. Guard the `aggregate()` call against zero rows and add an empty-network
   test (F1).
2. Name `color.categories` so each category keeps its documented colour
   (F2). **Consequence: categorised plot colours change for all users.**
3. Resolve the relationship/kinship scale: (a) halve `G` in the default
   path, or (b) relabel as relationship and document; plus custodian check
   of the `standardise` formula against Goudet et al. (F3).
   **Consequence of (a): links, categories, and the returned matrix change
   for every default call.**
4. Validate that `G` is square, named, and matches `indNames(x)` (F4).
   **Consequence: calls with a `G` smaller than `x` now error instead of
   plotting extra isolated nodes.**
5. Document and name the returned list elements (F5).
6. Drop `build =`, use `match.arg()` for `method`, fix the split message
   (F6). **Consequence: an invalid `method` now errors instead of falling
   back to "fr".**
7. Replace non-ASCII characters in the roxygen table (F7, proposed rule).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behaviour vs roxygen with constructed `G` matrices on
  `testset.gl` subsets — run; F1, F2, F4, F5 reproduced
- Scale check: `gl.grm()` on `possums.gl` (callrate-filtered, 500 loci) —
  run; supports F3. The Goudet comparison is from the published formula,
  not re-derived numerically: confidence medium
- `NA` in `G`: run — pairs with `NA` are dropped, no error
- DEP1 (`igraph` guard returns `-1`): not flagged, codebase-wide idiom
  (same decision as first review)
- FBM path (DAT6): not applicable — `x` supplies names and populations only
- Callers: grepped `dartR.*` siblings and dartr2shiny — only dartr2shiny
  calls it (named arguments; `G` from `gl.grm()` or EMIBD9 `$rel`)

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis Mijangos | |
| 2 | approved | Luis Mijangos | |
| 3 | approved | Luis Mijangos | option (a): halve G in the default path. The `standardise` formula question (Goudet et al.) is deferred to the custodian, Arthur Georges |
| 4 | approved | Luis Mijangos | |
| 5 | approved | Luis Mijangos | |
| 6 | approved | Luis Mijangos | |
| 7 | approved | Luis Mijangos | |

## Outcome

All seven approved changes applied to `R/gl.grm.network.r` on 2026-09-23;
`devtools::document()` run. `tests/testthat/test-gl.grm.network.R`: 25
passing. The four pre-change "baseline r2" snapshots were replaced by
tests of the approved behaviour; the first-round tests pass unchanged.

- Change 1 (F1): `aggregate()` runs only when `links_plot` has rows. An
  identity `G` now returns a plot with no links (was: `no rows to
  aggregate`).
- Change 2 (F2): colours named by category. With all three categories,
  Same Individual `#E63E94`, Full Siblings `#E5D44C`, Half Siblings
  `#3ED2E6` (was cyan / pink / yellow); removing the full-sib pair leaves
  the other two colours unchanged.
- Change 3 (F3, option a): default path uses `kinship = G / 2`. `G` = 0.4
  returns kinship 0.2; `G` = 0.2 (kinship 0.1) no longer draws a link at
  the default threshold 0.125. Caller grep: dartr2shiny passes `gl.grm()`
  and `gl.run.EMIBD9()$rel` (EMIBD9 `r(1,2)`, relatedness scale), both 2 x
  kinship, so halving is correct for both; arguments are passed by name.
  No callers in other `dartR.*` packages. NEWS entry added.
- Change 4 (F4): `G` must be square with row/column names equal to
  `indNames(x)` as a set. A 4 x 4 `G` for a 5-individual `x`, or a `G`
  without dimnames, now errors with a clear message. A `G` in a different
  row order than `x` still works.
- Change 5 (F5): return is `list(plot =, kinship =)`, documented. While
  checking this, the matrix's row order proved to follow `G`, not
  alphabetical order; the `@return` text states that.
- Change 6 (F6): `build =` removed, `method` checked with `match.arg()`,
  population message on one line.
- Change 7 (F7): en dashes, non-breaking hyphens, and non-breaking spaces
  replaced with ASCII; `grep -P '[^\x00-\x7F]'` on the file returns 0 lines.
- End-to-end: the `@examples` call on `possums.gl` at `verbose = 3`, all
  four layouts, `standardise = TRUE, categorise = TRUE`, and a row-shuffled
  `G` (same values once symmetrised) all run without error.
- PR #97 (`review-gl.grm.network` -> `dev`), commit 8b1ecbe.

```json
{
  "function": "gl.grm.network",
  "package": "dartR.captive",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "1b6a37e",
  "review_round": 2,
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "FS6", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "medium", "rule": "API1", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS3", "status": "approved", "change": 6},
    {"id": "F7", "severity": "INFO", "confidence": "high", "rule": "DOC6", "status": "approved", "change": 7}
  ],
  "coverage_skipped": ["DAT6: not applicable"],
  "status": "pr-open",
  "pr": 97
}
```
