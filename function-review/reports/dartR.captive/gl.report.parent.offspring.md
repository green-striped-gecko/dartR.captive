# Review: gl.report.parent.offspring (dartR.captive)
- Family mode: report
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 65b6c0f (origin/dev)
- Datasets: platypus.gl (81 x 1,000, has rdepth and RepAvg), testset.gl, testset.gs, a synthetic 250 x 10,000 genlight
- Baseline: tests/testthat/test-gl.report.parent.offspring.R (snapshot captured pre-review, 8 expectations)

## Verdict

**Standards: Needs work** — the input object is not modified and no
history is added, as a report function should; but the outlier set is
read out of a ggplot object, the plot always prints, and `gl.colors()`
prints three lines at `verbose = 0`.
**Spec: Rework** — the pair counts are correct (they match an independent
matrix computation), but the p-values use the wrong distribution, `range`
below 1.5 is ignored, and individuals with missing data are flagged as
parents or offspring of many others.

## Findings

**F1 [HIGH, confidence: high] — p-values are computed on the wrong
scale (DOC5)**
`R/gl.report.parent.offspring.r:240-249` — `zscore` is already
standardised, but `pnorm(q = zscore, mean = mean(count), sd = sd(count))`
treats it as a raw count.
Failure scenario: `platypus.gl`, pair T5-T3 (0 inconsistent loci, mean
18.5, sd 4.92): z = -3.76, correct one-sided p = `pnorm(-3.76)` =
8.4e-05; reported p = 2.98e-06, 28 times smaller. The error grows with
the mean count, so large datasets report vanishing p-values for any
outlier.
Proposed change: `p = pnorm(zscore)`; drop the rounding to 8 decimals
(which turns small p-values into 0).
**Consequence: every reported p-value changes (becomes larger).**

**F2 [HIGH, confidence: high] — `range` below 1.5 has no effect (DOC5,
PLT3)**
`R/gl.report.parent.offspring.r:201-212,299` — candidate pairs are the
outliers of a ggplot boxplot, which always uses 1.5 x IQR; `range` only
filters that set afterwards. The result therefore depends on a plot
object's internals.
Failure scenario: `platypus.gl` with `range = 0.5` or `1`: 10 pairs,
the same as `range = 1.5`, although 218 and 55 pairs lie below those
cutoffs.
Proposed change: compute the cutoff Q1 - range x IQR from the counts and
flag pairs strictly below it, without `ggplot_build()`. Default output
unchanged.
**Consequence: `range < 1.5` now flags more pairs, as documented.**

**F3 [HIGH, confidence: high] — missing data creates false
parent-offspring pairs (DOC5)**
`R/gl.report.parent.offspring.r:175-180` — the statistic is the raw count
of opposing homozygous loci, so a pair typed at fewer loci has fewer
chances to show an inconsistency.
Failure scenario: `platypus.gl` with 80% of T10's genotypes set to
missing: T10 appears in 44 of the 54 flagged pairs.
Proposed change (custodian's choice): (a) use the proportion of
pedigree-inconsistent loci among loci typed in both individuals, and add
columns `n.loci` and `prop` to the result; or (b) keep counts and warn
when individual call rates differ enough to bias them.
**Consequence of (a): flagged pairs, z-scores and p-values change for any
dataset with missing data; the result gains two columns.**

**F4 [MEDIUM, confidence: high] — pairwise counting is slow (STY2)**
`R/gl.report.parent.offspring.r:152-180` — one R function call per pair.
250 individuals x 10,000 loci takes 9.7 s; 1,000 individuals is 16 times
more pairs. The same counts come from two matrix products
(`H0 %*% t(H2) + H2 %*% t(H0)`, verified equal on `platypus.gl`).
Proposed change: matrix products. Same counts.

**F5 [MEDIUM, confidence: high] — SilicoDArT returns an empty result
(DAT1)**
`R/gl.report.parent.offspring.r:108,175-178` — presence/absence 0/1 never
produces the 0-vs-2 pattern, so every count is 0 and no pair is ever
reported, without a message.
Proposed change: stop for SilicoDArT with "only SNP data are supported".
**Consequence: SilicoDArT input errors instead of returning an empty
table.**

**F6 [LOW, confidence: high] — messages (VRB2, VRB3, FS3)**
The default `plot_colors = gl.colors(2)` prints "Starting gl.colors..."
at `verbose = 0`; progress messages contain line breaks and indentation
from the source; `verbose >= 3` prints all boxplot outliers, not the
returned pairs; `build = "Jody"`; `pop(x) <- x$ind.names` is unused.
Proposed change: `gl.colors(2, verbose = 0)` inside the function when
`plot_colors` is not given; one-line messages; print the returned table;
drop `build =` and the `pop` assignment.

**F7 [LOW, confidence: high] — documentation (DOC1, DOC5, DOC7)**
`@param plot.file` ends with "Creates a plot that shows the sex linked
markers" (copied from another function); `@return` does not say that
`Outlier` is the count of inconsistent loci or what `p` is; typos ("no
used", "miss-called"); `@author` has no `Author(s):` part.
Proposed change: fix those texts.

**F8 [INFO] — `gl.filter.parent.offspring` holds a copy of the same code**
`R/gl.filter.parent.offspring.r:162-222` repeats the pairwise count and
the boxplot-based outliers, so F2-F4 apply there too. Review next; it
could call this function instead of duplicating it. The argument names
`plot_theme`/`plot_colors` predate the `plot.theme`/`plot.colors`
convention (PLT1); renaming them is an API change and is not proposed.

## Proposed changes

1. `p = pnorm(zscore)`, no rounding (F1).
   **Consequence: every p-value changes.**
2. Outliers computed from the counts with the documented cutoff, strictly
   below Q1 - range x IQR (F2).
   **Consequence: `range < 1.5` flags more pairs; default unchanged.**
3. Missing data: (a) proportions with `n.loci` and `prop` columns, or (b)
   counts plus a warning (F3).
   **Consequence of (a): flagged pairs change whenever data are missing.**
4. Matrix-product counting (F4). Same counts.
5. Stop for SilicoDArT (F5).
   **Consequence: SilicoDArT input errors.**
6. Messages and small cleanups (F6).
7. Documentation (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Report-mode checks: input returned untouched (local copy only), no
  history appended — pass
- Spec: counts checked against an independent matrix computation on
  `platypus.gl` (T5-T3: 0 in both); p-value recomputed by hand; `range`
  swept 0.5-3; missing-data effect with 80% missing for T10 — run
- Performance: 250 x 10,000 synthetic — 9.7 s
- Issues: no GitHub issues found for "parent.offspring"; Google Group not
  searched (no access)
- Callers: dartr2shiny shows it as a table and plot
  (`slick_table_plot`); `gl.filter.parent.offspring` does not call it
  (it has its own copy, F8); no other `dartR.*` package calls it
- FBM path (DAT6): not tested; `as.matrix(x)` densifies once

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis Mijangos | |
| 2 | approved | Luis Mijangos | |
| 3 | approved | Luis Mijangos | option (a): proportions, with `n.loci` and `prop` columns |
| 4 | approved | Luis Mijangos | |
| 5 | approved | Luis Mijangos | |
| 6 | approved | Luis Mijangos | |
| 7 | approved | Luis Mijangos | |

## Outcome

All seven approved changes applied to `R/gl.report.parent.offspring.r` on
2026-09-23; `devtools::document()` run.
`tests/testthat/test-gl.report.parent.offspring.R`: 14 expectations pass;
the pre-change baseline tests were replaced by tests of the approved
behaviour.

- Default output on `platypus.gl`: the same 10 pairs as before the review
  (T5-T3, T42-T28 with 0 inconsistent loci, etc.), now ordered by `prop`
  and with `n.loci` and `prop` columns.
- Change 1 (F1): `p = pnorm(zscore)`; T5-T3 p = 1.06e-04 (was 2.98e-06).
  `zscore` is now computed on `prop` (-3.70; was -3.76 on counts).
- Change 2 (F2): cutoff computed from the proportions; `range` = 0.5, 1,
  1.5, 3 flags 218, 55, 10, 0 pairs (was 10, 10, 10, 0). The boxplot
  whiskers use the same `range`.
- Change 3 (F3, option a): with 80% of T10's genotypes missing, T10 is in
  5 of 15 flagged pairs (was 44 of 54). Remaining limitation: with only
  ~100 loci typed in both, a proportion is noisier than with ~480, so an
  individual with heavy missing data can still produce a few false pairs;
  `n.loci` in the output shows when this applies.
- Change 4 (F4): counts from `tcrossprod()`; equal to an independent
  matrix computation for every returned pair; 250 x 10,000 synthetic
  genlight 1.2 s (was 9.7 s).
- Change 5 (F5): `testset.gs` errors "Only SNP data are supported".
- Change 6 (F6): 0 lines printed at `verbose = 0` (`plot_colors` default
  is now `NULL`, resolved with `gl.colors(2, verbose = 0)`; dartr2shiny
  does not pass it); `verbose = 3` prints the returned table.
- Change 7 (F7): `plot.file` text, `@return` columns, `@details` on
  proportions and the cutoff, typos, `Author(s)`.
- The `@examples` block runs (`devtools::run_examples`).
- Not in this change: the NAMESPACE drift on `dev` (`dnorm`, `qnorm`).
- PR #102 (`review-gl.report.parent.offspring` -> `dev`).

```json
{
  "function": "gl.report.parent.offspring",
  "package": "dartR.captive",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "65b6c0f",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "PLT3", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "STY2", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DAT1", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 7},
    {"id": "F8", "severity": "INFO", "confidence": "high", "rule": "none", "status": "noted", "change": null}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "Google Group search: no access"],
  "status": "done",
  "pr": 102
}
```
