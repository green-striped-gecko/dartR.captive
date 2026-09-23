# Review: gl.filter.parent.offspring (dartR.captive)
- Family mode: modify
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 93275ab (branch `review-gl.report.parent.offspring`, PR #102, based on origin/dev 4de436e)
- Datasets: platypus.gl (81 x 1,000), testset.gl, testset.gs
- Baseline: tests/testthat/test-gl.filter.parent.offspring.R (snapshot captured pre-review, 8 expectations)

## Verdict

**Standards: Needs work** — metadata stay in sync (removal goes through
`gl.drop.ind`) and history is appended, but the function carries a copy
of the pair-detection code from `gl.report.parent.offspring`, prints at
`verbose = 0`, and adds two history entries per call.
**Spec: Rework** — `method = "best"` keeps the individual with the
lower call rate in 9 of 10 pairs on `platypus.gl`, the opposite of its
documentation; `range` makes "best" crash and is ignored by "random".

## Findings

**F1 [HIGH, confidence: high] — "best" keeps the worse individual
(DOC5; inverted-filter class)**
`R/gl.filter.parent.offspring.r:332-339` — missing-data counts are put
into a character column (`cbind`), sorted as text, and the first
(fewest missing, when the text order agrees) is removed.
Failure scenario: `platypus.gl`, default call: T5 (59 NAs) is removed and
T3 (68 NAs) kept; the same inversion in 9 of 10 pairs. The one exception,
T39 (140 NAs) vs T35 (64), is removed only because "140" sorts before
"64" as text.
Proposed change: compare numeric missing-data counts and remove the
member with more missing data, as documented.
**Consequence: a different individual is removed from most pairs; the
retained individuals are the better-genotyped ones.**

**F2 [HIGH, confidence: high] — pair detection copied from the report
function, with its defects (FS1, DOC5)**
`R/gl.filter.parent.offspring.r:129-291` — a copy of the code reviewed
in `gl.report.parent.offspring` (PR #102): raw counts that let missing
data create false pairs, the boxplot that ignores `range` below 1.5, the
per-pair loop, and silent SilicoDArT. Specific to the filter:
- `range = 3` with "best" crashes ("Subsetting resulted in zero
  individuals"): no pair passes the cutoff, and `for (i in 1:nrow(...))`
  runs over `1:0`.
- "random" does not apply the cutoff at all: with `range = 3` it still
  removes 10 individuals.
Proposed change: call `gl.report.parent.offspring()` to find the pairs
(same arguments, plots passed through) and keep only the selection step
here. Requires PR #102.
**Consequence: the filter removes individuals from the same pairs the
report lists (proportion-based; identical to before on `platypus.gl` at
the default `range`); "random" respects `range`; SilicoDArT errors.**

**F3 [MEDIUM, confidence: high] — both members of a pair can be removed
(DOC5)**
`R/gl.filter.parent.offspring.r:301-343` — each pair is resolved on its
own, so an individual in several pairs can lead to removing more
individuals than needed, and both members of one pair can go.
Failure scenario: `platypus.gl` pairs SUS19-SUS34 and SUS19-SUS22: SUS19
is removed for the first pair and SUS22 for the second, so the
SUS19-SUS22 pair loses both members; 10 individuals are removed for 10
pairs where 9 suffice.
Proposed change: resolve pairs in order of evidence (lowest `prop`
first) and skip a pair once either member has been removed.
**Consequence: fewer individuals removed when individuals appear in more
than one pair.**

**F4 [LOW, confidence: high] — two history entries per call (FS8)**
`R/gl.filter.parent.offspring.r:343,409` — the internal `gl.drop.ind()`
call appends its own entry before this function appends one (baseline:
history grows by 2).
Proposed change: restore the history after the internal call so only
this call is recorded.

**F5 [LOW, confidence: high] — messages, argument check, docs (VRB2,
VRB3, FS5, DOC1, DOC5, DOC7)**
The `plot_colors = gl.colors(2)` default prints three lines at
`verbose = 0`; messages carry source line breaks; any `method` other than
"best" (e.g. "Best") silently selects "random"; `@return` says "NULL if
no parent-offspring relationships were found" (the object is returned
unchanged); typos ("no used", "miss-called"); `@author` has no
`Author(s):`; `build = "Jody"`.
Proposed change: `plot_colors = NULL` as in PR #102;
`match.arg(method, c("best", "random"))`; one-line messages; doc fixes.

## Proposed changes

1. "best" removes the member with more missing data, compared as numbers
   (F1).
   **Consequence: a different individual is removed from most pairs.**
2. Find pairs with `gl.report.parent.offspring()` (F2).
   **Consequence: pair detection as in PR #102; "random" respects
   `range`; `range = 3` no longer crashes; SilicoDArT errors.**
3. Resolve pairs strongest first, skipping pairs already broken (F3).
   **Consequence: fewer individuals removed when an individual is in
   several pairs.**
4. One history entry per call (F4).
5. Messages, `match.arg(method)`, docs (F5).
   **Consequence: an unknown `method` value errors instead of meaning
   "random".**

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Modify-mode checks: metadata sync via `gl.drop.ind` (ind.metrics, pop,
  flags) — delegated, not re-verified here; history — F4
- Spec: removed vs kept individuals compared with per-individual missing
  counts on `platypus.gl`; `range` = 3 for both methods; SilicoDArT — run
- FBM path: the code has an FBM branch for missing counts (`hold@fbm`);
  not tested (no FBM fixture)
- Callers: dartr2shiny (named arguments, no `plot_colors`); no other
  `dartR.*` package calls it

## Addendum (found while applying)

**A1 [MEDIUM, confidence: high] — change 3 as proposed still lets a pair
lose both members**
Resolving the strongest pair first and skipping broken pairs removed
SUS34 (pair SUS19-SUS34) and then SUS19 (pair SUS19-SUS22, because SUS19
has 64 missing genotypes and SUS22 62), so pair SUS19-SUS34 lost both
members and 9 individuals were removed. Decision (Luis Mijangos): remove
individuals that appear in the most unresolved pairs first, then resolve
the rest in order of evidence by call rate ('best') or at random. On
`platypus.gl` this removes 8 and every pair keeps one member. In general
(e.g. a chain A-B, B-C, C-D) a pair can still lose both members; the
documentation does not claim otherwise.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis Mijangos | |
| 2 | approved | Luis Mijangos | stacked on PR #102 |
| 3 | approved | Luis Mijangos | amended by A1: most-pairs-first rule |
| 4 | approved | Luis Mijangos | |
| 5 | approved | Luis Mijangos | |

## Outcome

All five approved changes (change 3 as amended by A1) applied to
`R/gl.filter.parent.offspring.r` on 2026-09-23; `devtools::document()`
run. `tests/testthat/test-gl.filter.parent.offspring.R`: 15 expectations
pass; `test-gl.report.parent.offspring.R` still 14. The pre-change
baseline tests were replaced by tests of the approved behaviour.

- Change 1 (F1): in every pair whose members are in no other pair, the
  member with more missing genotypes is removed (T3 68 vs T5 59, T42 99
  vs T28 68, T36 76 vs T38 64, T13 80 vs T16 72, T39 140 vs T35 64, T22 78
  vs T20 77). Was: the reverse in 9 of 10 pairs.
- Change 2 (F2): pairs come from `gl.report.parent.offspring()`; the
  function body drops the copied counting, boxplot and plotting code.
  `range = 3` returns all 81 individuals for both methods (was: error for
  'best', 10 removed for 'random'); SilicoDArT errors "Only SNP data are
  supported" (from the report function).
- Change 3 + A1 (F3): `platypus.gl` default removes 8 (SUS36, SUS19 first,
  then one per remaining pair), every pair broken, none losing both
  members (was 10, with SUS19-SUS22 losing both).
- Change 4 (F4): history grows by 1 (was 2), also with
  `rm.monomorphs = TRUE`.
- Change 5 (F5): 0 lines printed at `verbose = 0`; `method = "Best"`
  errors; `@return`, `@details`, `@param method`, typos, `Author(s)`.
- Metadata: after filtering (with `rm.monomorphs = TRUE`), `ind.metrics`
  rows, `pop` length and `loc.metrics` rows match the object.
- The `@examples` block runs (`devtools::run_examples`).
- Depends on PR #102 (`gl.report.parent.offspring`); merge #102 first.

```json
{
  "function": "gl.filter.parent.offspring",
  "package": "dartR.captive",
  "family": "modify",
  "skill_version": "2.0.0",
  "commit": "93275ab",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "FS1", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "FS8", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 5},
    {"id": "A1", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3}
  ],
  "coverage_skipped": ["DAT6/FBM branch: no FBM fixture"],
  "status": "pr-open",
  "pr": null
}
```
