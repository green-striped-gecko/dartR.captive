# Review: gl.report.ind.remove (dartR.captive)
- Family mode: report
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: fdb8967 (origin/dev, after #121; the reference-kinship guard from #113 is in place)
- Datasets: testset2.gl (EmmacCaptBred, 24 individuals, with kinship from the full 274; subsets of 60/120/200 individuals for timing), the same after `gl.filter.callrate(method = "loc", threshold = 0.95)`, testset2.gs; dartR.data 1.2.6
- Baseline: tests/testthat/test-gl.report.ind.remove.R (snapshot captured pre-review, 16 expectations, all pass)

## Verdict

**Standards: Needs work**: the structure and read-only behaviour conform,
and the reference-kinship guard from #113 works. `dGD`, `MK` and the
greedy `gd.after` all match an independent recomputation exactly. After
two removals, the greedy set equals the best pair found by exhaustive
search.
**Spec: Needs work**: the method is sound, but one `NA` kinship crashes
the function with an opaque R error. An invalid `n.best`, including `0`,
silently means "unlimited". Low call rates change which animals are
removed. The greedy search slows steeply with population size.

## Findings

**F1 [MEDIUM, confidence: high] — one NA kinship crashes the function (DAT integrity)**
`R/gl.report.ind.remove.r:124-155` — `utils.kin.dgd` without `na.rm`
makes every `gain` `NA`, and `if (gains[best] <= gd.current)` stops with
"missing value where TRUE/FALSE needed".
Failure scenario: `kin` has an `NA` for one poorly genotyped animal, and
the user gets an R-internal error that names no cause.
Proposed change: `na.rm = TRUE` in the `utils.kin.dgd` calls (the
argument from #120), with a `verbose >= 1` warning giving the number of
`NA` pairs. The same treatment as #119, #120 and #122.

**F2 [MEDIUM, confidence: high] — low call rates change the removal set (DOC5, proposed rule)**
`R/gl.report.ind.remove.r:101` — the captive-bred animals have call rates
of 0.70-0.80, and mean imputation in `gl.kin` pulls their kinships toward
0. Unfiltered, the greedy set removes 7 animals. After filtering loci at
call rate 0.95 it removes 11, and 6 of those are not in the unfiltered
set. The rank correlation of `dGD` is 0.92.
Failure scenario: a manager culls or transfers the unfiltered removal
set, which differs substantially from the set that better-called loci
support.
Proposed change: a `verbose >= 1` call-rate warning and a `@details`
sentence, as in #115, #116, #119, #120 and #122.

**F3 [LOW, confidence: high] — invalid n.best becomes "unlimited" (FS5; API1, proposed rule)**
`R/gl.report.ind.remove.r:104-108` — `n.best = 0`, a negative value or a
character string gives a warning (at `verbose >= 1` only) and then
removes without limit: `n.best = 0` returns 7 removals.
Failure scenario: a user who wants the ranking only sets `n.best = 0`
and receives a 7-animal removal set, silently at `verbose = 0`.
Proposed change: invalid `n.best` is a fatal error (fail fast, FS5).
Values above `nInd - 1` keep their clamp and warning.

**F4 [LOW, confidence: medium] — greedy search scales as n^4 (STY2)**
`R/gl.report.ind.remove.r:150-153` — each greedy step recomputes the
full mean of an `n x n` matrix for every remaining candidate, so the cost
grows with the fourth power of population size. It takes 0.03 s at 60
individuals, 0.33 s at 120 and 1.6 s at 200. Extrapolated, that is about
1 minute at 500 and about 25 minutes at 1000.
Failure scenario: a large managed population (1000 animals in a breeding
program) blocks an interactive session for tens of minutes.
Proposed change: keep running sums (the total of the remaining block and
each candidate's row sum over it), and update them after each removal.
This gives the same results in `O(n^2)` per step. Verified by the
existing exact-value tests.

**F5 [LOW, confidence: high] — convention gaps (FS3, VRB2, DOC)**
`:90-92` uses the outdated `build =`, and `:175-186` prints the
`verbose >= 3` summary with raw `cat()`. The first `@details` paragraph
repeats `@param kin` in the pre-#113 wording ("should be estimated on",
"inflates estimates").
Failure scenario: none for results.
Proposed change: drop `build =`, use `report()`, and replace the stale
paragraph with a one-line pointer to `@param kin`.

Checked, nothing found: the input comes back identical, self-referenced
`kin` is rejected (#113), `n.best` caps the greedy set, and the
SilicoDArT path runs.

## Proposed changes

1. `NA`-tolerant GD with a warning (F1). **Consequence: a `kin` with
   `NA`, which errored before, now returns results.**
2. Call-rate warning and `@details` sentence (F2). No numerical change.
3. Invalid `n.best` is an error (F3). **Consequence: `n.best = 0`,
   negative or non-numeric values stop with an error instead of removing
   without limit.**
4. Incremental greedy update (F4). No numerical change (same removal set
   and `gd.after`), faster.
5. Standards cleanup and the stale `@details` paragraph (F5). No output
   change.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: behaviour against roxygen on the captive colony with full-dataset kinship, self-referenced kinship, `NA` kinship, `n.best` edge cases, testset2.gs — run
- Independent numerical check: `dGD`, `MK` and every `gd.after` recomputed from block means (maximum difference 0); greedy step 2 against exhaustive search over all pairs (equal)
- Missing data: removal set compared unfiltered against call-rate filtered
- Scaling: timed at 60, 120 and 200 individuals; 500 and 1000 extrapolated, not run
- Callers (API3): none outside `@seealso`
- FBM path (DAT6): SKIPPED — no FBM fixture
- Known complaints: not checked — added in the 2026 kinship series, no release history
- PLT: not applicable — no plot

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence (results instead of a crash for NA-containing kin) approved |
| 2 | approved | Luis | |
| 3 | approved | Luis | consequence (invalid n.best errors) approved |
| 4 | approved | Luis | |
| 5 | approved | Luis | |

## Outcome

- Change 1: `na.rm = TRUE` in the ranking's `utils.kin.dgd` calls and in `MK`; the greedy search counts only non-missing kinships; warning with the number of missing values. With one `NA` pair the function returns results, and its greedy set equals a brute-force greedy on the same matrix.
- Change 2: call-rate warning (24 individuals) and `@details` paragraph.
- Change 3: invalid `n.best` stops with `Fatal Error: n.best must be NULL or a single number >= 1`; values above `nInd - 1` still clamp with a warning.
- Change 4: running-sum greedy. Identical removal set and `gd.after` on the colony (baseline expectations unchanged). Timing: 0.021 / 0.057 / 0.081 / 0.171 s at 60 / 120 / 200 / 270 individuals (was 0.032 / 0.327 / 1.603 s at 60 / 120 / 200).
- Change 5: `build =` dropped, `report()` in the summary, stale pre-#113 paragraph replaced by a pointer to `@param kin`.
- Snapshot: 2 baseline expectations changed (NA crash, `n.best = 0`), mapped to changes 1 and 3; everything else unchanged. 21 expectations pass; example runs.
- NEWS entry added. PR: #123

```json
{
  "function": "gl.report.ind.remove",
  "package": "dartR.captive",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "fdb8967",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "DAT", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "FS5,API1", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "medium", "rule": "STY2", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS3,VRB2", "status": "approved", "change": 5}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "forum/issues: no release history", "timing at n=500/1000 extrapolated"],
  "status": "pr-open",
  "pr": 123
}
```
