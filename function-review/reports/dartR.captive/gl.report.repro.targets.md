# Review: gl.report.repro.targets (dartR.captive)
- Family mode: report
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: a2976e3 (origin/dev, after #125)
- Datasets: testset2.gl, EmmacCaptBred (24 individuals: 11 males, 13 females) with `kin` from the full testset2.gl; the same after `gl.filter.callrate(method = "loc", threshold = 0.95)`; testset2.gs; dartR.data 1.2.6
- Baseline: tests/testthat/test-gl.report.repro.targets.R (snapshot captured pre-review, 15 expectations, all pass)

## Verdict

**Standards: Needs work**: the structure, the read-only behaviour and the
#113 reference guard are in place, and each sex's targets sum exactly to
`n.target`. Three inputs crash with R's own error message, and there are
small convention gaps.
**Spec: Needs work**: the allocation does what the name promises. On the
captive colony, the offspring generation's expected mean kinship is
0.0582, against 0.0667 for equal contributions and 0.0560 for the optimum
(optimal contributions, non-negative). But targets for 9 of 24
individuals change after call-rate filtering. The description also claims
that founder representation is "equalised", which the method does not
guarantee.

## Findings

**F1 [MEDIUM, confidence: high] — low call rates reshuffle the targets (DOC5, proposed rule)**
`R/gl.report.repro.targets.r:140` — `MK` is a row mean of kinship that
`gl.kin` shrinks toward 0 for individuals with missing genotypes. On the
colony (call rates 0.70-0.80), filtering loci at call rate 0.95 changes
the target of 9 of 24 individuals: CB_X_05 goes from 2 to 0, CB_X_03 and
CB_X_04 from 1 to 0, and CB_Y_02 and CB_Y_04 from 3 to 4. `MK`
unfiltered and filtered correlate at 0.93.
Failure scenario: a manager breeds CB_X_05 twice on unfiltered data; on
filtered data its target is 0.
Proposed change: the series call-rate warning at `verbose >= 1` (any
individual below 0.8), a `@details` sentence with these numbers, and an
example that filters loci on call rate before `gl.kin`.

**F2 [MEDIUM, confidence: high] — one NA kinship crashes the function (DAT integrity)**
`R/gl.report.repro.targets.r:140, 146-150` — `rowMeans(kin)` gives `NA`
mean kinship for both individuals of an `NA` pair; `max()` of the sex
pool is then `NA` and the function stops with
`missing value where TRUE/FALSE needed`.
Failure scenario: `kin` from the dominant estimator has an `NA` for two
individuals that share no scored loci; the report fails without naming
the cause.
Proposed change: `rowMeans(kin, na.rm = TRUE)` with a `verbose >= 1`
warning giving the number of `NA` values.

**F3 [MEDIUM, confidence: high] — an NA sex crashes the function (DAT integrity)**
`R/gl.report.repro.targets.r:104-105` — `ids[sex == "Male"]` returns an
`NA` id for an `NA` sex; the function warns that 1 individual is excluded,
then stops with `missing value where TRUE/FALSE needed`. The same defect
was fixed in `gl.report.mate.suitability` (#126).
Failure scenario: an empty sex cell in the metadata makes the report
fail with an error that does not mention sex.
Proposed change: `sex %in% "Male"` and `sex %in% "Female"`, so `NA` counts
as unknown sex, as documented.

**F4 [LOW, confidence: high] — an invalid `n.target` is replaced, not rejected (FS5)**
`R/gl.report.repro.targets.r:117-132` — `n.target = 0` or a negative
value is reset to `nInd(x)` with a warning that is silent at
`verbose = 0`; `n.target = NA` stops with R's
`missing value where TRUE/FALSE needed`.
Failure scenario: `n.target = 0` in a quiet script returns targets for 24
offspring.
Proposed change: `stop(error(...))` for a missing, non-numeric or
below-1 `n.target`; a fractional value is still rounded down with a
warning. **Consequence: calls with an invalid `n.target` stop instead of
running with `nInd(x)`.**

**F5 [LOW, confidence: high] — the description claims more than the method delivers (DOC5, proposed rule)**
`R/gl.report.repro.targets.r:5-11` — "founder-genome representation is
equalised in the next generation" is not something a weighting of mean
kinship can guarantee; the method lowers the expected mean kinship of the
next generation. It uses only each individual's mean kinship, not the
kinship between the parents it favours, so it stops short of the
optimum: 0.0582 against 0.0560 on the colony (it achieves 79% of the
possible reduction from equal contributions, 0.0667).
Failure scenario: a user reads the targets as founder-equalising and does
not compare them with an optimal-contribution tool.
Proposed change: reword `@description` to "lower the expected mean
kinship of the next generation" and add the colony comparison to
`@details`. Docs only.

**F6 [LOW, confidence: high] — convention gaps (FS3, DOC2)**
`:84-86` passes the outdated `build =`; the `verbose` text differs from
DOC2.
Failure scenario: none for results.
Proposed change: drop `build =`, use the DOC2 text.

Checked, nothing found: the input comes back identical; each sex's
targets sum exactly to `n.target` (24 and 12 tested); `kin = NULL` stops
with the #113 message; SilicoDArT runs; ties are broken deterministically.

## Proposed changes

1. Call-rate warning, `@details` sentence and filtered example (F1). No
   numerical change.
2. `na.rm = TRUE` in the mean kinship, with a warning (F2).
   **Consequence: a `kin` with `NA` now returns targets instead of an
   error.**
3. `NA` sex treated as unknown (F3). **Consequence: data with an `NA` sex
   now returns targets instead of an error.**
4. Invalid `n.target` is a fatal error (F4). **Consequence: calls with an
   invalid `n.target` stop instead of running with `nInd(x)`.**
5. Reworded `@description` and the optimum comparison in `@details` (F5).
   Docs only.
6. Standards cleanup (F6). No output change.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: behaviour against roxygen on EmmacCaptBred (SNP and SilicoDArT), `n.target` 24, 12, 0 and `NA`, `NA` sex, `NA` kinship, `kin = NULL` — run
- Independent check: expected offspring mean kinship `c' K c` for the targets, equal contributions, unconstrained optimal contributions (0.0533, with negative contributions) and non-negative optimal contributions by `quadprog::solve.QP` (0.0560); 200 random allocations (median 0.0720)
- Missing data: targets compared unfiltered against call-rate filtered
- PMx Repro Goals (Manual pp. 90-91): SKIPPED — manual not available; the function documents that it replaces the PMx simulation with a closed-form allocation
- Callers (API3): none outside `@seealso`; grep of dartR.base, dartR.popgen, dartR.sim, dartR.sexlinked, dartR.spatial and dartr2shiny finds no code callers
- FBM path (DAT6): SKIPPED — no FBM fixture
- Known complaints: not checked — added in the 2026 kinship series, no release history
- PLT: not applicable — no plot

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | consequence (targets instead of an error for NA kinship) approved |
| 3 | approved | Luis | consequence (targets instead of an error for NA sex) approved |
| 4 | approved | Luis | consequence (invalid n.target stops) approved |
| 5 | approved | Luis | |
| 6 | approved | Luis | |

## Outcome

- Change 1: warning at `verbose >= 1` (24 individuals below 0.8 in the unfiltered colony, none after filtering); `@details` paragraph; the example filters loci at call rate 0.95 and passes the full-dataset `kin`.
- Change 2: `rowMeans(kin, na.rm = TRUE)`; warning with the number of `NA` values; every individual gets a target and each sex still sums to `n.target`.
- Change 3: `%in%` selection; an `NA` sex is excluded with the existing unknown-sex warning (23 rows, sums 24 and 24).
- Change 4: 0, -3, `NA`, `"12"` and `c(10, 12)` stop with `n.target must be a single number >= 1, or NULL for nInd(x)`; 12.7 is still rounded down to 12 with a warning.
- Change 5: `@description` says the targets lower the next generation's expected mean kinship; `@details` gives 0.0582 / 0.0667 / 0.0560.
- Change 6: `build =` dropped, DOC2 `verbose` text.
- Snapshot: targets, mean kinship, sums and the offspring mean kinship (0.0582) unchanged. Three baseline expectations changed, each mapped to an approved change: invalid `n.target` (4), `NA` sex (3), `NA` kinship (2). 25 expectations pass; full suite 510 pass, 0 fail; example runs.
- NEWS entry added. PR: pending

```json
{
  "function": "gl.report.repro.targets",
  "package": "dartR.captive",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "a2976e3",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS3,DOC2", "status": "approved", "change": 6}
  ],
  "coverage_skipped": ["PMx manual: not available", "DAT6: no FBM fixture", "forum/issues: no release history"],
  "status": "awaiting-approval",
  "pr": null
}
```
