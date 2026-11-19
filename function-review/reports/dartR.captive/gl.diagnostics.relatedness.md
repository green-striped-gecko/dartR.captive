# Review: gl.diagnostics.relatedness (dartR.captive)
- Family mode: analysis (with helpers in `utils.functions.diagnostics.relatedness.r` and `utils.classes.diagnostics.relatedness.r`)
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: c66ad81 (origin/dev)
- Datasets: testset.gl[1:15, ] and [1:30, ] (gl.filter.allna); a 3-generation Wright-Fisher simulation from `dartR.sim` (`ref_variables.csv`, `sim_variables.csv` shipped with dartR.sim); a hand-made 7-individual pedigree
- Dependencies: `related` (GitHub, installed here); EMIBD9 path not exercised in this review (reviewed in PR #100)
- Baseline: tests/testthat/test-gl.diagnostics.relatedness.R (snapshot captured pre-review, 7 expectations)
- Earlier work: informal review and fixes on 2026-08-17 (commit 3bd045f; PRs #83, #84; commit 9a3f94a). This is the first review under the campaign.

## Verdict

**Standards: Needs work** — the August fixes hold (no crashes on the
tested inputs; names restored), but the result table carries duplicate
ID columns, the default simulation inputs are unusable, and there were
no tests.
**Spec: Rework** — the accuracy statistics this function exists to
produce are biased: pairs below 0.05 are discarded before summarising,
founders in an attached pedigree are labelled half sibs, some pairs get
two "true" relationships, and "RMSE" is the mean absolute error.

## Findings

**F1 [HIGH, confidence: high] — missing parents make founders half
sibs (DOC5)**
`R/utils.functions.diagnostics.relatedness.r:224-231` — the half-sib
joins merge on `dad` and on `mom` without excluding missing values;
data.table matches missing to missing. With an attached pedigree
(`includedPed = TRUE`), 0 is converted to missing
(`generateRelatedTableBaseInput`), so every founder pairs with every
other founder.
Failure scenario: 4 unrelated founders and 3 offspring: all 6 founder
pairs are labelled half sibs; their unrelated offspring O1-O3 and O2-O3
become "half first cousins"; the full-sib pair O1-O2 is also labelled
half first cousins. RMSE and variance are then computed against these
false relationships.
Proposed change: exclude missing parents from the half-sib (and derived
cousin) joins.

**F2 [HIGH, confidence: high] — pairs with two "true" relationships
(DOC5)**
`R/utils.functions.diagnostics.relatedness.r:289-298,456` — relationship
classes are built independently and all kept, so inbred families give a
pair two labels.
Failure scenario: 3-generation simulation: 15 of 467 pairs are both
full sibs (0.25) and first cousins (0.0625); each appears twice in
`MergedDf` and is scored against both values, inflating the cousin
error (e.g. `wang` 0.44 counted as a first-cousin estimate).
Proposed change (member's choice): (a) compute the exact pedigree
kinship of each pair (standard recursive kinship from the pedigree)
and use it as the true value, labelling the pair with its closest class;
or (b) keep one label per pair, the one with the highest expected
kinship.
**Consequence: RMSE and variance change for inbred pedigrees; each pair
appears once.**

**F3 [HIGH, confidence: high] — pairs below 0.05 are discarded before
the statistics (DOC5)**
`R/utils.functions.diagnostics.relatedness.r:80-81` — `cleanup_rel`
keeps only pairs with rrBLUP and every estimator at or above 0.05.
Failure scenario: `testset.gl[1:15, ]`: 25 of 105 pairs remain. In the
simulation, 36 of 222 second-cousin pairs (true kinship 0.016) remain,
all with estimates of at least 0.05, so their RMSE and variance describe
a truncated, upward-biased sample; unrelated pairs never enter the
comparison, so estimator bias at zero cannot be seen.
Proposed change: keep all pairs; with a pedigree, pairs with no
relationship in it are kept as "unrelated" (true kinship 0).
**Consequence: `MergedDf` holds all pairs; RMSE and variance change,
especially for distant classes; a new "unrelated" class appears.**

**F4 [MEDIUM, confidence: high] — "RMSE" is the mean absolute error
(DOC5)**
`R/utils.functions.diagnostics.relatedness.r:537-539` — `rmse()` is
applied to each value separately (giving |error|) and then averaged.
Failure scenario: simulated full sibs, `wang`: reported 0.064 (mean
absolute error); RMSE 0.080.
Proposed change: `sqrt(mean((estimate - truth)^2))` per class and
estimator.
**Consequence: every RMSE value changes (larger).**

**F5 [MEDIUM, confidence: high] — `run_sim = TRUE` fails with the
default arguments (DOC5)**
`R/gl.diagnostics.relatedness.r:91-92,257-262` — `ref_variables` and
`sim_variables` default to `NULL` (documented "[optional]"), but the
`DartSim` class requires file paths: "invalid class 'DartSim' object:
... got class NULL". The `@examples` call (`run_sim = TRUE`, no files)
fails.
Proposed change: default to the files shipped with dartR.sim
(`system.file("extdata", "ref_variables.csv", package = "dartR.sim")`
and `sim_variables.csv`); document them.

**F6 [LOW, confidence: high] — output table and smaller defects (FS3,
FS5, DAT1, DOC1, DOC7)**
- `MergedDf` carries leftover join columns `i.ind1`, `i.ind2`.
- SilicoDArT is documented as accepted, but `related` and `gl.grm` need
  SNP genotypes; stop for SilicoDArT.
- `cleanup = TRUE` calls `gl.filter.heterozygosity()` and
  `gl.filter.allna()` at the default verbosity, which print.
- `build = "Jody"`; `@author` "Ethan, Luis" without `Author(s)`/
  `Custodian`; typos ("attache", "bewteen" in class descriptions);
  `@param verbose` non-standard; non-ASCII en dash.
- Not fixable here: `related`'s Fortran code prints its progress
  ("Dyad of individuals of ...") directly to the console at every
  verbosity; R cannot capture it.
Proposed change: fix each fixable item.

## Proposed changes

1. Exclude missing parents from sibling/cousin joins (F1).
   **Consequence: attached pedigrees no longer create false relatives;
   RMSE/variance change for `includedPed = TRUE`.**
2. One true value per pair: (a) exact pedigree kinship, or (b) closest
   class only (F2).
   **Consequence: each pair appears once; statistics change for inbred
   pedigrees.**
3. Keep all pairs; "unrelated" class for pairs not related in the
   pedigree (F3).
   **Consequence: larger `MergedDf`; RMSE/variance change; new class.**
4. True RMSE (F4).
   **Consequence: RMSE values change.**
5. Default simulation files from dartR.sim (F5).
6. Output columns, SilicoDArT check, quiet cleanup filters, docs (F6).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: no-simulation run on testset.gl; simulation run with dartR.sim's
  variable files; pedigree classifier on the simulation (full sibs 179
  = independent count) and on a hand-made pedigree with missing parents;
  RMSE computed by hand — run
- EMIBD9 path (`run.e9 = TRUE`): not run here; `runE9` uses
  `gl.run.EMIBD9()$rel` (kinship) without halving, consistent with the
  other estimators (halved to kinship)
- Callers: none in dartr2shiny or other `dartR.*` packages
- Tests: 5 warnings in the baseline come from the estimators on small
  data, not from the function

## Addendum (found while applying)

**A1 [LOW, confidence: high] — internal calls print at `verbose = 0`**
`coanct_clean()` calls `gl2related()` and `GRM_clean()` calls `gl.grm()`
without `verbose`, so "Starting gl2related / Completed: gl2related /
Starting gl.grm / Completed: gl.grm" print at `verbose = 0`. Change 6
covered the cleanup filters only. Proposed: pass `verbose = 0`.
Approved by Luis Mijangos and applied: 0 R-level lines at `verbose = 0`.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis Mijangos | |
| 2 | approved | Luis Mijangos | option (a): exact pedigree kinship; class = closest relationship the classifier finds |
| 3 | approved | Luis Mijangos | all pairs kept, "unrelated" class |
| 4 | approved | Luis Mijangos | |
| 5 | approved | Luis Mijangos | |
| 6 | approved | Luis Mijangos | `@author` keeps "Ethan" without a surname (not known); custodian set to Luis Mijangos (assumption) |

## Outcome

All six approved changes applied on 2026-09-23 to
`R/gl.diagnostics.relatedness.r`, `R/utils.functions.diagnostics.relatedness.r`
and `R/utils.classes.diagnostics.relatedness.r`; `devtools::document()`
run. `tests/testthat/test-gl.diagnostics.relatedness.R`: 18 expectations
pass (estimator runs skip without `related`). The pre-change baseline
tests were replaced by tests of the approved behaviour.

- Change 1 (F1): the half-sib joins exclude missing parents. Hand-made
  pedigree: 0 half-sib pairs (was 6 founder pairs) and 0 half first
  cousins (was 3).
- Change 2 (F2, option a): new `pedigreeKinship()` computes exact
  kinship recursively (parents ordered before offspring; parents not in
  the pedigree added as unrelated founders; loops and duplicate IDs
  stop with an error). Equal to `kinship2::kinship()` on the hand-made
  pedigree (including the inbred O4, 0.625) and on the simulated
  pedigree (186 IDs). Each pair keeps one class, the closest relationship
  the classifier finds; `rel` is the exact kinship, so inbred full sibs
  carry their true kinship, not 0.25. In the simulation, pairs appear
  once (was: 15 pairs twice).
- Change 3 (F3): the 0.05 filter is removed. `testset.gl[1:15, ]` keeps
  all 105 pairs (was 25); the simulation keeps all 11,175 pairs, 10,422
  of them "unrelated".
- Change 4 (F4): RMSE = `sqrt(mean((estimate - rel)^2))` per class;
  full sibs, `wang`: 0.0953, equal to the hand computation. Classes
  without pairs give NA (was a misleading 0 or NaN).
- Change 5 (F5): `ref_variables`/`sim_variables` default to the files
  shipped with dartR.sim; `run_sim = TRUE` runs without supplying them
  (4.3 s for 30 individuals, 3 generations).
- Change 6 (F6): `MergedDf` columns `ind1, ind2, RelDegree, rel,
  <estimators>` (no `i.ind1`/`i.ind2`); SilicoDArT errors; cleanup
  filters run at `verbose = 0`; docs, `@family`, `Author(s)`/`Custodian`,
  typos, ASCII.
- Defect caught in verification before commit: the first version of the
  attached-pedigree reader left founders' parents as "0"; the class
  labels then treated "0" as a shared parent (67 "full sibs" in a
  pedigree with one full-sib pair; `rel` was already correct). Fixed by
  converting 0 and "" to NA before classifying; the attached-pedigree
  test now checks 1 full-sib, 6 parent-offspring and 98 unrelated pairs.
- A1: `gl2related()` and `gl.grm()` called with `verbose = 0`.
- PR #107 (`review-gl.diagnostics.relatedness` -> `dev`).
- Not in this change: the EMIBD9 path was not run
  (no change to it beyond the shared merge).

```json
{
  "function": "gl.diagnostics.relatedness",
  "package": "dartR.captive",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "c66ad81",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved", "change": 6},
    {"id": "A1", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "approved", "change": "A1"}
  ],
  "coverage_skipped": ["run.e9 path: EMIBD9 run not exercised (reviewed in PR #100)"],
  "status": "done",
  "pr": 107
}
```
