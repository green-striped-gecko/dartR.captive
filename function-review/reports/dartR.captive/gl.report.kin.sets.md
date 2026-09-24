# Review: gl.report.kin.sets (dartR.captive)
- Family mode: report
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: fdb8967 (origin/dev, after #121; `meanMK` is block-based since #113)
- Datasets: testset2.gl (274 individuals, 31 populations), EmmacCaptBred alone and split into random halves, the same after `gl.filter.callrate(method = "loc", threshold = 0.95)`, testset2.gs; dartR.data 1.2.6
- Baseline: tests/testthat/test-gl.report.kin.sets.R (snapshot captured pre-review, 16 expectations, all pass)

## Verdict

**Standards: Needs work**: the structure and read-only behaviour conform,
and `meanMK`, `GD.w`, `meanF`, `mkb` and `fst` match their documented
formulas exactly. `@param kin` still carries the pre-#113 wording.
**Spec: Needs work**: with one population and `kin` estimated on it
(including the default `kin = NULL`), `meanMK` is 0 and `GD.w` is 1 by
construction. The only warning printed is that the between-set matrices
are trivial. One `NA` voids a population's statistics and all its Fst
values. The kinship Fst has an upward sampling term from self-kinship,
and the documentation does not mention it.

## Findings

**F1 [MEDIUM, confidence: high] — one population with self-referenced kinship gives GD.w = 1 (DOC5, proposed rule)**
`R/gl.report.kin.sets.r:100,132-133` — the block mean is informative
when `x` holds several populations. With one population and `kin`
estimated on it, the block is the whole reference, so `meanMK` is 0
(-5e-18) and `GD.w` is 1 exactly. That is the #113 degeneracy, reached
through the single-population case. `gl.report.kin.sets(colony)` returns
`GD.w = 1` with only the "between-set matrices are trivial" warning.
Failure scenario: a manager reports the captive colony's gene diversity
as 1.
Proposed change: when `nPop(x) == 1`, call `utils.kin.check` with
`need.reference = TRUE`. Multi-population calls keep accepting `kin`
estimated on `x`, where the blocks are informative.

**F2 [MEDIUM, confidence: high] — one NA kinship voids a population's statistics (DAT integrity)**
`R/gl.report.kin.sets.r:132-160` — `mean()` without `na.rm`. One `NA`
pair inside EmmacCaptBred makes its `meanMK` and `GD.w` `NA`, along with
all 60 Fst entries involving it.
Failure scenario: `kin` has an `NA` for one poorly genotyped animal, and
the colony drops out of the Fst matrix with no message.
Proposed change: block means with `na.rm = TRUE`, with a `verbose >= 1`
warning giving the count, as in #119-#124.

**F3 [MEDIUM, confidence: high] — low call rates bias GD.w and meanF (DOC5, proposed rule)**
`R/gl.report.kin.sets.r:100` — the captive colony's `GD.w` is 0.934
(0.905 after filtering loci at call rate 0.95), and its `meanF` is -0.33
(-0.07 after filtering).
Failure scenario: strongly negative `meanF` is read as outbreeding.
Proposed change: a `verbose >= 1` call-rate warning and a `@details`
sentence, as in the other kinship functions.

**F4 [LOW, confidence: medium] — Fst carries an upward self-kinship term, undocumented (DOC5, proposed rule)**
`R/gl.report.kin.sets.r:159-162` — `GDt.pair` averages over the pooled
block including the diagonal (self-kinship about 0.5), while `GDb.pair`
uses cross-pairs only. Even with no differentiation, `Fst` is therefore
positive, by roughly `1/(2(n_s + n_t))`. Two random halves of the
captive colony give `Fst = 0.0098`. Across the 465 population pairs of
testset2.gl, Fst runs from 0.032 to 0.30 (median 0.10). For pairs of
about 10 individuals the term is about 0.025. The formula follows the
PMx definition that `@details` cites, which I could not check against
the PMx manual. It differs from Nei's `1 - mean(GD.w)/GD.t`, which gives
0.084 against 0.055 for the colony and EmmacMaclGeor.
Failure scenario: a user reads Fst = 0.03 between two small sets as
differentiation, when it is mostly the sampling term.
Proposed change: document the term, and the random-halves value as a
reference point, in `@details`. The formula stays unchanged pending the
custodian's confirmation of the PMx definition.

**F5 [LOW, confidence: high] — convention and documentation gaps (FS3, VRB2, DOC5)**
- `:92-94` uses the outdated `build =`.
- `:171-177` prints the `verbose >= 3` summary with raw `cat()`.
- `@param kin` says kinship must have "row and column names identical to
  `indNames(x)`", but since #113 a larger reference matrix is accepted
  and restricted.
- `@details` says kinship "is centred on the individuals of x", which is
  true only for `kin` estimated on `x`.
- SilicoDArT `meanF` is 0 by construction, which is undocumented.

Failure scenario: none for results. Users are misled about which `kin`
to pass.
Proposed change: drop `build =`, use `report()`, and update `@param kin`
and `@details` to the #113 wording, including the SilicoDArT `meanF` note.

Checked, nothing found: the input comes back identical, the `mkb`
diagonal equals `meanMK`, the `fst` diagonal is 0, singleton populations
get a warning, and the SilicoDArT path runs.

## Proposed changes

1. Require reference kinship when `x` has one population (F1).
   **Consequence: single-population calls with `kin = NULL` or `kin`
   estimated on `x` now error instead of returning `GD.w = 1`.**
2. `NA`-tolerant block means with a warning (F2). **Consequence: `kin`
   with `NA` gives values instead of `NA` for the affected population and
   its Fst.**
3. Call-rate warning and `@details` sentence (F3). No numerical change.
4. Document the Fst self-kinship term (F4). Docs only.
5. Standards and documentation cleanup (F5). No output change.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: behaviour against roxygen on the full testset2.gl, a single population, random halves of one population, `NA` kinship, testset2.gs — run
- Independent numerical check: `mkb`, `fst`, `GD.w` recomputed for one pair; Nei's Fst computed for comparison
- PMx Fst definition: not checked against the PMx manual — reason for medium confidence on F4 and for proposing documentation only
- Callers (API3): none outside `@seealso`
- FBM path (DAT6): SKIPPED — no FBM fixture
- Known complaints: not checked — added in the 2026 kinship series, no release history
- PLT: not applicable — no plot

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence (single-population self-referenced kin errors) approved |
| 2 | approved | Luis | consequence (values instead of NA) approved |
| 3 | approved | Luis | |
| 4 | approved | Luis | formula kept pending the custodian's confirmation of the PMx definition |
| 5 | approved | Luis | |

## Outcome

- Change 1: `need.reference = (nPop(x) == 1)`; `gl.report.kin.sets(colony)` stops with the standard message; with full-dataset `kin` the colony gets `GD.w` 0.934.
- Change 2: `na.rm = TRUE` in every block mean and `meanF`; warning with the count. With one `NA` pair the colony's `GD.w` and all Fst values are computed.
- Change 3: call-rate warning (24 individuals) and a `@details` paragraph.
- Change 4: `@details` paragraph on the Fst self-kinship term with the random-halves reference value.
- Change 5: `@param kin` and `@details` in the #113 wording, SilicoDArT `meanF` note, `build =` dropped, `report()` in the summary.
- Snapshot: 3 baseline expectations changed (single-population GD.w, two NA counts), mapped to changes 1 and 2; values on the full dataset unchanged. 20 expectations pass; examples run.
- NEWS entry added. PR: #125

```json
{
  "function": "gl.report.kin.sets",
  "package": "dartR.captive",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "fdb8967",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "medium", "rule": "DOC5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS3,VRB2,DOC5", "status": "approved", "change": 5}
  ],
  "coverage_skipped": ["PMx manual not consulted for the Fst definition", "DAT6: no FBM fixture", "forum/issues: no release history"],
  "status": "pr-open",
  "pr": 125
}
```
