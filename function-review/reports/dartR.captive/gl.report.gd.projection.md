# Review: gl.report.gd.projection (dartR.captive)
- Family mode: report
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: fdb8967 (origin/dev, after #121; the reference-kinship guard from #113 is in place)
- Datasets: testset2.gl (EmmacCaptBred, 24 individuals, with kinship from the full 274), the same after `gl.filter.callrate(method = "loc", threshold = 0.95)`; direct `gd.now` inputs
- Baseline: tests/testthat/test-gl.report.gd.projection.R (snapshot captured pre-review, 24 expectations, all pass)

## Verdict

**Standards: Needs work**: the structure and read-only behaviour conform,
the #113 guard works, the plot is saved through `utils.plot.save`, and
every formula (projection, `years.to.target`, `ne.required`, generation
length) matches an independent computation. There are small gaps in
ordering and conventions.
**Spec: Needs work**: the retention target is measured against *current*
gene diversity. The retention outputs (`prop.retained`,
`years.to.target`, `ne.required`) are therefore identical whatever the
data: `gd.now = 0.5` and `0.95` give the same values. The PMx Goals
convention that `@details` cites is a proportion of the *source*
population's gene diversity, and `gd.now` computed from reference
kinship is already on that scale. For the captive colony at Ne = 50,
source-relative GD falls below 0.9 after 3.7 years, not 10.5, and it
needs Ne ≈ 1360, not 475. This is a design choice for the custodian.

## Findings

**F1 [HIGH, confidence: medium] — retention is relative to current, not source, gene diversity (DOC5, proposed rule)**
`R/gl.report.gd.projection.r:206-213` — `prop.retained = gd / gd.now`,
`years.to.target = L log(target) / log(lambda)` and `ne.required` depend
only on `ne`, `gen.length`, `years` and `gd.target`. The genomic input
(`x` and `kin`, or `gd.now`) changes only the `gd` column and the plot
title. `@details` presents the function as the analogue of the PMx Goals
screen and its convention of 90% over 100 years. That goal, following
Soulé et al. (1986), is 90% of the gene diversity of the *source*
population, and PMx reports current GD as a proportion of source GD.
With `kin` estimated on a reference that includes the source
populations, `gd.now = 1 - mean(kin)` is on that scale: the reference
itself has GD = 1 (`gl.kin`, grm method). Source-relative, the captive
colony (`gd.now` 0.934) falls to 0.9 after 3.7 years at Ne = 50, where
the current-relative figure is 10.5. Holding GD at 0.9 or above for
100 years needs Ne = 1360, where the current-relative figure is 475.
Failure scenario: a manager reads "Ne required for target: 474.8" as the
PMx goal and plans for a population that loses 90% of source GD within
four years.
Confidence is medium because PMx's exact definition is cited from the
literature, not checked against the PMx manual pages named in
`@details`.
Proposed change: keep the current-relative outputs unchanged, and add
source-relative ones: `years.to.target.source =
L log(target / gd.now) / log(lambda)` and `ne.required.source =
1 / (2 (1 - (target / gd.now)^(L / years)))`. Both are `NA`, with a
warning, when `gd.now <= gd.target` (the goal is already missed). Add a
`gd.target` line on the absolute scale to the plot, and state in
`@details` which output answers the PMx question.

**F2 [MEDIUM, confidence: high] — one NA kinship gives an opaque error (DAT integrity)**
`R/gl.report.gd.projection.r:161` — `1 - mean(kin)` without `na.rm`
makes `gd.now` `NA`. The range check then fails with "missing value where
TRUE/FALSE needed".
Failure scenario: `kin` has an `NA` for one poorly genotyped animal, and
the user gets an R-internal error.
Proposed change: `mean(kin, na.rm = TRUE)` with a `verbose >= 1` warning
giving the count, as in #119, #120, #122 and #123.

**F3 [MEDIUM, confidence: high] — low call rates inflate gd.now (DOC5, proposed rule)**
`R/gl.report.gd.projection.r:159-161` — the captive colony's `gd.now` is
0.934, and 0.905 after filtering loci at call rate 0.95, because mean
imputation pulls kinship toward 0. Under F1's source-relative reading
that difference decides the outcome: 0.905 leaves 0.6 years before the
90% goal is missed, and 0.934 leaves 3.7 years.
Failure scenario: an unfiltered `gd.now` overstates the diversity left
and the time available to act.
Proposed change: a `verbose >= 1` call-rate warning when `x` is supplied,
and a `@details` sentence, as in the other kinship functions.

**F4 [LOW, confidence: high] — convention gaps (FS5, FS3, VRB2)**
- `:153-167` validates `ne` only after `gl.kin`/`utils.kin.check` has
  run, so a missing `ne` is reported after the kinship work.
- `:132-134` uses the outdated `build =`.
- `:221-246` prints the `verbose >= 3` summary with raw `cat()`.

Failure scenario: none for results.
Proposed change: validate the scalar arguments before resolving `gd.now`,
drop `build =`, and use `report()`.

Checked, nothing found: the input comes back identical, the formulas
match, `gen.length` scales time correctly, self-referenced `kin` is
rejected, invalid arguments error, and the plot saves to `plot.dir`.

## Proposed changes

1. Add `years.to.target.source` and `ne.required.source` to `summary`,
   `NA` with a warning when `gd.now <= gd.target`; add an absolute-scale
   target line to the plot; `@details` explains both readings (F1).
   **Consequence: `summary` gains two elements; existing elements and
   `projection` are unchanged.**
2. `NA`-tolerant `gd.now` with a warning (F2). **Consequence: `kin` with
   `NA`, which errored, now returns results.**
3. Call-rate warning and `@details` sentence (F3). No numerical change.
4. Validate scalars first; standards cleanup (F4). No output change.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: behaviour against roxygen on the captive colony with reference kinship, direct `gd.now`, `gen.length`, `NA` kinship, invalid inputs, plot file — run
- Independent numerical check: projection, `years.to.target` and `ne.required` recomputed; source-relative alternatives computed
- PMx Goals definition: from the literature (Soulé et al. 1986; the PMx 90%/100-year convention), not from the PMx manual pages — reason for medium confidence on F1
- Callers (API3): none outside `@seealso`
- FBM path (DAT6): SKIPPED — no FBM fixture
- Known complaints: not checked — added in the 2026 kinship series, no release history
- PLT: `plot.theme`, `plot.colors` and `utils.plot.save` idioms present (PLT1, PLT2); results do not depend on plotting (PLT3)

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence (summary gains two elements, plot gains a line) approved; worth confirming the PMx reading with Arthur |
| 2 | approved | Luis | consequence (results instead of an error for NA-containing kin) approved |
| 3 | approved | Luis | |
| 4 | approved | Luis | |

## Outcome

- Change 1: `years.to.target.source` and `ne.required.source` added (3.66 years and Ne 1360.4 for the colony at Ne = 50), `NA` with a warning when `gd.now <= gd.target`; the `verbose >= 3` summary has current-relative and source-relative sections; the plot adds a dot-dash line at `gd.target / gd.now`; `@details` and `@return` explain both readings.
- Change 2: `gd.now = 1 - mean(kin, na.rm = TRUE)` with a warning giving the count.
- Change 3: call-rate warning when `x` is supplied (24 individuals) and a `@details` paragraph.
- Change 4: `ne`, `gen.length`, `years`, `gd.target` and `n` are validated before `gd.now` is resolved (a missing `ne` now errors before `utils.kin.check`); `build =` dropped; `report()` in the summary.
- Snapshot: 2 baseline expectations changed (`summary` names, the NA error), mapped to changes 1 and 2; existing values unchanged. 35 expectations pass; examples run.
- NEWS entry added. PR: #124

```json
{
  "function": "gl.report.gd.projection",
  "package": "dartR.captive",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "fdb8967",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "medium", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "FS5,FS3,VRB2", "status": "approved", "change": 4}
  ],
  "coverage_skipped": ["PMx manual not consulted directly", "DAT6: no FBM fixture", "forum/issues: no release history"],
  "status": "pr-open",
  "pr": 124
}
```
