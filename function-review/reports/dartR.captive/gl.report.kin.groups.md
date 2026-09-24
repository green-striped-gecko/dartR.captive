# Review: gl.report.kin.groups (dartR.captive)
- Family mode: report
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: e5b5708 (origin/dev, after #113-#116)
- Datasets: testset2.gl (274 x 755 SNPs; EmmacCaptBred subset of 24 in 5 cohorts), the same after `gl.filter.callrate(method = "loc", threshold = 0.95)` (438 loci), testset2.gs; dartR.data 1.2.6
- Baseline: tests/testthat/test-gl.report.kin.groups.R (snapshot captured pre-review, 18 expectations, all pass)

## Verdict

**Standards: Needs work**: the structure is sound, the function is
read-only, and the block-mean algebra is correct: `GD` equals
`1 - mean(kin)` to machine precision, as `@details` states. It also
accepts a larger reference matrix through `utils.kin.check` (#113). The
remaining gaps are small.
**Spec: Needs work**: `MK` and `GD` carry no information when `kin` is
self-referenced, which includes the default `kin = NULL`. Every group
gets `MK = 0` and `GD = 1`. The #113 audit wrongly listed this function
as unaffected. Missing data also biases `meanF`, `MK` and `GD`
substantially.

## Findings

**F1 [HIGH, confidence: high] — self-referenced kinship makes MK 0 and GD 1 (DOC5, proposed rule; corrects the #113 audit)**
`R/gl.report.kin.groups.r:98,150-151` — `MK_g = sum_h n_h f[g,h] / N`
averages group g's members over every individual, so it equals the mean
of their dataset-wide row means. `gl.kin` centres kinship on the
individuals it is estimated on, so each of those row means is 0.
`gl.report.kin.groups(testset2.gl)` returns `max |MK| = 1e-17` and
`gd = 1`. Run on the captive colony alone, it returns `MK = 0` for all
five cohorts and `gd = 1`. The group matrix `kin.groups` stays
informative. `kinship-suite-rowsum-audit.md` lists this function under
"Not affected (use population blocks)". That is wrong: the blocks feed
`MK` through a weighted sum over all groups, which recovers the row means.
Failure scenario: a manager calls
`gl.report.kin.groups(colony, group.col = "enclosure")`, gets identical
`MK = 0` for every enclosure, and concludes that no group is more related
than another.
Proposed change: call `utils.kin.check(..., need.reference = TRUE)`, as
the five functions in #113 do, so a self-referenced `kin` stops with the
standard message. Correct the audit entry.

**F2 [MEDIUM, confidence: high] — missing data biases meanF, MK and GD (DOC5, proposed rule)**
`R/gl.report.kin.groups.r:154` — `gl.kin` fills missing genotypes with
the locus mean, which shrinks kinship and self-kinship toward 0 for
individuals with low call rates (all 24 captive-bred animals have
0.70-0.80). Unfiltered, the cohorts' `meanF` is -0.25 to -0.38 and
`GD` is 0.934. After filtering loci at call rate 0.95, `meanF` is -0.15
to +0.01, `MK` rises by 30-70% (`F2_X` 0.078 to 0.132), and `GD` is
0.905.
Failure scenario: strongly negative `meanF` is read as outbreeding, and
`GD` overstates the diversity retained.
Proposed change: warn at `verbose >= 1` when any individual's call rate is
below 0.8 (the same warning as #115 and #116), and document the effect
in `@details`.

**F3 [MEDIUM, confidence: high] — one NA kinship makes MK and GD NA (DAT integrity)**
`R/gl.report.kin.groups.r:145` — `mean()` without `na.rm`. One `NA` pair
makes the `MK` of both groups involved `NA`, and `gd` `NA`.
Failure scenario: an estimator returns `NA` for one poorly genotyped
individual, and the whole-population `GD` is lost without a message.
Proposed change: block means with `na.rm = TRUE`, plus a `verbose >= 1`
warning with the number of `NA` pairs. This follows #112 and #114.

**F4 [LOW, confidence: high] — SilicoDArT meanF is 0 by construction, undocumented (DOC5, proposed rule)**
`R/gl.report.kin.groups.r:154` — the dominant estimator fixes the
diagonal at 0.5, so `meanF` is 0 for every group. `gl.report.kinship`
documents this; this function does not.
Failure scenario: a user reads `meanF = 0` as "no inbreeding".
Proposed change: one sentence in `@details` and `@return`, matching
`gl.report.kinship`.

**F5 [LOW, confidence: high] — convention gaps (FS3, VRB2)**
`:90-92` uses the outdated `build =`, and `:165-173` prints the
`verbose >= 3` summary with raw `cat()`. Failure scenario: none for
results.
Proposed change: drop `build =`, and use `report()`.

Checked, nothing found: the input comes back identical, the block means
and `GD` identity hold, missing-group individuals are dropped with a
warning, and an absent column is a fatal error. The example's "F1_AB vs
F1_AE elevated (shared sire)" holds: 0.082, the highest F1-F1 value.

## Proposed changes

1. Require reference kinship (`need.reference = TRUE`) and correct the
   rowsum audit (F1). **Consequence: `kin = NULL`, or a `kin` estimated on
   `x` alone, now stops with an error instead of returning `MK = 0`,
   `GD = 1`.**
2. Warn on individual call rate below 0.8 and document the missing-data
   effect (F2). No numerical change.
3. `NA`-tolerant block means with a warning (F3). **Consequence:
   numerical output changes for a `kin` containing `NA` (`MK` and `GD`
   are computed from the remaining pairs instead of being `NA`).**
4. Document that SilicoDArT `meanF` is 0 by construction (F4). Docs
   only.
5. Standards cleanup (F5). No output change.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: behaviour against roxygen on testset2.gl (populations as groups, captive cohorts with reference and self-referenced kinship, a larger reference matrix), testset2.gs, `NA` kinship, missing group values, absent column — run
- Independent numerical check: `GD` against `1 - mean(kin)` on the subset; filtered against unfiltered data
- Callers (API3): none in dartR.captive other than `@seealso`
- FBM path (DAT6): SKIPPED — no FBM fixture
- Known complaints: not checked — added in the 2026 kinship series, no release history
- PLT: not applicable — no plot

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence (self-referenced kin errors) approved |
| 2 | approved | Luis | |
| 3 | approved | Luis | consequence (MK/GD from remaining pairs instead of NA) approved |
| 4 | approved | Luis | |
| 5 | approved | Luis | |

## Outcome

- Change 1: `utils.kin.check(..., need.reference = TRUE)`; `@param kin` and `@details` follow the #113 wording. `gl.report.kin.groups(testset2.gl)` and a colony-only `kin` now stop with the standard message. Audit entry corrected in `kinship-suite-rowsum-audit.md`.
- Change 2: call-rate warning (24 individuals on the captive subset) and `@details` paragraph with the filtered comparison.
- Change 3: `na.rm = TRUE` in block means and `meanF`; warning with the number of `NA` pairs. With one `NA` pair, `MK` and `GD` are values, not `NA`.
- Change 4: `@details` and `@return` state that SilicoDArT `meanF` is 0 by construction.
- Change 5: `build =` dropped; `verbose >= 3` summary through `report()`.
- Snapshot: 3 baseline tests errored, all from change 1 (tests passing a self-referenced `kin`); rewritten with a reference `kin`. Reference-kinship values (`MK`, `meanF`, `GD`, `kin.groups`) unchanged, which confirms changes 2, 4 and 5 leave output alone. 20 expectations pass; example runs.
- NEWS entry added. PR: (pending)

```json
{
  "function": "gl.report.kin.groups",
  "package": "dartR.captive",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "e5b5708",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS3,VRB2", "status": "approved", "change": 5}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "forum/issues: no release history"],
  "status": "in-apply",
  "pr": null
}
```
