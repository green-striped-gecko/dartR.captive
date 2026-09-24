# Review: gl.report.mate.suitability (dartR.captive)
- Family mode: report
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 15d100d (origin/dev, after #124)
- Datasets: testset2.gl, EmmacCaptBred (24 individuals: 11 males, 13 females, 143 pairings) with `kin` from the full testset2.gl; the same after `gl.filter.callrate(method = "loc", threshold = 0.95)` (755 to 438 loci); testset2.gs; 200 random individuals of testset2.gl for timing; dartR.data 1.2.6
- Baseline: tests/testthat/test-gl.report.mate.suitability.R (snapshot captured pre-review, 17 expectations, all pass)

## Verdict

**Standards: Needs work**: the structure, the read-only behaviour and the
`need.reference` guard from #113 are in place, and `dgd` and `mkdiff` match
an independent recomputation. Two inputs crash with R's own error message,
and there are small convention gaps.
**Spec: Needs work**: the documented PMx rule set is implemented as
described, but the do-not-breed screen compares shrunken kinship with a
fixed threshold. With call rates of 0.70-0.80, 17 of the 48 pairings that
are 'NoWay' after call-rate filtering receive a rating instead, and the
MSI changes for 61 of 143 pairings.

## Findings

**F1 [HIGH, confidence: high] — low call rates hide do-not-breed pairings (DOC5, proposed rule)**
`R/gl.report.mate.suitability.r:211, 258` — `f.off` is the pairwise
kinship, and `gl.kin` fills missing genotypes with the locus mean, which
pulls kinship toward 0. The No Way point (0.125) is absolute, so shrinkage
moves pairs below it. On the example colony (call rates 0.70-0.80), full
sibs CB_AB_01 x CB_AB_02 have `f.off` 0.214 unfiltered and 0.252 after
filtering loci at call rate 0.95 (expected 0.25). Across all pairings,
filtered `f.off` is 1.36 times unfiltered (regression slope). Unfiltered,
17 of the 48 filtered 'NoWay' pairings are rated instead: 4 as MSI 4,
5 as MSI 5 and 8 as MSI 6. Their filtered kinship is 0.126-0.191. The
MSI differs for 61 of 143 pairings overall.
Failure scenario: a manager screens the colony on unfiltered data and
breeds a pair of half sibs that the function rated MSI 4 ("slightly
detrimental") instead of 'NoWay'.
Proposed change: the series call-rate warning at `verbose >= 1` (any
individual below 0.8, as in #115, #116, #119, #120, #122), a `@details`
paragraph stating that the No Way point is compared with kinship that
missing data shrinks, with these numbers, and an example that filters
loci on call rate before `gl.kin`.

**F2 [MEDIUM, confidence: high] — one NA kinship crashes the function (DAT integrity)**
`R/gl.report.mate.suitability.r:218-245` — `utils.kin.dgd` and
`rowMeans`/`mean` run without `na.rm`. With one `NA` pair, `f.mean` is
`NA` and the function stops at line 245 with
`missing value where TRUE/FALSE needed`.
Failure scenario: `kin` from the dominant estimator has an `NA` for two
individuals that share no scored loci; the user gets R's error with no
hint of the cause.
Proposed change: ignore `NA` in the gene diversity, mean kinship and
mean-F calculations (`na.rm = TRUE`, added to `utils.kin.dgd` in #120).
A pairing whose own `f.off` is `NA` cannot be screened, so its MSI is
`NA` (neither rated nor 'NoWay'). A `verbose >= 1` warning gives the
number of `NA` pairs.

**F3 [MEDIUM, confidence: high] — an NA sex crashes the function (DAT integrity)**
`R/gl.report.mate.suitability.r:170-172` — `ids[sex == "Male"]` returns
`NA` ids when `sex` is `NA`. The function warns that 1 individual is
excluded, then stops at line 211 with `subscript out of bounds`.
Failure scenario: an individual with an empty sex cell in the metadata
file makes the whole report fail with an error that does not mention sex.
Proposed change: select with `sex %in% "Male"` and `sex %in% "Female"`,
so `NA` counts as unknown sex, as documented. `gl.report.repro.targets`
(line 104) and `gl.select.pairs` (line 127) have the same pattern; they
are next in the campaign and are left to their own reviews.

**F4 [LOW, confidence: high] — an invalid `f.noway` is replaced, not rejected (FS5)**
`R/gl.report.mate.suitability.r:185-191` — a value outside (0, 0.5] is
reset to 0.125 with a warning, and the warning is silent at
`verbose = 0`. An invalid `unknown.breaks` is a fatal error.
Failure scenario: `f.noway = 25` (a percentage) runs at 0.125 without a
message in a quiet script.
Proposed change: `stop(error(...))` for an invalid `f.noway`.
**Consequence: calls with an out-of-range `f.noway` stop instead of
running at 0.125.**

**F5 [LOW, confidence: high] — per-pair matrix rebuild and densification (STY2, DAT6 proposed rule)**
`R/gl.report.mate.suitability.r:221-226, 236` — each pairing calls
`utils.kin.dgd`, which copies the `(n + 1) x (n + 1)` matrix. Work grows
with pairings times `n^2`: 1.7 s for 95 x 94 pairings with `n = 200`; a
1,000-animal colony with 500 of each sex takes about 700 times as long
(estimated, not run). The call rate uses `as.matrix(x)`.
Failure scenario: a large colony takes tens of minutes for a report
that has an exact closed form.
Proposed change: compute `dgd` in closed form,
`S / n^2 - (S + rs[m] + rs[f] + 0.5 * (1 + kin[m, f])) / (n + 1)^2`,
where `S` is the sum of `kin` and `rs` its row sums (verified equal to the current output within 1.1e-16), and the call rate
from `NA.posi`, as in #122.

**F6 [LOW, confidence: high] — convention and example gaps (FS3, DOC2, DOC5)**
`:150-152` passes the outdated `build =`. The `verbose` text differs from
DOC2. The example comments say full sibs have `f.off` "~0.23" (it is
0.214) and the cross-family pair "~0.01" (it is -0.009).
Failure scenario: none for results.
Proposed change: drop `build =`, use the DOC2 text, correct the example
comments (they change again with F1's filtered example).

Checked, nothing found: the input comes back identical; `kin = NULL` and
self-referenced `kin` stop with the #113 message; `dgd` and `mkdiff`
equal an independent recomputation; bin edges give ranks 1-6 at the
range limits and break points; SilicoDArT runs; unknown sex warns at
`verbose >= 1`.

## Proposed changes

1. Call-rate warning, `@details` paragraph and filtered example (F1). No
   numerical change.
2. `na.rm = TRUE` in the GD, MK and mean-F calculations, MSI `NA` for
   pairings with `NA` kinship, with a warning (F2). **Consequence: a `kin`
   with `NA` now returns results instead of an error.**
3. `NA` sex treated as unknown (F3). **Consequence: data with an `NA` sex
   now returns results instead of an error.**
4. Invalid `f.noway` is a fatal error (F4). **Consequence: calls with an
   out-of-range `f.noway` stop instead of running at 0.125.**
5. Closed-form `dgd` and `NA.posi` call rate (F5). Output identical
   within 1.1e-16.
6. Standards and example-comment cleanup (F6). No output change.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: behaviour against roxygen on EmmacCaptBred (SNP and SilicoDArT), `unknown.breaks` at the PMx defaults, invalid `f.noway`, `NA` sex, `NA` kinship, self-referenced `kin` — run
- Independent numerical check: `dgd` by closed form and by an explicit augmented matrix, `mkdiff` from row means (maximum difference 1.1e-16)
- Missing data: MSI compared unfiltered against call-rate filtered
- Tulsa rules and bin definitions against the PMx manual (pp. 108-111): SKIPPED — manual not available in this session; the code was checked against the roxygen description only
- Callers (API3): none outside `@seealso`; grep of dartR.base, dartR.popgen, dartR.sim, dartR.sexlinked, dartR.spatial and dartr2shiny finds no code callers
- FBM path (DAT6): SKIPPED — no FBM fixture
- Known complaints: not checked — added in the 2026 kinship series, no release history
- PLT: not applicable — no plot

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | consequence (results instead of an error for NA kinship) approved |
| 3 | approved | Luis | consequence (results instead of an error for NA sex) approved |
| 4 | approved | Luis | consequence (out-of-range f.noway stops) approved |
| 5 | approved | Luis | |
| 6 | approved | Luis | |

## Outcome

- Change 1: warning at `verbose >= 1` (24 individuals below 0.8 in the unfiltered colony, none after filtering); `@details` paragraph; the example filters loci at call rate 0.95 and passes the full-dataset `kin` (NoWay 48, full sibs 0.252, CB_AB_01 x CB_CD_02 `f.off` 0.016 and MSI 5).
- Change 2: `na.rm = TRUE` in the gene diversity, mean kinship and mean-F calculations; the pairing with `NA` kinship gets MSI `NA`, every other pairing is rated; warning with the number of `NA` values.
- Change 3: `%in%` selection; an `NA` sex is excluded with the existing unknown-sex warning.
- Change 4: out-of-range `f.noway` stops with `f.noway must be a single value in (0, 0.5]`.
- Change 5: closed-form `dgd` when `kin` has no `NA` (the `utils.kin.dgd` loop with `na.rm = TRUE` otherwise); call rate from `NA.posi`. 200 individuals: 0.04 s (was 1.7 s).
- Change 6: `build =` dropped, DOC2 `verbose` text, example comments corrected.
- Snapshot: MSI table, break points, `f.off` and the closed-form `dgd` check unchanged. Three baseline expectations changed, each mapped to an approved change: invalid `f.noway` (4), `NA` sex (3), `NA` kinship (2). 22 expectations pass; full suite 487 pass, 0 fail; example runs.
- NEWS entry added. PR: pending

```json
{
  "function": "gl.report.mate.suitability",
  "package": "dartR.captive",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "15d100d",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "STY2,DAT6", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS3,DOC2,DOC5", "status": "approved", "change": 6}
  ],
  "coverage_skipped": ["PMx manual: not available", "DAT6: no FBM fixture", "forum/issues: no release history"],
  "status": "awaiting-approval",
  "pr": null
}
```
