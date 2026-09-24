# Review: gl.report.ind.add (dartR.captive)
- Family mode: report
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: d6e62f2 (origin/dev, after #117)
- Datasets: testset2.gl (274 x 755 SNPs; target EmmacCaptBred, 24 individuals; 250 wild candidates), the same after `gl.filter.callrate(method = "loc", threshold = 0.95)`, testset2.gs; dartR.data 1.2.6
- Baseline: tests/testthat/test-gl.report.ind.add.R (snapshot captured pre-review, 15 expectations, all pass)

## Verdict

**Standards: Needs work**: the structure, input validation and read-only
behaviour conform. `dgd` matches an independent recomputation exactly
(maximum difference 0). There are small convention gaps.
**Spec: Needs work**: the method is right. The recorded founders of the
captive colony (`AA019158`, `AA000307`) rank last among 250 wild
candidates, as they should. One `NA` in `kin` silently erases every
result, and the target's low call rates shift the ranking. The #113 audit
listed this function as unaffected by self-referenced kinship, and that
holds: it uses only the target's block, and `x` always contains more than
the target.

## Findings

**F1 [MEDIUM, confidence: high] — one NA in the target block makes every dgd NA (DAT integrity)**
`R/gl.report.ind.add.r:134-139` — `utils.kin.dgd` takes `mean(kin)`
without `na.rm`, so one `NA` among the target's pairs makes `gd.target`
`NA`. Every candidate's `dgd` then becomes `NA` too, and the ranking is
arbitrary.
Failure scenario: an estimator returns `NA` for one poorly genotyped
colony member, and the user gets a ranked table of `NA` with no message.
Proposed change: add an `na.rm` argument (default `FALSE`) to
`utils.kin.dgd`, so the other callers keep their behaviour. Pass
`na.rm = TRUE` from this function, with a `verbose >= 1` warning giving
the number of `NA` pairs involved. This follows #112, #114 and #119.

**F2 [MEDIUM, confidence: high] — low call rates in the target shift the ranking (DOC5, proposed rule)**
`R/gl.report.ind.add.r:125` — the 24 captive-bred individuals have call
rates of 0.70-0.80, and mean imputation in `gl.kin` pulls their kinships
toward 0. The target's GD is 0.934 unfiltered and 0.905 after filtering
loci at call rate 0.95. The candidate ranking changes too: rank
correlation 0.90, and only 5 of the top 10 candidates are the same.
Failure scenario: a manager imports the unfiltered top candidates, half
of which drop out of the top 10 on better-called loci.
Proposed change: a `verbose >= 1` warning when any target or candidate
individual has call rate below 0.8, and a sentence in `@details`. This is
the same warning as #115, #116 and #119.

**F3 [LOW, confidence: high] — "the reported gain is exact" overstates it (DOC5, proposed rule)**
`R/gl.report.ind.add.r:32` — `dgd` is exact given `kin`, but `kin` is
relative to the individuals it was estimated on. With `kin` estimated on
the target and one candidate population only, instead of the full
dataset, the 11 `dgd` values correlate at 0.72 with the full-reference
values.
Failure scenario: a user estimates `kin` on a small subset, believes the
gain is exact, and gets a different ranking.
Proposed change: `@details` says the gain is exact *given `kin`*, and that
`kin` should be estimated on the widest dataset available. Docs only.

**F4 [LOW, confidence: high] — duplicate candidates give duplicate rows**
`R/gl.report.ind.add.r:136` — `candidates = c("A", "A")` returns two
identical rows with ranks 1 and 2.
Failure scenario: a candidate list assembled from two sources lists an
animal twice, which inflates the table and misleads the ranks.
Proposed change: drop duplicates with a `verbose >= 2` note.

**F5 [LOW, confidence: high] — convention gaps (FS3, VRB2)**
`:84-86` uses the outdated `build =`, and `:153-157` prints the
`verbose >= 3` summary with raw `cat()`.
Failure scenario: none for results.
Proposed change: drop `build =`, and use `report()`.

Checked, nothing found: the population-name expansion, the error for a
candidate already in the target, missing ids, an invalid target, the
SilicoDArT path, and the input coming back identical.

## Proposed changes

1. `NA`-tolerant GD via a new `na.rm` argument to `utils.kin.dgd`
   (default `FALSE`), with a warning (F1). **Consequence: numerical output
   changes for a `kin` with `NA` in the target block or candidate rows
   (values instead of `NA`). Other `utils.kin.dgd` callers are
   unchanged.**
2. Call-rate warning and a `@details` sentence (F2). No numerical change.
3. `@details` wording on exactness and the reference (F3). Docs only.
4. Drop duplicate candidates (F4). **Consequence: a candidate listed
   twice now appears once.**
5. Standards cleanup (F5). No output change.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: behaviour against roxygen on testset2.gl (one candidate population; all 250 wild candidates), testset2.gs, `NA` kinship, duplicates, invalid inputs — run
- Independent numerical check: `dgd` recomputed from the block-mean formula (maximum difference 0)
- Self-reference (#113 audit claim): verified as unaffected; a subset reference changes values (correlation 0.72) without degenerating them
- Missing data: ranking compared unfiltered against call-rate filtered
- Callers (API3): none outside `@seealso`; `utils.kin.dgd` callers (`gl.report.ind.remove`, `gl.report.ind.move`, `gl.report.mate.suitability`, `gl.select.pairs`) are unaffected by a default-`FALSE` argument
- FBM path (DAT6): SKIPPED — no FBM fixture
- Known complaints: not checked — added in the 2026 kinship series, no release history
- PLT: not applicable — no plot

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence (dgd values instead of NA) approved |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | consequence (duplicate rows removed) approved |
| 5 | approved | Luis | |

## Outcome

- Change 1: `utils.kin.dgd(..., na.rm = FALSE)` added; `gl.report.ind.add` passes `na.rm = TRUE` and warns with the number of `NA` pairs among target and candidates. With one `NA` target pair, `dgd` is computed (was all `NA`); the engine still returns `NA` by default.
- Change 2: call-rate warning over target and candidates (24 on testset2.gl) and `@details` paragraph.
- Change 3: `@details` says the gain is exact given `kin`, and gives the subset-reference comparison.
- Change 4: duplicated candidates dropped with a `verbose >= 2` note.
- Change 5: `build =` dropped; `verbose >= 3` summary through `report()`.
- Snapshot: 2 baseline expectations changed (NA test, duplicate test), mapped to changes 1 and 4; values, order and founders-last unchanged. 18 expectations pass; example runs.
- Full suite: 6 failures in `test-gl.grm.R` and `test-gl.assign.grm.R`, identical on clean `origin/dev` (d6e62f2), so pre-existing and not caused by this change; all `utils.kin.dgd` callers' tests pass.
- NEWS entry added. PR: (pending)

```json
{
  "function": "gl.report.ind.add",
  "package": "dartR.captive",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "d6e62f2",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "DAT", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "DAT", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS3,VRB2", "status": "approved", "change": 5}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "forum/issues: no release history"],
  "status": "in-apply",
  "pr": null
}
```
