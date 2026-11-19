# Review: gl.report.ind.move (dartR.captive)
- Family mode: report
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: fdb8967 (origin/dev, after #121)
- Datasets: testset2.gl (EmmacCaptBred, EmmacMaclGeor and EmmacBurnBara: 46 individuals, 92 moves; all 274 individuals, 8220 moves), the same after `gl.filter.callrate(method = "loc", threshold = 0.95)`, testset2.gs; dartR.data 1.2.6
- Baseline: tests/testthat/test-gl.report.ind.move.R (snapshot captured pre-review, 30 expectations, all pass)

## Verdict

**Standards: Needs work**: the structure and read-only behaviour conform,
and `dgd.source`/`dgd.dest` match an independent recomputation exactly.
All 8220 moves on the full dataset take 0.2 s. There are small convention
gaps.
**Spec: Needs work**: the calculation is right, and the #113 audit's
"not affected by self-referenced kinship" holds (the blocks are
populations, and at least two are required). One `NA` kinship silently
voids every move into or out of a population. Low call rates reshuffle
the ranking. The documentation does not say what `net` rewards: moves of
genetically divergent individuals into small populations. On testset2.gl,
9 of the top 10 moves send captive-bred animals of mixed ancestry into
wild drainages.

## Findings

**F1 [MEDIUM, confidence: high] — one NA inside a population voids all its moves (DAT integrity)**
`R/gl.report.ind.move.r:116-139` — `utils.kin.dgd` without `na.rm`. With
one `NA` pair inside EmmacCaptBred, 70 of 92 moves have `net = NA`: every
move out of the population and every move into it. They are sorted to
the bottom without a message.
Failure scenario: `kin` has an `NA` for one poorly genotyped animal, and
the whole colony drops out of the transfer table.
Proposed change: pass `na.rm = TRUE` to `utils.kin.dgd` (the argument
added in #120), with a `verbose >= 1` warning giving the number of `NA`
pairs.

**F2 [MEDIUM, confidence: high] — low call rates reshuffle the ranking (DOC5, proposed rule)**
`R/gl.report.ind.move.r:95` — mean imputation in `gl.kin` pulls kinship
and self-kinship toward 0. The captive-bred animals (call rates
0.70-0.80) have mean self-kinship 0.34, against 0.46 after filtering loci
at call rate 0.95. `dgd.dest` correlates with the mover's self-kinship
(r = -0.59). After filtering, only 5 of the top 10 moves stay in the top
10.
Failure scenario: a manager acts on the top moves from unfiltered data,
and half of them are artefacts of missing genotypes.
Proposed change: a `verbose >= 1` call-rate warning and a `@details`
sentence, as in #115, #116, #119 and #120.

**F3 [LOW, confidence: high] — what `net` rewards is undocumented (DOC5, proposed rule)**
`R/gl.report.ind.move.r:153` — `net` is the unweighted sum of two
populations' GD changes, and it is not the change in GD across all
populations combined. Adding one individual to a population of `n`
changes its GD by roughly `1/n`, so moves into small populations score
high. An individual with low kinship to the destination scores high too.
With pooled allele frequencies, that means an individual from a different
lineage, or an admixed one. On the example, 9 of the top 10 moves (all 10
after call-rate filtering) send captive-bred animals descended from
several drainages into a single wild drainage.
Failure scenario: a user takes the top-ranked move as a management
recommendation and translocates admixed captive stock into a wild
population. GD rises, but genetic integrity, which the metric does not
measure, is lost.
Proposed change: a `@details` paragraph stating that `net` is a sum of
per-population changes, that it favours small destinations and divergent
movers, and that it measures gene diversity only. The example should also
estimate `kin` on the full dataset, which the function family recommends
(the effect here is small: `net` correlates at 0.996).

**F4 [LOW, confidence: high] — convention gaps (FS3, VRB2)**
`:84-86` uses the outdated `build =`, and `:159-163` prints the
`verbose >= 3` summary with raw `cat()`.
Failure scenario: none for results.
Proposed change: drop `build =`, and use `report()`.

Checked, nothing found: the input comes back identical, singleton
populations get `dgd.source = NA` with a warning (as documented), fewer
than two populations is an error, and the SilicoDArT path runs.

## Proposed changes

1. `na.rm = TRUE` in the GD calls, with a warning (F1). **Consequence:
   numerical output changes for a `kin` with `NA` (moves into and out of
   the affected population get values instead of `NA`).**
2. Call-rate warning and `@details` sentence (F2). No numerical change.
3. `@details` on what `net` measures and rewards; the example uses
   full-dataset `kin` (F3). Docs only.
4. Standards cleanup (F4). No output change.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: behaviour against roxygen on a 3-population subset and the full testset2.gl, testset2.gs, `NA` kinship, a single population — run
- Independent numerical check: `dgd.source` and `dgd.dest` recomputed from block means (maximum difference 0)
- Self-reference (#113 audit claim): verified as unaffected; subset-reference and full-reference `net` correlate at 0.996
- Missing data: ranking compared unfiltered against call-rate filtered
- Callers (API3): none outside `@seealso`
- FBM path (DAT6): SKIPPED — no FBM fixture
- Known complaints: not checked — added in the 2026 kinship series, no release history
- PLT: not applicable — no plot

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence (values instead of NA for NA-containing kin) approved |
| 2 | approved | Luis | |
| 3 | approved | Luis | |
| 4 | approved | Luis | |

## Outcome

- Change 1: `na.rm = TRUE` in all three `utils.kin.dgd` calls; warning with the number of `NA` pairs. With one `NA` pair in EmmacCaptBred, 0 of 92 moves are `NA` (was 70).
- Change 2: call-rate warning (24 individuals in the example) and `@details` paragraph.
- Change 3: `@details` paragraph on what `net` measures and rewards; example estimates `kin` on the full testset2.gl.
- Change 4: `build =` dropped; `verbose >= 3` summary through `report()`.
- Snapshot: 1 baseline expectation changed (the `NA` count), mapped to change 1; values, order and "9 of top 10 from EmmacCaptBred" unchanged. 32 expectations pass; example runs.
- NEWS entry added. PR: #122

```json
{
  "function": "gl.report.ind.move",
  "package": "dartR.captive",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "fdb8967",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "MEDIUM", "confidence": "high", "rule": "DAT", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "FS3,VRB2", "status": "approved", "change": 4}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "forum/issues: no release history"],
  "status": "pr-open",
  "pr": 122
}
```
