# Review: gl.select.pairs (dartR.captive)
- Family mode: report
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 103b47a (origin/dev, after #127; gl.select.pairs unchanged since 73e7cc1)
- Datasets: testset2.gl, EmmacCaptBred (24 individuals: 11 males, 13 females) with `kin` from the full testset2.gl; the same after `gl.filter.callrate(method = "loc", threshold = 0.95)`; testset2.gs; 200 random individuals of testset2.gl for timing; dartR.data 1.2.6
- Baseline: tests/testthat/test-gl.select.pairs.R (snapshot captured pre-review, 26 expectations, all pass)

## Verdict

**Standards: Needs work**: the structure, the read-only behaviour and the
#113 reference guard are in place; `dgd.cum` matches an independent
recomputation and every selected pair respects `f.max` and the
capacities. `NA` inputs crash or give silent `NA` output, invalid
arguments are replaced instead of rejected, and the `dynamic` scheme is
slow for colonies of a few hundred.
**Spec: Needs work**: the three schemes do what the documentation says.
But `f.max` is compared with kinship that missing data shrinks: one of
the 10 pairs selected on unfiltered data has kinship 0.185 after
call-rate filtering, above the 0.125 ceiling, and only 3 of the 10 pairs
survive filtering. The output also does not say that later pairs can
lower gene diversity: on the example, `gd.projected` ends below
`gd.start`.

## Findings

**F1 [HIGH, confidence: high] — low call rates let related pairs through the `f.max` ceiling (DOC5, proposed rule)**
`R/gl.select.pairs.r:231, 245, 268` — `gl.kin` fills missing genotypes
with the locus mean, which pulls kinship toward 0; `f.max` is a fixed
value. On the colony (call rates 0.70-0.80), the default `dynamic` run
selects CB_X_01 x CB_X_02 at kinship 0.110; after filtering loci at call
rate 0.95 their kinship is 0.185. The `static` and `ranked` schemes also
select one such pair. After filtering, the `dynamic` scheme selects 9
pairs, and only 3 of them are among the 10 selected unfiltered.
Failure scenario: a manager breeds the selected pairs on unfiltered data,
including a pair related at about the half-sib level that the function
promised to exclude.
Proposed change: the series call-rate warning at `verbose >= 1` (any
individual below 0.8), a `@details` paragraph with these numbers, and an
example that filters loci on call rate before `gl.kin`.

**F2 [MEDIUM, confidence: high] — one NA kinship crashes one scheme and silently corrupts the others (DAT integrity)**
`R/gl.select.pairs.r:187-188, 231, 245, 268-281` — with one `NA` pair
(CB_CD_01 x CB_Y_03), `dynamic` stops with the misleading
`parent id(s) in add.pairs row 1 not found ...: NA`. `static` and
`ranked` return 11 pairs with `dgd.cum`, `gd.start` and `gd.projected`
all `NA`, and the two individuals of the `NA` pair get `NA` mean kinship,
which `order()` sorts last, so the individual with the lowest mean
kinship in the colony (CB_CD_01, 0.040) is paired last instead of first.
Failure scenario: `kin` from the dominant estimator has one `NA`; the
user either gets an error that names no cause or a pair list whose
ranking is wrong without a message.
Proposed change: ignore `NA` in the mean kinships and gene diversities
(`na.rm = TRUE`), treat a pair with `NA` kinship as infeasible (its
kinship cannot be checked against `f.max`), and warn at `verbose >= 1`
with the number of `NA` values.

**F3 [MEDIUM, confidence: high] — an NA sex crashes the function (DAT integrity)**
`R/gl.select.pairs.r:127-128` — `ids[sex == "Male"]` returns an `NA` id
for an `NA` sex; `dynamic` stops with
`parent id(s) ... not found ...: NA`. The same defect was fixed in #126
and #127.
Failure scenario: an empty sex cell in the metadata stops pair selection
with an error that does not mention sex.
Proposed change: `sex %in% "Male"` and `sex %in% "Female"`, so `NA`
counts as unknown sex, as documented.

**F4 [LOW, confidence: high] — invalid arguments are replaced, not rejected (FS5)**
`R/gl.select.pairs.r:143-184` — an unknown `scheme` (for example
`"Static"`), `max.per.sire`/`max.per.dam` below 1, `f.max` outside
(0, 1] and `n.pairs` below 1 are each replaced by a default with a
warning that is silent at `verbose = 0`. An `NA` in any of them stops
with R's `missing value where TRUE/FALSE needed`.
Failure scenario: `scheme = "Static"` in a quiet script runs the
`dynamic` scheme; the user reports results as static.
Proposed change: `stop(error(...))` for each invalid argument, naming
what is allowed; `NULL` still means the default `n.pairs`, and a
fractional `n.pairs` is still rounded down.
**Consequence: calls with an invalid `scheme`, `max.per.sire`,
`max.per.dam`, `f.max` or `n.pairs` stop instead of running with a
default.**

**F5 [LOW, confidence: high] — the dynamic scheme rebuilds the whole matrix for every candidate (STY2)**
`R/gl.select.pairs.r:276-280` — each round calls `utils.kin.dgd` once per
feasible pair, and each call rebuilds the kinship matrix with every
virtual offspring so far. For 200 individuals (95 males, 94 females), 10
pairs take 22 s; the default `n.pairs` (94 pairs) takes 164 s.
Failure scenario: a manager runs the default scheme on a colony of a few
hundred and waits many minutes, or abandons it for `static`.
Proposed change: keep running sums (total, and each individual's row sum
over the growing matrix) and score each candidate in closed form, as in
#126; selections unchanged.

**F6 [LOW, confidence: high] — later pairs can lower gene diversity without comment (DOC5, proposed rule)**
`R/gl.select.pairs.r:58-63, 84-89` — `dgd.cum` peaks at pair 4 (0.93608)
and falls to 0.93355 at pair 10, below `gd.start` (0.93370): each
further offspring of well-represented parents adds more kinship than it
removes. The function fills `n.pairs` regardless, and the documentation
does not say that `dgd.cum` can fall. On this colony all three schemes
end with the same 10 pairs.
Failure scenario: a user breeds all 10 pairs believing each adds gene
diversity; the last six lower it.
Proposed change: a `@details` sentence saying that `dgd.cum` can fall and
that the row where it peaks is the size beyond which more pairs lower
projected gene diversity. Docs only.

**F7 [LOW, confidence: high] — convention gaps (FS3, DOC2)**
`:107-109` passes the outdated `build =`; the `verbose` text differs from
DOC2.
Failure scenario: none for results.
Proposed change: drop `build =`, use the DOC2 text.

Checked, nothing found: the input comes back identical; `dgd.cum`
equals `utils.kin.dgd` on the selected pairs; the first `dynamic` pair
equals a brute-force search over all feasible pairs; capacities hold with
`max.per.sire = max.per.dam = 2`; `kin = NULL` stops with the #113
message; SilicoDArT runs.

## Proposed changes

1. Call-rate warning, `@details` paragraph and filtered example (F1). No
   numerical change.
2. `na.rm = TRUE` in mean kinship and gene diversity, `NA`-kinship pairs
   infeasible, warning (F2). **Consequence: a `kin` with `NA` now returns
   pairs with gene diversity values in all schemes (it crashed
   `dynamic` and gave `NA` values and a wrong order in `static` and
   `ranked`).**
3. `NA` sex treated as unknown (F3). **Consequence: data with an `NA` sex
   now returns pairs instead of an error.**
4. Invalid arguments are fatal errors (F4). **Consequence: calls with an
   invalid `scheme`, `max.per.sire`, `max.per.dam`, `f.max` or `n.pairs`
   stop instead of running with a default.**
5. Closed-form scoring in the `dynamic` scheme (F5). Selections
   unchanged.
6. `@details` on falling `dgd.cum` (F6). Docs only.
7. Standards cleanup (F7). No output change.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: behaviour against roxygen on EmmacCaptBred (SNP and SilicoDArT), three schemes, capacities 1 and 2, invalid and `NA` arguments, `NA` sex, `NA` kinship, `kin = NULL` — run
- Independent checks: `dgd.cum` recomputed with `utils.kin.dgd`; first `dynamic` pair against brute force over all feasible pairs
- Missing data: selections compared unfiltered against call-rate filtered
- PMx Auto Pair (Manual pp. 87-92): SKIPPED — manual not available
- Callers (API3): none outside `@seealso`; grep of dartR.base, dartR.popgen, dartR.sim, dartR.sexlinked, dartR.spatial and dartr2shiny finds no code callers
- FBM path (DAT6): SKIPPED — no FBM fixture
- Known complaints: not checked — added in the 2026 kinship series, no release history
- PLT: not applicable — no plot
- A concurrent observer session sent overlapping findings (NA sex, NA arguments, `build =`); they are covered by F3, F4 and F7

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | consequence (values in all schemes for NA kinship; NA pair never selected) approved |
| 3 | approved | Luis | consequence (pairs instead of an error for NA sex) approved |
| 4 | approved | Luis | consequence (invalid arguments stop) approved |
| 5 | approved | Luis | |
| 6 | approved | Luis | |
| 7 | approved | Luis | |

## Outcome

- Change 1: warning at `verbose >= 1` (24 individuals below 0.8 unfiltered); `@details` paragraph; the example filters loci at call rate 0.95 and passes the full-dataset `kin`.
- Change 2: `na.rm = TRUE` in mean kinship and every `utils.kin.dgd` call; pairs with `NA` kinship are infeasible in all three schemes; warning with the number of `NA` values. With one `NA` pair, all schemes return gene diversity values, the `NA` pair is not selected, and `static` pairs CB_CD_01 first again.
- Change 3: `%in%` selection; an `NA` sex is excluded with the existing unknown-sex warning.
- Change 4: 11 invalid values across `scheme`, `max.per.sire`, `max.per.dam`, `f.max` and `n.pairs` each stop with `<argument> must be ...`; `n.pairs = 3.9` still gives 3 pairs, `Inf` capacities still accepted.
- Change 5: `dynamic` scores candidates from running sums (total and row sums of the growing matrix) when `kin` has no `NA`, otherwise through `utils.kin.dgd` with `na.rm = TRUE`. 200 individuals, default `n.pairs`: same 94 pairs, `dgd.cum` difference 0, 0.68 s (was 164 s). A new test checks the third pick against brute force.
- Change 6: `@details` on falling `dgd.cum`.
- Change 7: `build =` dropped, DOC2 `verbose` text.
- Snapshot: selections, `f.off`, `dgd.cum`, `gd.start` and `gd.projected` unchanged for all schemes. Baseline expectations changed only for invalid arguments (4), `NA` sex (3) and `NA` kinship (2). 51 expectations pass; full suite 583 pass, 0 fail; example runs.
- NEWS entry added. PR: #128

```json
{
  "function": "gl.select.pairs",
  "package": "dartR.captive",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "103b47a",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "MEDIUM", "confidence": "high", "rule": "DAT", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "STY2", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "FS3,DOC2", "status": "approved", "change": 7}
  ],
  "coverage_skipped": ["PMx manual: not available", "DAT6: no FBM fixture", "forum/issues: no release history"],
  "status": "pr-open",
  "pr": 128
}
```
