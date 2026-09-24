# Review: gl.report.kin.classes (dartR.captive)
- Family mode: report
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 463361e (origin/dev, after #112; PR #113 open, does not touch this function)
- Datasets: testset2.gl (274 x 755 SNPs, 31 populations, 24 captive-bred offspring in EmmacCaptBred with sire/dam columns), testset2.gs; dartR.data 1.2.5
- Baseline: tests/testthat/test-gl.report.kin.classes.R (snapshot captured pre-review, 20 expectations, all pass)

## Verdict

**Standards: Needs work**: the structure, the read-only behaviour (input
comes back identical, no history append) and the scale handling through
`utils.kin.check` all conform. There are three small convention gaps.
**Spec: Needs work**: the single median baseline cannot absorb population
structure, and the documentation does not say so. On the packaged example,
all 127 full-sib calls fall within populations, and 19 of the 48 recorded
parent-offspring links are called second-degree. A single missing kinship
value also turns every pair's class into `NA`.

## Findings

**F1 [HIGH, confidence: high] — population structure drives the classes, undocumented (DOC5, proposed rule)**
`R/gl.report.kin.classes.r:130,141` — kinship is computed with allele
frequencies pooled over the whole sample, then shifted by one global
median. In a structured sample, pairs within a population sit above the
pooled baseline and pairs between populations sit below it. One shift
cannot correct both. On `testset2.gl` (31 drainages, about 10 turtles each):
- all 127 full-sib calls and 977 of the 1000 second-degree calls are
  within-population pairs. The pedigree records only 46 full-sib pairs, and
  23 of those are called second-degree.
- 19 of 48 recorded parent-offspring links are called second-degree
  (`kinship.adj` 0.11-0.19). Every one is a cross-drainage link: a wild
  founder from one drainage with its captive offspring, or offspring of a
  cross between captive lines founded from different drainages.
  Within-lineage links classify correctly.

The `@details` warns only about a sample dominated by one family. The
example comment calls the conflicts "a few near the class boundary", but
they are 40% of the recorded links, and structure, not the boundary,
explains them.
Failure scenario: a manager runs the function on a multi-population
dataset, reads the pedigree conflicts as pedigree or sample errors, and
"corrects" a pedigree that is right.
Proposed change: document the structure bias in `@details`. Print a
`verbose >= 1` warning when `x` has more than one population, recommending
classification within one population or a panmictic sample. Rewrite the
example comment so it names the cause.

**F2 [HIGH, confidence: high] — one NA kinship makes every class NA (DAT integrity, no rule fits exactly)**
`R/gl.report.kin.classes.r:130` — `stats::median()` is called without
`na.rm`, so a single `NA` off-diagonal value makes `baseline` `NA`. Every
`kinship.adj` and every `class` then becomes `NA`, and the
`class != "unrelated"` filter keeps all 37,401 rows as `NA` rows.
Failure scenario: `kin` from an estimator that returns `NA` for an
individual with too few called loci. The user gets a full-size table with
no classes and no message.
Proposed change: take the median with `na.rm = TRUE`. Pairs with `NA`
kinship get class `NA` and are excluded from `pairs` unless
`all.pairs = TRUE`. Print a `verbose >= 1` warning with the count. This
follows the approach approved for `gl.report.kinship` in #112.

**F3 [LOW, confidence: medium] — few individuals make the baseline meaningless (DOC5, proposed rule)**
`R/gl.report.kin.classes.r:121,130` — with two individuals, the median
off-diagonal value is the pair's own kinship, so `kinship.adj` is 0 and the
pair is always "unrelated". With a handful of individuals, the median is set
by whichever pairs happen to be related.
Failure scenario: `gl.report.kin.classes(x[1:2, ])` on a known
parent-offspring pair returns zero related pairs and no message.
Proposed change: require at least three individuals. Warn at
`verbose >= 1` when there are fewer than 10 that the median baseline is
unreliable. The cutoff of 10 is a judgement call, not a derived value.

**F4 [LOW, confidence: high] — convention gaps (FS3, VRB2, DAT6 proposed rule)**
- `:106-108` `utils.flag.start(build = "v.2026.1")` uses the outdated
  `build =` argument. The functions reviewed earlier in this series have
  dropped it.
- `:230-240` the `verbose >= 3` summary prints with raw `cat()`, without
  `report()`.
- `:161` `as.matrix(x)` densifies every individual to split what may be a
  handful of first-degree pairs. It only needs the rows of the individuals
  in those pairs.

Failure scenario: none for results. The first two are inconsistencies
with the house style. The third costs memory on large datasets.
Proposed change: drop `build =`, wrap the summary in `report()`, and
densify only the rows needed. Output is unchanged.

**F5 [INFO, confidence: high] — `conflicts` is NULL rather than an empty data frame (consistency)**
`R/gl.report.kin.classes.r:191,258` — the documentation says `NULL`, and
the code does that. #96 changed `gl.report.parent.offspring` to return an
empty data frame so callers can use `nrow()` without a `NULL` check.
Failure scenario: `nrow(res$conflicts)` returns `NULL`, not 0, and
`if (nrow(res$conflicts) > 0)` errors.
Proposed change: return a zero-row data frame with the documented columns
when there are no conflicts or no sire/dam columns.

Checked, nothing found: the PO/FS split by opposite homozygotes separates
cleanly on this dataset (parent-offspring maximum rate 0.0034, full-sib
minimum 0.0050, threshold 0.005). `gl.grm` input is halved through the
scale tag and gives the same classes as `gl.kin`. Factor `sire`/`dam`
columns and a self-listed parent do not break the cross-check.

## Proposed changes

1. Document the population-structure bias in `@details`, add a
   `verbose >= 1` warning when `nPop(x) > 1`, and rewrite the example
   comment on conflicts (F1). No numerical change.
2. Make `NA` kinship values local: median with `na.rm = TRUE`, `NA` pairs
   classed `NA`, excluded unless `all.pairs = TRUE`, and a warning with the
   count (F2). **Consequence: numerical output changes for any `kin`
   containing `NA`. Previously every pair was `NA`; now only the affected
   pairs are.**
3. Require at least three individuals, and warn at `verbose >= 1` when
   there are fewer than 10 (F3). **Consequence: a two-individual call that
   returned an empty table now errors.**
4. Standards cleanup: drop `build =`, use `report()` in the `verbose >= 3`
   summary, and densify only the needed rows (F4). No output change.
5. Return a zero-row data frame instead of `NULL` for `conflicts` (F5).
   **Consequence: code testing `is.null(res$conflicts)` changes behaviour.**

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: behaviour against roxygen on testset2.gl (SNP), testset2.gs (SilicoDArT), gl.grm input, NA kinship, two-individual subset, factor and self-listed sire/dam — run
- Pedigree accuracy against recorded sire/dam links in testset2.gl — run
- Callers (API3): none in dartR.captive; none in dartR.base, dartR.popgen or dartr2shiny sources (dartR.base mentions it only in its manifest)
- Known complaints (Google Group / GitHub issues): not checked — the function was added in the 2026 kinship series and has no release history
- FBM path (DAT6): SKIPPED — no FBM fixture for dartR.captive
- PLT: not applicable — no plot
- Independent numerical check: the class thresholds and the baseline were recomputed by hand in the baseline test; the PO/FS rates were recomputed independently

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | consequence (output changes for NA-containing kin) approved |
| 3 | approved | Luis | consequence (two-individual call errors) approved |
| 4 | approved | Luis | |
| 5 | approved | Luis | consequence (`conflicts` no longer NULL) approved |

## Outcome

- Change 1: `@details` paragraph on structure, `verbose >= 1` warning when `nPop(x) > 1`, example comments rewritten. Test "population-structure warning" passes (warns on testset2.gl, silent on the EmmacCaptBred subset). Class counts unchanged (127 / 30 / 1000 / 3233).
- Change 2: median with `na.rm = TRUE`, `NA` pairs excluded unless `all.pairs = TRUE`, warning with the count. With one NA pair: 1 pair classed NA (was all 37,401), the rest classified.
- Change 3: fewer than 3 individuals is an error; fewer than 10 warns.
- Change 4: `build =` dropped, `verbose >= 3` summary through `report()`, only the individuals in first-degree pairs densified. Output identical (baseline class counts and 19 conflicts unchanged).
- Change 5: `conflicts` is a zero-row data frame with the documented columns.
- Snapshot: 3 baseline expectations changed, all mapped to changes 2 and 3; tests updated, 28 expectations pass. Example runs; `verbose = 3` run end to end on testset2.gl.
- `devtools::document()` also dropped two stale `importFrom(stats, dnorm/qnorm)` lines from NAMESPACE that no roxygen tag produces; reverted to keep this PR to one function.
- NEWS entry added. PR: (pending)

```json
{
  "function": "gl.report.kin.classes",
  "package": "dartR.captive",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "463361e",
  "verdict_standards": "needs_work",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "medium", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "FS3,VRB2,DAT6", "status": "approved", "change": 4},
    {"id": "F5", "severity": "INFO", "confidence": "high", "rule": "API1", "status": "approved", "change": 5}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "forum/issues: function has no release history"],
  "status": "in-apply",
  "pr": null
}
```
