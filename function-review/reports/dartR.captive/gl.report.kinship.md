# Review: gl.report.kinship (dartR.captive)
- Family mode: report
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 4510d67 (origin/dev, includes the gl.kin fixes #108 and the scale contract #110)
- Datasets: testset2.gl (274 x 755 SNPs, 30 populations, 24 captive-bred offspring in EmmacCaptBred), the same after `gl.filter.callrate(threshold = 0.95)`, testset2.gs (242 x 755); dartR.data 1.2.5
- Baseline: tests/testthat/test-gl.report.kinship.R (snapshot captured pre-review, 24 expectations)

## Verdict

**Standards: Ready** — house structure, read-only behaviour (no history
append, input untouched), results identical at every verbosity and
independent of the plot (PLT3); two small input-handling gaps.
**Spec: Needs work** — the per-individual MK and its rank are taken over
the whole dataset, so the "most valuable" captive animal is chosen by its
kinship to 250 wild individuals rather than to the captive group, and the
"overall" row reports GD = 1 and FGE = 2.35e17, both fixed by how
genomic kinship is centred rather than by the data.

## Findings

**F1 [HIGH, confidence: high] — MK and MKrank are not relative to the
managed population (DOC5, proposed rule)**
`R/gl.report.kinship.r:144,155-156` — `MK <- rowMeans(kin)` averages
each individual's kinship over every individual in the dataset, and
`MKrank` ranks within sex across the whole dataset. PMx defines MK as the
mean kinship of an individual with the living managed population,
including itself, and the rank is used to choose breeders within that
population. The per-population table already uses within-population
blocks, which is why `@details` has to warn that the two tables do not
agree.
Failure scenario: testset2.gl, the 24 captive-bred individuals. Their MK
over the whole dataset and their MK within EmmacCaptBred have Spearman
correlation 0.27; the rank-1 (most valuable) male and female both differ
between the two, and 9 of 11 male and 12 of 13 female ranks change. The
example in `@examples` reads the captive-bred rows from the whole-dataset
table.
Proposed change: compute `MK` within each individual's population (row
means of its population block, including self) and rank within
population and sex; keep the dataset-wide value as a new column
`MK.all`. Update `@details`, `@return` and the example.
**Consequence: `MK` and `MKrank` values change for every call with more
than one population; `$ind` gains a column `MK.all`.**

**F2 [HIGH, confidence: high] — the "overall" row reports GD and FGE
fixed by the centring (DOC5, proposed rule)**
`R/gl.report.kinship.r:175-181,189` — genomic kinship is centred on the
supplied individuals, so the mean of the full matrix is 0 (for `grm`,
exactly, up to rounding). The overall row then reports `GD = 1` and
`FGE = 1 / (2 * meanMK)` of a rounding residual.
Failure scenario: testset2.gl: overall `meanMK = 2.1e-18`, `GD = 1`,
`FGE = 2.35e17`; after `gl.filter.callrate(0.95)`: `FGE = 1.86e17`.
testset2.gs (dominant): overall `meanMK = -0.00018`, `GD = 1.0002`,
`FGE = NA`, and every call prints a warning blaming "one or more
populations" for what is the overall row. Per-population rows are not
affected: their blocks are subsets, with mean kinship clearly above 0.
Proposed change: report `GD` and `FGE` as `NA` in the overall row
(keep `n`, `meanMK` and `meanF`), and document why; treat block means
within 1e-8 of 0 as 0 so FGE is NA rather than astronomically large;
warn only for population rows.
**Consequence: overall `GD` and `FGE` become NA; the spurious warning
for dominant data disappears.**

**F3 [LOW, confidence: high] — two "unknown" sex groups (DOC5, proposed
rule)**
`R/gl.report.kinship.r:154` — NA and empty sexes are grouped as
`"unknown"`, but the literal value `"Unknown"` (as in testset2.gl) is a
group of its own.
Failure scenario: testset2.gl with 3 sexes set to NA: those 3 are ranked
1-3 in one group and the 16 "Unknown" individuals 1-16 in another, so
two individuals share rank 1 among animals of unknown sex.
Proposed change: treat NA, empty and case-insensitive "unknown" as one
group.

**F4 [LOW, confidence: medium] — one missing kinship value empties MK
and the overall row (DOC5, proposed rule)**
`R/gl.report.kinship.r:144,175` — `rowMeans(kin)` and `mean(sub)` have
no `na.rm`. `utils.kin.dominant` produces NA for pairs sharing no scored
loci.
Failure scenario: a supplied kinship matrix with one NA pair: MK is NA
for both individuals and the overall `meanMK` is NA, without a message.
Proposed change: average with `na.rm = TRUE` and warn at
`verbose >= 1` with the number of NA pairs.

## Addendum (found while applying)

**Correction to F1.** F1 understated the defect. For method `grm` (the
SNP default) every row of the kinship matrix sums to 0 by construction:
`rrBLUP::A.mat` centres each locus on the sample allele frequency, so
`G %*% 1 = 0`. The dataset-wide MK was therefore 0 for every individual
up to rounding (testset2.gl: -1.7e-17 to 2.6e-17), and the pre-review
`MKrank` ranked floating-point noise. The Spearman correlation of 0.27
quoted in F1 is a correlation with that noise. With `dominant` kinship
the row means are near 0 (-0.0061 to 0.0004); with `emibd9`, which is
not centred on the sample, they vary (0.064 to 0.137 on a 37-individual
subset).

**A1 [LOW, confidence: high] — the approved `MK.all` column is 0 by
construction under the default method (DOC5, proposed rule)**
Change 1 kept the dataset-wide value as `MK.all` "for reference". Under
`grm` it holds rounding residuals (about 1e-17) for every individual,
and near-zero values under `dominant`, so it invites the same misreading
that F1 fixes.
Proposed change: drop `MK.all`; document in `@details` why a
dataset-wide MK carries no information for sample-centred kinship.
**Consequence: `$ind` keeps its original six columns; no new column.**

## Proposed changes

1. MK and MKrank within each individual's population; dataset-wide MK
   kept as `MK.all`; docs and example updated (F1).
   **Consequence: `MK` and `MKrank` values change for every call with
   more than one population; `$ind` gains a column `MK.all`.**
2. Overall row `GD` and `FGE` as NA, documented; block means within 1e-8
   of 0 give FGE NA; warning only for population rows (F2).
   **Consequence: overall `GD` and `FGE` become NA; the spurious
   warning for dominant data disappears.**
3. One "unknown" sex group for NA, empty and "unknown" in any case (F3).
4. `na.rm = TRUE` for MK and block means, with a warning counting NA
   pairs (F4).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. Notes
  without findings: `utils.flag.start(build = ...)` outdated argument
  (FS3), still accepted; `plot.dir` resolved before the start flag (FS
  order, harmless); `exists("p1")` is always TRUE.
- Report family: input object untouched and no history append (FS8) —
  run; results identical at `verbose = 0` and `verbose = 3`, plot
  suppressed or not (PLT3) — run.
- Spec: behaviour vs roxygen on testset2.gl, testset2.gl filtered,
  testset2.gs — run. `gl.grm()` output passed as `kin` gives the same
  MK as the default (scale contract, #110) — run.
- Noted, not a finding: per-population GD depends on the reference set.
  EmmacCaptBred (full sibs from crosses between localities) has GD 0.93
  (0.90 after filtering) against 0.77-0.82 for wild localities, because
  kinship is relative to the pooled allele frequencies of 30 localities.
  `gl.kin` documents this caveat.
- FBM path (DAT6): SKIPPED — the function works on the kinship matrix,
  not genotypes; `gl.kin` densifies when kin is computed internally.
- Callers: none in dartr2shiny or elsewhere in the package.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis Mijangos | |
| 2 | approved | Luis Mijangos | |
| 3 | approved | Luis Mijangos | |
| 4 | approved | Luis Mijangos | |
| A1 | approved | Luis Mijangos | drop MK.all |

## Outcome

- Change 1 (F1): `MK` equals the row means of each population's block;
  ranks run 1..n within every population x sex group; the captive-group
  rank-1 male is its lowest within-group MK; population `meanMK` equals
  the mean of that population's `MK`.
- A1: `MK.all` not added; `@details` explains that dataset-wide row means
  are 0 by construction (max |row mean| < 1e-12 on testset2.gl).
- Change 2 (F2): overall `GD` and `FGE` NA; population rows unchanged
  (EmmacCaptBred GD 0.9336865, FGE 7.539946); no FGE warning for
  testset2.gs.
- Change 3 (F3): NA and "Unknown" sexes ranked as one group.
- Change 4 (F4): one NA pair gives no NA MK or overall meanMK, and a
  warning at `verbose = 1`.
- Snapshot: the pre-review baseline (24 expectations) run on changes 1-4
  gave 11 failures, each mapped: 4 to change 1 (new column, MK values,
  rank range), 3 to change 2, 2 to change 3, 2 to change 4. Per-population
  GD/FGE and F unchanged.
- Specification tests: tests/testthat/test-gl.report.kinship.R, 30
  expectations, 0 failures. `gl.report.kinship(testset2.gl, verbose = 3)`
  runs end to end; the example runs.
- PR: pending.

```json
{
  "function": "gl.report.kinship",
  "package": "dartR.captive",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "4510d67",
  "verdict_standards": "ready",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "medium", "rule": "DOC5", "status": "approved", "change": 4},
    {"id": "A1", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": "A1"}
  ],
  "coverage_skipped": ["DAT6: operates on the kinship matrix"],
  "status": "approved-applied",
  "pr": null
}
```
