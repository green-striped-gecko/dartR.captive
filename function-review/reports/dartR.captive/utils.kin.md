# Review: utils.kin.check, utils.kin.dgd, utils.kin.dominant (dartR.captive)

- Family mode: analysis (internal helpers of the captive management series), reviewed together: `utils.kin.check` validates the kinship matrix every consumer receives, `utils.kin.dominant` estimates it for presence/absence data, and `utils.kin.dgd` scores gene diversity from it
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 586d93e (origin/dev)
- Datasets: testset.gl, testset.gs (subsets); a 6 × 6 synthetic kinship matrix for `utils.kin.dgd`
- Baseline: `tests/testthat/test-utils.kin.dgd.R` (10 expectations) and `test-utils.kin.dominant.R` (5), captured pre-review; existing `test-utils.kin.check.R` (27). All pass.
- Author and custodian: Arthur Georges
- Callers: `utils.kin.check` in 11 exported functions; `utils.kin.dgd` in `gl.report.ind.add`, `gl.report.ind.move`, `gl.report.ind.remove`, `gl.report.mate.suitability`, `gl.select.pairs`; `utils.kin.dominant` in `gl.kin`

## Verdict

**Standards: Ready** — documented contracts, argument validation with clear fatal errors, matrix algebra instead of loops; two cosmetic fixes.
**Spec: Ready** — each helper reproduces an independent computation exactly; one latent edge case no caller reaches today.

## Findings

**F1 [LOW, confidence: high] — `utils.kin.dgd` errors on an empty `add.pairs` (spec)**
`R/utils.kin.dgd.r:116-119` — with zero rows, `paste0("offspring_", seq_len(0))` returns `"offspring_"` (length 1, not 0), so the id vector is one longer than the matrix and `matrix(..., dimnames = )` fails: "length of 'dimnames' [1] not equal to array extent".
Failure scenario: none today; `gl.select.pairs` and `gl.report.mate.suitability` always pass at least one pair. A caller that scores an empty selection (for example, GD before the first pairing) would stop instead of getting the GD of the unmodified matrix.
Proposed change: treat a zero-row `add.pairs` like `NULL` (no virtual offspring).

**F2 [LOW, confidence: high] — `utils.kin.dominant` prints its own start and end lines inside `gl.kin` (VRB1)**
`R/utils.kin.dominant.r:71-75, 141-144` — unlike the other two helpers, it runs `utils.flag.start()` and prints "Completed:", so `gl.kin(<SilicoDArT>, verbose = 1)` shows "Starting utils.kin.dominant" / "Completed: utils.kin.dominant" nested inside its own start and end lines.
Failure scenario: a user sees an internal helper announced as if they had called it.
Proposed change: remove the flag start and end lines; keep the verbosity resolution, the datatype check and the verbose >= 2/3 progress and summary lines.

## Proposed changes

1. `utils.kin.dgd`: an empty `add.pairs` adds no offspring (F1).
2. `utils.kin.dominant`: drop the nested Starting/Completed lines (F2).

## Coverage

- `utils.kin.dominant` vs a pairwise loop over shared informative loci (40 individuals of testset.gs): maximum difference 2.8e-16; symmetric, diagonal 0.5, dimnames = `indNames(x)`; SNP input refused
- `utils.kin.dgd` vs explicit construction of the extended matrix, including an offspring of two virtual offspring and a selfed offspring: identical; removals, unknown ids, dropping all, a dropped parent, `data.frame` pairs, numeric-looking ids, NA with `na.rm` both ways — run
- `utils.kin.check`: pass-through, restriction of a larger matrix, reordered matrix of the same individuals, relatedness halving, `need.reference` with superset (accepted), self (rejected) and reordered self (rejected), missing individual, `data.frame` input (rejected with message), `kin = NULL` on SilicoDArT (computed with method `dominant`) — run; plus the existing 27-expectation test file
- Standards walk: DOC, VRB, DAT, STY — run. FS structure: `utils.kin.check` and `utils.kin.dgd` are deliberately lightweight (documented).
- FBM path: not applicable to `utils.kin.dominant` (`gl.gen2fbm` accepts SNP data only); the other two take matrices.

## Report notes (not findings)

- `gl.report.kin.confidence` repeats the dominant estimator inline for its locus bootstrap instead of calling `utils.kin.dominant`, because it resamples loci with band frequencies fixed from the original data. The two copies must be kept in step if the estimator changes.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | |

## Outcome

- Changes 1 and 2 applied on branch `review-utils.kin` (from `origin/dev`).
- Characterization tests: `utils.kin.dgd` 10, `utils.kin.dominant` 6 expectations pass. Diffs from baseline map to approved changes only: an empty `add.pairs` returns the GD of the unmodified matrix (1); `gl.kin(verbose = 1)` on SilicoDArT no longer mentions `utils.kin.dominant` (2).
- Full dartR.captive test suite: 29 files, 0 failures, 0 skipped (including `gl.kin` 29, `gl.select.pairs` 51, `gl.report.mate.suitability` 22, `utils.kin.check` 27).
- Roxygen unchanged; NEWS entry added.
- Full `R CMD check` not run.
- PR: pending.

## Machine block

```json
{
  "function": "utils.kin.check, utils.kin.dgd, utils.kin.dominant",
  "package": "dartR.captive",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "586d93e",
  "verdict_standards": "ready",
  "verdict_spec": "ready",
  "findings": [
    {"id": "F1", "severity": "LOW", "confidence": "high", "rule": "spec", "status": "approved", "change": 1},
    {"id": "F2", "severity": "LOW", "confidence": "high", "rule": "VRB1", "status": "approved", "change": 2}
  ],
  "coverage_skipped": ["FBM: not applicable"],
  "status": "pr-open",
  "pr": null
}
```
