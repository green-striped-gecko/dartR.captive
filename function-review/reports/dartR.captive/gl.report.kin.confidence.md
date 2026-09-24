# Review: gl.report.kin.confidence (dartR.captive)
- Family mode: report
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: cb8d955 (origin/dev, after #113)
- Datasets: testset2.gl (274 x 755 SNPs, individual call rate 0.70-0.89), testset2.gs; dartR.data 1.2.5
- Baseline: tests/testthat/test-gl.report.kin.confidence.R (snapshot captured pre-review, seeded, 12 expectations, all pass)

## Verdict

**Standards: Needs work**: the structure and read-only behaviour conform,
seeded runs are reproducible, and the internal SNP relationship matrix
closely matches `gl.kin` on testset2.gl (off-diagonal difference at most
8e-5 against an SD of 0.037; diagonal at most 0.007). There are small
convention gaps.
**Spec: Rework**: the standard errors do not describe the estimator users
get from `gl.kin`. The SNP bootstrap keeps the centring that #108 removed
from `gl.kin`, which inflates the standard errors by 17%. The SilicoDArT
bootstrap uses a different estimator altogether. Individuals with more
missing data get smaller standard errors, which is the reverse of the
trust the function claims to measure.

## Findings

**F1 [HIGH, confidence: high] — SNP bootstrap uses centring that gl.kin no longer applies (DOC5, proposed rule)**
`R/gl.report.kin.confidence.r:158-164` — each resample computes
`G/2 - MS` with `MS = mean(diag(G) - 1)`, and the comment calls it
"copied exactly from gl.kin's grm method". Since #108, `gl.kin` returns
`G/2` with no `MS` shift. On testset2.gl the bootstrap distribution sits
0.017 above `gl.kin`'s values, and because `MS` varies between resamples it
adds variance. The median off-diagonal SE is 0.0258 with the shift and
0.0220 for `G/2`, so the reported SE is 17% too large.
Failure scenario: a user compares the kinship difference between two
candidate pairs against 2 SE and concludes they are indistinguishable when
the correct SE would separate them.
Proposed change: bootstrap `G/2`, the estimator `gl.kin` returns. Update
`@details` to match.

**F2 [HIGH, confidence: high] — SilicoDArT bootstrap is a different estimator from gl.kin's (DOC5, proposed rule)**
`R/gl.report.kin.confidence.r:172-178` — each resample computes the
Pearson correlation of two individuals' band profiles, divided by 2.
`gl.kin` (method `dominant`, via `utils.kin.dominant`) uses band
covariance standardised by pooled band frequencies, pairwise-complete,
with monomorphic loci removed. On testset2.gs the two estimators correlate
at 0.75. The off-diagonal means are 0.144 and -0.002, and their spreads
differ by a factor of 1.22.
Failure scenario: every SilicoDArT SE describes a quantity the user never
sees, so the "pairs resolvable" figure and any interval built from it are
wrong.
Proposed change: bootstrap the `utils.kin.dominant` estimator, with band
frequencies from the original data, monomorphic loci dropped, and
pairwise-complete denominators.

**F3 [HIGH, confidence: high] — low call rate gives smaller, not larger, SE (DOC5, proposed rule)**
`R/gl.report.kin.confidence.r:147` — missing genotypes are set to 0 after
centring, the mean imputation `gl.grm` also uses. An individual with many
missing calls contributes 0 at those loci in every resample, so its
kinships shrink toward 0 and vary less. On testset2.gl the correlation
between individual call rate and mean SE is +0.76. With one individual
made 90% missing, its mean SE drops to 0.014, against 0.025 for everyone
else.
Failure scenario: the least reliable individual in a colony looks like
the most precisely estimated one. The function exists to weight kinships
by how much they can be trusted (the PMx rationale in `@details`), so this
works against its purpose.
Proposed change: keep the estimator matched to `gl.kin`, and document that
the bootstrap cannot see the bias that mean imputation introduces. Warn
at `verbose >= 1` when any individual's call rate is below 0.8,
recommending `gl.filter.callrate(method = "ind")` first. The 0.8 cutoff is
a judgement call.

**F4 [MEDIUM, confidence: high] — SEs are always for grm/dominant, whatever kin is (DOC5, proposed rule)**
`R/gl.report.kin.confidence.r:217-218` — with `kin` from
`gl.kin(method = "emibd9")`, the `verbose >= 3` "pairs resolvable" line
compares EMIBD9 estimates with SEs of the grm estimator. The two
estimators differ in scale and dispersion.
Failure scenario: an EMIBD9 user reads the resolvable percentage as
applying to their matrix.
Proposed change: when `attr(kin, "method")` is not `grm` (SNP) or
`dominant` (SilicoDArT), warn at `verbose >= 1` that the SEs describe that
estimator, not the supplied one. Say so in `@param kin`.

**F5 [LOW, confidence: high] — convention gaps and a fragile diagnostic (FS3, VRB2, STY1)**
- `:98-100` uses the outdated `utils.flag.start(build = ...)`.
- `:210-220` the `verbose >= 3` block prints with raw `cat()`.
- `:139-140` the warning string has a literal line break inside it.
- `:217` `median(kin)` has no `na.rm`, so one `NA` in `kin` makes the
  resolvable percentage `NaN` (the same defect as `gl.report.kin.classes`
  F2).
- `:106` when `kin` is `NULL`, `gl.kin` runs even at `verbose < 3`, where
  `kin` is never used. That is cheap for grm, but it is wasted work.

Failure scenario: cosmetic, except the `NaN` percentage.
Proposed change: drop `build =`, route the summary through `report()`,
fix the string, use `na.rm = TRUE`, and compute `kin` only at
`verbose >= 3`. Returned values are unchanged.

Checked, nothing found: the input comes back identical, there is no
history append, seeded runs are identical, and `nboots`/`conf` validation
works. The memory cost (`nInd^2 x nboots` doubles) is documented.

## Proposed changes

1. SNP bootstrap on `G/2`, matching `gl.kin`; update `@details` (F1).
   **Consequence: numerical output changes; SNP standard errors drop by
   about 15% (median 0.0258 to 0.0220 on testset2.gl).**
2. SilicoDArT bootstrap with the `utils.kin.dominant` estimator (F2).
   **Consequence: numerical output changes for all SilicoDArT standard
   errors and interval widths.**
3. Document the missing-data effect and warn when any individual's call
   rate is below 0.8 (F3). No numerical change.
4. Warn when `kin` comes from an estimator other than the one bootstrapped;
   document it in `@param kin` (F4). No numerical change.
5. Standards cleanup and an `na.rm` fix to the printed diagnostic;
   compute `kin` only when it is printed (F5). Returned values unchanged.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- Spec: behaviour against roxygen on testset2.gl and testset2.gs — run
- Independent numerical check: the internal SNP matrix was recomputed and compared with `gl.kin` (grm); the SilicoDArT estimator was compared with `gl.kin` (dominant); the SE with and without the `MS` shift was recomputed with the same seed; call rate against SE, with a 90%-missing individual constructed as a check
- Callers (API3): none in dartR.captive other than `@seealso` in `gl.report.kin.classes`
- EMIBD9 input: SKIPPED as a run — the EMIBD9 binary is not used in tests; F4 follows from reading the code
- FBM path (DAT6): SKIPPED — no FBM fixture
- Known complaints: not checked — the function was added in the 2026 kinship series and has no release history
- PLT: not applicable — no plot

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence (SNP SEs drop ~15%) approved |
| 2 | approved | Luis | consequence (all SilicoDArT SEs change) approved |
| 3 | approved | Luis | |
| 4 | approved | Luis | |
| 5 | approved | Luis | |

## Outcome

- Change 1: resamples are `G/2`. Seeded median SE 0.0258 -> 0.02203, matching the independent recomputation in F1; `CB_AB_01`/`CB_AB_02` SE 0.0236 -> 0.0198.
- Change 2: resamples use the `utils.kin.dominant` estimator; on all loci it equals `gl.kin(testset2.gs)` exactly (max difference 0). Seeded median SE 0.0182 -> 0.0179.
- Change 3: `@details` paragraph, `verbose >= 1` warning; on testset2.gl it fires for 24 individuals, all of them EmmacCaptBred (call rate 0.70-0.80), which are the pairs the example reports. Example comment updated.
- Change 4: warning when `attr(kin, "method")` differs; `@param kin` updated.
- Change 5: `build =` dropped, `report()` in the summary, warning string fixed, `na.rm` in the resolvable diagnostic (no `NaN` with an `NA` in `kin`), `kin` computed only at `verbose >= 3`.
- Addendum correction to the verdict: the internal SNP matrix matches `gl.kin` closely, not exactly (off-diagonal at most 8e-5, diagonal at most 0.007); `@details` and the test state the measured bound.
- Snapshot: 3 baseline values changed (SNP median SE, pair SE, SilicoDArT median SE), all mapped to changes 1 and 2; 17 expectations pass. Examples run.
- `devtools::document()` again dropped stale `importFrom(stats, dnorm/qnorm)` from NAMESPACE; reverted to keep scope.
- NEWS entry added. PR: (pending)

```json
{
  "function": "gl.report.kin.confidence",
  "package": "dartR.captive",
  "family": "report",
  "skill_version": "2.0.0",
  "commit": "cb8d955",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS3,VRB2,STY1", "status": "approved", "change": 5}
  ],
  "coverage_skipped": ["EMIBD9 run: binary not used in tests", "DAT6: no FBM fixture", "forum/issues: no release history"],
  "status": "in-apply",
  "pr": null
}
```
