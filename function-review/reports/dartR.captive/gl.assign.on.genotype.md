# Review: gl.assign.on.genotype (dartR.captive)
- Family mode: analysis
- Date: 2026-09-04
- Reviewer: Claude (claude-sonnet-5), dartr-function-review v1.0.0
- Package commit: ad61194
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl.assign.on.genotype.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — the FS1-FS9 anatomy is followed correctly, but
a results table prints unconditionally regardless of `verbose` (VRB3), a
deprecated argument to `gl.join` produces an unconditional warning, and the
returned (filtered) genlight object gets no history entry (FS8).
**Spec: Rework** — the documented "unknown-only, with a warning" fallback
for zero or too-few qualifying populations does not exist in the code;
both plausible parameter values (`aic.threshold = 1`, `n.best` larger than
the candidate count) crash with opaque internal errors instead. Separately,
the genotype-likelihood model is diploid-HWE-only, silently misapplied to
SilicoDArT (ploidy = 1) data with no dispatch guard or warning.

## Findings

**F1 [HIGH, confidence: high] — crashes instead of the documented "unknown-only" fallback when no population meets `aic.threshold` (DOC5)**
`R/gl.assign.on.genotype.r:293-299` — `@return` states: "If no such
populations, the genlight object contains only data for the unknown
individual with a warning." When every candidate population's `aic.wt`
falls below `aic.threshold`, `pop.keep <- result$pop[result$assign ==
"yes"]` is `character(0)`, and `gl.keep.pop(knowns, pop.list =
character(0), verbose = 0)` throws `Fatal Error: no populations listed to
keep!` — an uncaught error, not the documented warning-and-continue path.
Failure scenario: `gl.assign.on.genotype(sub, unknown = "UC_00126", nmin =
5, aic.threshold = 1, verbose = 0)` (reproduced) stops with that internal
error instead of returning a genlight containing only `UC_00126`.
`aic.threshold = 1` is a valid value within the documented `[0,1]` range,
and the same crash occurs with any threshold high enough that no
population clears it, which happens routinely for a genuine outlier.
Proposed change: before calling `gl.keep.pop`, check
`length(pop.keep) == 0` and, if so, skip straight to `gl.out <-
unknown.ind` with the documented warning, matching the `@return` contract.

**F2 [HIGH, confidence: high] — crashes when `n.best` exceeds the number of qualifying populations (DOC5)**
`R/gl.assign.on.genotype.r:293-294` — `pop.keep <- result$pop[1:n.best]`
has no bound check against `nrow(result)`. When `n.best` is larger than
the number of candidate populations, `1:n.best` indexes past the end of
`result$pop`, producing `NA` entries in `pop.keep`.
Failure scenario: `gl.assign.on.genotype(sub, unknown = "UC_00126", nmin =
5, n.best = 5, verbose = 0)` with only 2 candidate populations available
(reproduced) fails downstream with `Subsetting resulted in zero
individuals` — again an opaque internal error, not the documented
"unknown-only, with a warning" fallback, and not an error that names
`n.best` as the cause.
Proposed change: cap `n.best` at `nrow(result)` before indexing (e.g.
`pop.keep <- result$pop[seq_len(min(n.best, nrow(result)))]`), optionally
with a `verbose >= 2` note that fewer than `n.best` populations were
available.

**F3 [HIGH, confidence: high] — no SNP/SilicoDArT dispatch: diploid HWE genotype model silently applied to presence/absence data (DAT1)**
`R/gl.assign.on.genotype.r:76,209-246` — `datatype <-
utils.check.datatype(x, verbose = verbose)` is computed but never read
again. `utils.gen.prob()` always treats each locus value as a diploid
genotype (`g == 0/1/2` → `pAA`/`pAB`/`pBB` under Hardy-Weinberg), with no
branch for SilicoDArT (ploidy 1, presence/absence coded 0/1, where `g ==
2` never occurs and `g == 1` is not a meaningful "heterozygote").
Failure scenario: `gl.assign.on.genotype(testset.gs[1:30,1:200], unknown =
indNames(testset.gs)[2], nmin = 5, verbose = 0)` (reproduced) runs to
completion, prints a results table, and returns a genlight object with no
error or warning that the underlying statistic is invalid for this data
type. A user working with SilicoDArT data — one of the package's two
primary supported types — gets a confident-looking but statistically
meaningless AIC ranking.
Proposed change: branch on `datatype` — either implement a presence/absence
appropriate likelihood (e.g. binomial on the dominant-marker allele
frequency) or `stop(error(...))` for SilicoDArT input until that model
exists, matching the pattern other `gl.assign.*` siblings use to gate
input type.

**F4 [MEDIUM, confidence: high] — results table prints unconditionally regardless of `verbose` (VRB3)**
`R/gl.assign.on.genotype.r:284` — `print(result)` has no `verbose` guard,
unlike the `verbose >= 3` block immediately after it (lines 286-288).
Failure scenario: `gl.assign.on.genotype(sub, unknown = "UC_00126", nmin =
5, verbose = 0)` (reproduced) prints the full population/AIC/assign table
to the console even though `verbose = 0` is documented as "silent or fatal
errors" — a script capturing only fatal-error output gets an unexpected
table on stdout every call.
Proposed change: gate the `print(result)` call behind `verbose >= 3`,
matching the "results summary" verbosity level.

**F5 [LOW, confidence: high] — stale deprecated argument to `gl.join` produces an unconditional warning (STY3)**
`R/gl.assign.on.genotype.r:300` — `gl.join(gl.out, unknown.ind,
method = "join.by.loc", verbose = 0)` still passes `method`, which
`gl.join`'s current signature (`function(x1, x2, method = NULL, verbose =
NULL)`) treats as deprecated: any non-`NULL` value triggers `cat(warn("
Warning: The parameter method is deprecated, no longer required"))`
unconditionally inside `gl.join`, independent of the `verbose` this
function passes.
Failure scenario: every call — including at `verbose = 0` (reproduced) —
prints that deprecation warning; and `gl.join` also has a
`stop(error(...))` branch for unrecognized `method` values, so this call
site is one `gl.join` release away from a hard error for a value that
already does nothing.
Proposed change: drop the `method = "join.by.loc"` argument from the call.

**F6 [LOW, confidence: high] — per-individual population log-likelihoods are computed and discarded (efficiency, no rule fits)**
`R/gl.assign.on.genotype.r:232-242` — inside `utils.gen.prob()`,
`logLs_pop <- apply(popmat, 1, function(ind) {...})` and `logL_pop_totals
<- colSums(logLs_pop, na.rm = TRUE)` compute a per-locus log-likelihood
for every individual in the candidate population, then sum it — and
neither `logLs_pop` nor `logL_pop_totals` is read anywhere afterward; only
`logL_focal_total` (the unknown's likelihood) is returned.
Failure scenario: none incorrect — this is `O(n_loci x n_ind)` work
repeated for every candidate population on every call, discarded
immediately, roughly doubling the loop's runtime for no effect.
Proposed change: delete lines 231-242 (`logLs_pop`/`logL_pop_totals`).

**F7 [LOW, confidence: high] — no history entry on the returned (filtered) genlight object (FS8)**
`R/gl.assign.on.genotype.r:299-300` — `gl.out` is built by filtering
`knowns` to the retained populations and joining in the unknown individual
— a modified genlight object handed back to the caller — but no
`x@other$history` entry is appended, unlike the FS8 convention for
functions that return a modified object.
Failure scenario: none crash-worthy; a caller inspecting
`gl.out@other$history` to audit how the object was derived sees no trace
of this call. Note: the sibling functions `gl.assign.mahal` and
`gl.assign.mahalanobis` have the same gap, so this may reflect a
family-wide, not function-specific, omission.
Proposed change: append
`nh <- length(gl.out@other$history); gl.out@other$history[[nh + 1]] <-
match.call()` before `return(gl.out)`.

**F8 [INFO, confidence: high] — local closure named `utils.gen.prob` shadows the package's reserved `utils.*` naming convention (FS1)**
`R/gl.assign.on.genotype.r:209` — `utils.gen.prob` is defined as a local
closure inside `gl.assign.on.genotype`, not as its own `R/utils.gen.prob.r`
file. FS1 reserves the `utils.*` prefix for internal helper files. No
conflicting exported function exists today, but the name reads as a
package-level utility to anyone grepping for one.
Proposed change: rename the local closure (e.g. `calc.genotype.logL`) to
avoid the `utils.*` naming collision.

**F9 [INFO, confidence: high] — typos in `@param` text**
`R/gl.assign.on.genotype.r:31,34` — "for which their is considered some
support" and "or more if their are ties" should read "there is" and
"there are".
Proposed change: fix both typos.

**F10 [INFO, confidence: high] — `@param` roxygen order doesn't match the function signature order**
`R/gl.assign.on.genotype.r:29-39` vs `:59-64` — roxygen documents
`nmin`, `aic.threshold`, `n.best`, `verbose`; the signature is `nmin`,
`n.best`, `aic.threshold`, `verbose`.
Proposed change: reorder either the roxygen block or the signature so they
match (house order per DOC1 follows the signature).

## Proposed changes

1. Handle the zero-qualifying-population case per the documented contract:
   return `unknown.ind` alone with a warning instead of crashing inside
   `gl.keep.pop` (F1).
2. Bound `n.best` to the number of available candidate populations before
   indexing `result$pop[1:n.best]` (F2).
3. Add a datatype dispatch guard for `utils.gen.prob`: implement or reject
   SilicoDArT input explicitly instead of silently applying the diploid
   HWE model (F3).
   **Consequence: calls on SilicoDArT (ploidy = 1) data that currently run
   to completion will either use a different statistic or stop with an
   error, depending on the option chosen.**
4. Gate `print(result)` behind `verbose >= 3` (F4).
5. Drop the deprecated `method = "join.by.loc"` argument from the
   `gl.join()` call (F5).
6. Delete the unused `logLs_pop`/`logL_pop_totals` computation inside
   `utils.gen.prob` (F6).
7. Append an FS8 history entry to `gl.out` before returning it (F7).
8. Rename the local `utils.gen.prob` closure to avoid the reserved
   `utils.*` naming convention (F8).
9. Fix the two `@param` typos ("their is"/"their are" → "there is"/"there
   are") (F9).
10. Reorder the `@param` roxygen block to match the function signature
    (F10).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behavior vs roxygen on testset.gl (SNP) and testset.gs
  (SilicoDArT) — run
- Edge cases exercised: unknown not in object, `pop(x)` NULL, duplicate
  individual names, zero qualifying populations (`aic.threshold = 1`),
  `n.best` exceeding candidate count, SilicoDArT input — run
- FBM path (DAT6): SKIPPED — no FBM fixture available for this function;
  the function densifies only per-candidate-population matrices via
  `seppop()`/`as.matrix()`, not the whole object, so the risk is lower
  than a full `as.matrix(x)` call, but this was not verified against an
  FBM-backed object
- Plot output (PLT): not applicable — function has no plotting path

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Bernd Gruber | |
| 2 | approved | Bernd Gruber | |
| 3 | approved | Bernd Gruber | reject SilicoDArT explicitly (minimal guard, not a new model) |
| 4 | approved | Bernd Gruber | |
| 5 | approved | Bernd Gruber | |
| 6 | approved | Bernd Gruber | |
| 7 | approved | Bernd Gruber | |
| 8 | approved | Bernd Gruber | |
| 9 | approved | Bernd Gruber | |
| 10 | approved | Bernd Gruber | |

## Outcome

All 10 approved changes applied in `R/gl.assign.on.genotype.r`.
Characterization test (`tests/testthat/test-gl.assign.on.genotype.R`, 7
tests) re-run after the change: every diff from the pre-review baseline
maps to an approved change, no unexplained diffs.

- F1/F2: `gl.assign.on.genotype(sub, unknown = "UC_00126", nmin = 5,
  aic.threshold = 1)` now returns a genlight with only `UC_00126` and
  prints the documented warning, instead of crashing in `gl.keep.pop`.
  `n.best = 5` (2 candidate populations available) now returns `nInd = 20`
  (both populations retained, capped) instead of crashing.
- F3: `gl.assign.on.genotype(testset.gs[1:30,1:200], unknown =
  indNames(testset.gs)[2])` now stops with a clear error naming the
  diploid-HWE/SilicoDArT mismatch, instead of returning a genlight object
  with statistically invalid results.
- F4/F5: at `verbose = 3` on `testset.gl[1:30,1:200]`, the results table
  prints once (gated), and the previously unconditional `gl.join`
  deprecation warning no longer appears.
- F7: `gl.out@other$history` carries a `match.call()` entry for this call
  (verified: last history entry is
  `gl.assign.on.genotype(x = testset.gl[1:30, 1:200], unknown =
  "UC_00126", nmin = 5, verbose = 3)`).
- F6/F8/F9/F10: verified by reading the diff; no behavioral effect.
- `devtools::document()` run; `man/gl.assign.on.genotype.Rd` regenerated.
- PR: #91 (`review/gl.assign.on.genotype` -> `dev`).

```json
{
  "function": "gl.assign.on.genotype",
  "package": "dartR.captive",
  "family": "analysis",
  "skill_version": "1.0.0",
  "commit": "ad61194",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DAT1", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "STY3", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "principle: dead computation", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "FS8", "status": "approved", "change": 7},
    {"id": "F8", "severity": "INFO", "confidence": "high", "rule": "FS1", "status": "approved", "change": 8},
    {"id": "F9", "severity": "INFO", "confidence": "high", "rule": "DOC6", "status": "approved", "change": 9},
    {"id": "F10", "severity": "INFO", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 10}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture available for this function"],
  "status": "pr-open",
  "pr": 91
}
```
