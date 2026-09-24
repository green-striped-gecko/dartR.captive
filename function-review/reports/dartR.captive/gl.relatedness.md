# Review: gl.relatedness (dartR.captive)
- Family mode: analysis
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 103b47a (origin/dev) plus `R/gl.relatedness.R`, `man/gl.relatedness.Rd` and `tests/testthat/test-gl.relatedness.R` from `dev_luis` (b310230) including uncommitted working-copy edits (install instructions, pre-filter message through `cat(warn())`), at the custodian's request. The function is not on `dev`; this branch adds it.
- Engine: dartR.coancestry 0.1.0 (closed-source binary); independent reference: related (Coancestry Fortran)
- Datasets: platypus.gl (subsets of 12-81 individuals); simulated biallelic SNP data with known relationships (1,500 loci; 20, 50 and 120 individuals; unrelated pairs and full sibs); dartR.data 1.2.6
- Baseline: tests/testthat/test-gl.relatedness.R (10 existing tests plus 2 characterization tests captured pre-review, 42 expectations, all pass)

## Verdict

**Standards: Rework**: adding the function to `dev` as it stands makes
R CMD check warn that `dartR.coancestry` is undeclared, and the CI
workflow fails on warnings. Arguments are not validated, and the
plotting path has small defects.
**Spec: Needs work**: `lynchrd`, `ritland` and `quellergt` match the
Coancestry reference (`related`) to its 4-decimal rounding. `wang` and
`lynchli` sit above it by an offset that shrinks with sample size
(+0.047 and +0.018 on 30 platypus). An individual with no genotype calls
is reported as unrelated (0) to everyone by the moment estimators. The
matrices carry no scale tag, so the kinship series reads them as kinship,
twice too large.

## Findings

**F1 [BLOCKER, confidence: high] — the engine is an undeclared dependency; CI fails (DEP1)**
`R/gl.relatedness.R:118, 151` — `requireNamespace("dartR.coancestry")`
and `dartR.coancestry::relatedness_cpp()` name a package that is not in
DESCRIPTION. `R CMD check` (run on this branch) gives
`WARNING: '::' or ':::' import not declared from: 'dartR.coancestry'`.
`.github/workflows/check-standard.yml` uses `check-r-package@v2`, which
fails on warnings by default.
Failure scenario: merging this function into `dev` turns the `dev` CI red
for every later PR.
Proposed change: the idiom `gl.diagnostics.relatedness` uses for
`related`: hold the package name in a variable for `requireNamespace()`
and fetch the engine with `getExportedValue(pkg, "relatedness_cpp")`.
Declaring it in Suggests is not an option, because it is not in a
CRAN-like repository.

**F2 [HIGH, confidence: high] — an individual with no calls is reported as unrelated (DAT integrity)**
`R/gl.relatedness.R:144-186` — for a pair that shares no called loci,
the engine returns 0 for the moment estimators and `NaN` for `dyadml`;
the wrapper passes both through. With `n.boots > 0` the same individual
is dropped instead, so the matrices are 11 x 11 for a 12-individual
object, contrary to `@return` ("rownames/colnames = indNames(x)").
Failure scenario: a sample that failed genotyping appears unrelated
(r = 0) to every animal and is paired or released as unrelated.
Proposed change: count shared called loci per pair
(`crossprod(!is.na(dosage))`); set every estimator to `NA` for pairs with
none, with a `verbose >= 1` warning giving the count; with `n.boots > 0`,
return matrices over all of `indNames(x)` with `NA` rows for the dropped
individuals.
**Consequence: numerical output changes for pairs with no shared called
loci (0 or NaN becomes NA), and bootstrap matrices keep dropped
individuals as NA rows.**

**F3 [MEDIUM, confidence: medium] — `wang` and `lynchli` depart from Coancestry by a sample-size term (DOC5, proposed rule)**
Engine behaviour, reached through `R/gl.relatedness.R:151`. On 30
platypus with no missing data, `lynchrd`, `ritland` and `quellergt`
differ from `related::coancestry` by at most 7e-5 (its rounding);
`wang` is on average 0.047 higher (maximum 0.104) and `lynchli` 0.018
higher (correlation 1: an offset). In simulations the gap for unrelated
pairs is 0.050 at n = 20, 0.023 at n = 50 and 0.013 at n = 120, so it
behaves as a small-sample correction applied differently. Neither
implementation recovers r = 0 for unrelated pairs on biallelic SNPs
(both about -0.16), while `quellergt` does (-0.009 at n = 120); for full
sibs all give 0.57-0.58.
Failure scenario: a user compares `gl.relatedness` Wang estimates with
published Coancestry or `related` values from the same data and finds a
systematic 0.05 difference that the documentation does not mention.
Proposed change: document the difference in `@details` with these
numbers, and pin the agreement and the offset in the tests. Whether the
engine should match Coancestry is a decision for the engine's author and
lies outside this package.

**F4 [MEDIUM, confidence: high] — relatedness matrices are read as kinship by the kinship series (#110 contract)**
`R/gl.relatedness.R:166-170` — the matrices carry no `scale` attribute.
`utils.kin.as.kinship` treats an untagged matrix as kinship, so passing
`res$wang` to `gl.report.kin.classes` runs without a message on values
twice the kinship scale (relatedness = 2 x kinship). `gl.grm` output is
tagged `"relatedness"` since #110.
Failure scenario: a user passes Wang relatedness as `kin`; half-sib pairs
(r = 0.25) are classed as full-sib-level kinship.
Proposed change: `attr(M, "scale") <- "relatedness"` on every estimator
matrix, so the series halves it.
**Consequence: kinship-series functions given these matrices now halve
them.**

**F5 [LOW, confidence: high] — arguments are not validated; `plot.stat` is checked after the run (FS5)**
`R/gl.relatedness.R:108-117, 130, 189-194` — `n.boots = NA` stops with
`missing value where TRUE/FALSE needed`; `n.boots = -1`, `num.trios = 0`,
`n.threads = 0` and `rng.seed = NA` run without a message. An invalid
`plot.stat` is only detected after the engine has run, and only when
`plot.out = TRUE`.
Failure scenario: a `dyadml` run on a few hundred individuals finishes and
then stops on a mistyped `plot.stat`, discarding the results.
Proposed change: validate `n.boots` (>= 0), `num.trios` and `n.threads`
(>= 1), `rng.seed` (a number), `allow.inbreeding` and `plot.out`
(TRUE/FALSE) and `plot.stat` before the engine runs, each with
`stop(error(...))`.
**Consequence: calls with invalid values of these arguments stop instead
of running.**

**F6 [LOW, confidence: high] — plot saving depends on display; `gl.colors` prints at `verbose = 0` (PLT3, VRB1)**
`R/gl.relatedness.R:189-199` — with `plot.out = FALSE` and a `plot.file`,
nothing is saved. `gl.colors("div")` is called without `verbose`, so
"Starting gl.colors" prints at `verbose = 0`.
Failure scenario: a batch script sets `plot.out = FALSE, plot.file = "h"`
and finds no file.
Proposed change: build the heatmap when either `plot.out` or `plot.file`
asks for it, display it only when `plot.out = TRUE`, and pass
`verbose = 0` to `gl.colors`.

**F7 [LOW, confidence: high] — convention gaps (FS3, FS10, DOC2, DOC7 proposed rule)**
`:100` passes the outdated `build =`; `:204` returns the list visibly, so
an unassigned call prints every matrix; the `verbose` text is not the
DOC2 text; `@author` has no `Author(s):` part.
Failure scenario: `gl.relatedness(gl)` without assignment prints
thousands of lines.
Proposed change: drop `build =`, return `invisible(out)`, use the DOC2
text, and add `Author(s): Luis Mijangos` to `@author`.

Checked, nothing found: the matrices are symmetric and equal the long
table; names map back through all optional frames; the bootstrap
pre-filter leaves point estimates unchanged (maximum difference 0); the
SilicoDArT and non-diploid guards stop; unknown estimators stop.

## Proposed changes

1. Undeclared-engine idiom (`getExportedValue`) (F1). No output change.
2. `NA` for pairs with no shared called loci, full-size bootstrap
   matrices, warning (F2). **Consequence: numerical output changes for
   pairs with no shared called loci (0 or NaN becomes NA); bootstrap
   matrices keep dropped individuals as NA rows.**
3. `@details` on the `wang`/`lynchli` offset from Coancestry and on
   estimator behaviour for unrelated SNP pairs; tests pin the comparison
   (F3). Docs and tests only.
4. `scale = "relatedness"` tag on the matrices (F4). **Consequence:
   kinship-series functions given these matrices now halve them.**
5. Argument validation before the engine runs (F5). **Consequence: calls
   with invalid argument values stop instead of running.**
6. Plot saved whenever `plot.file` is given; `gl.colors` silenced (F6).
7. Convention cleanup, `invisible()` return (F7). **Consequence: an
   unassigned call no longer prints the result.**

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run
- R CMD check (no tests, no examples) on this branch: the dependency WARNING in F1; the other WARNING (packages built under R 4.4.3) and NOTEs are local or pre-existing
- Spec: behaviour against roxygen on platypus.gl subsets; all estimators; bootstrap; all-NA individual with and without bootstrap; invalid arguments; plotting and saving — run
- Independent numerical check: `related::coancestry` on 30 platypus (420 complete loci; 566 loci with 4.1% missing) and on simulated data with known relationships
- Kinship-series interoperability: `res$wang` passed to `gl.report.kin.classes`
- Engine internals: SKIPPED — closed-source binary; findings about estimator values are observations of the engine's output
- `trioml` against `related`: SKIPPED — its reference trios are drawn at random, so values are not comparable pair by pair
- Callers (API3): none — the function is new to `dev`
- FBM path (DAT6): SKIPPED — no FBM fixture
- Known complaints: not checked — unreleased
- Note: `related::coancestry` writes and reads files in the working directory; the characterization test calls it from the test directory

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |
| 2 | approved | Luis | consequence (NA instead of 0/NaN; NA rows for dropped individuals) approved |
| 3 | approved | Luis | the engine itself is out of scope |
| 4 | approved | Luis | consequence (kinship series halves these matrices) approved |
| 5 | approved | Luis | consequence (invalid arguments stop) approved |
| 6 | approved | Luis | |
| 7 | approved | Luis | consequence (no printing on an unassigned call) approved |

## Outcome

- Change 1: the engine name is held in a variable for `requireNamespace()`, and `relatedness_cpp` is fetched with `getExportedValue()`. R CMD check on the branch: "checking dependencies in R code ... OK" (was the WARNING in F1).
- Change 2: pairs sharing no called locus (`tcrossprod` of the call matrix) get `NA` in `$dyads`, `$delta19`, `$trio_delta` and the matrices; an individual with no calls gets `NA` in `$inbreeding`; warning with the count. Matrices span all of `indNames(x)`; with `n.boots > 0` a dropped individual keeps an `NA` row (12 x 12 for 12 individuals, was 11 x 11).
- Change 3: `@details` paragraph with the comparison numbers; a test pins agreement of `lynchrd`, `ritland` and `quellergt` with `related` (< 1e-4) and the `wang`/`lynchli` offsets. The test runs `related` in a temporary directory, because it writes files to the working directory.
- Change 4: every estimator matrix carries `attr(, "scale") = "relatedness"`; `utils.kin.as.kinship()` halves it.
- Change 5: `n.boots`, `num.trios`, `n.threads`, `rng.seed`, `allow.inbreeding`, `plot.out` and `plot.stat` are validated before the engine runs; 8 invalid values tested, each stops with `<argument> must be ...`.
- Change 6: the heatmap is built when `plot.out` or `plot.file` asks for it, drawn on a null device when not displayed; `gl.colors(verbose = 0)`. With `plot.out = FALSE, plot.file = "h"` one file is saved and no device is left open.
- Change 7: `build =` dropped, `invisible(out)`, DOC2 `verbose` text, `Author(s):` in `@author`.
- Snapshot: the 10 pre-existing tests pass unchanged; the characterization expectations changed only for approved changes (2: NA and matrix size, 4: scale tag, 5: argument errors, 6: saved file, 7: visibility). 56 expectations pass; 11 warnings come from `dartR.base::utils.heatmap` (a `dendextend` message, also seen in the `gl.run.EMIBD9` tests). Full suite 599 pass, 0 fail. The documented example runs at `verbose = 3`.
- NEWS entry added. PR: #129

## Addendum

**A1 [LOW, confidence: high] — `verbose = 3` prints no results summary (VRB1)**
Found while running the example in Phase C. The verbosity scale promises a
results summary at level 3; the function prints only the start and end
lines. Proposed change: at `verbose >= 3`, print the number of pairs, the
number set to `NA`, and the mean and range of each estimator. Not applied;
awaiting approval.

```json
{
  "function": "gl.relatedness",
  "package": "dartR.captive",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "103b47a+dev_luis:b310230+uncommitted",
  "verdict_standards": "rework",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "DEP1", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "medium", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "API1", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "PLT3,VRB1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "FS3,FS10,DOC2,DOC7", "status": "approved", "change": 7},
    {"id": "A1", "severity": "LOW", "confidence": "high", "rule": "VRB1", "status": "pending", "change": null}
  ],
  "coverage_skipped": ["engine internals: closed-source", "trioml vs related: random reference trios", "DAT6: no FBM fixture", "forum/issues: unreleased"],
  "status": "pr-open",
  "pr": 129
}
```
