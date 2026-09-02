# Review: gl.assign.grm (dartR.captive)
- Family mode: analysis
- Date: 2026-09-01
- Reviewer: Claude (claude-sonnet-5), dartr-function-review v1.0.0
- Package commit: 571d6a6
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl.assign.grm.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — the function follows the canonical FS1-FS8
anatomy, but has no completion message at any verbosity (FS9), and lacks
the input guards its sibling assignment functions carry, so malformed
input crashes with an internal R error instead of a clear one.
**Spec: Rework** — `@description` and `@return` promise assignment
probabilities returned as a `data.frame`; the code returns unnormalized
mean pairwise relatedness scores as a named numeric vector. No probability
is ever computed.

## Findings

**F1 [HIGH, confidence: high] — crashes with an internal R error when `pop(x)` is NULL (DAT5)**
`R/gl.assign.grm.r:57-60` — `vec <- as.vector(pop(x))` is never checked for
`NULL`. Downstream, `split()`/`order()` on the resulting `rel_list` fails
with `argument 1 is not a vector` — an error that names no dartR function
and gives the user no indication their genlight object lacks population
assignments.
Failure scenario: `pop(sub) <- NULL; gl.assign.grm(sub, unknown = "AA019237")`
throws `Error in order(rel_list, decreasing = TRUE) : argument 1 is not a
vector` with no hint of the cause (reproduced). The sibling `gl.assign.pa`
guards the same precondition explicitly:
`if (is.null(pop(x))) stop(error("Fatal Error: Population assignments
(pop(x)) are NULL.\n"))` (`R/gl.assign.pa.r:109-111`); `gl.grm` instead
auto-assigns a default population. `gl.assign.grm` does neither.
Proposed change: add the same `is.null(pop(x))` guard as `gl.assign.pa`,
in "FUNCTION SPECIFIC ERROR CHECKING", before "DO THE JOB".

**F2 [HIGH, confidence: high] — documented output ("assignment probabilities", `data.frame`) doesn't match what the function returns (DOC5)**
`R/gl.assign.grm.r:4-6,19-20` vs `:83-91` — `@description` says the
function "estimates their probability of coming from individual
populations"; `@return` says "A `data.frame` consisting of assignment
probabilities for each population." The code computes, per candidate
population, the mean of the unknown individual's pairwise GRM values
against that population's members (`rel_list <- unlist(lapply(x_split,
function(y) mean(y$Freq)))`) and returns that as a sorted named `numeric`
vector. No normalization to `[0,1]`, no probability model, and no
`data.frame` construction happen anywhere in the function.
Failure scenario: `class(gl.assign.grm(testset.gl[1:30,1:200],
unknown = "UC_00126"))` is `"numeric"`, not `"data.frame"` (reproduced).
A caller following the documented contract — e.g. indexing a `pop` column
or expecting values that sum to 1 — gets a type error or silently
misreads a raw relatedness score as a probability.
Proposed change: rewrite `@description`/`@return` to describe the actual
output (mean pairwise relatedness score per candidate population, as a
named numeric vector sorted from most to least related) — this is the
minimal, non-behaviour-changing fix. Reworking the code to genuinely
return probabilities in a `data.frame` is a larger, separate change that
needs the custodian's design input; not proposed here.

**F3 [MEDIUM, confidence: high] — no progress or completion output at any verbosity (FS9, VRB1, DOC5)**
`R/gl.assign.grm.r` — the function contains no `cat(report(...))` or
`cat(important(...))` call anywhere in its body, and no FS9 "FLAG SCRIPT
END" block. `@param verbose` promises "3, progress and results summary; 5,
full report", but nothing beyond `utils.flag.start`'s own start banner and
`utils.check.datatype`'s generic message ever prints, at any verbose
level.
Failure scenario: `gl.assign.grm(sub, unknown = "AA019237", verbose = 3)`
prints only "Starting gl.assign.grm" and a generic NA-loci warning, then
returns silently — no results summary, no "Completed" line (reproduced).
A user scripting at `verbose = 1` (the documented "begin and end" level)
gets a begin with no matching end.
Proposed change: add the standard FS9 completion line
(`if (verbose >= 1) cat(report("Completed:", funname, "\n"))`) before
`return(rel_list)`, and a `verbose >= 3` summary (e.g., print the top
candidate population and its score).

**F4 [MEDIUM, confidence: high] — duplicate individual names are silently merged into "unknown" (DAT5)**
`R/gl.assign.grm.r:57-58` — `vec[indNames(x) == unknown] <- "unknown"`
matches on name equality with no check that `indNames(x)` is unique. If
two individuals happen to share the value passed as `unknown`, both are
reassigned to population `"unknown"` and both are removed from their true
populations, with no warning.
Failure scenario: duplicating `indNames(sub)[1]` onto `sub`'s second
individual, then calling `gl.assign.grm(sub, unknown = indNames(sub)[1])`,
runs to completion and returns a plausible-looking numeric vector, but two
individuals — not one — have been pulled out of their populations and
folded together as a single "unknown" (reproduced; no error or message of
any kind). `gl.assign.pa` guards this exact precondition:
`if (any(duplicated(indNames(x)))) stop(error("Fatal Error: Duplicate
individual names in genlight object.\n"))` (`R/gl.assign.pa.r:91-93`).
Proposed change: add the same duplicate-name guard as `gl.assign.pa`.
**Consequence: a genlight object with duplicate individual names, which
currently runs and returns a (silently wrong) numeric result, will error
instead.**

**F5 [LOW, confidence: high] — dead code: the best-match population is computed and discarded**
`R/gl.assign.grm.r:89` — `res <- names(rel_list)[1]` computes the
top-scoring population name (since `rel_list` was just sorted decreasing)
but `res` is never read again; only `rel_list` is returned.
Failure scenario: none (no incorrect behaviour) — this is unused work that
also reads as an abandoned attempt to surface a "best guess" answer,
which the function's name and documentation both suggest should exist.
Proposed change: either delete the unused assignment, or use it — e.g.
print "Best-matching population: `<res>`" at `verbose >= 3`.

**F6 [LOW, confidence: high] — `@examples` requires `gplots` even though this function's call path never uses it (DOC3)**
`R/gl.assign.grm.r:25` — the example gates on
`requireNamespace("gplots", quietly = TRUE)` alongside `rrBLUP`. Internally,
`gl.assign.grm` calls `gl.grm(x, plotheatmap = FALSE, verbose = 0)`, and
`gl.grm`'s `gplots::heatmap.2()` call — and its own `gplots` guard — only
run when `plotheatmap == TRUE`. `gplots` is Suggests-only and unrelated to
any code path `gl.assign.grm` can reach.
Failure scenario: a machine with `rrBLUP` installed but not `gplots` skips
the example during `R CMD check`/CRAN checks even though the example would
run fine.
Proposed change: drop the `gplots` clause from the example's
`requireNamespace` guard, keeping only `rrBLUP`.

**F7 [INFO, confidence: high] — stale, non-roxygen `@param` line for an unimplemented argument**
`R/gl.assign.grm.r:10` — `# @param inbreeding_par The inbreeding parameter
[default 0].` uses a single `#`, so roxygen never renders it, and no
`inbreeding_par` argument exists in the function signature. It reads as an
abandoned feature stub.
Proposed change: remove the line, or implement `inbreeding_par` and
document it properly if the feature is still wanted.

**F8 [INFO, confidence: high] — `@author` has no separate `Author(s):` line (proposed rule DOC7)**
`R/gl.assign.grm.r:21-22` — `@author Custodian: Luis Mijangos -- Post to
...` names only a custodian. DOC7 (proposed) expects a distinct
`Author(s):` label even when author and custodian are the same person.
Proposed change: prefix with `Author(s): Luis Mijangos. ` before the
existing `Custodian:` text.

## Proposed changes

1. Add an `is.null(pop(x))` guard, matching `gl.assign.pa`, so missing
   population assignments fail with a clear message instead of an opaque
   internal error (F1).
2. Rewrite `@description`/`@return` to describe the actual output (a named
   numeric vector of mean pairwise relatedness scores, most to least
   related), instead of "assignment probabilities" / `data.frame` (F2).
3. Add an FS9 completion message and a `verbose >= 3` results summary (F3).
4. Add a duplicated-`indNames` guard, matching `gl.assign.pa` (F4).
   **Consequence: a genlight object with duplicate individual names, which
   currently runs and returns a (silently wrong) numeric result, will
   error instead.**
5. Remove the unused `res <- names(rel_list)[1]` line, or use it to print
   the best-matching population at `verbose >= 3` (F5).
6. Drop the unnecessary `gplots` clause from the `@examples`
   `requireNamespace` guard (F6).
7. Remove the stale commented-out `inbreeding_par` `@param` line (F7).
8. Add the `Author(s):` label to the `@author` block (F8, proposed rule).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behavior vs roxygen on testset.gl, verified empirically via
  `devtools::load_all()` — run; confirmed the return type/semantics
  mismatch (F2) and cross-checked the relatedness averages against an
  independent manual computation (`tapply` on the `gl.grm` matrix
  directly) — numbers match exactly, so the core arithmetic is correct
  even though its documented framing is not
- SNP vs SilicoDArT dispatch (analysis family check): `gl.assign.grm` has
  no local gate, but inherits `gl.grm`'s `datatype == "SilicoDArT"` stop
  (fixed in the `gl.grm` review, PR #85) — verified
  `gl.assign.grm(testset.gs[1:10,1:100], unknown = indNames(testset.gs)[1])`
  errors with "SilicoDArT" — run
- Edge cases: NULL `pop(x)` (F1, reproduced), duplicate `indNames` (F4,
  reproduced), unknown individual not in `x` (already guarded, verified) —
  run
- `@family` tag: absent, per DOC1. Not raised as a finding — only 4 of 13
  `gl.*` files in this package carry `@family` at all, including no other
  `gl.assign.*` sibling, so this is a repo-wide gap predating this
  function, not a defect specific to it.
- PLT/PLT3 (plot/results coupling): not applicable — this function has no
  plotting parameters of its own.
- FBM path (DAT6): not run — `gl.assign.grm`'s own body never densifies;
  densification, if any, happens inside `gl.grm` (covered by that review).

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Bernd Gruber | |
| 2 | approved | Bernd Gruber | |
| 3 | approved | Bernd Gruber | |
| 4 | approved | Bernd Gruber | consequence (duplicate names now error) explicitly approved |
| 5 | approved | Bernd Gruber | implemented as "use it" (print best match at verbose >= 3), not deletion |
| 6 | approved | Bernd Gruber | |
| 7 | approved | Bernd Gruber | |
| 8 | approved | Bernd Gruber | |

## Outcome

All eight approved changes applied to `R/gl.assign.grm.r` on 2026-09-01.

- Change 1 (F1): added an `is.null(pop(x))` guard, matching `gl.assign.pa`.
  Verified: `pop(sub) <- NULL; gl.assign.grm(sub, unknown = "AA019237")`
  now errors with "Fatal Error: Population assignments (pop(x)) are NULL."
  instead of the opaque `order()` error.
- Change 2 (F2): rewrote `@description` and `@return` to describe the
  actual output (mean pairwise GRM relatedness score per candidate
  population, as a named numeric vector). No code/behaviour change.
  `devtools::document()` run; `man/gl.assign.grm.Rd` updated.
- Change 3 (F3): added a `verbose >= 3` results summary (best-matching
  population and its score) and the FS9 `verbose >= 1` completion message.
  Verified: `gl.assign.grm(sub, unknown = "AA019237", verbose = 3)` now
  prints both.
- Change 4 (F4): added an `any(duplicated(indNames(x)))` guard, matching
  `gl.assign.pa`. Verified: duplicating an individual's name then calling
  the function now errors with "Fatal Error: Duplicate individual names in
  genlight object." instead of silently merging two individuals into
  `'unknown'`. Grepped `dartR.base`, `dartR.data`, `dartR.popgen`,
  `dartR.sexlinked`, `dartR.sim`, `dartR.spatial` (local sibling checkouts
  under `D:\Bernd\R\dartRs`) for `gl.assign.grm(` — no external callers
  found. NEWS.md entry added.
- Change 5 (F5): the previously-unused `res <- names(rel_list)[1]` is now
  used in the `verbose >= 3` summary added for change 3.
- Change 6 (F6): dropped the unnecessary `gplots` clause from the
  `@examples` `requireNamespace` guard.
- Change 7 (F7): removed the stale, non-functional `inbreeding_par`
  `@param` line.
- Change 8 (F8): `@author` now states `Author(s): Luis Mijangos.
  Custodian: Luis Mijangos -- ...`.

Snapshot verification: `tests/testthat/test-gl.assign.grm.R` re-run against
the patched function. One diff from baseline: the `pop(x) = NULL` test's
expected error message changed from the opaque `"argument 1 is not a
vector"` to `"Population assignments"` — this maps directly to approved
change 1, and the test was updated to assert the new, approved behaviour.
A new test for change 4 (duplicate `indNames` now errors) was added. All
11 expectations pass. The numeric-output snapshot test (mean relatedness
values on `testset.gl[1:30, 1:200]`) is unchanged, confirming changes 2,
3, 5, 6, 7, 8 didn't alter computed output for the default path.
`devtools::document()` run end to end; ran
`gl.assign.grm(sub, unknown = "AA019237", verbose = 3)` on reference data
and confirmed begin/end messages, the results summary, and the returned
vector all behave as expected.
PR: #88 (https://github.com/green-striped-gecko/dartR.captive/pull/88),
merged into `dev` 2026-09-01.

```json
{
  "function": "gl.assign.grm",
  "package": "dartR.captive",
  "family": "analysis",
  "skill_version": "1.0.0",
  "commit": "571d6a6",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS9/VRB1", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "principle:dead-code", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC3", "status": "approved", "change": 6},
    {"id": "F7", "severity": "INFO", "confidence": "high", "rule": "STY1", "status": "approved", "change": 7},
    {"id": "F8", "severity": "INFO", "confidence": "high", "rule": "DOC7 (proposed)", "status": "approved", "change": 8}
  ],
  "coverage_skipped": ["DAT6: not applicable, no densification in this function's own body"],
  "status": "done",
  "pr": 88
}
```
