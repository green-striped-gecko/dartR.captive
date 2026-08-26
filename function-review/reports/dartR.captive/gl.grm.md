# Review: gl.grm (dartR.captive)
- Family mode: analysis
- Date: 2026-08-26
- Reviewer: Claude (claude-sonnet-5), dartr-function-review v1.0.0
- Package commit: bdb623d
- Datasets: testset.gl, testset.gs
- Baseline: tests/testthat/test-gl.grm.R (snapshot captured pre-review)

## Verdict

**Standards: Needs work** — structure follows the canonical anatomy, but a
variable used outside the branch that defines it crashes a plausible
argument combination, and one documented default doesn't match the code.
**Spec: Rework** — the function computes an SNP-only additive relationship
matrix but never gates on marker type, so calling it on SilicoDArT data
returns numbers that contradict the function's own documented range,
without any warning.

## Findings

**F1 [HIGH, confidence: high] — no SNP/SilicoDArT dispatch guard**
`R/gl.grm.r:87,146` — `datatype <- utils.check.datatype(x, verbose = verbose)`
is computed and never read again. The `@details` section states diagonal
elements "range from 1 to 2" (mean `1+f`), a claim that only holds for
diploid SNP dosage (0/1/2) centered by `A.mat`. Presence/absence
(SilicoDArT, ploidy 1) data is 0/1, and `as.matrix(x) - 1` maps it to -1/0,
which `A.mat` treats as valid SNP dosage input without complaint.
Failure scenario: `gl.grm(testset.gs[1:10, 1:100], plotheatmap = FALSE)`
runs to completion and returns a 10x10 matrix with diagonal range
`0.182–0.438` and off-diagonal values as low as `-0.165` — nowhere near
the documented 1–2 range, and not a meaningful IBD/kinship estimate. No
error, no warning; a user would take the output at face value.
Proposed change: gate on `datatype` from `utils.check.datatype()`
(`accept = c("SNP")`, or an explicit `stop(error(...))` when
`datatype == "SilicoDArT"`) before calling `A.mat`.

**F2 [HIGH, confidence: high] — `p3` used outside the branch that defines it**
`R/gl.grm.r:172,196` — `p3` is assigned only inside
`if (plotheatmap == TRUE) { ... }` (line 172), but the plot-saving block
(line 196, `tmp <- utils.plot.save(p3, ...)`) runs whenever
`plot.file` is non-`NULL`, unconditional on `plotheatmap`.
Failure scenario:
`gl.grm(testset.gl[1:10, 1:100], plotheatmap = FALSE, plot.file = "grm_test", plot.dir = tempdir())`
throws `Error in ... : object 'p3' not found`. This is a reasonable call —
compute the matrix, skip the interactive heatmap, still keep a saved plot
object — and it currently cannot succeed.
Proposed change: either compute and save the heatmap object regardless of
`plotheatmap` (only gating the on-screen draw), or make the save block
conditional on `plotheatmap == TRUE` and state that constraint in the
`plot.file` parameter doc.

**F3 [MEDIUM, confidence: high] — documented default for `legendy` doesn't match the code (DOC1)**
`R/gl.grm.r:15,68` — roxygen states
`@param legendy y coordinates for the legend[default 1].` but the function
signature default is `legendy = 0.5`. `man/gl.grm.Rd` inherits the same
wrong value, since it's generated from this header.
Failure scenario: a user reading `?gl.grm` places the legend expecting
`legendy = 1` behavior and gets `0.5` instead; harmless to results but
directly contradicts the documentation.
Proposed change: change the roxygen default to `[default 0.5]` (or change
the code default to `1`, whichever value the custodian intends), then run
`devtools::document()`.

**F4 [MEDIUM, confidence: high] — `gplots` dependency guard fires even when the plot is disabled (DEP1)**
`R/gl.grm.r:102-110` — the `requireNamespace("gplots", ...)` guard runs
unconditionally in "FUNCTION SPECIFIC ERROR CHECKING", before
`plotheatmap` is consulted.
Failure scenario: a user without `gplots` installed who calls
`gl.grm(x, plotheatmap = FALSE)` — explicitly asking to skip the plot —
still gets `Package gplots needed for this function to work` and no
matrix, even though nothing plot-related is needed for that call.
Proposed change: move the `gplots` guard inside
`if (plotheatmap == TRUE) { ... }`, ahead of the `gplots::heatmap.2()`
call.

**F5 [LOW, confidence: high] — plot-only work runs unconditionally**
`R/gl.grm.r:129-159` — `colors_pops`, `df_colors_temp_1/2`, and the
`merge()`-built `df_colors` are computed regardless of `plotheatmap`, but
are only consumed inside the `if (plotheatmap == TRUE)` block (line
172-190, `ColSideColors`/`RowSideColors`/`legend()`).
Failure scenario: no incorrect output, but every `plotheatmap = FALSE`
call pays for a `gl.select.colors()` call and a data-frame `merge()` it
never uses — wasted work that scales with individual count.
Proposed change: move this block inside `if (plotheatmap == TRUE)`, ahead
of the `heatmap.2()` call it feeds.

**F6 [LOW, confidence: medium] — `plot.dir` isn't resolved via `gl.check.wd()` (FS7)**
`R/gl.grm.r:195-201` — `plot.dir` is passed straight to
`utils.plot.save(dir = plot.dir, ...)`. The sibling function
`gl.grm.network` resolves it first:
`plot.dir <- gl.check.wd(plot.dir, verbose = 0)` (`R/gl.grm.network.r:160`).
`utils.plot.save()` does have its own `NULL`-and-missing-directory
fallback (defaults to `getwd()`, then `tempfile()` if that path doesn't
exist), so this doesn't currently produce a wrong or unchecked write —
it's an inconsistency with the house pattern and with the sibling file,
not a reproduced failure.
Proposed change: add `plot.dir <- gl.check.wd(plot.dir, verbose = 0)`
near the top of the function, matching `gl.grm.network`.

**F7 [INFO, confidence: high] — `@author` has no separate `Author(s):` line (proposed rule DOC7)**
`R/gl.grm.r:47-48` — `@author Custodian: Arthur Georges -- Post to ...`
names only a custodian. DOC7 (proposed) expects a distinct `Author(s):`
label even when author and custodian are the same person.
Proposed change: prefix with `Author(s): Arthur Georges. ` before the
existing `Custodian:` text.

**F8 [INFO, confidence: medium] — full densification for FBM-backed input (proposed rule DAT6)**
`R/gl.grm.r:146` — `rrBLUP::A.mat(as.matrix(x) - 1, ...)` densifies the
whole genotype matrix. The function's own `@examples` demonstrate calling
it on an FBM-converted object, but only after subsetting to 10
individuals x 100 loci. On a genuinely large FBM-backed object (the
scenario FBM support exists for) this call would materialize the entire
matrix in memory, which is what FBM backing is meant to avoid.
`rrBLUP::A.mat` requires a dense matrix, so this is likely a hard
constraint of the third-party algorithm rather than something fixable
inside `gl.grm` — the more realistic fix is a documented size/FBM caveat
rather than an implementation change.
Proposed change: no code change proposed; add a note to `@details`
that `gl.grm` densifies its input and is not suited to full-size
FBM-backed objects.

## Proposed changes

1. Gate `gl.grm` on `datatype` and refuse (or clearly warn on)
   SilicoDArT/ploidy-1 input before calling `A.mat` (F1).
   **Consequence: SilicoDArT input that currently returns a matrix will
   instead error (or the returned values will change if a warning-only
   path is chosen) — this is an API1-class behavior change and needs
   explicit custodian sign-off.**
2. Fix the `p3`-undefined crash for `plotheatmap = FALSE` +
   `plot.file` set, by scoping the heatmap object construction to cover
   both consumers or gating the save on `plotheatmap` (F2).
3. Correct the `legendy` documented default to match the code, then run
   `devtools::document()` (F3).
4. Move the `gplots` dependency guard inside `if (plotheatmap == TRUE)`
   (F4).
5. Move the color/`df_colors` computation inside
   `if (plotheatmap == TRUE)` (F5).
6. Add `plot.dir <- gl.check.wd(plot.dir, verbose = 0)` near the top of
   the function, matching `gl.grm.network` (F6).
7. Add the `Author(s):` label to the `@author` block (F7, proposed rule).
8. Add an FBM/size caveat to `@details`; no code change (F8, proposed
   rule).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: behavior vs roxygen/`@details` on testset.gl and testset.gs,
  verified empirically via `devtools::load_all()` — run
- SNP vs SilicoDArT dispatch (analysis family check) — run, defect found
  (F1)
- FBM path (DAT6): partially run — `gl.gen2fbm()` + subset + `as.matrix()`
  verified to work at small scale; full-scale FBM densification behavior
  not run (no large FBM fixture available), reasoned from the code and
  `rrBLUP::A.mat`'s known interface instead
- `gl.grm.network` (sibling, same `@family`) read for comparison only; not
  in scope for this review's findings

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Bernd Gruber | error on SilicoDArT input (not warn-only) |
| 2 | approved | Bernd Gruber | |
| 3 | approved | Bernd Gruber | |
| 4 | approved | Bernd Gruber | |
| 5 | approved | Bernd Gruber | |
| 6 | approved | Bernd Gruber | |
| 7 | approved | Bernd Gruber | |
| 8 | approved | Bernd Gruber | |

## Outcome

All eight approved changes applied to `R/gl.grm.r` on 2026-08-26.

- Change 1 (F1): added a `datatype == "SilicoDArT"` gate that
  `stop(error(...))`s before `A.mat` is called. Verified:
  `gl.grm(testset.gs[1:10, 1:100], plotheatmap = FALSE)` now errors with
  "gl.grm computes an additive relationship matrix for SNP (diploid) data;
  it is not valid for SilicoDArT (presence/absence) data." Grepped
  `dartR.base`, `dartR.data`, `dartR.popgen`, `dartR.sexlinked`,
  `dartR.sim`, `dartR.spatial` (local sibling checkouts under
  `D:\Bernd\R\dartRs`) for `gl.grm(` — no external callers found. The two
  internal callers (`R/gl.assign.grm.r:62`,
  `R/utils.functions.diagnostics.relatedness.r:62`) both feed SNP-workflow
  genlight objects, unaffected. NEWS.md entry added.
- Change 2 (F2): moved the plot-save block inside
  `if (plotheatmap == TRUE)`; added an `else if (!is.null(plot.file))`
  branch that warns nothing was saved. Verified:
  `gl.grm(testset.gl[1:10, 1:100], plotheatmap = FALSE, plot.file = "grm_test", plot.dir = tempdir())`
  now returns the matrix (no error) and prints the warning; no file is
  written. NEWS.md entry added.
- Change 3 (F3): `@param legendy` roxygen default corrected to `0.5`;
  `devtools::document()` run, `man/gl.grm.Rd` updated.
- Change 4 (F4): the `gplots` `requireNamespace` guard moved inside
  `if (plotheatmap == TRUE)`.
- Change 5 (F5): color/`df_colors` computation moved inside
  `if (plotheatmap == TRUE)`.
- Change 6 (F6): added `plot.dir <- gl.check.wd(plot.dir, verbose = 0)`
  near the top of the function.
- Change 7 (F7): `@author` now states `Author(s): Arthur Georges.
  Custodian: Arthur Georges -- ...`.
- Change 8 (F8): added an FBM/densification caveat paragraph to
  `@details`.

Snapshot verification: `tests/testthat/test-gl.grm.R` re-run against the
patched function — all 12 expectations pass. The two tests that
characterized F1 and F2's pre-fix behavior were updated to assert the new,
approved behavior (SilicoDArT now errors; `plotheatmap = FALSE` +
`plot.file` now succeeds with a warning instead of crashing) — both diffs
map to approved changes 1 and 2. The baseline SNP-path numeric snapshot
(diagonal/off-diagonal values on `testset.gl[1:10, 1:100]`) is unchanged,
confirming changes 3-8 didn't alter computed output for the default path.
PR: not yet opened.

```json
{
  "function": "gl.grm",
  "package": "dartR.captive",
  "family": "analysis",
  "skill_version": "1.0.0",
  "commit": "bdb623d",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "principle:analysis-dispatch/DAT1", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "principle:scope-leak", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "STY3", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "medium", "rule": "FS7", "status": "approved", "change": 6},
    {"id": "F7", "severity": "INFO", "confidence": "high", "rule": "DOC7 (proposed)", "status": "approved", "change": 7},
    {"id": "F8", "severity": "INFO", "confidence": "medium", "rule": "DAT6 (proposed)", "status": "approved", "change": 8}
  ],
  "coverage_skipped": ["DAT6 full-scale FBM run: no large FBM fixture available"],
  "status": "pr-open",
  "pr": 85
}
```
