# Review: gl.run.EMIBD9 (dartR.captive)
- Family mode: io (wrapper for an external program)
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 3f54843 (origin/dev)
- Datasets: testset.gl[1:10, 1:120] and [1:10, 121:250] (gl.filter.allna), testset.gs[1:10, 1:100]
- External program: EMIBD9 1.0.0.0 build 20211115 (`~/programs/emibd9-v1.0/EM_IBD_P`, macOS, serial)
- Baseline: tests/testthat/test-gl.run.EMIBD9.R (snapshot captured pre-review, 10 expectations; skips without EMIBD9)

## Correction to earlier work (gl.grm.network, PR #97, merged 2026-09-23)

The `gl.grm.network` review (round 2) stated that EMIBD9's `$rel` is on
the relatedness scale (2 x kinship) and approved halving `G`. That is
wrong. EMIBD9's `r(1,2)` column equals the kinship coefficient:
self-comparisons are 0.5 (not 1), and for every pair
`r(1,2) = D1 + (D3 + D5 + D7)/2 + D8/4` (e.g. D7 = 0.1145, D8 = 0.4059
gives 0.1587, the reported value). The baseline test checks this for all
45 pairs. With PR #97 as it stands, the Shiny app, which passes
`gl.run.EMIBD9()$rel` to `gl.grm.network()`, would plot EMIBD9 kinship
at half its value. `gl.grm()` output is unaffected (it is on the
relatedness scale). See change 8.

## Verdict

**Standards: Needs work** — the name sanitising added in August works,
but the exit status is ignored, EMIBD9 prints about 60 lines at every
verbosity, and the error for a missing executable is empty.
**Spec: Rework** — a failed run returns the previous run's results; a
documented argument value (`OutAlleleFre = TRUE`) crashes EMIBD9; and
`rel` is documented as relatedness but holds kinship.

## Findings

**F1 [HIGH, confidence: high] — a failed run returns the previous run's
results (FS6)**
`R/gl.run.EMIBD9.r:155,273,286,290` — every call runs in the shared
session `tempdir()`, `system(cmd)` is not checked, and the output file is
read whether or not this run wrote it. EMIBD9 exits with status 0 even
when it stops ("Parameter file ... does not exist!", exit 0).
Failure scenario: run on `testset.gl[1:10, 1:120]`, then on loci
`121:250` with a setting that crashes EMIBD9 (F2). The second call
returns normally with a `rel` matrix identical to the first run's
(`identical(r1$rel, r2$rel)` is TRUE). If the individuals differ, the
call fails later with an unrelated `dimnames` error instead.
Proposed change: run each call in its own new folder
(`tempfile("EMIBD9_")`); check the exit status; stop with EMIBD9's last
console lines if the output file was not written.

**F2 [HIGH, confidence: high] — `OutAlleleFre = TRUE` crashes EMIBD9
(DOC5)**
`R/gl.run.EMIBD9.r:127,254` — documented as "A boolean", but the value
is pasted into `MyData.par` as `TRUE`; EMIBD9 expects 0/1 and stops with
"SIGSEGV: Segmentation fault". Combined with F1, the call then returns
old results.
Proposed change: write `as.integer(isTRUE(OutAlleleFre))`, as already
done for `Inbreed`.

**F3 [MEDIUM, confidence: high] — `rel` holds kinship, documented as
relatedness (DOC5)**
`R/gl.run.EMIBD9.r:94,359-360` — `@return` and the verbose message call
`rel` "pairwise relatedness"; the same message also says "kinship". The
values are the kinship coefficient (see Correction above); the diagonal
is 0.5. The Shiny app labels it "EMIBD9 relatedness".
Failure scenario: a user compares `rel` with a relatedness threshold
(0.5 for full sibs) and misses every full-sib pair, whose kinship is
0.25.
Proposed change: document `rel` as the kinship coefficient (theta), with
the diagonal meaning, and fix the message. No change to the numbers.

**F4 [MEDIUM, confidence: high] — SilicoDArT accepted (DAT1)**
`R/gl.run.EMIBD9.r:151,260-265` — the datatype is detected but not used;
presence/absence 0/1 is written as SNP genotypes 0/1 and EMIBD9 returns
kinship values that have no meaning for dominant markers.
Proposed change: stop for SilicoDArT with "only SNP data are supported".
**Consequence: SilicoDArT input errors instead of returning numbers.**

**F5 [MEDIUM, confidence: high] — numbers returned as text (DOC5)**
`R/gl.run.EMIBD9.r:299-318` — every column of `raw` and `processed` is
character (`"0.1587"`), and `raw$#IIS1` is a list column.
Failure scenario: `mean(res$processed$Delta9)` returns `NA` with a
warning; `res$processed[order(-"r(1,2)")]` sorts as text.
Proposed change: convert all columns except the individual IDs to
numeric.
**Consequence: column types change from character to numeric; no known
caller reads them (dartr2shiny uses only `$rel`).**

**F6 [MEDIUM, confidence: high] — genotype file build is quadratic in
loci (STY2)**
`R/gl.run.EMIBD9.r:263-265` — `Reduce(paste0, y)` rebuilds the string
at every locus: 0.84 s per individual at 20,000 loci against 0.001 s
for `paste(y, collapse = "")`; about 5 s per individual at 50,000 loci,
so 1,000 individuals spend over an hour before EMIBD9 starts.
Proposed change: `paste(y, collapse = "")`. Same file content.

**F7 [LOW, confidence: high] — `plot.file` without `plot.out` errors
(PLT3)**
`R/gl.run.EMIBD9.r:370-396` — the plot is built only when `plot.out`
is TRUE, but saved whenever `plot.file` is set.
Failure scenario: `plot.out = FALSE, plot.file = "p"` runs EMIBD9, then
fails with `object 'p1' not found`, losing the results.
Proposed change: build the plot when either is requested.

**F8 [LOW, confidence: high] — messages, preconditions, docs (FS3, FS5,
VRB3, DOC1, DOC6, DOC7)**
- EMIBD9's console output (~60 lines per run) prints at `verbose = 0`.
- A missing executable prints the reason with `message()` and then calls
  `stop()` with an empty message (the error text is blank).
- The summary message says the list contains "the input gl object" (it
  does not) and prints at `verbose >= 1`, not 3.
- `inbreeding` is always returned (EMIBD9 writes the table even with
  `Inbreed = FALSE`); `@return` says "if requested".
- `outfile` is documented as "path and name", but a path breaks the run
  (it is resolved inside the run folder).
- `build = "Jody"`; no `@family`; `@author` has no `Author(s):` part;
  non-ASCII table characters; typos ("vakue", "emidb").
Proposed change: gate EMIBD9 output at `verbose >= 2`; `stop(error(...))`
naming the missing file; fix the summary message and gate it at 3;
document `inbreeding` as always present; document `outfile` as a file
name; the doc items.

**F9 [INFO] — not tested: parallel (MPI) path and Windows**
`parallel = TRUE` runs `mpirun ... --use-hwthread-cpus`, an Open MPI
flag; the Shiny app always sets `parallel = TRUE`. No `EM_IBD_P_mpi` or
Windows binary was available here. The app also stores `rel` in its
distance-matrix list (`MyDis`), although kinship is a similarity; worth
a look on the dartr2shiny side.

## Proposed changes

1. Isolated run folder, exit-status check, stop with EMIBD9's message
   when no output is written (F1).
2. Write `OutAlleleFre` as 0/1 (F2).
3. Document `rel` as kinship; fix the message (F3). Docs only.
4. Stop for SilicoDArT (F4).
   **Consequence: SilicoDArT input errors instead of returning numbers.**
5. Numeric columns in `raw` and `processed` (F5).
   **Consequence: column types change from character to numeric.**
6. `paste(collapse = "")` for the genotype file (F6).
7. Plot/`plot.file` coupling, messages, preconditions, docs (F7, F8).
8. Correction for PR #97 (gl.grm.network): stop halving EMIBD9 kinship.
   Options: (a) `gl.run.EMIBD9()` tags `rel` with
   `attr(rel, "scale") <- "kinship"` and `gl.grm.network()` skips the
   halving for tagged matrices (no Shiny change); (b) a new
   `gl.grm.network()` argument `G.scale = c("relatedness", "kinship")`,
   set by the Shiny app for the EMIBD9 container.
   **Consequence: `gl.grm.network()` on EMIBD9 input plots kinship at
   its true value; `gl.grm()` input is unchanged.**

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: real EMIBD9 runs on macOS (serial) — run; F1, F2, F4, F5, F7
  reproduced; kinship identity checked on 45 pairs
- Performance: `Reduce(paste0)` vs `paste(collapse)` timed on 20,000 loci
- Parallel/MPI and Windows: SKIPPED — no binaries (F9)
- Callers: dartr2shiny (`$rel` only, `parallel = TRUE`, `outpath =
  global$temp.dir`); no other `dartR.*` package calls it
- FBM path (DAT6): not tested; `as.matrix(x2)` densifies once, as the
  text genotype file requires

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis Mijangos | |
| 2 | approved | Luis Mijangos | |
| 3 | approved | Luis Mijangos | |
| 4 | approved | Luis Mijangos | |
| 5 | approved | Luis Mijangos | |
| 6 | approved | Luis Mijangos | |
| 7 | approved | Luis Mijangos | `@author` lists Luis Mijangos as author (assumption: git history does not show the original author) |
| 8 | approved | Luis Mijangos | option (a): tag `rel` with `attr(, "scale") = "kinship"`. PR #97 was merged before this correction, so the `gl.grm.network` half goes in its own PR from `dev` |

## Outcome

All eight approved changes applied on 2026-09-23 (change 8 split across
this PR and a `gl.grm.network` PR). `tests/testthat/test-gl.run.EMIBD9.R`:
7 tests pass with EMIBD9 1.0.0.0 on macOS (they skip without the binary);
the pre-change baseline tests were replaced by tests of the approved
behaviour.

- Numbers unchanged: on `testset.gl[1:10, 1:200]`, `rel` (ignoring the new
  attribute), all `processed` values and `inbreeding` are identical to the
  pre-change run, so change 6 (`paste(collapse = "")`) writes the same
  genotype file.
- Change 1 (F1): runs in `tempfile("EMIBD9_")`; a run where EMIBD9 fails
  (`EM_Method = "bad"`) stops with "EMIBD9 did not write its results" and
  EMIBD9's last console lines (was: returned the previous call's results).
  No executables are left in the session temp folder.
- Change 2 (F2): `OutAlleleFre = TRUE` is written as 1; EMIBD9 runs and
  the results differ from the previous call's.
- Change 3 (F3): `@return` documents `rel` as kinship with the 0.5 (1 + F)
  diagonal; the summary message says kinship.
- Change 4 (F4): `testset.gs` errors "Only SNP data are supported".
- Change 5 (F5): all `raw`/`processed` columns except the IDs are numeric;
  the list column is gone.
- Change 7 (F7, F8): `plot.file` with `plot.out = FALSE` saves the plot
  (drawn on a null device); EMIBD9 console output only at `verbose >= 2`
  (otherwise kept in `EMIBD9_console.txt` in the run folder and shown on
  failure); a missing executable errors "Cannot find EM_IBD_P in the
  folder given by emibd9.path"; the output file is copied to `outpath`
  from the run folder; docs (`@family`, author, ASCII, typos, `inbreeding`
  always returned, `outfile` is a file name).
- Change 8 (a): `rel` carries `attr(rel, "scale") = "kinship"`. Checked
  end to end with the `gl.grm.network` change: EMIBD9 kinship stored in a
  list as the Shiny app does and passed to `gl.grm.network()` comes out
  equal to `r(1,2)` (was: halved on `dev` since PR #97).
- Not in this change: the NAMESPACE drift on `dev` (`dnorm`, `qnorm`),
  reverted as in PR #98.

```json
{
  "function": "gl.run.EMIBD9",
  "package": "dartR.captive",
  "family": "io",
  "skill_version": "2.0.0",
  "commit": "3f54843",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "FS6", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT1", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "STY2", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "PLT3", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 7},
    {"id": "F9", "severity": "INFO", "confidence": "medium", "rule": "none", "status": "noted", "change": null}
  ],
  "coverage_skipped": ["parallel/MPI and Windows: no binaries", "DAT6: no FBM fixture"],
  "status": "pr-open",
  "pr": null
}
```
