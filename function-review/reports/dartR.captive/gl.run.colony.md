# Review: gl.run.colony (dartR.captive)
- Family mode: io (wrapper for an external program)
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 3f54843 (origin/dev, includes the gl2colony fixes from PR #98)
- Datasets: testset.gl[1:20, 1:30] with constructed offspring/father/mother columns
- External program: COLONY 2.0.7.1 (`~/programs/colony2s.out`, macOS)
- Baseline: tests/testthat/test-gl.run.colony.R (snapshot captured pre-review; skips when COLONY is not installed)

## Verdict

**Standards: Needs work** — the exit-status check (added for OP 8821)
is in place, but COLONY's full console output prints at every verbosity,
the command is built without quoting, and the parameter docs are a stale
copy of the ones corrected in `gl2colony`.
**Spec: Rework** — the default call cannot work (`outpath = NULL` sends
COLONY an empty file name), and the success check accepts a failed run
whenever output files from an earlier run sit in the working directory.

## Findings

**F1 [HIGH, confidence: high] — the default `outpath = NULL` always fails
(DOC5)**
`R/gl.run.colony.r:146-147` — `gl.check.wd()` is commented out, so
`outfilespec <- file.path(NULL, outfile)` is `character(0)`; `gl2colony()`
writes the file to `tempdir()`, but the command becomes `colony2s.out IFN:`.
Failure scenario: `gl.run.colony(x, colony.path = "~/programs")` — COLONY
prints "Input file does not exist!" and the function stops with
"COLONY finished without writing any output files", which points the
user at the wrong cause.
Proposed change: use the file path that `gl2colony()` returns (it already
resolves `outpath` with `gl.check.wd()`).

**F2 [HIGH, confidence: high] — failed runs pass the success check**
`R/gl.run.colony.r:221-235` — COLONY exits with status 0 when it rejects
its input, and the only other check is that some file starting with
`output.name` exists in `getwd()`, from any date.
Failure scenario: a user runs COLONY once, then again with a changed
object or path. The second run fails ("Program stopped in subroutine
ReadData") but the first run's `my_project.*` files satisfy the check, so
the function returns normally and the user reads the old results.
Reproduced with `outpath = NULL` and with an `outpath` containing a space.
Proposed change: record the start time; stop with the text of
`Colony2.ErrorMessage` if COLONY wrote one during the run; otherwise
require `<output.name>.BestConfig` modified after the start time.

**F3 [HIGH, confidence: high] — paths containing spaces break the run**
`R/gl.run.colony.r:209-219` — the command is pasted unquoted, so the
shell splits `colony.path` and `outfilespec` at spaces.
Failure scenario: `outpath = "~/Google Drive/colony"` — COLONY reports
"Input file .../Google does not exist!"; combined with F2 the call can
return as if it succeeded. Quoting the argument
(`"IFN:/a b/colony2.dat"`) was verified to work with COLONY 2.0.7.1.
Proposed change: build the command with `shQuote()` on the executable and
the `IFN:` argument.

**F4 [MEDIUM, confidence: high] — COLONY output goes to the R working
directory, not `outpath`**
`R/gl.run.colony.r:209-229` — COLONY writes about 36 files into its
working directory, which is the user's R working directory; `outpath`
holds only the input file, yet the function returns `outpath`.
Failure scenario: running from a project folder leaves 36 `my_project.*`
files next to the user's scripts; the returned path contains no results.
Proposed change: run COLONY with its working directory set to `outpath`
(restored on exit).
**Consequence: output files move from `getwd()` to `outpath`. dartr2shiny
zips COLONY output from `getwd()` (`gl.run.colony` `render_output` slot)
and passes `outpath = global$temp.dir`, so its download handler needs
`list.files(path = global$temp.dir, ...)` in the same release.**

**F5 [MEDIUM, confidence: high] — return value does not match the docs
(DOC5)**
`R/gl.run.colony.r:83-84,249` — `@return` says "the output filename";
the function returns `outpath` (`NULL` by default), and no COLONY result
is read back into R.
Failure scenario: `res <- gl.run.colony(...)` gives `NULL` or a folder
path; the user must find and parse COLONY's files by hand.
Proposed change (member's choice): (a) return, invisibly, the full paths
of the output files COLONY wrote in this run, and document it; or (b)
return a list with those paths plus `best.config`, the
`<output.name>.BestConfig` table (offspring, assigned father, mother, and
clusters) read as a data frame.
**Consequence: the return value changes from a folder path to a vector or
list; no known caller uses it (dartr2shiny ignores it).**

**F6 [LOW, confidence: high] — messages and preconditions (VRB3, FS3,
FS5)**
Lines 209-219 let COLONY print hundreds of progress lines at every
verbosity, including `verbose = 0`. A missing executable surfaces as
"exit status 127" after a shell error, not as a message naming the
expected file. Line 246 prints "Completed:, "; line 152 passes the
outdated `build = "Jody"`.
Proposed change: check that the executable exists before running and name
the expected file; send COLONY's console output to the console only at
`verbose >= 2`; fix the message; drop `build =`.

**F7 [LOW, confidence: high] — documentation (DOC1, DOC5, DOC6, DOC7)**
The 33 shared `@param` entries repeat the pre-#98 `gl2colony` text, so
they now contradict it (e.g. `sibship.prior` "0-4" here, "only 0" there).
The example calls `gl2colony()`, not `gl.run.colony()`. Citation
"Wang, J. (2011)" should be Jones & Wang (2010). Backticks around
`Colony2.DAT` print literally (roxygen markdown is off). No `@family`;
`@author` has no custodian. Non-ASCII hyphens and dashes (lines 4, 5, 38,
43, 50, 52, 80, 98).
Proposed change: `@inheritParams gl2colony` for the shared parameters,
correct example (kept in `\dontrun{}`, it needs the COLONY binary),
citation, `\code{}`, `@family captive management`, custodian, ASCII.

**F8 [INFO] — Windows and Linux executable names not verified**
The function expects `Colony2p.exe` (Windows) and `colony2s.ifort.out`
(Linux). Only the macOS binary was available here. If a COLONY release
ships different names, the F6 existence check will at least name the
file it looked for.

## Proposed changes

1. Use the path returned by `gl2colony()` so the default `outpath` works
   (F1).
2. Detect failed runs: `Colony2.ErrorMessage` from this run, or no fresh
   `<output.name>.BestConfig`, stops with COLONY's message (F2).
3. Quote the executable and `IFN:` argument (F3).
4. Run COLONY inside `outpath` so its output lands there (F4).
   **Consequence: output moves from `getwd()` to `outpath`; dartr2shiny's
   download handler must read from `global$temp.dir` in the same
   release.**
5. Return value: (a) output file paths, or (b) paths plus the
   `BestConfig` table (F5).
   **Consequence: return changes from a folder path to a vector/list.**
6. Executable check, verbose-gated COLONY output, message fix, drop
   `build =` (F6).
7. Documentation via `@inheritParams gl2colony` and the fixes listed in
   F7.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (PLT, DAT
  genotype checks delegated to `gl2colony`, reviewed in PR #98)
- Spec: behaviour vs roxygen with COLONY 2.0.7.1 on macOS — run; F1-F4
  reproduced
- Windows/Linux: SKIPPED — no binaries on this machine (F8)
- Callers: dartr2shiny `gl.run.colony` slots (named arguments,
  `outpath = global$temp.dir`, output collected from `getwd()`); no other
  `dartR.*` package calls it
- Issues: OP 8821 (silent failure) is the origin of the current exit-status
  check; F2 is the remaining gap in that fix

## Addendum (found while applying)

**A1 [MEDIUM, confidence: high] — COLONY truncates IDs to 20 characters
(DAT5)**
A 51-character individual name comes back in `BestConfig` as its first 20
characters. `best.config` IDs then do not match `indNames(x)`, and two
individuals whose names share the first 20 characters are merged by
COLONY.
Decision (Luis Mijangos): keep the names by replacing them for the run.
Applied in `gl.run.colony()`: if any name is longer than 20 characters or
contains whitespace, COLONY runs on `ind1..indN`; the original names are
restored in `best.config` (inferred parents `*k`/`#k` keep their labels)
and `<output.name>.IDmap.csv` in `outpath` maps the IDs in COLONY's own
files. Names within the limit are used as-is and no map is written.
Evidence: two offspring sharing their first 20 characters plus a father
named "father eleven" run and come back with their full names; a
30-character name is restored at `verbose = 2`; 18 tests pass.
`gl2colony()` used on its own still writes the names unchanged, so a user
running COLONY by hand still meets the limit.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis Mijangos | |
| 2 | approved | Luis Mijangos | |
| 3 | approved | Luis Mijangos | |
| 4 | approved | Luis Mijangos | "Apply + fix Shiny": dartr2shiny slot updated in a separate MR |
| 5 | approved | Luis Mijangos | option (b): paths + BestConfig table |
| 6 | approved | Luis Mijangos | |
| 7 | approved | Luis Mijangos | |

## Outcome

All seven approved changes applied to `R/gl.run.colony.r` on 2026-09-23;
`devtools::document()` run. `tests/testthat/test-gl.run.colony.R`: 18
passing (including A1) (COLONY 2.0.7.1, macOS); `test-gl2colony.R` still 21 passing.
The pre-change baseline tests were replaced by tests of the approved
behaviour.

- Change 1 (F1): the default call runs; the wrapper uses the path
  returned by `gl2colony()`.
- Change 2 (F2): success requires no new `Colony2.ErrorMessage` and a
  `<output.name>.BestConfig` written by this run. A broken input file with
  an hour-old `BestConfig` in `outpath` now errors (was: returned
  normally). COLONY run directly on that broken file exits with status 0,
  confirming the exit status alone cannot be trusted.
- Change 3 (F3): `outpath = ".../out dir"` runs end to end. `system2()`
  quotes the command itself in R 4.4, so only the `IFN:` argument is
  passed through `shQuote()` (quoting both gave exit status 127).
- Change 4 (F4): COLONY runs inside `outpath` (working directory restored
  on exit, checked in the test); 36 output files in `outpath`, none in
  the caller's working directory. dartr2shiny: `gl.run.colony`
  `render_output` slot reads from `global$temp.dir` with the pattern
  anchored to `^<output.name>\.`, so the zip cache is not re-zipped;
  `generate/main.R` changes only `Fun_gl.run.colony.R`; `tools/verify.R`
  321/321 identical; two simulated download requests each serve 36 files.
- Change 5 (F5, option b): returns `list(files, best.config)`;
  `best.config` read with `comment.char = ""` because COLONY labels
  inferred mothers `#1`, `#2`.
- Change 6 (F6): a missing executable errors "COLONY executable not
  found: <path>"; 0 lines printed at `verbose = 0`; `build =` removed.
- Change 7 (F7): `@inheritParams gl2colony`, example calls
  `gl.run.colony()`, citation, `\code{}`, `@family captive management`,
  custodian, ASCII (author names keep their accents).
- End to end at `verbose = 3` on the fixture: completes, 36 files,
  `best.config` 10 rows.
- Not in this change: the NAMESPACE drift on `dev` (`dnorm`, `qnorm`),
  reverted as in PR #98.
- PR #99 (`review-gl.run.colony` -> `dev`), commit 1f0ced0; dartr2shiny
  MR !54 (`fix/colony-output-dir` -> `fix/sexlinked-filter-migration`).

```json
{
  "function": "gl.run.colony",
  "package": "dartR.captive",
  "family": "io",
  "skill_version": "2.0.0",
  "commit": "3f54843",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "FS6", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "FS6", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "FS7", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 7},
    {"id": "A1", "severity": "MEDIUM", "confidence": "high", "rule": "DAT5", "status": "approved", "change": "A1"},
    {"id": "F8", "severity": "INFO", "confidence": "low", "rule": "none", "status": "noted", "change": null}
  ],
  "coverage_skipped": ["Windows/Linux executables: no binaries available"],
  "status": "done",
  "pr": 99
}
```
