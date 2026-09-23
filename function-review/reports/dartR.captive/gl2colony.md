# Review: gl2colony (dartR.captive)
- Family mode: io
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: d5dad54 (origin/dev)
- Datasets: testset.gl[1:20, 1:30] with constructed offspring/father/mother columns, testset.gs, platypus.gl
- External validation: COLONY2 binary `~/programs/colony2s.out` run on each exported file
- Baseline: tests/testthat/test-gl2colony.R (snapshot captured pre-review, 12 expectations)

## Verdict

**Standards: Needs work** — the file layout is correct for the default
settings (COLONY reads it and completes a run), but progress messages
ignore `verbose`, the helper functions shadow `dartR.base::gl2structure`,
and there are no input checks.
**Spec: Rework** — six documented arguments produce files COLONY
rejects, including the function's own `@examples` call; a space in an
individual name silently corrupts that individual's genotypes; and one
missing role column erases the other two.

## Findings

**F1 [HIGH, confidence: high] — names with whitespace corrupt genotypes
silently (DAT5)**
`R/gl2colony.r:301-344` — `write.table()` writes row names unquoted and
space-separated; COLONY reads the first token as the ID.
Failure scenario: `indNames(x)[1] <- "ind 1"`. COLONY runs to completion,
reports the ID as `ind`, and reads the `1` from the name as the first
allele, so every allele of that individual shifts by one position
(`OffGenotype`: `ind,mk1,1,2` where the true genotype is `2 2`). No
error or warning at any step.
Proposed change: stop with an error naming the individuals whose names
contain whitespace.
**Consequence: objects with such names now error instead of exporting.**

**F2 [HIGH, confidence: high] — a missing role column resets all roles
(DOC5, DAT5)**
`R/gl2colony.r:170-179` — if any of `offspring`, `mother`, `father` is
absent, all three are overwritten (`offspring = "yes"`, parents `"no"`).
Failure scenario: a `mother` column marking 5 dams but no `father`
column: output is "20 offspring, 0 fathers, 0 mothers" — the dams become
offspring and COLONY runs with no candidate mothers. Related:
column names are matched case-sensitively here but case-insensitively in
`parental.ids()`, so `Offspring/Father/Mother` prints the "not found,
setting all as offspring" warning and then uses the user's columns
anyway; with no `ind.metrics` at all the call fails with `incorrect number
of dimensions`.
Proposed change: match column names case-insensitively; add only the
missing columns (offspring `"yes"`, parents `"no"`) and name them in the
warning; create `ind.metrics` when it is `NULL`; trim whitespace in
values before matching `"yes"`.
**Consequence: exported roles change for objects that have one or two of
the three columns.**

**F3 [HIGH, confidence: high] — six arguments produce files COLONY
rejects (DOC5)**
`R/gl2colony.r:237-378` — the function writes the single value of each
argument, but COLONY expects extra data blocks for these settings:
- `sibship.prior > 0` needs a mean sibship size line: "Error in reading
  mean paternal & maternal sibship size!"
- `known.allele.freq = 1` needs allele frequencies: "Error in reading
  population allele number per locus"
- `paternity.exclusion.threshold`, `maternity.exclusion.threshold` are
  "number of offspring with known father/mother, threshold"; `"1 0"`
  needs the list: "Errors in DATA"
- `paternal.sibship`, `maternal.sibship`, `excluded.*` > 0 need lists:
  "Error in reading known excluded paternity" (and similar)
- `allelic.dropout`, `other.typ.err`, `marker.id`, `marker.type` without
  a trailing `@` are read as the value of locus 1 only: "Error in reading
  Dropout Rate per locus". The `@examples` call passes
  `allelic.dropout = '0.01'`, so the documented example fails.
COLONY exits with status 0 on these errors, so `gl.run.colony()` does
not detect them (see F8).
Proposed change: stop with a clear error for values that need data the
function cannot write (`sibship.prior > 0`, `known.allele.freq = 1`,
non-zero known/excluded counts); append `@` to a per-locus string that is
a single value with no `@`; fix the example; correct the `@param` text of
the two `*.exclusion.threshold` arguments.
**Consequence: calls with those values error in R instead of writing a
file COLONY rejects.**

**F4 [MEDIUM, confidence: high] — SilicoDArT exported as codominant SNPs
(DAT1)**
`R/gl2colony.r:162,212` — the datatype is detected but never used.
`gl2structure()` maps presence/absence 0/1 to genotypes `1 1`/`1 2`, and
the header declares codominant markers (`0@`).
Failure scenario: `gl2colony(testset.gs, ...)` exports without warning;
COLONY treats every presence as a heterozygote.
Proposed change: stop for SilicoDArT with a message that only SNP data
are supported.
**Consequence: SilicoDArT input errors instead of exporting.**

**F5 [LOW, confidence: high] — messages ignore `verbose` (VRB1, VRB2,
VRB3, FS3)**
`R/gl2colony.r:168,176,196,211,324,336,381` print at `verbose = 0`
(baseline test confirms). The warning has typos ("colums",
"genligth"); line 391 prints "Completed:, "; line 158 passes the outdated
`build = "Jody"`.
Failure scenario: `gl.run.colony(verbose = 0)` or a Shiny session still
prints six progress lines.
Proposed change: gate progress at `verbose >= 2`, results at `>= 3`,
warnings that change the result at `>= 1`; fix the typos; drop `build =`.

**F6 [LOW, confidence: high] — helpers shadow a dartR.base export; sink()
writes (FS1, STY1)**
`R/gl2colony.r:401,431` define `parental.ids()` and `gl2structure()` in
the package namespace. `gl2structure` has the same name as the exported
`dartR.base::gl2structure()` (different signature and output), so inside
dartR.captive the call resolves to the local copy; a reader or a future
edit can mistake one for the other. Lines 292-298, 310-320, 373-378
write through `sink()`; an interrupt inside a sink block leaves console
output diverted into the file.
Proposed change: rename to `utils.colony.parental.ids()` /
`utils.colony.genotypes()` (internal, not exported); write with
`cat(..., file = , append = TRUE)`. No change to the file content.

**F7 [LOW, confidence: high] — documentation gaps (DOC1, DOC6, DOC7)**
`R/gl2colony.r:1-108` — the reference is attributed to "Wang, J. (2011)";
the paper is Jones, O. R. & Wang, J. (2010), Mol Ecol Resour 10:551-555.
No `@family`. `@author` has no `Custodian:` part (proposed rule DOC7).
Non-ASCII en dashes (lines 35, 40, 47, 49, 105) (proposed rule DOC6).
The example is in `\dontrun{}` and fails as written (F3); it can run on
`testset.gl` into `tempdir()`.
Proposed change: fix the citation, add `@family` and the custodian (name
to be confirmed by the team), ASCII dashes, and a runnable example.

**F8 [INFO] — `gl.run.colony()` cannot detect a rejected file**
COLONY prints "Program stopped in subroutine StopOnDataError" and exits
with status 0, so a failing run looks like a successful one. Outside this
function; record for the `gl.run.colony` review. Also seen while verifying: COLONY writes its
output files into the R working directory, not `outpath`.

## Proposed changes

1. Stop with an error when individual names contain whitespace (F1).
   **Consequence: such objects now error instead of exporting.**
2. Role columns: case-insensitive match, add only missing columns, create
   `ind.metrics` if absent, trim values (F2).
   **Consequence: exported roles change for objects with one or two of
   the three columns.**
3. Validate COLONY settings: error for values needing unsupported data
   blocks, append `@` to single-value per-locus strings, correct the
   `*.exclusion.threshold` docs (F3).
   **Consequence: those calls error in R instead of writing a rejected
   file.**
4. Stop for SilicoDArT input (F4).
   **Consequence: SilicoDArT input errors instead of exporting.**
5. Gate messages by `verbose`, fix typos, drop `build =` (F5).
6. Rename the internal helpers and replace `sink()` with `cat(file =)`
   (F6). No change to file content.
7. Documentation: citation, `@family`, custodian, ASCII, runnable example
   (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run (PLT not
  applicable: no plot)
- Spec: behaviour vs roxygen on constructed metadata — run
- External validation: 9 exported files run through COLONY2
  (`colony2s.out`): default settings pass; F1 and F3 cases confirmed from
  COLONY's own output
- Performance: `platypus.gl` (81 individuals, 1,000 loci) exports in
  0.09 s; larger datasets not timed
- FBM path (DAT6): not tested — no FBM fixture; `as.matrix(x)` densifies
  once, which COLONY's text format requires anyway
- Issues: no GitHub issues on dartR.captive or dartRverse mention COLONY;
  Google Group not searched (no access from this session)
- Callers: `gl.run.colony()` (named arguments), dartr2shiny `gl2colony`
  slots (named arguments); no other `dartR.*` package calls it

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis Mijangos | |
| 2 | approved | Luis Mijangos | |
| 3 | approved | Luis Mijangos | including the '@' append for single values |
| 4 | approved | Luis Mijangos | |
| 5 | approved | Luis Mijangos | |
| 6 | approved | Luis Mijangos | |
| 7 | approved | Luis Mijangos | custodian not named in the approval; set to Luis Mijangos (assumption) |

## Outcome

All seven approved changes applied to `R/gl2colony.r` on 2026-09-23;
`devtools::document()` run. `tests/testthat/test-gl2colony.R`: 21
passing. Baseline tests that recorded the defects were replaced by tests
of the approved behaviour; the default-layout baseline passes unchanged.

- File content: the default export of the fixture and of `platypus.gl`
  (81 individuals, 1,000 loci) is byte-identical to the pre-change output
  (`cmp`), so changes 5 and 6 did not alter the file.
- Change 1 (F1): `"ind 1"` now errors naming the individual.
- Change 2 (F2): with no `father` column, the export has 10 offspring and
  "0 5" candidates (was 20 offspring, "0 0"); `Offspring/Father/Mother`
  and a `" YES "` value are matched with no warning; no `ind.metrics`
  exports all individuals as offspring (was `incorrect number of
  dimensions`). IDs now come from `indNames(x)`, not an `id` column.
- Change 3 (F3): `sibship.prior = 2`, `known.allele.freq = 1`,
  `paternity.exclusion.threshold = "1 0"`, and `excluded.paternity = 1`
  each error naming the argument; `allelic.dropout = '0.01'` is written
  as `0.01@`, and COLONY now reads that file with no data error (was
  "Error in reading Dropout Rate per locus").
- Change 4 (F4): `testset.gs` errors "Only SNP data are supported".
- Change 5 (F5): 0 lines printed at `verbose = 0` (was 6).
- Change 6 (F6): helpers renamed `utils.colony.parental.ids()` and
  `utils.colony.genotypes()` (not exported); `sink()` replaced by
  `cat(file =, append = TRUE)`.
- Change 7 (F7): citation, `@family captive management` (adds the
  cross-link to 15 sibling Rd files), `Author(s)`/`Custodian`, ASCII
  dashes, runnable example (runs via `example("gl2colony")`).
- End-to-end: `gl.run.colony()` with `~/programs/colony2s.out` on the
  fixture at `verbose = 3` completes and writes COLONY output.
- Not in this change: `devtools::document()` on `origin/dev` also drops
  `importFrom(stats, dnorm)` and `importFrom(stats, qnorm)` from
  NAMESPACE (no source declares them); reverted here to keep the PR to
  one function.

```json
{
  "function": "gl2colony",
  "package": "dartR.captive",
  "family": "io",
  "skill_version": "2.0.0",
  "commit": "d5dad54",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DAT5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DAT1", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "FS1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 7},
    {"id": "F8", "severity": "INFO", "confidence": "high", "rule": "none (gl.run.colony)", "status": "noted", "change": null}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "Google Group search: no access"],
  "status": "pr-open",
  "pr": null
}
```
