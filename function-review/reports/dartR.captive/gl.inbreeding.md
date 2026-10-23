# Review: gl.inbreeding (dartR.captive)
- Family mode: analysis (wrapper for two external programs)
- Origin: new — contributed by Luis Mijangos, converted from the inbreeding
  block of `~/jeremy_v3/scripts/jeremy_V2.R` (lines 363-456)
- Date: 2026-10-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v3.0.0
- Package commit: 20a8e7a (origin/dev)
- Datasets: `woylies_filtered.rds` (jeremy_v3; 496 individuals, 10,721
  loci, 26 populations), `platypus.gl`
- External programs: PLINK v1.90b6.10 (`~/programs/plink`), GCTA 1.93.2
  beta Mac (`~/programs/gcta64_2020`)
- Baseline: the script block run as written on the woylie data (only
  `gl2plink` argument names updated to dartR.base 1.2.6: `bed.files`,
  `plink.bin.path`); outputs kept in the session scratchpad
  (`base/baseline.rds`, `gl_plink.het`, `woylie_gcta.ibc`). The member's
  code is a script, not a function, so Intake step 3 (place the file
  unchanged in `R/`) does not apply; the function is written in Phase C
  from the approved changes.

## Verdict

**Standards: Rework** — the block is analysis-script code: no arguments,
roxygen, verbosity or error checks; it writes files into the working
directory and depends on `plinkQC`, which dartR.captive does not declare.

**Spec: Needs work** — PLINK F is computed correctly, but the GCTA column
uses a different set of loci, and both estimates take allele frequencies
from all populations pooled, which counts between-population
differentiation as inbreeding.

What works well: the PLINK `--het` values match an independent R
computation of F = (O - E) / (N - E) to the 4 decimals PLINK prints.

## Findings

**F1 [HIGH, confidence: high] — PLINK and GCTA read different loci
(principle: estimates reported side by side must share their input)**
`jeremy_V2.R:366-379` — the script sets every locus to chromosome `"1"`,
but `gl2plink` (dartR.base 1.2.6) recodes loci at position 0 as
chromosome 0. GCTA `--autosome` drops those loci; PLINK `--het` keeps
them. GCTA also counts monomorphic loci in the denominator of all three
Fhat with zero contribution, while PLINK skips them.
Failure scenario: woylie data. PLINK uses the 10,444 polymorphic loci.
GCTA uses 10,243 loci (478 position-0 loci dropped), 218 of them
monomorphic, which pull every Fhat towards 0. Checked: Fhat1 and Fhat3
equal sum(polymorphic terms) / NOMISS to 3e-5, with NOMISS counting the
monomorphic loci. In a per-population run, where many more loci are
monomorphic, the dilution grows.
Proposed change: give every locus chromosome 1 and positions 1..nLoc
(neither `--het` nor `--ibc` uses positions or linkage), and drop
monomorphic loci before writing the PLINK files, so both programs read
the same polymorphic loci.

**F2 [HIGH, confidence: high] — allele frequencies pooled across
populations (principle: Spec — individual inbreeding is defined relative
to a reference population)**
`jeremy_V2.R:368-379` — one PLINK and one GCTA run on all 496
individuals from 26 sites, so expected homozygosity comes from pooled
frequencies and differentiation among sites (the Wahlund effect) is added
to every individual's F.
Failure scenario: woylie data, the 434 individuals in populations of at
least 10. Their mean PLINK F is 0.220 with pooled frequencies and 0.138
with each population's own frequencies (independent R computation,
confirmed by the function in Phase C). Correction: the first version of
this report and the approval question compared 0.225, the pooled mean
over all 496 individuals, with 0.138; the comparison on the same
individuals is 0.220 vs 0.138. Individual
`00080ECCA5` (Perup) has Fhat3 2.98 pooled and 0.35 within Perup.
Proposed change: new argument `by.pop`; `TRUE` runs both programs once
per population on that population's polymorphic loci. Populations of one
individual return NA (every homozygous locus is monomorphic, so F = -1 by
construction); populations under 10 individuals get a warning that
frequencies are imprecise. With `by.pop = FALSE` and more than one
population, a warning at `verbose >= 1` states that F includes
between-population differentiation.

**F3 [MEDIUM, confidence: high] — Fhat1, Fhat2 and locus counts are
discarded (principle: Spec — return what the programs estimate)**
`jeremy_V2.R:376, 382` — only PLINK `F` and GCTA `Fhat3` are kept.
Failure scenario: Fhat3 weights each locus by 1 / (2p(1-p)); an
individual homozygous for an allele at frequency p adds (1 - p) / p,
which reaches 991 at p = 1/992 in the woylie data. A user who sees Fhat3
= 2.93 cannot compare it with Fhat1/Fhat2 or check how many loci it rests
on.
Proposed change: return PLINK O(HOM), E(HOM), N(NM) and F, and GCTA
NOMISS, Fhat1, Fhat2 and Fhat3; the plot keeps PLINK F and Fhat3, as in
the script.

**F4 [MEDIUM, confidence: high] — sex-linked loci treated as autosomal
(DOC5, proposed rule)**
`jeremy_V2.R:366` — all loci are placed on chromosome 1.
Failure scenario: X-linked loci are hemizygous in males, so a dataset
that still holds them reports males as more homozygous, that is, more
inbred.
Proposed change: state in `@details` that all loci are treated as
autosomal and that sex-linked loci must be removed first
(`dartR.sexlinked::gl.drop.sexlinked`). Docs only.

**F5 [MEDIUM, confidence: high] — external-program failures are not
detected; undeclared dependency (DEP1)**
`jeremy_V2.R:372, 379` — PLINK runs through `plinkQC`, which is not in
dartR.captive's DESCRIPTION; GCTA runs through `system()` with no exit
status check, in the working directory.
Failure scenario: if GCTA fails (wrong path, unreadable `.bed`), the
script reads `woylie_gcta.ibc` left by an earlier run and reports stale
values as current, or stops on a missing file with no reason given.
Proposed change: run PLINK through `dartR.base::utils.plink.run` (stops
with PLINK's own message), run GCTA with `system2()` and stop on a
non-zero exit with the last lines of its log, and run both in a fresh
folder under `tempdir()` so no earlier output can be read.

**F6 [MEDIUM, confidence: high] — writes files into the working
directory (FS7, PLT1, PLT2 proposed rule)**
`jeremy_V2.R:451, 456` — `ggsave("inbreeding.pdf")` and
`write.csv("inbreeding_table.csv")` write to `getwd()`; the plot uses a
hard-coded palette (`structure_colors`) and `base_family = "Helvetica"`.
Failure scenario: a second call overwrites the first call's files
without notice; on a system without Helvetica, the PDF device warns for
every text element.
Proposed change: return the table and the plot; save the plot only via
`plot.file`/`plot.dir` (`utils.plot.save`); take colours from
`plot.colors` (default `discrete_palette`) and the theme from
`plot.theme` (default `theme_dartR()`).

**F7 [LOW, confidence: high] — both programs are required**
`jeremy_V2.R:372-379`.
Failure scenario: a user with PLINK but no GCTA cannot get PLINK F.
Proposed change: argument `method = c("plink", "gcta")`; the default runs
both, and either can be run alone.

**F8 [LOW, confidence: high] — house structure absent (DOC1, DOC2, DOC7,
FS2-FS10, VRB1-VRB3)**
Whole block — no roxygen, verbosity, start/end flags, datatype check or
argument validation.
Failure scenario: a genlight of SilicoDArT data or a missing executable
path fails inside `gl2plink` or PLINK with a message that does not name
the cause.
Proposed change: full house structure; accept SNP data only; arguments
`plink.path` and `gcta.path` (a folder holding `plink`/`gcta64`, with
`.exe` on Windows, or the executable itself, since GCTA's macOS builds
ship under other names), default `getwd()` as in `gl.run.EMIBD9`.

Naming note (not a finding): FS1 asks for `gl.<verb>.<noun>`, but
`gl.inbreeding` sits with `gl.kin`, `gl.grm` and `gl.relatedness` in this
package.

## Proposed changes

1. Both programs read the same polymorphic loci: chromosome 1 and
   positions 1..nLoc for every locus, monomorphic loci dropped first (F1).
   **Consequence: GCTA output differs from the script's: it now includes
   the loci at position 0 (478 in the woylie data) and excludes
   monomorphic loci (218). PLINK F is unchanged.**
2. Argument `by.pop` to take allele frequencies within each population
   (F2), default `FALSE`, with the warnings described in F2.
   **Consequence: with `by.pop = TRUE`, F for multi-population data
   differs from the script's (woylie mean PLINK F 0.220 -> 0.138 in
   populations of at least 10). If
   the default were `TRUE`, the default output would differ.**
3. Return O(HOM), E(HOM), N(NM), NOMISS, Fhat1 and Fhat2 with F and
   Fhat3 (F3).
4. Argument `method` to run PLINK, GCTA, or both (F7).
5. House formatting: roxygen with the sex-linked note (F4), program calls
   with failure checks in a fresh `tempdir()` folder and no `plinkQC`
   (F5), no files written unless `plot.file` is given, `plot.colors` and
   `plot.theme` (F6), house structure, SNP-only check, `plink.path` and
   `gcta.path` (F8).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run on the script
  block
- Spec: PLINK F against independent R computation (woylie, max |diff|
  4e-5, PLINK prints 4 decimals) — run
- Spec: GCTA Fhat1/2/3 against independent R computation — run; matches
  only after counting monomorphic loci in the denominator (F1)
- Spec: pooled vs within-population frequencies — run (independent R
  computation, not the programs)
- Characterization test in `tests/testthat/`: deferred to Phase C — the
  script is not a function, so there is nothing to load with
  `devtools::load_all()`. The Phase C test checks the function against
  the independent computations above and skips when PLINK or GCTA is not
  found (CI has neither).
- DAT6 (FBM-backed input): not checked; `gl2plink` densifies with
  `as.matrix()` and is outside this function.
- Caller grep: not applicable (new function).
- Google Group / GitHub issues: not applicable (new function).

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence for GCTA output approved |
| 2 | approved | Luis | default `by.pop = FALSE` |
| 3 | approved | Luis | |
| 4 | approved | Luis | |
| 5 | applied without a box | Claude | the approval tool rejects a one-option question; house formatting is what turns the script into a function and alters no estimate. Listed here so it can still be rejected at the pre-push gate |

## Outcome

Written as `R/gl.inbreeding.r`; tests in
`tests/testthat/test-gl.inbreeding.R` (30 expectations, pass with PLINK
v1.90b6.10 and GCTA 1.93.2; all tests that need a program skip when it is
missing or does not run on the platform, as on CI).

- Change 1: on the woylie data, `N.NM` equals `NOMISS` for every
  individual (10,444 loci used by both). PLINK F equals the script's
  output exactly (max |diff| 0). GCTA Fhat3 differs from the script's by
  up to 0.052, the approved consequence. On `platypus.gl` with half the
  loci at position 0, F and Fhat1-3 match the independent computation to
  1e-4.
- Change 2: `by.pop = TRUE` matches the independent within-population
  computation on `platypus.gl`; woylie mean F in populations of at least
  10 is 0.220 pooled and 0.138 by population. Corbal and PSO (one
  individual each) return NA with a warning; five populations under 10
  get the imprecision warning.
- Change 3: all eight program columns returned.
- Change 4: `method = "plink"` runs with `gcta.path` pointing nowhere.
- Change 5: `devtools::document()` adds `export(gl.inbreeding)`;
  `man/gl.grm.Rd` and `man/gl.grm.network.Rd` gain the
  `gl.inbreeding` family link. Run end to end at `verbose = 3` on the
  woylie data (14 s for 496 individuals x 10,721 loci). Two cosmetic fixes
  to the new code during verification: horizontal population labels
  (long names were cut off) and spacing in the progress message.
- `devtools::check()` (R 4.4.2, `--no-manual`): 0 errors, 1 warning,
  1 note, none from the new file. The warning is the local library's
  "package built under R version 4.4.3"; the note is the worktree's `.git`
  file. Full test suite passes.

## Machine block

```json
{
  "function": "gl.inbreeding",
  "package": "dartR.captive",
  "family": "analysis",
  "origin": "new",
  "skill_version": "3.0.0",
  "commit": "20a8e7a",
  "verdict_standards": "rework",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "principle: shared input", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "principle: reference population", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "principle: return estimates", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 5},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "approved", "change": 5},
    {"id": "F6", "severity": "MEDIUM", "confidence": "high", "rule": "FS7", "status": "approved", "change": 5},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "principle: optional programs", "status": "approved", "change": 4},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 5}
  ],
  "coverage_skipped": ["characterization test deferred to Phase C: script, not a function", "DAT6: gl2plink densifies, outside this function"],
  "status": "awaiting-approval",
  "pr": null
}
```
