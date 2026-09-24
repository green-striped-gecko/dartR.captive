# Review: gl.kin (dartR.captive)
- Family mode: analysis
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 573a30f (origin/dev, includes the gl.run.EMIBD9 fixes from PR #100 and the gl.grm.network kinship fix from PR #101)
- Datasets: testset2.gl (274 individuals x 755 SNPs, 15.2% missing; 24 captive-bred offspring with `sire`/`dam` in ind.metrics), testset2.gl filtered with `gl.filter.callrate(threshold = 0.95)` (438 loci), testset2.gs (242 x 755), platypus.gl (81 individuals); dartR.data 1.2.5
- External program: EMIBD9 1.0 (macOS, serial), run on 37 individuals (captive-bred offspring, their sampled sires, EmmacBrisWive) x 306 loci
- Baseline: tests/testthat/test-gl.kin.R (snapshot captured pre-review, 25 expectations; the EMIBD9 test skips without the binary, the testset2 tests skip with dartR.data < 1.2.5)

## Verdict

**Standards: Needs work** — the preamble, datatype dispatch and error
checks follow the house structure; the gaps are a missing path argument
for EMIBD9, a verbose summary line that prints a constant, and an
untagged output scale.
**Spec: Rework** — both SNP methods return kinship on the wrong scale.
`method = "emibd9"` halves values that are already kinship, and the
default `method = "grm"` subtracts the mean inbreeding coefficient from
every pairwise kinship, which after a routine call-rate filter moves
every known parent-offspring pair out of the first-degree class. Plain
`G / 2` recovers the pedigree values (parent-offspring 0.255, full sibs
0.243, half sibs 0.118).

## Findings

**F1 [BLOCKER, confidence: high] — EMIBD9 kinship halved a second time
(spec: documented contract)**
`R/gl.kin.r:162` — `kin <- kin / 2` assumes `gl.run.EMIBD9()$rel` is
relatedness (2 x kinship). Since PR #100 it is kinship (EMIBD9 `r(1,2)`)
and carries `attr(rel, "scale") = "kinship"`; the `@details` text of
`gl.kin` predates that correction.
Failure scenario: 37 individuals including 24 offspring whose sire was
sampled: mean parent-offspring value in `$rel` 0.206, in
`gl.kin(method = "emibd9")` 0.103. Every diagonal element is 0.25,
contradicting the function's own contract (self-kinship 0.5 x (1 + F),
at least 0.5). Every downstream function of the series then sees
first-degree relatives as second-degree.
Proposed change: use `$rel` as kinship without rescaling (symmetrise
only); keep the diagonal fallback of 0.5; update `@details`.
**Consequence: numerical output doubles for `method = "emibd9"`.**

**F2 [BLOCKER, confidence: high] — GRM conversion subtracts mean
inbreeding from pairwise kinship (analysis: numerical correctness)**
`R/gl.kin.r:142-145` — `MS <- mean(diag(G) - 1); kin <- G / 2 - MS`
shifts every off-diagonal element by the mean inbreeding coefficient of
the sample, while the diagonal stays `diag(G) / 2`. `MS` is not a
property of pairs: it absorbs the Wahlund effect of pooling populations
(positive) and the shrinkage of the diagonal by mean imputation of
missing genotypes (negative). Goudet et al. (2018) standardise by the
mean allele sharing among pairs, applied to both diagonal and
off-diagonal, which is a different quantity.
Failure scenario (independent check against the pedigree in
testset2.gl; 24 parent-offspring, 46 full-sib, 20 half-sib pairs):

| Data | `MS` | P-O current / `G/2` | Full sibs current / `G/2` | Half sibs current / `G/2` | Median of all pairs current / `G/2` |
|---|---|---|---|---|---|
| testset2.gl, raw (15% NA) | -0.017 | 0.204 / 0.187 | 0.196 / 0.179 | 0.100 / 0.082 | 0.011 / -0.006 |
| after `gl.filter.callrate(0.95)` | 0.158 | 0.097 / 0.255 | 0.085 / 0.243 | -0.040 / 0.118 | -0.166 / -0.007 |

With the filtered data, the midpoint classes used by
`gl.report.kin.classes` (0.1875 / 0.09375 / 0.03125) place 11 of the 24
parent-offspring pairs in third degree and 13 in second degree; none in
first degree. `G / 2` places all 24 in first degree. Unrelated wild
individuals from different populations get kinship -0.215, so their
virtual offspring (`utils.kin.dgd`) has self-kinship 0.393, below the
0.5 minimum of a non-inbred individual. The same formula is the
`standardise = TRUE` option of `gl.grm.network`; that function is out of
scope here.
Proposed change: return `kin <- G / 2` (diagonal 0.5 x (1 + F),
off-diagonal kinship, both relative to the allele frequencies of the
supplied data, as the existing caveat states); update `@details`. The
choice of reference for the GRM is the custodian's; `G / 2` is proposed
because it is the scale `gl.grm.network` uses by default and it matches
the pedigree above.
**Consequence: numerical output changes for every `method = "grm"` call
(the default for SNP data), and with it every value computed from it by
the kinship series.**

**F3 [MEDIUM, confidence: high] — no way to tell gl.kin where EMIBD9 is
(spec: usability)**
`R/gl.kin.r:152` — `gl.run.EMIBD9(x, plot.out = FALSE, verbose = 0)` is
called with its default `emibd9.path = getwd()`, and `gl.kin` has no
argument to change it.
Failure scenario: EMIBD9 installed in `~/programs/emibd9-v1.0`, R working
directory is the project folder: `gl.kin(x, method = "emibd9")` stops
with "Cannot find EM_IBD_P in the folder given by emibd9.path: <project
folder>". The only workaround is `setwd()` into the binary folder, which
the review had to use.
Proposed change: add `emibd9.path = getwd()` to the signature (the same
default as `gl.run.EMIBD9`) and pass it through; document it.
**Consequence: new argument (additive); existing calls unchanged.**

**F4 [MEDIUM, confidence: high] — the verbose "gene diversity" is a
constant of the method (DOC5, proposed rule)**
`R/gl.kin.r:188-191` — at `verbose >= 3` the function prints
`GD = 1 - mean(kin)` over the full matrix. `rrBLUP::A.mat` centres on the
sample allele frequencies, so `sum(G) = 0` and the printed value is
exactly `1 + MS (n - 1) / n` whatever the diversity of the data.
Failure scenario: testset2.gl prints 0.9828 (= 1 + (-0.0172) x 273/274);
testset2.gs (dominant) prints 1.0002, a gene diversity above 1. After
change 2 the `grm` value would always be exactly 1. Gene diversity of a
subset relative to the full set is informative (0.80-0.81 for single
wild populations of testset2.gl); the full-set value is not. The
`@details` claim that the conventions "match those used by PMx" omits
that PMx measures kinship from founders, whereas here the base is the
supplied sample.
Proposed change: replace the GD line with the mean off-diagonal kinship,
and add to `@details` that GD of the full reference set is about 1 by
construction, so GD values are meaningful as comparisons between subsets
or changes (removals, additions), not as absolute diversity.
Docs and one message line.

**F5 [MEDIUM, confidence: high] — examples need dartR.data >= 1.2.5,
which DESCRIPTION does not require (DOC3, TST1)**
`R/gl.kin.r:70-76`, `DESCRIPTION:27` — `testset2.gl`/`testset2.gs` were
added in dartR.data 1.2.5; `Depends` lists `dartR.data` with no version.
Failure scenario: R CMD check with dartR.data 1.2.2 (the version the
PR #97 CI installed) fails on the example with "object 'testset2.gl' not
found". The same applies to every function of the kinship series.
Proposed change: set `dartR.data (>= 1.2.5)` in `Depends`. This is a
package-level decision pending with the kinship-suite author; it can be
deferred from this PR.

**F6 [LOW, confidence: high] — output not tagged as kinship
(API3, proposed rule)**
`R/gl.kin.r:179-181` — `grm` and `dominant` output carries no
`attr(kin, "scale")`; `emibd9` output inherits `"kinship"` from
`gl.run.EMIBD9` by accident of `pmax()`.
Failure scenario: `gl.grm.network(gl.kin(x), x)` treats an untagged
matrix as relatedness and halves it, so a user plotting the series'
kinship matrix sees half the values and different categories.
Proposed change: set `attr(kin, "scale") <- "kinship"` for every method.

**F7 [LOW, confidence: medium] — pairs with kinship above 0.5 pass
without warning (VRB4, proposed rule)**
`R/gl.kin.r:174` (via `utils.kin.dominant`) — the dominant diagonal is
fixed at 0.5, but off-diagonal values are not bounded by it.
Failure scenario: testset2.gs returns 16 pairs above 0.5 (maximum 1.05,
UC_206 with its neighbour, 698 shared loci), which are most likely
duplicate samples. A virtual offspring of such a pair gets self-kinship
above 1 in `utils.kin.dgd`. SNP methods can produce the same for
duplicates.
Proposed change: after the matrix is built, if any off-diagonal value
exceeds 0.5, warn at `verbose >= 1` with the number of pairs and a hint
that they may be duplicate samples (see `gl.report.replicates` in
dartR.base).

**F8 [LOW, confidence: medium] — missing data attenuate GRM kinship; the
documentation is silent (DOC5, proposed rule)**
`R/gl.kin.r:28-65` — `gl.grm` passes mean-imputed genotypes to
`rrBLUP::A.mat`, which shrinks kinship toward 0.
Failure scenario: raw testset2.gl (15% missing), `G / 2`:
parent-offspring 0.187; after `gl.filter.callrate(0.95)`: 0.255 against
the expected 0.25.
Proposed change: add a sentence to `@details` recommending call-rate
filtering before `gl.kin` for kinship thresholds and classes. Docs only.

## Proposed changes

1. Stop halving EMIBD9 kinship; symmetrise only; update `@details` (F1).
   **Consequence: numerical output doubles for `method = "emibd9"`.**
2. Return `G / 2` for `method = "grm"` instead of `G / 2 - MS`; update
   `@details` (F2).
   **Consequence: numerical output changes for every `method = "grm"`
   call, the default for SNP data, and for every kinship-series result
   built on it.**
3. Add `emibd9.path = getwd()` and pass it to `gl.run.EMIBD9` (F3).
   **Consequence: new argument (additive); existing calls unchanged.**
4. Replace the verbose GD line with mean off-diagonal kinship; document
   that full-set GD is about 1 by construction and qualify the PMx
   statement (F4).
5. Add `dartR.data (>= 1.2.5)` to `Depends` (F5).
6. Tag output with `attr(kin, "scale") <- "kinship"` (F6).
7. Warn at `verbose >= 1` when any off-diagonal kinship exceeds 0.5 (F7).
8. Add a missing-data sentence to `@details` (F8).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run. Notes
  without findings: `utils.flag.start(build = "v.2026.1")` uses the
  outdated `build=` argument (FS3), still accepted by dartR.base 1.2.3;
  `invisible(kin)` is appropriate for a matrix (FS10); no history append
  (FS8 not applicable, no genlight returned); rrBLUP is reached through
  `gl.grm`, which carries its own guard (DEP1).
- Spec: behaviour vs roxygen on testset2.gl, testset2.gs, platypus.gl —
  run; independent check against the pedigree columns of testset2.gl —
  run.
- EMIBD9 path: run on a 37-individual subset (EMIBD9 1.0, serial).
- Downstream effect on `gl.report.kin.classes`: checked by applying its
  documented midpoint breaks to the known parent-offspring pairs, not by
  running the function (it is reviewed separately).
- `utils.kin.dominant`, `utils.kin.check`, `utils.kin.dgd`: read as
  dependencies; not reviewed as functions in their own right.
- FBM path (DAT6): SKIPPED — no FBM fixture; `gl.grm` densifies with
  `as.matrix(x)` and `utils.kin.dominant` does the same.
- Callers in dartr2shiny: none found (`grep gl.kin` in ~/dartr2shiny).
  Callers in the package: `utils.kin.check`, `gl.report.gd.projection`.

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis Mijangos | |
| 2 | approved | Luis Mijangos | |
| 3 | approved | Luis Mijangos | |
| 4 | approved | Luis Mijangos | |
| 5 | approved | Luis Mijangos | |
| 6 | approved | Luis Mijangos | |
| 7 | approved | Luis Mijangos | |
| 8 | approved | Luis Mijangos | approved in chat (left out of the approval boxes) |

## Outcome

- Change 1 (F1): `gl.kin(method = "emibd9")` equals `gl.run.EMIBD9()$rel`;
  diagonal 0.5 (was 0.25). Test "emibd9: returns $rel unscaled".
- Change 2 (F2): `gl.kin()` equals `gl.grm() / 2`; on call-rate-filtered
  testset2.gl the 24 parent-offspring pairs have mean 0.255, all above
  0.1875 (was 0.097, none). Tests "kinship is G / 2" and "known
  parent-offspring pairs".
- Change 3 (F3): EMIBD9 found through `emibd9.path` with the working
  directory set to `tempdir()`.
- Change 4 (F4): verbose summary prints mean pairwise kinship (-0.0018 on
  testset2.gl); `@details` explains full-set GD.
- Change 5 (F5): `Depends: dartR.data (>= 1.2.5)`.
- Change 6 (F6): `attr(kin, "scale") = "kinship"` for all methods.
- Change 7 (F7): testset2.gs at `verbose = 1` warns about 16 pairs above
  0.5; silent at `verbose = 0`.
- Change 8 (F8): `@details` sentence on missing data and call-rate
  filtering (docs only).
- Snapshot: the pre-review baseline (25 expectations) run on the changed
  code gave 7 failures, each mapped: scale attribute NULL (change 6);
  `G/2 - MS` identity, mean and maximum of testset2.gl (change 2);
  "Gene diversity" line (change 4); EMIBD9 halving and 0.25 diagonal
  (change 1). Diagonal range and all dominant values unchanged.
- Specification tests: tests/testthat/test-gl.kin.R, 29 expectations,
  0 failures (EMIBD9 1.0 installed).
- `gl.kin(testset2.gl, verbose = 3)` runs end to end; the examples of all
  13 kinship-series functions run without error on the changed code
  (values not re-checked; those functions are reviewed separately).
- PR: #108.

```json
{
  "function": "gl.kin",
  "package": "dartR.captive",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "573a30f",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "BLOCKER", "confidence": "high", "rule": "spec: documented contract", "status": "approved", "change": 1},
    {"id": "F2", "severity": "BLOCKER", "confidence": "high", "rule": "analysis: numerical correctness", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "spec: usability", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DOC3", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "API3", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "medium", "rule": "VRB4", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "medium", "rule": "DOC5", "status": "approved", "change": 8}
  ],
  "coverage_skipped": ["DAT6: no FBM fixture", "gl.report.kin.classes not run; breaks applied directly"],
  "status": "done",
  "pr": 108
}
```
