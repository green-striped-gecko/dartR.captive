# Audit: self-referenced kinship in the kinship series (dartR.captive)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), cross-function audit following the gl.report.kinship review (#112)
- Package commit: 537c8b2 (origin/dev)
- Datasets: testset2.gl; its 24 captive-bred individuals (EmmacCaptBred) analysed alone

## The defect

`gl.kin(method = "grm")` centres each locus on the allele frequencies of
the individuals it is given, so every row of the kinship matrix sums to 0.
A mean kinship (MK) taken over all of those individuals is therefore 0 for
everyone (rounding residuals of about 1e-17), and the gene diversity of the
whole set is exactly 1. Removing individual i changes the sum of the matrix
by `kin[i, i]` only, so the change in GD depends on i's inbreeding and on
nothing else. The dominant estimator is centred the same way; with missing
data its row means are near 0 rather than exactly 0.

MK and whole-group GD are informative only when kinship is estimated on a
wider reference (for example, the wild source populations) and the matrix
is then restricted to the managed group. The series computes `kin` from
`x` itself whenever `kin = NULL`, which is the default.

## Evidence (24 captive-bred individuals; "self" = default `kin`, "ref" = `gl.kin(testset2.gl)` restricted to them)

| Function | Self-referenced output | With reference kinship |
|---|---|---|
| all | row means -2.7e-17 to 1.3e-17 | 0.040 to 0.107 |
| `gl.report.ind.remove` | MK column ~1e-17; dGD correlates 1.000 with -kin[i,i] (ranks by inbreeding only) | dGD vs -kin[i,i]: -0.27; rank correlation with self 0.65 |
| `gl.report.repro.targets` | MK ~1e-17; offspring targets allocated from rounding noise | targets differ for 16 of 24 individuals |
| `gl.select.pairs` | `gd.start = 1`; `mk.sire`, `mk.dam` ~1e-17 drive pair order | `gd.start = 0.934`; different pairs chosen |
| `gl.report.mate.suitability` | `mkdiff` 0 to 4e-17 ranked as a criterion; MSI | MSI differs for 113 of 143 pairings |
| `gl.report.gd.projection` (from `x`) | `gd.now = 1` | 0.934 |
| `gl.report.kin.sets` | `meanMK` ~1e-18 for every set (dataset-wide row means averaged per set); `GD.w` correct (block-based) | — |

Not affected: `gl.report.ind.add`, `gl.report.ind.move`,
`gl.report.kin.groups` (use population blocks), `gl.report.kin.classes`,
`gl.report.kin.confidence` (pairwise values or genotypes). Pairwise
quantities (offspring inbreeding `f.off`, kinship classes) are valid in
every function.

`gl.report.kinship` (#112) computes MK per population block, which is
informative when `x` holds more than one population; with a single
population it is self-referenced and has the same defect.

## Decision (Luis Mijangos, 2026-09-24)

Stop without reference kinship and accept a larger `kin` (option 1 of 4:
stop, warn, managed-group argument, hold for the suite author).
`gl.report.kin.sets` `meanMK` switched to the within-set block mean.

## Outcome

- `gl.kin()` output records `attr(, "ref.ids")`, the individuals kinship
  was estimated on.
- `utils.kin.check()` restricts a `kin` covering more individuals than `x`
  to `indNames(x)` (attributes kept), and with `need.reference = TRUE`
  rejects self-referenced kinship: `kin = NULL`, `ref.ids` equal to
  `indNames(x)`, or row means over `x` all within 1e-10 of 0.
- `need.reference = TRUE` in `gl.report.ind.remove`,
  `gl.report.repro.targets`, `gl.select.pairs`,
  `gl.report.mate.suitability`, and `gl.report.gd.projection` (x path,
  which now also uses `utils.kin.check`).
- `gl.report.kin.sets`: `meanMK` = block mean (= 1 - `GD.w`; captive set
  0.0663, was ~1e-18).
- `gl.report.kinship`: warning for one population with self-referenced
  kinship.
- Examples: the two that computed kinship internally
  (`gl.report.ind.remove` SilicoDArT, `gl.report.gd.projection`) now pass
  full-dataset kinship.
- Closed gap: a dominant (SilicoDArT) matrix estimated on `x` and subset by
  hand in `indNames(x)` order loses `ref.ids`, and its row means are only
  near 0 (-0.010 to 0.006 on the 24 captive-bred individuals), so the
  row-mean test misses it. `utils.kin.check()` now re-estimates dominant
  kinship on `x` and rejects a matching matrix (off-diagonals equal within
  1e-8). The dominant estimator is the only SilicoDArT method, so the match
  is exact; the full-dataset matrix restricted by hand differs by up to
  0.13 and is accepted. Before: accepted; after: rejected. (A reordered
  matrix was already caught, as its own rows become the reference.)
- Evidence: tests/testthat/test-utils.kin.check.R 27 pass; full suite 288
  expectations, 6 failures (the pre-existing gl.grm/gl.assign.grm
  snapshots); examples of all 17 kinship and GRM functions run. No callers
  in dartr2shiny or other `dartR.*` packages.
