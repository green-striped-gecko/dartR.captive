# Review: gl.sim.relatedness (dartR.captive)
- Family mode: analysis (simulation)
- Date: 2026-09-23
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 3e87a84 (origin/dev, includes the gl.run.EMIBD9 fixes from PR #100)
- Datasets: testset.gl[1:30, ] (gl.filter.allna); a simulated unrelated population of 60 individuals x 641 polymorphic loci (`dartR.sim::gl.sim.ind` from platypus.gl allele frequencies, no missing data)
- External program: EMIBD9 1.0.0.0 (macOS, serial)
- Baseline: tests/testthat/test-gl.sim.relatedness.R (snapshot captured pre-review, 5 expectations; EMIBD9 test skips without the binary)

## Verdict

**Standards: Needs work** — the simulated values are never returned, the
summary prints at every verbosity, each replicate draws an EMIBD9 heatmap,
and the function refuses to run without a package (`related`) it never
uses.
**Spec: Rework** — "full.sib" measures parent-offspring kinship, and the
reported "CI" is the confidence interval of the mean, not the range of
values the documentation invites users to use as thresholds. Half-sib and
first-cousin means are correct on clean data (0.126-0.157 and
0.053-0.077 intervals around the expected 0.125 and 0.0625).

## Findings

**F1 [HIGH, confidence: high] — "full.sib" returns parent-offspring
kinship (DOC5)**
`R/gl.sim.relatedness.R:138-163` — one offspring is simulated and the
value returned is the mean of its kinship with each parent.
Failure scenario: simulated population, 15 families: the function's
value has mean 0.259, sd 0.008; the kinship between the two full-sib
offspring of the same families has mean 0.251, sd 0.020. The means agree
(both 0.25 in expectation) but the spread of true full sibs is 2.5 times
wider, so any interval derived from "full.sib" is too narrow for full
sibs. Parents are also sampled with replacement, so the same individual
can be both parents (a selfed "sib").
Proposed change: simulate two offspring of the same parents (sampled
without replacement) and return their kinship.
**Consequence: "full.sib" values and their spread change.**

**F2 [HIGH, confidence: high] — the interval is the CI of the mean, not
the range of values (DOC5)**
`R/gl.sim.relatedness.R:253-256` — `confint(lm(relatedness ~ 1))` is the
confidence interval of the mean kinship; it shrinks with `nboots`. The
description presents the Speed & Balding ranges (e.g. full sibs 0.204-
0.296) "to guide the choosing of the relatedness threshold", which are
ranges of individual values.
Failure scenario: true full sibs, 15 replicates: 95% of values between
0.220 and 0.287; the function's interval 0.240-0.262. With `nboots = 100`
the interval narrows further, while the spread of individual pairs does
not.
Proposed change: report the empirical `(1 - conf)/2` and `1 - (1 - conf)/2`
quantiles of the simulated values as the interval, and keep the CI of the
mean as a separate, labelled output.
**Consequence: the reported interval widens and stops shrinking with
`nboots`.**

**F3 [MEDIUM, confidence: high] — results are printed, not returned
(DOC5, FS10)**
`R/gl.sim.relatedness.R:312-313` — the function ends with `print(sum)`
and `print(CI)`, so it returns the CI matrix invisibly; the simulated
values, the mean and the plot are lost. `@return` promises summary
statistics and a histogram. The Shiny app shows the return as a table.
Proposed change: return a named list: `values` (simulated kinship),
`mean`, `interval` (F2), `ci.mean`, `plot`; print the summary only at
`verbose >= 3`.
**Consequence: the return value changes from a 1 x 2 matrix to a list.**

**F4 [MEDIUM, confidence: high] — each replicate draws an EMIBD9 heatmap
and prints EMIBD9 output (VRB3, PLT3)**
`R/gl.sim.relatedness.R:158,186,224` — `plot.out` and `verbose` are
passed to every `gl.run.EMIBD9()` call, so `nboots = 100` draws 100
heatmaps of the whole dataset and, at the default verbosity, prints
EMIBD9's console output 100 times; the summary also prints at
`verbose = 0`.
Proposed change: call `gl.run.EMIBD9(..., plot.out = FALSE, verbose = 0)`
and report progress (replicate i of n) at `verbose >= 2`.

**F5 [MEDIUM, confidence: high] — requires `related`, which it never uses
(DEP1)**
`R/gl.sim.relatedness.R:120-131` — the function returns -1 unless the
GitHub-only package `related` is installed; nothing in the function calls
it, and it is not declared in DESCRIPTION. The `requireNamespace("dartR.captive")`
check inside dartR.captive itself is also redundant.
Failure scenario: a user without `related` gets "Package related needed"
and cannot run the simulation.
Proposed change: remove both checks.

**F6 [LOW, confidence: high] — the Shiny app's call fails (API2)**
dartr2shiny calls `gl.sim.relatedness(..., iseed = input$iseed, ...)`;
the argument is `ISeed`, so R stops with "unused argument (iseed = 42)"
(baseline test). Fix on the dartr2shiny side (`ISeed`); noted here
because the function's signature is the reference.
Proposed change: none in dartR.captive; a dartr2shiny merge request.

**F7 [LOW, confidence: high] — smaller defects (FS5, STY1, DOC1, DOC6,
DOC7)**
- `rel` is not validated: `rel = "fullsib"` fails late with "object 'rr'
  not found"; `match.arg()`.
- `first.cousin` draws the second pair of parents independently, so it
  can reuse the first pair (the cousins then share more ancestry); draw
  them from the remaining individuals.
- Simulations use R's random number generator, which `ISeed` (EMIBD9's
  seed) does not set, so results are not reproducible from the
  arguments; document that `set.seed()` controls it.
- Two plots lose their x/y labels (`labs()` after a missing `+`); three
  near-identical plot blocks.
- SilicoDArT fails only inside the first EMIBD9 run; check up front.
- Docs: title/labels say "relatedness" but values are kinship (EMIBD9
  `r(1,2)`, see the gl.run.EMIBD9 review); `@param EM_Method` says 0/1
  while gl.run.EMIBD9 documents 1/2/3; `@author` has no `Author(s):`;
  non-ASCII table characters; `build = "Jody"`.
Proposed change: fix each as listed.

## Proposed changes

1. "full.sib" = kinship between two full-sib offspring; parents without
   replacement (F1).
   **Consequence: full-sib values and spread change.**
2. Interval = empirical quantiles of the simulated values; CI of the mean
   kept as a separate output (F2).
   **Consequence: the interval widens and no longer shrinks with
   `nboots`.**
3. Return a named list (values, mean, interval, ci.mean, plot); print
   only at `verbose >= 3` (F3).
   **Consequence: the return value changes from a matrix to a list.**
4. EMIBD9 runs without plots and quietly; progress at `verbose >= 2`
   (F4).
5. Remove the `related` and `dartR.captive` checks (F5).
6. Smaller defects and docs (F7).
F6 is fixed on the dartr2shiny side (separate MR).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY — run
- Spec: each `rel` run on a simulated unrelated population (unrelated
  kinship mean 0.022, sd 0.017) and on `testset.gl`; true full-sib vs
  returned values compared on 15 families; interval vs quantiles — run
- `testset.gl` note: on `testset.gl[1:30, ]` the half-sib and first-cousin
  means were 0.245 and 0.211, far above expectation; that subset mixes
  populations and has few informative loci, so the simulated population
  is the reference here
- Callers: dartr2shiny (`slick_table_plot`; broken call, F6); no other
  `dartR.*` package calls it
- Windows / parallel EMIBD9: not tested (no binaries)

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis Mijangos | |
| 2 | approved | Luis Mijangos | |
| 3 | approved | Luis Mijangos | |
| 4 | approved | Luis Mijangos | |
| 5 | approved | Luis Mijangos | |
| 6 | approved | Luis Mijangos | `@author` lists Sam Amini as author (assumption: only the custodian was named) |
| F6 | deferred | Luis Mijangos | dartr2shiny MR not opened now |

## Outcome

All six approved changes applied on 2026-09-23; the function body was
rewritten (signature unchanged) and `devtools::document()` run.
`tests/testthat/test-gl.sim.relatedness.R`: 13 expectations pass (EMIBD9
tests skip without the binary); the tests use a simulated population, so
they do not depend on the dartR.data version of `testset.gl`.

On the simulated population (60 unrelated individuals, 641 loci), 15
replicates each:

| rel | mean | sd | interval (95% of values) | CI of the mean | expected |
|---|---|---|---|---|---|
| full.sib | 0.247 | 0.022 | 0.215-0.274 | 0.234-0.259 | 0.25 |
| half.sib | 0.122 | 0.024 | 0.086-0.156 | 0.109-0.136 | 0.125 |
| first.cousin | 0.075 | 0.024 | 0.042-0.117 | 0.062-0.089 | 0.0625 |

- Change 1 (F1): full sibs are two offspring of the same two parents
  (sampled without replacement); sd 0.022, as for true full sibs (0.020),
  instead of 0.008 for the old parent-offspring average.
- Change 2 (F2): `interval` holds the empirical quantiles; `ci.mean` the
  CI of the mean.
- Change 3 (F3): returns (invisibly) `list(values, mean, interval,
  ci.mean, plot)`; the summary prints at `verbose >= 3`.
- Change 4 (F4): EMIBD9 runs with `plot.out = FALSE, verbose = 0`; the
  missing-data message of `gl.sim.offspring` is captured; 0 lines printed
  at `verbose = 0`; "Simulating ... pair i of n" at `verbose >= 2`.
- Change 5 (F5): the `related` and `dartR.captive` checks are removed.
- Change 6 (F7): `rel` checked with `match.arg`; SilicoDArT and fewer than
  4 individuals error up front; first-cousin parents are four distinct
  individuals; one plot with axis labels; docs describe kinship, the
  interval, the simulation scheme and `set.seed()`; `EM_Method` doc
  matches gl.run.EMIBD9; `@family`, `Author(s)`, ASCII.
- Not in this change: the dartr2shiny call (`iseed =`, F6) still fails;
  the app also needs to show one element of the new return list.

```json
{
  "function": "gl.sim.relatedness",
  "package": "dartR.captive",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "3e87a84",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "FS10", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "VRB3", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "DEP1", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "API2", "status": "deferred", "change": null},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved", "change": 6}
  ],
  "coverage_skipped": ["Windows/parallel EMIBD9: no binaries"],
  "status": "pr-open",
  "pr": null
}
```
