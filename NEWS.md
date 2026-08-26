# dartR.captive NEWS

## Unreleased

* `gl.grm()`: now errors on SilicoDArT (presence/absence) input instead of
  silently returning a numerically meaningless matrix — the additive
  relationship algorithm and the documented diagonal range (1 to 2) only
  hold for SNP dosage data. Callers passing SilicoDArT data will need to
  filter to SNP data before calling `gl.grm()`.
* `gl.grm()`: fixed a crash (`object 'p3' not found`) when
  `plotheatmap = FALSE` was combined with a non-`NULL` `plot.file`; that
  combination now computes and returns the matrix and warns that nothing
  was saved, since no plot is generated.
* `gl.grm()`: corrected the documented default for `legendy` (was `1`,
  actual default is `0.5`).
