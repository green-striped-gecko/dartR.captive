# Tests of the reviewed gl.select.pairs
# (function-review/reports/dartR.captive/gl.select.pairs.md).

skip_without_testset2 <- function() {
  skip_if_not(exists("testset2.gl") && exists("testset2.gs"),
              "testset2.gl/gs need dartR.data >= 1.2.5")
}

cb.setup <- function() {
  kin <- gl.kin(testset2.gl, verbose = 0)
  cb <- gl.keep.pop(testset2.gl, pop.list = "EmmacCaptBred", verbose = 0)
  list(cb = cb, kin = kin[indNames(cb), indNames(cb)])
}

test_that("structure, read-only, three schemes", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  s <- cb.setup()
  cb0 <- s$cb
  p <- gl.select.pairs(s$cb, kin = s$kin, verbose = 0)
  expect_identical(s$cb, cb0)
  expect_named(p, c("sire", "dam", "f.off", "mk.sire", "mk.dam", "dgd.cum"))
  expect_equal(nrow(p), 10)
  expect_equal(p$sire[1:2], c("CB_CD_01", "CB_Y_02"))
  expect_equal(p$dam[1:2], c("CB_Y_03", "CB_CD_05"))
  expect_equal(attr(p, "gd.start"), 0.9337007, tolerance = 1e-6)
  expect_equal(attr(p, "gd.projected"), 0.9335477, tolerance = 1e-6)
  expect_true(all(p$f.off < 0.125))
  expect_equal(p$dgd.cum[3],
               utils.kin.dgd(s$kin, add.pairs = as.matrix(p[1:3, 1:2])))
  # the closed-form dynamic scoring picks the brute-force best pair
  sx <- s$cb@other$ind.metrics$sex
  cand <- expand.grid(sire = indNames(s$cb)[sx == "Male"],
                      dam = indNames(s$cb)[sx == "Female"],
                      stringsAsFactors = FALSE)
  cand <- cand[s$kin[cbind(cand$sire, cand$dam)] < 0.125 &
                 !cand$sire %in% p$sire[1:2] & !cand$dam %in% p$dam[1:2], ]
  g <- vapply(seq_len(nrow(cand)), function(i) {
    utils.kin.dgd(s$kin, add.pairs = rbind(as.matrix(p[1:2, 1:2]),
                                           c(cand$sire[i], cand$dam[i])))
  }, numeric(1))
  expect_equal(max(g), p$dgd.cum[3])
  for (sch in c("static", "ranked")) {
    q <- gl.select.pairs(s$cb, kin = s$kin, scheme = sch, verbose = 0)
    expect_equal(nrow(q), 10)
    expect_equal(attr(q, "gd.projected"), 0.9335477, tolerance = 1e-6)
  }
  p2 <- gl.select.pairs(s$cb, kin = s$kin, max.per.sire = 2,
                        max.per.dam = 2, n.pairs = 15, verbose = 0)
  expect_equal(nrow(p2), 15)
  expect_true(all(table(p2$sire) <= 2))
  expect_equal(attr(p2, "gd.projected"), 0.936508, tolerance = 1e-6)
})

test_that("argument handling, SilicoDArT", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  s <- cb.setup()
  # invalid arguments are errors (were reset to defaults; change 4)
  bad <- list(list(scheme = "Static"), list(scheme = NA),
              list(max.per.sire = 0), list(max.per.sire = NA_real_),
              list(max.per.dam = 0.5), list(f.max = 0), list(f.max = 2),
              list(f.max = NA_real_), list(n.pairs = 0),
              list(n.pairs = NA_real_), list(n.pairs = "5"))
  for (b in bad) {
    expect_error(do.call(gl.select.pairs,
                         c(list(s$cb, kin = s$kin, verbose = 0), b)),
                 paste(names(b), "must be"))
  }
  expect_equal(nrow(gl.select.pairs(s$cb, kin = s$kin, n.pairs = 3.9,
                                    max.per.sire = Inf, verbose = 0)), 3)
  kg <- gl.kin(testset2.gs, verbose = 0)
  cbs <- gl.keep.pop(testset2.gs, pop.list = "EmmacCaptBred", verbose = 0)
  expect_equal(nrow(gl.select.pairs(cbs, kin = kg, verbose = 0)), 10)
})

test_that("NA sex, NA kinship, low call rate, self-reference", {
  skip_if_not_installed("rrBLUP")
  skip_without_testset2()
  s <- cb.setup()
  # NA sex counts as unknown (was a crash; change 3)
  cb2 <- s$cb
  cb2@other$ind.metrics$sex[1] <- NA
  expect_output(p <- gl.select.pairs(cb2, kin = s$kin, verbose = 1),
                "1 individuals of unknown sex")
  expect_false(indNames(cb2)[1] %in% c(p$sire, p$dam))
  expect_gt(nrow(p), 0)
  # NA kinship: values in every scheme, the NA pair never selected
  # (dynamic crashed, static and ranked gave NA; change 2)
  k2 <- s$kin
  k2["CB_CD_01", "CB_Y_03"] <- k2["CB_Y_03", "CB_CD_01"] <- NA
  for (sch in c("dynamic", "static", "ranked")) {
    expect_output(q <- gl.select.pairs(s$cb, kin = k2, scheme = sch,
                                       verbose = 1),
                  "1 missing \\(NA\\) kinship")
    expect_false(anyNA(q$dgd.cum))
    expect_false(anyNA(q$mk.sire))
    expect_false(is.na(attr(q, "gd.start")))
    expect_false(any(q$sire == "CB_CD_01" & q$dam == "CB_Y_03"))
  }
  # CB_CD_01 (lowest mean kinship) is no longer sorted last in static
  q <- gl.select.pairs(s$cb, kin = k2, scheme = "static", verbose = 0)
  expect_equal(q$sire[1], "CB_CD_01")
  # low call rate warns (change 1)
  expect_output(gl.select.pairs(s$cb, kin = s$kin, verbose = 1),
                "24 individuals have call rate below 0.8")
  expect_error(gl.select.pairs(s$cb, verbose = 0), "wider")
})
