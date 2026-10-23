# Paper B ruling Q9: display layer (models/B/R/display.R). Synthetic entries only.

local({
  root <- if (exists("MODELS_ROOT")) MODELS_ROOT else normalizePath(file.path(getwd(), "..", "..", ".."), winslash = "/")
  source(file.path(root, "B", "R", "display.R"), local = globalenv())
})

test_that("fmt_p_sig gives 1-2 significant digits, <0.001 and >0.99", {
  expect_equal(fmt_p_sig(c(0.27, 0.049, 0.0021, 0.5, 0.0999, 0.00099, 0.001, 0.996, 1, 0.721)),
               c("0.27", "0.049", "0.0021", "0.50", "0.10", "<0.001", "0.0010", ">0.99", ">0.99", "0.72"))
  expect_true(is.na(fmt_p_sig(NA_real_)))
  expect_error(fmt_p_sig(1.2), "\\[0, 1\\]")
})

test_that("apply_b_display rewrites displays only; values are unchanged", {
  mk <- function(v, disp) list(value = v, display = disp, script = "s", description = "d")
  e <- list(B.h4.bw.T2_dry.p = mk(0.27443, "0.274"),
            B.h6.anaem.recent12_tmax.p = mk(0.00049, "<0.001"),
            B.supp.x.p = mk(NA, "not estimable"),
            B.med.oni.T2.flag = mk("ok", "ok"),
            B.med.oni.T2.prop_flag = mk("total_effect_ci_includes_0", "total_effect_ci_includes_0"),
            B.med.oni.T2.prop = mk(-0.1567, "not estimated (total_effect_ci_includes_0)"),
            B.med.oni.T3.prop = mk(0.2, "21% (95% CI 3-40%)"),
            B.sens.differ = mk("nonmover/xsbw", "nonmover/xsbw"),
            B.h4.verdict_terms = mk("neither", "neither"),
            B.h4.verdict = mk("not consistent with", "not consistent with"),
            B.h4.bw.T2_dry = mk(-0.4, "−0.4 (95% CI −3.9 to 3.0)"))
  out <- apply_b_display(e)
  expect_identical(lapply(out, `[[`, "value"), lapply(e, `[[`, "value"))       # values untouched
  expect_identical(out$B.h4.bw.T2_dry.p$display, "0.27")
  expect_identical(out$B.h6.anaem.recent12_tmax.p$display, "<0.001")
  expect_identical(out$B.supp.x.p$display, "not estimable")
  expect_identical(out$B.med.oni.T2.flag$display, "bootstrap completed (at most 10% of resamples failed)")
  expect_identical(out$B.med.oni.T2.prop_flag$display, "not estimated (total-association CI includes 0)")
  expect_identical(out$B.med.oni.T2.prop$display, "not estimated (total-association CI includes 0)")
  expect_identical(out$B.med.oni.T3.prop$display, e$B.med.oni.T3.prop$display)  # an estimated proportion is kept
  expect_identical(out$B.sens.differ$display, "non-movers, El Niño × sensitivity, birth weight")
  expect_identical(out$B.h4.verdict_terms$display, "neither the second- nor the third-trimester term met the rule")
  expect_identical(out$B.h4.verdict$display, "not consistent with")              # verdicts untouched
  expect_identical(out$B.h4.bw.T2_dry$display, e$B.h4.bw.T2_dry$display)
  # no raw code survives in a rewritten display
  expect_false(any(grepl("_", vapply(out[c("B.med.oni.T2.flag", "B.med.oni.T2.prop_flag", "B.sens.differ")], `[[`, "", "display"))))
})

test_that("differ_prose handles several pairs and none; unknown codes fail loudly", {
  expect_identical(differ_prose("none"), "none")
  expect_identical(differ_prose("card/h4bw, cluster1/xshaz"),
                   "card-recorded birth weight, trimester dryness, birth weight; district-only clustering, El Niño × sensitivity, HAZ 24–59 months")
  expect_error(differ_prose("foo/h4bw"), "no prose")
  expect_error(verdict_terms_prose("T1"), "no prose")
})
