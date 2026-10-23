# Paper B Task 7 helpers (models/B/R/prep.R): exposure construction, bands, verdict rules,
# contract key expansion and small-cell suppression. Synthetic values only.

local({
  root <- if (exists("MODELS_ROOT")) MODELS_ROOT else normalizePath(file.path(getwd(), "..", "..", ".."), winslash = "/")
  source(file.path(root, "B", "R", "prep.R"), local = globalenv())
})

test_that("ym_add moves YYYYMM months across year boundaries", {
  expect_equal(ym_add(c(201501L, 201512L, 201006L), c(-5L, 1L, -17L)), c(201408L, 201601L, 200901L))
})

test_that("enso_year is the June-May year of the T2 midpoint (birth month - 5)", {
  # birth Nov 2015 -> midpoint Jun 2015 -> 2015; birth Oct 2015 -> May 2015 -> 2014
  # birth Feb 2010 -> Sep 2009 -> 2009; birth May 2016 -> Dec 2015 -> 2015; birth Jun 2016 -> Jan 2016 -> 2015
  expect_equal(enso_year_of(c(201511L, 201510L, 201002L, 201605L, 201606L, 201610L, 201611L)),
               c(2015L, 2014L, 2009L, 2015L, 2015L, 2015L, 2016L))
  expect_true(is.na(enso_year_of(NA_integer_)))
})

test_that("mw_mean is the month-weighted mean and NaN if any component is missing", {
  r <- mw_mean(list(c(1, 1, NA), c(2, 3, 1)), c(1, 3))
  expect_equal(r[1:2], c((1 + 6) / 4, (1 + 9) / 4))
  expect_true(is.na(r[3]))
  expect_error(mw_mean(list(1, 2), c(1, 2, 3)), "weights")
})

test_that("spi3_1000d uses weights (3,3,3,12,12)/33 for 24-59 and (3,3,3,12)/21 for 0-23", {
  T1 <- c(1, 1); T2 <- c(0, 0); T3 <- c(-1, -1); i0 <- c(0.5, 0.5); i1 <- c(2, NA)
  e2459 <- (3 * 1 + 0 - 3 + 12 * 0.5 + 12 * 2) / 33
  r <- spi3_1000d(T1, T2, T3, i0, i1, "2459")
  expect_equal(r[1], e2459)
  expect_true(is.na(r[2]))                       # inf12_23 NaN -> NaN, never a partial mean
  r0 <- spi3_1000d(T1, T2, T3, i0, i1, "0023")
  expect_equal(r0, rep((3 * 1 + 0 - 3 + 12 * 0.5) / 21, 2))   # inf12_23 ignored
  expect_true(is.na(spi3_1000d(1, 1, 1, NA, 1, "0023")))
  # not the equal-weight mean of the five windows
  expect_false(isTRUE(all.equal(e2459, mean(c(1, 0, -1, 0.5, 2)))))
})

test_that("bands and factor recodes follow the contract", {
  expect_equal(as.character(band_bord(c(1, 2, 3, 4, 9))), c("1", "2", "3", "4+", "4+"))
  expect_equal(as.character(band_v012(c(15, 19, 20, 29, 30, 34, 35, 49))),
               c("15-19", "15-19", "20-24", "25-29", "30-34", "30-34", "35+", "35+"))
  expect_equal(as.character(caste_f(c(1, 2, 3, 4, 8, NA))), c("SC", "ST", "OBC", "none", "missing", "missing"))
  expect_equal(as.character(religion_f(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 96))),
               c("Hindu", "Muslim", "Christian", "Sikh", "Buddhist", rep("other", 5)))
  expect_error(band_v012(c(14)), "v012")
})

test_that("nonmover: always (95), or years in place >= child's age in years + 1; visitors out", {
  # age_m 30 -> 2 years -> need v104 >= 3
  expect_equal(nonmover(c(95, 3, 2, 96, 97, 0, 1), c(30, 30, 30, 30, 30, 5, 5)),
               c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE))
})

test_that("dry-orientation rule: negative estimate with CI upper bound below 0", {
  expect_equal(rule_dry_negative(c(-10, -10, 5, NA, -1), c(-1, 0.5, 9, NA, NA)),
               c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_equal(verdict_text(c(TRUE, FALSE)), c("consistent with", "not consistent with"))
  expect_equal(verdict_terms(c(TRUE, FALSE, TRUE, FALSE), c(FALSE, TRUE, TRUE, FALSE)),
               c("T2", "T3", "T2 and T3", "neither"))
})

test_that("flip_ci negates the estimate and swaps the bounds", {
  f <- flip_ci(2, 1, 3)
  expect_equal(c(f$est, f$lo, f$hi), c(-2, -3, -1))
})

test_that("expand_keys expands brace sets shell-style, including empty alternatives", {
  expect_equal(sort(expand_keys("B.h4.bw.{T2,T3}_dry{,.p}")),
               sort(c("B.h4.bw.T2_dry", "B.h4.bw.T2_dry.p", "B.h4.bw.T3_dry", "B.h4.bw.T3_dry.p")))
  expect_equal(expand_keys("B.n.children"), "B.n.children")
  expect_equal(length(expand_keys("B.{a,b,c}.{x,y}")), 6L)
})

test_that("contract_keys reads first-column backticked B. keys from markdown tables only", {
  f <- tempfile(fileext = ".md")
  writeLines(c("# t", "| Key | Meaning |", "|---|---|",
               "| `B.x.{a,b}` | mentions `B.not.a.key` in prose |",
               "| `B.y`, `B.z{,.p}` | two |", "| `ci1` | format row |",
               "Prose `B.prose.key` outside a table."), f)
  expect_equal(sort(contract_keys(f)), sort(c("B.x.a", "B.x.b", "B.y", "B.z", "B.z.p")))
})

test_that("suppress_small blanks values and n where n < 30", {
  dt <- data.table::data.table(g = c("a", "b"), n = c(29L, 30L), v = c(1.5, 2.5))
  s <- suppress_small(dt, "n", "v")
  expect_equal(s$suppressed, c(TRUE, FALSE))
  expect_true(is.na(s$v[1]) && is.na(s$n[1]))
  expect_equal(s$v[2], 2.5)
  expect_equal(s$n[2], 30L)
})

test_that("multi_birth_mothers keeps mothers with at least 2 rows", {
  expect_equal(multi_birth_mothers(c(1, 1, 2, 3, 3, 3)), c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE))
})

# ---------------------------------------------------------------- ruling P14 polish

test_that("with_vcov_fix_flag records and muffles fixest's non-PD vcov fix only", {
  r <- with_vcov_fix_flag({ warning("The VCOV matrix is not positive definite and was 'fixed' (see ?vcov)."); 1 })
  expect_equal(r$value, 1)
  expect_true(r$vcov_fixed)
  r2 <- with_vcov_fix_flag(2)
  expect_false(r2$vcov_fixed)
  expect_warning(with_vcov_fix_flag({ warning("something else"); 3 }), "something else")
})

test_that("unit_of labels scale and per-unit of every coefficient", {
  expect_equal(unit_of("bw_g", "gaussian", "T2_dry"), "grams per +1 SPI unit drier")
  expect_equal(unit_of("lbw", "logit", "T2_tmax_anom"), "odds ratio per +1 deg C")
  expect_true(startsWith(unit_of("lbw", "gaussian", "T2_dry"), "probability (x100 = percentage points) per +1 SPI unit drier"))
  expect_true(startsWith(unit_of("haz", "gaussian", "inter"), "HAZ (z) per unit of preg_elnino_dev x sens_dry"))
  expect_match(unit_of("bw_g", "gaussian", "preg_oni_neg"), "min(x, 0)", fixed = TRUE)
  expect_match(unit_of("haz", "gaussian", "dry"), "per +1 SPI unit drier", fixed = TRUE)
})

test_that("fingerprint_matches compares md5 and size with the conception_shift metadata", {
  meta <- list(b_source_md5 = "abc", b_source_size = "10", b_source_sha256 = "f00")
  expect_true(fingerprint_matches(meta, list(md5 = "abc", size = 10)))
  expect_false(fingerprint_matches(meta, list(md5 = "abd", size = 10)))
  expect_false(fingerprint_matches(meta, list(md5 = "abc", size = 11)))
  expect_false(fingerprint_matches(list(), list(md5 = "abc", size = 10)))
})

test_that("check_aggregate_table: guard names (P12b) and every count-like column", {
  ok <- data.table::data.table(state = "a", n_children = 30L, count_x = 31L, N = 40L, share = 0.1)
  expect_true(check_aggregate_table("ok", ok))
  for (bad in c("latitude", "LONGNUM", "dhsclust", "dhsid", "lng", "anaemia_any", "caseid", "mother_id"))
    expect_error(check_aggregate_table("t", data.table::data.table(x = 1, tmp = 1)[, (bad) := 1]), "forbidden")
  expect_error(check_aggregate_table("t", data.table::data.table(n_children = 29L)), "n < 30")
  expect_error(check_aggregate_table("t", data.table::data.table(count_births = 5L)), "n < 30")
  expect_error(check_aggregate_table("t", data.table::data.table(N = 2L)), "n < 30")
  expect_true(check_aggregate_table("t", data.table::data.table(n_children = NA_integer_, suppressed = TRUE)))
  expect_true(check_aggregate_table("t", data.table::data.table(note = "n<30 rows suppressed")))  # text, not a count
})
