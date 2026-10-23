# Paper B Phase 3b Task 1: supplementary models (09_supplementary.R), figure aggregates
# (10_figure_data.R) and the generated text keys (08_results.R). Synthetic data only.

local({
  root <- if (exists("MODELS_ROOT")) MODELS_ROOT else normalizePath(file.path(getwd(), "..", "..", ".."), winslash = "/")
  for (f in c("R/fe.R", "R/prep.R", "R/io.R", "R/display.R", "05_enso.R", "08_results.R", "09_supplementary.R", "10_figure_data.R"))
    source(file.path(root, "B", f), local = globalenv())
})

CAUSAL_WORDS <- "(?i)\\b(cause[sd]?|causal(ly)?|effects?|impacts?|led to|leads? to|reduc(e|es|ed|ing))\\b"

test_that("test_map_aggregates_suppress_small_n: a district with n = 29 keeps its row with NaN values", {
  set.seed(3)
  mk <- function(reg, n, ym) data.table::data.table(
    region_id = reg, state_gadm = "S1", round = "nfhs5", birth_ym = ym,
    T2_dry = stats::rnorm(n), weight = stats::runif(n, 0.5, 1.5),
    s_haz2459 = TRUE, stunted = as.numeric(stats::runif(n) < 0.3))
  d <- rbind(mk("IND.1.1_1", 29L, 201801L), mk("IND.1.2_1", 40L, 201801L))
  m <- district_map_data(d)
  small <- m[district_gid == "IND.1.1_1" & panel %in% c("T2_dry", "stunting_24_59") & round == "nfhs5"]
  expect_equal(nrow(small), 2L)                                   # rows kept, not dropped
  expect_true(all(is.na(small$value)))
  expect_equal(small$n_band, rep("<30", 2))
  expect_true(all(small$suppressed))
  big <- m[district_gid == "IND.1.2_1" & panel == "T2_dry" & round == "nfhs5"]
  expect_equal(big$n_band, "30\u201399")
  x <- d[region_id == "IND.1.2_1"]
  expect_equal(big$value, sum(x$weight * x$T2_dry) / sum(x$weight))
  expect_false(big$suppressed)
  st <- m[district_gid == "IND.1.2_1" & panel == "stunting_24_59" & round == "nfhs5"]
  expect_equal(st$value, sum(x$weight * x$stunted) / sum(x$weight))
  # full grid: a district with no 2015-16 cohort children still appears, suppressed
  coh <- m[panel == "T2_dry_cohort_2015_16"]
  expect_setequal(coh$district_gid, c("IND.1.1_1", "IND.1.2_1"))
  expect_true(all(coh$suppressed) && all(is.na(coh$value)))
  # the result passes the aggregate guard (no region_id / birth_ym column, n >= 30 or NA)
  expect_silent(check_aggregate_table("fig2_district_map", m))
  expect_false(any(c("region_id", "birth_ym", "n_children") %in% names(m)))
})

# Feasible values of a count from its published band.
band_range <- function(b) switch(b, "<30" = 0:29, "30\u201399" = 30:99, "100\u2013499" = 100:499, "500+" = 500:5000)

test_that("I-1 (Q2): a suppressed fig2 cohort cell cannot be recovered by differencing", {
  # district D: 10 children of the 2015-16 cohort (all NFHS-5), 50 other NFHS-4, 60 other NFHS-5.
  # Identity: n(cohort) + n(other) = n(T2_dry, nfhs4) + n(T2_dry, nfhs5).
  set.seed(8)
  mk <- function(n, rnd, ym) data.table::data.table(region_id = "IND.9.9_1", state_gadm = "S9", round = rnd,
    birth_ym = ym, T2_dry = stats::rnorm(n), weight = stats::runif(n, 0.5, 1.5), s_haz2459 = TRUE,
    stunted = as.numeric(stats::runif(n) < 0.3))
  d <- rbind(mk(10L, "nfhs5", 201510L), mk(50L, "nfhs4", 201301L), mk(60L, "nfhs5", 201801L))
  m <- district_map_data(d)
  # no exact count of any kind is published: the only numeric column is the (weighted) value
  num <- names(m)[vapply(m, is.numeric, logical(1))]
  expect_equal(num, "value")
  b <- function(p, r) m[panel == p & round == r, n_band]
  expect_equal(b("T2_dry_cohort_2015_16", "pooled"), "<30")
  cand <- unique(as.vector(outer(outer(band_range(b("T2_dry", "nfhs4")), band_range(b("T2_dry", "nfhs5")), "+"),
                                 band_range(b("T2_dry_other_cohorts", "pooled")), "-")))
  cand <- intersect(cand, band_range("<30"))
  expect_gt(length(cand), 1L)                      # the published columns admit many cohort counts
  expect_true(10L %in% cand)
})

test_that("I-1 (Q2): fig1 publishes weighted shares, and a suppressed month is not recoverable", {
  set.seed(9)
  months <- c(201001:201012)
  n_m <- c(9L, rep(100L, 11))                      # one month of 9 births, suppressed
  d <- data.table::data.table(round = "nfhs4", birth_ym = rep(months, n_m))
  d[, weight := stats::runif(.N, 0.2, 3)]
  d[, weight := weight * .N / sum(weight)]          # rescaled to the round size, as in Task 5
  enso <- data.table::CJ(year = 2010L, month = 1:12)[, `:=`(oni_anom = 0.1, roni = 0.1, phase = "neutral")]
  f <- fig1_timeline(d, enso = enso, years = c(2010L, 2010L))
  expect_false(any(grepl("^(n_|weighted_)", names(f))))
  expect_true(is.na(f$birth_share_pct_nfhs4[1]) && f$suppressed_nfhs4[1])
  w <- d[, .(w = sum(weight)), by = birth_ym]
  expect_equal(f$birth_share_pct_nfhs4[-1], round(100 * w$w[-1] / sum(d$weight), 1))
  # differencing: round total (published elsewhere) x (100 - shown shares) does not give 9
  rec <- round(nrow(d) * (100 - sum(f$birth_share_pct_nfhs4, na.rm = TRUE)) / 100)
  expect_false(rec == 9L)
})

test_that("band_n: <30 (NA = 0), 30-99, 100-499, 500+", {
  expect_equal(band_n(c(NA, 0, 29, 30, 99, 100, 499, 500, 1e5)),
               c("<30", "<30", "<30", "30\u201399", "30\u201399", "100\u2013499", "100\u2013499", "500+", "500+"))
})

test_that("unit_of labels 1,000-day dryness terms per +1 SPI unit drier (M-6)", {
  for (t in c("dry1000_2459", "b10r_dry1000_2459", "dist_dry1000_2459", "dry", "T2_dry"))
    expect_match(unit_of("haz", "gaussian", t), "per +1 SPI unit drier", fixed = TRUE)
  expect_match(unit_of("bw_g", "gaussian", "T2_oni"), "trimester mean ONI", fixed = TRUE)
})

test_that("Q2b: ENSO x sensitivity terms carry the interaction unit, not 'per +1 SPI unit drier'", {
  inter <- "per unit of preg_elnino_dev x sens_dry (sens_dry = -sens_beta, % JJAS rain loss per deg C ONI)"
  for (t in c("inter", "elnino_x_sens_dry", "elnino_x_sens_dry_xw")) {
    expect_equal(unit_of("bw_g", "gaussian", t), paste("grams", inter))
    expect_equal(unit_of("haz", "gaussian", t), paste("HAZ (z)", inter))
  }
  # every term name passed to run_model by 02-09 gets its intended per-unit text
  want <- list(
    "per +1 SPI unit drier" = c("T1_dry", "T2_dry", "T3_dry", "dry", "recent1_dry", "recent12_dry", "placebo_dry",
                                "T2_dry_sm1", "T3_dry_sp1", "dist_T2_dry", "T2_dry_b10r", "dry1000_2459",
                                "dist_dry1000_2459", "b10r_dry1000_2459"),
    "per +1 deg C" = c("T1_tmax_anom", "T3_tmax_anom_sm1", "dist_T2_tmax_anom", "recent1_tmax_anom",
                       "recent12_tmax_anom", "placebo_tmax_anom"),
    "per +1 deg C of max(x, 0) (El Nino term)" = c("preg_oni_pos", "preg_roni_pos", "preg_nino34_pos", "post_oni_pos"),
    "per +1 deg C of the trimester mean ONI" = c("T1_oni", "T2_oni", "T3_oni"))
  for (per in names(want)) for (t in want[[per]])
    expect_equal(unit_of("bw_g", "gaussian", t), paste("grams", per), info = t)
  for (t in c("preg_oni_neg", "post_oni_neg")) expect_match(unit_of("bw_g", "gaussian", t), "min(x, 0)", fixed = TRUE)
})

test_that("test_conclusion_sentence_matches_verdicts: null template, sign clause and alternative", {
  nc <- "not consistent with"
  s <- conclusion_sentence(h4 = nc, h4_terms = "neither", h5 = nc, xs = nc, enso_haz = 0.031, 0.012, 0.050)
  expect_match(s, "did not", fixed = TRUE)
  expect_match(s, "positive (95% CI excluding 0), opposite in sign to the hypothesis", fixed = TRUE)
  expect_false(grepl(CAUSAL_WORDS, s, perl = TRUE))
  # M-1: "opposite in sign" / "in the direction" only when the CI excludes 0
  expect_match(conclusion_sentence(nc, "neither", nc, nc, -0.02, -0.03, -0.01),
               "negative (95% CI excluding 0), in the direction of the hypothesis", fixed = TRUE)
  wide_pos <- conclusion_sentence(nc, "neither", nc, nc, 0.02, -0.01, 0.05)
  expect_match(wide_pos, "positive, with a 95% CI including 0.", fixed = TRUE)
  expect_false(grepl("opposite", wide_pos, fixed = TRUE))
  expect_match(conclusion_sentence(nc, "neither", nc, nc, -0.02, -0.05, 0.01), "negative, with a 95% CI including 0", fixed = TRUE)
  expect_match(conclusion_sentence(nc, "neither", nc, nc, 0.02), "not estimable", fixed = TRUE)   # no CI given
  expect_match(conclusion_sentence(nc, "neither", nc, nc, NA), "not estimable", fixed = TRUE)

  a <- conclusion_sentence(h4 = "consistent with", h4_terms = "T2", h5 = nc, xs = nc, enso_haz = 0.031, 0.012, 0.05)
  expect_false(grepl("did not show", a, fixed = TRUE))
  expect_match(a, "drier conditions in the second trimester of pregnancy were associated with lower birth weight", fixed = TRUE)
  expect_false(grepl("T2", a, fixed = TRUE))                                 # M-3: ordinal, not "T2"
  expect_match(a, "consistent with the pre-specified hypothesis;", fixed = TRUE)   # M-3: singular
  expect_match(a, "the height-for-age hypothesis (H5) and the El Ni\u00f1o \u00d7 district-sensitivity hypothesis did not meet the pre-specified rule", fixed = TRUE)
  expect_false(grepl(CAUSAL_WORDS, a, perl = TRUE))
  expect_match(conclusion_sentence("consistent with", "T3", nc, nc, 0), "in the third trimester of pregnancy", fixed = TRUE)
  all_ok <- conclusion_sentence("consistent with", "T2 and T3", "consistent with", "consistent with", -0.01, -0.02, -0.005)
  expect_match(all_ok, "in the second and third trimesters of pregnancy", fixed = TRUE)
  expect_match(all_ok, "consistent with the pre-specified hypotheses;", fixed = TRUE)
  expect_false(grepl("did not", all_ok, fixed = TRUE))

  v <- verdict_summary(h4 = nc, h4_terms = "neither", h5 = nc, xs = nc)
  expect_match(v, "the birth-weight results were not consistent with H4 (neither the second- nor the third-trimester term met the rule)", fixed = TRUE)
  expect_match(v, "and the results did not support the pre-specified El Ni\u00f1o \u00d7 district-sensitivity hypothesis.", fixed = TRUE)
  expect_false(grepl("pathway", v, fixed = TRUE))                            # M-2
  expect_match(verdict_summary(nc, "neither", nc, "consistent with"), "and the results were consistent with the pre-specified El Niño", fixed = TRUE)
  expect_error(conclusion_sentence("maybe", "neither", nc, nc, 0), "verdict")
})

test_that("build_text_keys reads the verdict keys and the ENSO-total HAZ value", {
  mk <- function(v, desc = "x") list(value = v, display = as.character(v), script = "x", description = desc)
  e <- list(B.h4.verdict = mk("not consistent with"), B.h4.verdict_terms = mk("neither"),
            B.h5.verdict = mk("not consistent with"), B.ensoxs.verdict = mk("not consistent with"),
            B.enso.oni.haz2459.elnino = mk(0.0307, "ENSO-total ...; per +1 deg C of the El Nino term max(x,0); lo=0.0115 hi=0.0499"))
  expect_equal(ci_from_description(e$B.enso.oni.haz2459.elnino$description), c(lo = 0.0115, hi = 0.0499))
  k <- build_text_keys(e)
  expect_setequal(names(k), c("B.abstract.conclusion_sentence", "B.discussion.verdict_summary"))
  expect_identical(k$B.abstract.conclusion_sentence$value, k$B.abstract.conclusion_sentence$display)
  expect_match(k$B.abstract.conclusion_sentence$value, "did not show", fixed = TRUE)
  expect_match(k$B.abstract.conclusion_sentence$value, "positive (95% CI excluding 0)", fixed = TRUE)
  e$B.enso.oni.haz2459.elnino$description <- "x; lo=-0.002 hi=0.05"
  expect_match(build_text_keys(e)$B.abstract.conclusion_sentence$value, "positive, with a 95% CI including 0", fixed = TRUE)
  e$B.h5.verdict <- NULL
  expect_error(build_text_keys(e), "B.h5.verdict")
})

test_that("test_trim_oni_terms_present: ENSO-trim returns rows T1_oni, T2_oni, T3_oni", {
  set.seed(11)
  n <- 3000L
  d <- data.table::data.table(
    region_id = factor(sample(sprintf("R%02d", 1:40), n, TRUE)),
    birth_ym = sample(c(201501:201512, 201601:201612, 201701:201712), n, TRUE),
    round = factor(sample(c("nfhs4", "nfhs5"), n, TRUE)),
    T1_oni = stats::rnorm(n), T2_oni = stats::rnorm(n), T3_oni = stats::rnorm(n),
    cmc_ns1 = stats::runif(n), weight = stats::runif(n, 0.5, 2), s_bw = TRUE)
  d[, birth_month_of_year := birth_ym %% 100L]
  d[, bw_g := 2800 + 10 * T2_oni + stats::rnorm(n, 0, 400)]
  t <- fit_trim_oni(d, "bw", "ENSO-trim-test", covars = "cmc_ns1")
  expect_equal(t$term, c("T1_oni", "T2_oni", "T3_oni"))
  expect_true(all(is.finite(t$est) & t$lo < t$est & t$est < t$hi))
  expect_equal(unique(t$fe), "region_id+birth_month_of_year+round")
  expect_equal(unique(t$cluster), "region_id+birth_ym")
})

test_that("test_hw57_kappa_bounds: kappa in [-1, 1]; identical vectors give 1", {
  set.seed(5)
  for (i in 1:20) {
    x <- stats::runif(200) < 0.4; y <- stats::runif(200) < stats::runif(1)
    k <- cohen_kappa(x, y)
    expect_true(k >= -1 && k <= 1)
  }
  x <- stats::runif(100) < 0.5
  expect_equal(cohen_kappa(x, x), 1)
  expect_equal(cohen_kappa(c(TRUE, FALSE, TRUE, FALSE), c(FALSE, TRUE, FALSE, TRUE)), -1)
  expect_equal(cohen_kappa(rep(TRUE, 5), rep(TRUE, 5)), 1)       # pe = 1, po = 1
  expect_error(cohen_kappa(c(TRUE, NA), c(TRUE, FALSE)), "missing")
})

test_that("hw57 classification: levels 1-3 anaemic, 4 not, other codes excluded", {
  expect_equal(hw57_anaemic(c(1, 2, 3, 4, 9, NA)), c(TRUE, TRUE, TRUE, FALSE, NA, NA))
})
