# Paper B fixest layer (models/B/R/fe.R): weighted, two-way-clustered FE models and
# fixest mediation. Simulated data only.

local({
  root <- if (exists("MODELS_ROOT")) MODELS_ROOT else normalizePath(file.path(getwd(), "..", "..", ".."), winslash = "/")
  source(file.path(root, "B", "R", "fe.R"), local = globalenv())
})

# 40 districts x 60 birth months, several children per cell. Exposure x varies by child and
# carries district and month components that are correlated with the district and month
# effects (so omitting the FE biases beta) plus a district-month shock; the
# error has its own, independent district-month shock (both cluster-correlated, so clustering
# matters, without biasing beta). y = alpha_d + gamma_t + beta * x + e.
sim_fe <- function(n_dist = 40, n_month = 60, per_cell = 3, beta = -25, seed = 11) {
  set.seed(seed)
  dist <- data.table::data.table(region_id = sprintf("D%02d", seq_len(n_dist)), xd = stats::rnorm(n_dist))
  dist[, alpha := 3000 + 150 * xd + stats::rnorm(.N, 0, 100)]
  mon <- data.table::data.table(birth_ym = sprintf("2015-%03d", seq_len(n_month)), xm = stats::rnorm(n_month))
  mon[, gamma := 80 * xm + stats::rnorm(.N, 0, 50)]
  d <- data.table::CJ(region_id = dist$region_id, birth_ym = mon$birth_ym, k = seq_len(per_cell))
  d <- dist[d, on = "region_id"]
  d <- mon[d, on = "birth_ym"]
  d[, `:=`(ux = stats::rnorm(1), ue = stats::rnorm(1)), by = .(region_id, birth_ym)]
  d[, x := xd + xm + 0.8 * ux + stats::rnorm(.N)]
  d[, y := alpha + gamma + beta * x + 60 * ue + stats::rnorm(.N, 0, 200)]
  d[, weight := stats::runif(.N, 0.5, 1.5)]
  d[, weight := weight / mean(weight)]
  d[]
}

test_that("fit_fe recovers the known beta within 3 SE (weighted, two-way FE)", {
  d <- sim_fe()
  fit <- fit_fe(y ~ x | region_id + birth_ym, d)
  expect_s3_class(fit, "fixest")
  est <- unname(stats::coef(fit)["x"]); se <- unname(fixest::se(fit)["x"])
  expect_lt(abs(est - (-25)), 3 * se)
  expect_equal(stats::nobs(fit), nrow(d))
  # the FE matter: without them beta is biased by more than 3 SE
  nofe <- fit_fe(y ~ x, d)
  expect_gt(abs(unname(stats::coef(nofe)["x"]) - (-25)), 3 * unname(fixest::se(nofe)["x"]))
})

test_that("SEs are two-way clustered by region_id and birth_ym, not iid", {
  d <- sim_fe()
  fit <- fit_fe(y ~ x | region_id + birth_ym, d)
  se_fit <- unname(fixest::se(fit)["x"])
  se_iid <- unname(fixest::se(summary(fit, vcov = "iid"))["x"])
  se_cl <- unname(fixest::se(summary(fit, cluster = ~region_id + birth_ym))["x"])
  expect_gt(abs(se_fit - se_iid) / se_iid, 0.05)
  expect_equal(se_fit, se_cl, tolerance = 1e-10)
  # independently specified reference: weighted feols, two-way clustered
  ref <- fixest::feols(y ~ x | region_id + birth_ym, d, weights = ~weight, cluster = ~region_id + birth_ym)
  expect_equal(se_fit, unname(fixest::se(ref)["x"]), tolerance = 1e-10)
  expect_equal(stats::coef(fit), stats::coef(ref), tolerance = 1e-10)
  expect_match(attr(fit$cov.scaled, "vcov_type"), "region_id")
  expect_match(attr(fit$cov.scaled, "vcov_type"), "birth_ym")
  # fe_table uses the clustered SE
  tab <- fe_table(fit, "x")
  expect_equal(tab$est, unname(stats::coef(fit)["x"]))
  expect_equal(c(tab$lo, tab$hi), unname(unlist(stats::confint(fit)["x", ])), tolerance = 1e-10)
  expect_equal(tab$p, unname(fixest::pvalue(fit)["x"]), tolerance = 1e-10)
  expect_identical(tab$n, as.integer(nrow(d)))
  expect_identical(names(tab), c("term", "est", "lo", "hi", "p", "n", "vcov", "df"))
  expect_identical(tab$vcov, "Clustered (region_id & birth_ym)")
  expect_equal(tab$df, 39)                                           # t(min(G) - 1), G = 40 districts
  expect_equal(tab$hi - tab$est, stats::qt(0.975, 39) * se_fit, tolerance = 1e-8)
})

test_that("fe_table refuses a fit whose vcov is not two-way clustered, unless cluster = NULL", {
  d <- sim_fe(n_dist = 20, n_month = 24)
  fit <- fit_fe(y ~ x | region_id + birth_ym, d)
  expect_error(fe_table(summary(fit, vcov = "iid"), "x"), "not the vcov clustered")
  dflt <- fixest::feols(y ~ x | region_id + birth_ym, d, weights = ~weight)     # no vcov argument (0.14: iid)
  expect_error(fe_table(dflt, "x"), "not the vcov clustered")
  one <- fixest::feols(y ~ x | region_id + birth_ym, d, weights = ~weight, cluster = ~region_id)
  expect_error(fe_table(one, "x"), "not the vcov clustered")
  two <- fixest::feols(y ~ x | region_id + birth_ym, d, weights = ~weight, cluster = ~region_id + birth_ym)
  expect_equal(fe_table(two, "x")$lo, fe_table(fit, "x")$lo, tolerance = 1e-10)
  opt <- fe_table(summary(fit, vcov = "iid"), "x", cluster = NULL)               # explicit opt-out
  expect_identical(opt$vcov, "IID")
  expect_error(fe_table(one, "x", cluster = ~region_id), NA)                     # matching one-way spec
})

test_that("weights are applied: they change the estimate when they correlate with the effect", {
  d <- sim_fe(seed = 12)
  d[, grp := as.integer(substr(region_id, 2, 3)) <= 20]             # half the districts
  d[, y := alpha + gamma + ifelse(grp, -40, -5) * x + stats::rnorm(.N, 0, 200)]
  d[, weight := ifelse(grp, 4, 0.25)]                               # heavy where the effect is large
  w <- fit_fe(y ~ x | region_id + birth_ym, d)
  u <- fit_fe(y ~ x | region_id + birth_ym, d, weights = NULL)
  bw <- unname(stats::coef(w)["x"]); bu <- unname(stats::coef(u)["x"])
  expect_lt(bw, bu - 5)                                              # weighted pulled towards -40
  expect_equal(bw, unname(stats::coef(fixest::feols(y ~ x | region_id + birth_ym, d, weights = ~weight))["x"]),
               tolerance = 1e-10)
  expect_equal(bu, unname(stats::coef(fixest::feols(y ~ x | region_id + birth_ym, d))["x"]), tolerance = 1e-10)
})

test_that("fit_fe refuses missing, non-positive or non-finite weights", {
  d <- sim_fe(n_dist = 10, n_month = 12)
  expect_error(fit_fe(y ~ x | region_id + birth_ym, d, weights = ~nope), "weight")
  d2 <- data.table::copy(d); d2[1, weight := 0]
  expect_error(fit_fe(y ~ x | region_id + birth_ym, d2), "weight")
  d3 <- data.table::copy(d); d3[1, weight := NA]
  expect_error(fit_fe(y ~ x | region_id + birth_ym, d3), "weight")
})

test_that("logit fit_fe and fe_table return odds ratios", {
  d <- sim_fe(n_dist = 40, n_month = 24, per_cell = 8, seed = 13)
  d[, eta := -0.5 + 0.4 * x + (alpha - 3000) / 300 + gamma / 160]
  d[, z := stats::rbinom(.N, 1, stats::plogis(eta))]
  fit <- fit_fe(z ~ x | region_id + birth_ym, d, family = "logit")
  expect_identical(fit$family$link, "logit")
  tab <- fe_table(fit, "x")
  b <- unname(stats::coef(fit)["x"]); ci <- unname(unlist(stats::confint(fit)["x", ]))
  expect_equal(tab$est, exp(b), tolerance = 1e-10)
  expect_equal(c(tab$lo, tab$hi), exp(ci), tolerance = 1e-10)
  expect_true(tab$lo < exp(0.4) && exp(0.4) < tab$hi)
  expect_true(tab$lo > 1)
  expect_identical(tab$df, Inf)                                       # feglm: normal reference
  expect_equal(log(tab$hi) - log(tab$est), stats::qnorm(0.975) * unname(fixest::se(fit)["x"]), tolerance = 1e-8)
  expect_match(attr(fit$cov.scaled, "vcov_type"), "birth_ym")
  # gaussian table stays on the coefficient scale
  g <- fe_table(fit_fe(y ~ x | region_id + birth_ym, d), "x")
  expect_lt(g$est, 0)
})

test_that("fe_table errors on an unknown term", {
  d <- sim_fe(n_dist = 10, n_month = 12)
  fit <- fit_fe(y ~ x | region_id + birth_ym, d)
  expect_error(fe_table(fit, c("x", "nope")), "nope")
})

# ---- mediation ----

sim_med <- function(seed = 21) {
  d <- sim_fe(seed = seed)
  set.seed(seed)
  d[, m := 0.5 * x + 0.3 * xm + stats::rnorm(.N)]                    # a = 0.5
  d[, cov1 := stats::rnorm(.N)]
  d[, y := alpha + gamma + 10 * x + 40 * m + 5 * cov1 + stats::rnorm(.N, 0, 200)]   # c' = 10, b = 40
  d[]
}

test_that("mediate_fe recovers a*b and its cluster-bootstrap CI covers the truth (B = 40)", {
  d <- sim_med()
  r <- mediate_fe(d, "y", "x", "m", fe = c("region_id", "birth_ym"), covars = "cov1", B = 40, cores = 1)
  expect_lt(abs(r$a - 0.5), 0.05)
  expect_lt(abs(r$b - 40), 4)
  expect_lt(abs(r$indirect - 20), 3)
  expect_equal(r$indirect, r$a * r$b)
  expect_true(r$ci$indirect[1] < 20 && 20 < r$ci$indirect[2])
  expect_identical(r$flag, "ok")
  expect_identical(r$prop_flag, "ok")
  expect_equal(r$prop, r$indirect / r$c_total)
  # linear: c_total = c_direct + a*b exactly on the same rows and weights
  expect_equal(r$c_total, r$c_direct + r$indirect, tolerance = 1e-8)
  expect_identical(r$B_ok, 40L)
  expect_equal(nrow(r$boot), 40L)
  q <- function(v) stats::quantile(v, c(0.025, 0.975), names = FALSE)
  expect_equal(r$ci$indirect, q(r$boot$indirect))
  expect_equal(r$ci$prop, q(r$boot$prop))
})

test_that("mediate_fe point estimates equal the separate weighted fixest fits", {
  d <- sim_med()
  r <- mediate_fe(d, "y", "x", "m", fe = c("region_id", "birth_ym"), covars = "cov1", B = 2, cores = 1)
  fa <- fixest::feols(m ~ x + cov1 | region_id + birth_ym, d, weights = ~weight)
  fb <- fixest::feols(y ~ x + m + cov1 | region_id + birth_ym, d, weights = ~weight)
  ft <- fixest::feols(y ~ x + cov1 | region_id + birth_ym, d, weights = ~weight)
  expect_equal(r$a, unname(stats::coef(fa)["x"]), tolerance = 1e-10)
  expect_equal(r$b, unname(stats::coef(fb)["m"]), tolerance = 1e-10)
  expect_equal(r$c_direct, unname(stats::coef(fb)["x"]), tolerance = 1e-10)
  expect_equal(r$c_total, unname(stats::coef(ft)["x"]), tolerance = 1e-10)
})

test_that("mediate_fe fits all three models on the same complete rows", {
  d <- sim_med()
  d[1, m := NA]; d[2, y := NA]; d[3, cov1 := NA]
  r <- mediate_fe(d, "y", "x", "m", fe = c("region_id", "birth_ym"), covars = "cov1", B = 2, cores = 1)
  expect_identical(r$n, nrow(d) - 3L)
  expect_equal(unname(r$n_by_model), rep(nrow(d) - 3L, 3))
  expect_setequal(r$dropped$variable, c("y", "m", "cov1"))
  expect_equal(sum(r$dropped$n), 3L)
})

test_that("mediate_fe resamples whole districts as distinct clusters, seeded seed + i - 1", {
  d <- sim_med()
  b1 <- .fe_bootstrap_units(d, "region_id", seed = 3)
  expect_equal(length(unique(b1$region_id)), length(unique(d$region_id)))
  expect_true(all(grepl("_[0-9]+$", b1$region_id)))
  k <- as.character(b1$region_id)[1]
  expect_equal(b1[region_id == k]$y, d[region_id == sub("_[0-9]+$", "", k)]$y)
  set.seed(99); before <- .Random.seed
  r1 <- mediate_fe(d, "y", "x", "m", fe = c("region_id", "birth_ym"), covars = "cov1", B = 3, cores = 1, seed = 5)
  expect_identical(.Random.seed, before)
  r2 <- mediate_fe(d, "y", "x", "m", fe = c("region_id", "birth_ym"), covars = "cov1", B = 3, cores = 1, seed = 5)
  expect_identical(r1$boot, r2$boot)
  # resample i is bootstrap draw with seed + i - 1
  r3 <- mediate_fe(d, "y", "x", "m", fe = c("region_id", "birth_ym"), covars = "cov1", B = 1, cores = 1, seed = 6)
  expect_equal(r3$boot$a, r1$boot$a[2], tolerance = 1e-12)
  # boot draw 2 is the weighted a-path fit on the resample with seed 6
  di <- .fe_bootstrap_units(d, "region_id", seed = 6)
  aw <- unname(stats::coef(fixest::feols(m ~ x + cov1 | region_id + birth_ym, di, weights = ~weight))["x"])
  au <- unname(stats::coef(fixest::feols(m ~ x + cov1 | region_id + birth_ym, di))["x"])
  expect_equal(r1$boot$a[2], aw, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(aw, au, tolerance = 1e-8)))
})

test_that("mediate_fe: more than 10% failed resamples gives NA CIs and the flags", {
  d <- sim_med()
  n <- 0
  flaky <- function(formula, data, weights) {
    n <<- n + 1; if (n > 3) stop("fail") else fe_ols(formula, data, weights)
  }
  r <- mediate_fe(d, "y", "x", "m", fe = c("region_id", "birth_ym"), covars = "cov1", B = 5, cores = 1,
                  fitter = flaky)
  expect_identical(r$flag, "too_many_failed_fits")
  expect_identical(r$prop_flag, "not_assessed")
  expect_true(all(is.na(unlist(r$ci))))
  expect_false(is.na(r$prop))
  expect_identical(r$B_failed, 5L)
})

test_that("mediate_fe: proportion CI withheld when the total-effect CI includes 0", {
  d <- sim_med()
  d[, y := alpha + gamma + 0 * x + 0.5 * m + stats::rnorm(.N, 0, 400)]   # tiny, noisy total effect
  r <- mediate_fe(d, "y", "x", "m", fe = c("region_id", "birth_ym"), covars = character(), B = 20, cores = 1)
  expect_identical(r$flag, "ok")
  expect_true(r$prop_flag %in% c("total_effect_ci_includes_0", "opposite_signs"))
  expect_true(all(is.na(r$ci$prop)))
})

test_that("mediate_fe parallel bootstrap (PSOCK, cores = 2) reproduces the serial result", {
  d <- sim_med()
  s <- mediate_fe(d, "y", "x", "m", fe = c("region_id", "birth_ym"), covars = "cov1", B = 4, cores = 1, seed = 9)
  p <- mediate_fe(d, "y", "x", "m", fe = c("region_id", "birth_ym"), covars = "cov1", B = 4, cores = 2, seed = 9)
  expect_identical(p$boot, s$boot)
  expect_identical(p$ci, s$ci)
  expect_identical(p$cores, 2L)
  closure <- local({ k <- 0; function(formula, data, weights) { k <<- k + 1; fe_ols(formula, data, weights) } })
  expect_error(mediate_fe(d, "y", "x", "m", fe = "region_id", B = 2, cores = 2, fitter = closure), "closure")
})

test_that("mediate_fe validates its inputs", {
  d <- sim_med()
  expect_error(mediate_fe(d, "y", "x", "m", fe = "region_id", B = 0, cores = 1), "B")
  expect_error(mediate_fe(d, "y", "x", "m", fe = "region_id", B = 2, cores = 0), "cores")
  d2 <- data.table::copy(d)[, region_id := NULL]
  expect_error(mediate_fe(d2, "y", "x", "m", fe = "birth_ym", B = 2, cores = 1), "region_id")
})

test_that("the fit guard rejects an iid vcov, a one-way cluster and unapplied weights", {
  d <- sim_fe(n_dist = 10, n_month = 12)
  iid <- fixest::feols(y ~ x | region_id + birth_ym, d, weights = ~weight, vcov = "iid")
  expect_error(.fe_check_fit(iid, ~region_id + birth_ym, ~weight), "not the vcov clustered")
  one <- fixest::feols(y ~ x | region_id + birth_ym, d, weights = ~weight, cluster = ~region_id)
  expect_error(.fe_check_fit(one, ~region_id + birth_ym, ~weight), "not the vcov clustered")
  unw <- fixest::feols(y ~ x | region_id + birth_ym, d, cluster = ~region_id + birth_ym)
  expect_error(.fe_check_fit(unw, ~region_id + birth_ym, ~weight), "weights were not applied")
  expect_silent(.fe_check_fit(unw, ~region_id + birth_ym, NULL))
})

test_that("mediate_fe keeps singleton FE levels (fixef.rm = 'none'); fit_fe drops them (fixest default)", {
  d <- sim_med()
  extra <- data.table::copy(d[1])[, region_id := "D99"]                 # a singleton district
  d2 <- rbind(d, extra)
  r <- mediate_fe(d2, "y", "x", "m", fe = c("region_id", "birth_ym"), covars = "cov1", B = 2, cores = 1)
  expect_identical(r$n, nrow(d2))
  expect_equal(unname(r$n_by_model), rep(nrow(d2), 3))
  expect_equal(stats::nobs(fit_fe(y ~ x | region_id + birth_ym, d2)), nrow(d2) - 1L)
  dropper <- function(formula, data, weights) fixest::feols(formula, data = data, weights = weights, notes = FALSE)
  expect_error(mediate_fe(d2, "y", "x", "m", fe = c("region_id", "birth_ym"), B = 2, cores = 1, fitter = dropper),
               "a model dropped rows")
})

test_that("mediate_fe 10% boundary: 4 of 40 failed resamples keep the CI, 5 of 40 do not", {
  d <- sim_med()
  failing <- function(k) {
    n_a <- 0
    function(formula, data, weights) {
      if (all.vars(formula)[1] == "m") {                       # a-path: first fit of each resample
        n_a <<- n_a + 1
        if (n_a >= 2 && n_a <= k + 1) stop("planned failure")  # call 1 is the point estimate
      }
      fe_ols(formula, data, weights)
    }
  }
  r4 <- mediate_fe(d, "y", "x", "m", fe = c("region_id", "birth_ym"), B = 40, cores = 1, fitter = failing(4))
  expect_identical(r4$B_failed, 4L)
  expect_identical(r4$flag, "ok")
  expect_false(anyNA(r4$ci$indirect))
  r5 <- mediate_fe(d, "y", "x", "m", fe = c("region_id", "birth_ym"), B = 40, cores = 1, fitter = failing(5))
  expect_identical(r5$B_failed, 5L)
  expect_identical(r5$flag, "too_many_failed_fits")
  expect_true(all(is.na(unlist(r5$ci))))
})

# Time-only exposure (ENSO-like): x and a common shock s vary by month only, identical across
# districts; months come in blocks ("ENSO years") of 4. True a = 0.5.
sim_time <- function(n_dist = 20, n_year = 30, per_year = 4, seed = 31) {
  set.seed(seed)
  n_t <- n_year * per_year
  mon <- data.table::data.table(birth_ym = sprintf("M%03d", seq_len(n_t)),
                                enso_year = sprintf("Y%02d", rep(seq_len(n_year), each = per_year)),
                                x = stats::rnorm(n_t), s = stats::rnorm(n_t))
  d <- data.table::CJ(region_id = sprintf("D%02d", seq_len(n_dist)), birth_ym = mon$birth_ym)
  d <- mon[d, on = "birth_ym"]
  d[, alpha := stats::rnorm(1), by = region_id]
  d[, m := 0.5 * x + s + stats::rnorm(.N, 0, 0.5)]
  d[, y := alpha + 1 * x + 2 * m + stats::rnorm(.N)]
  d[, weight := 1]
  d[]
}

test_that("time-only exposure: the time-block bootstrap is calibrated, the district bootstrap is not", {
  d <- sim_time()
  xt <- unique(d[, .(birth_ym, x)])$x
  sd_true <- sqrt((1 + 0.5^2 / 20) / sum((xt - mean(xt))^2))      # analytic SD of a-hat given x
  blk <- mediate_fe(d, "y", "x", "m", fe = "region_id", B = 100, cores = 1, boot_unit = "enso_year")
  dis <- mediate_fe(d, "y", "x", "m", fe = "region_id", B = 100, cores = 1, boot_unit = "region_id")
  expect_identical(blk$boot_unit, "enso_year")
  expect_equal(blk$a, dis$a)                                       # same point estimate
  sd_blk <- stats::sd(blk$boot$a); sd_dis <- stats::sd(dis$boot$a)
  expect_true(sd_blk > sd_true / 2 && sd_blk < sd_true * 2)
  expect_true(blk$ci$a[1] < 0.5 && 0.5 < blk$ci$a[2])
  expect_lt(sd_dis, sd_true / 3)                                   # district bootstrap overconfident
})

test_that("boot_unit relabels the unit and the FEs nested in it, not crossing FEs", {
  d <- sim_time()
  expect_identical(.fe_nested_fe(d, "enso_year", c("region_id", "birth_ym")), "birth_ym")
  expect_identical(.fe_nested_fe(d, "region_id", c("region_id", "birth_ym")), "region_id")
  b <- .fe_bootstrap_units(d, "enso_year", seed = 4, relabel = c("enso_year", "birth_ym"))
  expect_equal(data.table::uniqueN(b$enso_year), 30L)
  expect_equal(data.table::uniqueN(b$birth_ym), 120L)              # each copy has its own months
  expect_identical(sub(".*_", "", as.character(b$birth_ym)), sub(".*_", "", as.character(b$enso_year)))
  expect_setequal(as.character(b$region_id), unique(d$region_id))  # crossing FE untouched
  d[, x2 := x + stats::rnorm(.N)]                                  # varies within month (x alone is collinear with birth_ym)
  r <- mediate_fe(d, "y", "x2", "m", fe = c("region_id", "birth_ym"), B = 2, cores = 1, boot_unit = "enso_year")
  expect_identical(r$relabelled, "birth_ym")
  expect_error(mediate_fe(d, "y", "x", "m", fe = "region_id", B = 2, cores = 1, boot_unit = "nope"), "boot_unit")
})
