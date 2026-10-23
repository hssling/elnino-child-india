# Paper B ruling Q3: covariate distribution and missingness table (models/B/11_covariates.R).
# Synthetic data only.

local({
  root <- if (exists("MODELS_ROOT")) MODELS_ROOT else normalizePath(file.path(getwd(), "..", "..", ".."), winslash = "/")
  for (f in c("R/fe.R", "R/prep.R", "R/io.R", "11_covariates.R")) source(file.path(root, "B", f), local = globalenv())
})

# cells of one sample x round x variable group
grp <- function(n, sample = "cov", variable = "caste", levels = paste0("l", seq_along(n))) {
  data.table::data.table(sample = sample, round = "nfhs4", variable = variable, level = levels,
                         n = as.integer(n), w = as.numeric(n))
}

test_that("Q3: one small cell triggers complementary suppression; the zero cell stays published", {
  cells <- grp(c(10, 40, 500, 1000, 0), levels = c("sc", "st", "obc", "none", "missing"))
  s <- suppress_cov(cells)
  expect_equal(s$suppressed, c(TRUE, TRUE, FALSE, FALSE, FALSE))     # 40 = smallest cell >= 30
  expect_equal(s$zero, c(FALSE, FALSE, FALSE, FALSE, TRUE))
  # differencing against the published group total gives only the sum of two cells, not 10
  total <- sum(cells$n)
  rec <- total - sum(s$n[!s$suppressed & !s$zero])
  expect_equal(rec, 50L)
  expect_false(rec == 10L)
  # no small cell, no suppression
  expect_false(any(suppress_cov(grp(c(30, 31, 0)))$suppressed))
  # two small cells and no cell >= 30 among them -> a complement is still added
  s2 <- suppress_cov(grp(c(5, 7, 300, 400)))
  expect_equal(s2$suppressed, c(TRUE, TRUE, TRUE, FALSE))
})

test_that("Q3: a nested covbw cell whose difference from cov is 1-29 is suppressed (with a complement)", {
  cov <- grp(c(100, 900, 2000))
  bw <- grp(c(95, 800, 1500), sample = "covbw")                     # 100 - 95 = 5 non-BW children
  s <- suppress_cov(rbind(cov, bw))
  sb <- s[sample == "covbw"]
  expect_true(sb$suppressed[1])
  expect_true(sum(sb$suppressed) >= 2L)                              # complement in the covbw group
  expect_false(any(s[sample == "cov"]$suppressed))
  # the published columns cannot give the covbw cell: every published covbw n differs from cov by 0 or >= 30
  pub <- merge(s[sample == "cov" & !suppressed, .(level, n_cov = n)], s[sample == "covbw" & !suppressed, .(level, n_bw = n)], by = "level")
  d <- pub$n_cov - pub$n_bw
  expect_true(all(d == 0 | d >= 30))
  rec <- sum(bw$n) - sum(sb$n[!sb$suppressed])
  expect_false(rec == 95L)
})

test_that("Q3: cov_level maps codes to the modelled levels and missing", {
  d <- data.table::data.table(b4 = c(1, 2, NA), bord_band = factor(c("1", "4+", NA), levels = c("1", "2", "3", "4+")),
                              v012_band = factor(c("15-19", "35+", "20-24")), v106 = c(0, 3, 9), v190 = c(1, 5, 3),
                              v025 = c(1, 2, 2), s116 = c(1, 8, NA), v130 = c(1, 96, 7))
  expect_equal(cov_level(d, "sex"), c("male", "female", "missing"))
  expect_equal(cov_level(d, "birth_order"), c("1", "4plus", "missing"))
  expect_equal(cov_level(d, "mother_age"), c("15_19", "35plus", "20_24"))
  expect_equal(cov_level(d, "education"), c("none", "higher", "missing"))
  expect_equal(cov_level(d, "wealth"), c("poorest", "richest", "middle"))
  expect_equal(cov_level(d, "caste"), c("sc", "missing", "missing"))
  expect_equal(cov_level(d, "religion"), c("hindu", "other", "other"))
})

test_that("Q3: covariate table shares sum to 100 and keys follow the contract pattern", {
  set.seed(4)
  n <- 4000L
  d <- data.table::data.table(round = factor(rep(c("nfhs4", "nfhs5"), each = n / 2)),
    b4 = sample(1:2, n, TRUE), bord_band = factor(sample(c("1", "2", "3", "4+"), n, TRUE), levels = c("1", "2", "3", "4+")),
    v012_band = factor(sample(c("15-19", "20-24", "25-29", "30-34", "35+"), n, TRUE)),
    v106 = sample(0:3, n, TRUE), v190 = sample(1:5, n, TRUE), v025 = sample(1:2, n, TRUE),
    s116 = sample(c(1:4, 8), n, TRUE, prob = c(.25, .25, .3, .195, .005)), v130 = sample(c(1, 2, 96), n, TRUE),
    weight = stats::runif(n, 0.5, 1.5), s_bw = stats::runif(n) < 0.8)
  r <- covariate_table(d)
  tab <- r$table
  expect_setequal(names(tab), c("sample", "round", "variable", "level", "n", "zero", "weighted_pct", "unweighted_pct", "suppressed"))
  ok <- tab[variable == "wealth" & sample == "cov" & round == "nfhs4"]
  expect_equal(sum(ok$weighted_pct), 100)
  x <- d[round == "nfhs4"]
  expect_equal(ok[level == "poorest", weighted_pct], 100 * sum(x$weight[x$v190 == 1]) / sum(x$weight))
  expect_true(all(is.na(tab$n) | tab$n >= 30))                       # check_aggregate_table rule
  expect_silent(check_aggregate_table("supp_covariates", tab))
  e <- r$entries
  expect_true("B.supp.cov.wealth.poorest.nfhs4" %in% names(e))
  expect_true("B.supp.covbw.caste.missing.nfhs5" %in% names(e))
  expect_true("B.supp.cov.round.nfhs4.all" %in% names(e))
  z <- e[["B.supp.cov.sex.missing.nfhs4"]]
  expect_identical(z$value, 0L)
  expect_identical(z$display, "0 (0.0%)")
  # caste missing (~10 per round) is suppressed with a complement; its key is null with the fixed display
  sup <- tab[suppressed == TRUE]
  expect_true(all(tab[sample == "cov" & variable == "caste" & level == "missing", suppressed]))
  expect_true(all(tab[sample == "cov" & variable == "caste", sum(suppressed), by = round]$V1 >= 2L))
  if (nrow(sup)) {
    k <- sprintf("B.supp.%s.%s.%s.%s", sup$sample[1], sup$variable[1], sup$level[1], sup$round[1])
    expect_true(is.na(e[[k]]$value)); expect_identical(e[[k]]$display, "suppressed (n < 30)")
  }
})
