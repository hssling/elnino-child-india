# Paper B ruling Q5: exposure distributions and rescaled precision (models/B/12_precision.R).
# Synthetic data only.

local({
  root <- if (exists("MODELS_ROOT")) MODELS_ROOT else normalizePath(file.path(getwd(), "..", "..", ".."), winslash = "/")
  for (f in c("R/fe.R", "R/prep.R", "R/io.R", "02_birthweight.R", "12_precision.R")) source(file.path(root, "B", f), local = globalenv())
})

test_that("wquantile is the inverse weighted ECDF", {
  expect_equal(wquantile(1:10, rep(1, 10), c(0.1, 0.5, 0.55, 1)), c(1, 5, 6, 10))
  expect_equal(wquantile(1:10, c(rep(1, 9), 100), 0.5), 10)        # heavy weight on the maximum
  expect_equal(wquantile(c(3, NA, 1, 2), c(1, 1, 1, NA), 0.5), 1)  # non-finite pairs dropped
})

test_that("exp_stats: percentiles, unweighted min/max, weighted SD, IQR and n", {
  x <- c(5, 1, 4, 2, 3); w <- c(1, 1, 1, 1, 1)
  s <- exp_stats(x, w)
  expect_equal(unname(s[c("p25", "p50", "p75", "min", "max")]), c(2, 3, 4, 1, 5))
  expect_equal(s[["iqr"]], 2)
  expect_equal(s[["sd"]], sqrt(mean((x - 3)^2)))
  expect_equal(s[["n"]], 5)
})

test_that("rescale_ci scales the estimate and bounds and swaps them for a negative scale", {
  expect_equal(rescale_ci(-2, -5, 1, 2), list(est = -4, lo = -10, hi = 2))
  expect_equal(rescale_ci(-2, -5, 1, -1), list(est = 2, lo = -1, hi = 5))
})

test_that("refit_sample returns exactly the rows fixest used (singletons dropped)", {
  set.seed(2)
  n <- 600L
  d <- data.table::data.table(region_id = factor(c(sample(sprintf("R%02d", 1:20), n - 1L, TRUE), "R99")),
                              birth_ym = sample(201501:201512, n, TRUE), weight = stats::runif(n, 0.5, 2),
                              s_bw = TRUE, x1 = stats::rnorm(n))
  d[, bw_g := 2800 + 20 * x1 + stats::rnorm(n, 0, 300)]
  m <- list(id = "T", lhs = "bw_g", rhs = "x1", fe = c("region_id", "birth_ym"), sample = "s_bw", terms = "x1")
  r <- refit_sample(d, m)
  expect_equal(nrow(r$data), r$tab$n)
  expect_false("R99" %in% as.character(r$data$region_id))             # the singleton region is not in the sample
  expect_equal(r$tab$term, "x1")
})
