# Paper B ruling Q5: exposure distributions in each main model's estimation sample and the
# primary null estimates rescaled to interpretable exposure contrasts (results_keys.md
# "Amendment 2026-10-23 (ruling Q5)", section 11). Descriptive, no verdict.
# Each main model is refitted exactly as in 02-05 (fit_fe / fe_table, same data, weights and
# cluster); fixest::obs() gives the rows it used. check_precision() verifies, on every run, that
# the refits reproduce primary_models.csv / enso_total.csv.

EXPDIST_STATS <- c("p5", "p10", "p25", "p50", "p75", "p90", "p95", "min", "max", "sd")
EXPDIST_P <- c(p5 = 0.05, p10 = 0.10, p25 = 0.25, p50 = 0.50, p75 = 0.75, p90 = 0.90, p95 = 0.95)

# Weighted percentile: smallest x whose cumulative normalised weight is >= p (inverse weighted ECDF).
wquantile <- function(x, w, p) {
  ok <- is.finite(x) & is.finite(w)
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; cw <- cumsum(w[o]) / sum(w)
  vapply(p, function(pp) x[which(cw >= pp - 1e-12)[1]], numeric(1))
}

exp_stats <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  q <- wquantile(x[ok], w[ok], EXPDIST_P)
  out <- c(q, min = min(x[ok]), max = max(x[ok]), sd = unname(wmean_sd(x[ok], w[ok])["sd"]))
  c(out, iqr = out[["p75"]] - out[["p25"]], n = sum(ok))
}

# s * (est, lo, hi); bounds swapped when s < 0.
rescale_ci <- function(est, lo, hi, s) {
  b <- c(s * lo, s * hi)
  list(est = s * est, lo = min(b), hi = max(b))
}

PREC_MODELS <- list(
  h4bw = list(id = "H4-bw", lhs = "bw_g", rhs = c(H4_RHS, COV), fe = c("region_id", "birth_ym"),
              sample = "s_bw", terms = c("T2_dry", "T3_dry")),
  h5haz2459 = list(id = "H5-haz2459", lhs = "haz", rhs = c("dry", COV), fe = c("region_id", "birth_ym", "age_fe"),
                   sample = "s_haz2459", dry = "dry1000_2459", terms = "dry"),
  h5haz0023 = list(id = "H5-haz0023", lhs = "haz", rhs = c("dry", COV), fe = c("region_id", "birth_ym", "age_fe"),
                   sample = "s_haz0023", dry = "dry1000_0023", terms = "dry"),
  h6diarr = list(id = "H6-diarr", lhs = "diarrhoea", rhs = c("recent1_dry", "recent1_tmax_anom", COV),
                 fe = c("region_id", "interview_ym", "age_fe"), sample = "s_diarr", family = "logit", terms = "recent1_dry"),
  h6anaem = list(id = "H6-anaem", lhs = "anaemia_any", rhs = c("recent12_dry", "recent12_tmax_anom", COV),
                 fe = c("region_id", "interview_ym", "age_fe"), sample = "s_anaem", family = "logit", terms = "recent12_dry"),
  ensobw = list(id = "ENSO-total-oni-bw", lhs = "bw_g", rhs = c("preg_oni_pos", "preg_oni_neg", COV, CMC_NS),
                fe = c("region_id", "birth_month_of_year", "round"), sample = "s_bw", terms = "preg_oni_pos"),
  xsbw = list(id = "ENSOxS-bw", lhs = "bw_g", rhs = c("elnino_x_sens_dry", COV), fe = c("region_id", "birth_ym"),
              sample = "s_bw", terms = "elnino_x_sens_dry", table_term = "inter"),
  xshaz = list(id = "ENSOxS-haz2459", lhs = "haz", rhs = c("elnino_x_sens_dry", COV), fe = c("region_id", "birth_ym", "age_fe"),
               sample = "s_haz2459", terms = "elnino_x_sens_dry", table_term = "inter"))

# exposure key -> (model, column in that model's data)
EXPDIST_MAP <- list(
  T1_dry = c("h4bw", "T1_dry"), T2_dry = c("h4bw", "T2_dry"), T3_dry = c("h4bw", "T3_dry"),
  T1_tmax = c("h4bw", "T1_tmax_anom"), T2_tmax = c("h4bw", "T2_tmax_anom"), T3_tmax = c("h4bw", "T3_tmax_anom"),
  dry1000_2459 = c("h5haz2459", "dry"), dry1000_0023 = c("h5haz0023", "dry"),
  recent1_dry = c("h6diarr", "recent1_dry"), recent1_tmax = c("h6diarr", "recent1_tmax_anom"),
  recent12_dry = c("h6anaem", "recent12_dry"), recent12_tmax = c("h6anaem", "recent12_tmax_anom"),
  preg_oni = c("ensobw", "preg_oni"), sens_dry_bw = c("xsbw", "sens_dry"), sens_dry_haz2459 = c("xshaz", "sens_dry"))

# Refit one main model; return its estimation sample (rows of its data) and its table rows.
refit_sample <- function(d, m) {
  x <- d[get(m$sample) == TRUE]
  if (!is.null(m$dry)) x[, dry := get(m$dry)]
  fam <- if (is.null(m$family)) "gaussian" else m$family
  r <- with_vcov_fix_flag({
    fit <- fit_fe(fml(m$lhs, m$rhs, m$fe), x, weights = ~weight, cluster = B_CLUSTER, family = fam)
    list(fit = fit, tab = fe_table(fit, m$terms, cluster = B_CLUSTER))
  })$value
  idx <- fixest::obs(r$fit)
  if (length(idx) != stats::nobs(r$fit)) stop("refit_sample: obs() and nobs() disagree for ", m$id, call. = FALSE)
  tab <- r$tab
  if (!is.null(m$table_term)) tab[, term := m$table_term]
  list(data = x[idx], tab = tab[, .(model = m$id, term, est, lo, hi, n)])
}

run_precision <- function() {
  S <- "models/B/12_precision.R"
  d <- load_b()
  fits <- lapply(PREC_MODELS, function(m) { log_msg("  precision refit", m$id); refit_sample(d, m) })
  e <- list(); rows <- list()
  dg <- function(x) if (startsWith(x, "sens_dry")) 1 else 2
  st <- list()
  for (k in names(EXPDIST_MAP)) {
    mm <- EXPDIST_MAP[[k]]; f <- fits[[mm[1]]]
    s <- exp_stats(as.numeric(f$data[[mm[2]]]), f$data$weight)
    st[[k]] <- s
    rows[[k]] <- data.table::data.table(exposure = k, model = PREC_MODELS[[mm[1]]]$id, n = as.integer(s[["n"]]),
                                        t(s[c(names(EXPDIST_P), "min", "max", "sd", "iqr")]))
    desc <- sprintf("%s (%s) in the %s estimation sample (fixest::obs), weighted", k, mm[2], PREC_MODELS[[mm[1]]]$id)
    for (stat in EXPDIST_STATS)
      e[[sprintf("B.supp.expdist.%s.%s", k, stat)]] <- .ent(s[[stat]], fmt_num(s[[stat]], dg(k)), S,
        paste0(stat, " of ", desc, if (stat %in% c("min", "max")) " (unweighted)" else ""))
    e[[sprintf("B.supp.expdist.%s.n", k)]] <- ent_int(s[["n"]], S, paste("n of", desc))
  }
  # rescaled precision
  pr <- list()
  add <- function(key, model, term, scale, sv, r, digits, unit) {
    z <- rescale_ci(r$est, r$lo, r$hi, sv)
    e[[key]] <<- ent_ci(z$est, z$lo, z$hi, digits, S, sprintf("%s %s rescaled: coefficient x %s (s = %s); %s",
                                                            model, term, scale, signif(sv, 6), unit))
    pr[[key]] <<- data.table::data.table(model = model, term = term, scale = scale, s = sv,
                                         est = z$est, lo = z$lo, hi = z$hi, unit = unit)
  }
  t4 <- fits$h4bw$tab
  for (tm in c("T2_dry", "T3_dry")) {
    r <- t4[term == tm]
    add(sprintf("B.precision.h4bw.%s.per_iqr", tm), "H4-bw", tm, "IQR", st[[tm]][["iqr"]], r, 1, "grams per IQR of trimester dryness")
    add(sprintf("B.precision.h4bw.%s.per_sd", tm), "H4-bw", tm, "SD", st[[tm]][["sd"]], r, 1, "grams per SD of trimester dryness")
  }
  r <- fits$h5haz2459$tab
  add("B.precision.h5haz2459.dry1000.per_iqr", "H5-haz2459", "dry1000_2459", "IQR", st$dry1000_2459[["iqr"]], r, 3, "HAZ per IQR of 1,000-day dryness")
  add("B.precision.h5haz2459.dry1000.per_sd", "H5-haz2459", "dry1000_2459", "SD", st$dry1000_2459[["sd"]], r, 3, "HAZ per SD of 1,000-day dryness")
  for (o in c("bw", "haz2459")) {
    f <- fits[[if (o == "bw") "xsbw" else "xshaz"]]$tab; s <- st[[paste0("sens_dry_", o)]]
    u <- if (o == "bw") "grams" else "HAZ"; dgt <- if (o == "bw") 1 else 3; id <- paste0("ENSOxS-", o)
    add(sprintf("B.precision.ensoxs.%s.per_iqr", o), id, "elnino_x_sens_dry", "IQR of sens_dry", s[["iqr"]], f, dgt,
        paste(u, "change in the El Nino-pregnancy difference across the IQR of sens_dry"))
    add(sprintf("B.precision.ensoxs.%s.per_sd", o), id, "elnino_x_sens_dry", "SD of sens_dry", s[["sd"]], f, dgt,
        paste(u, "change in the El Nino-pregnancy difference per SD of sens_dry"))
    add(sprintf("B.precision.ensoxs.%s.at_p50", o), id, "elnino_x_sens_dry", "p50 of sens_dry", s[["p50"]], f, dgt,
        paste(u, "implied El Nino-pregnancy difference at the p50 of sens_dry"))
    add(sprintf("B.precision.ensoxs.%s.at_p90", o), id, "elnino_x_sens_dry", "p90 of sens_dry", s[["p90"]], f, dgt,
        paste(u, "implied El Nino-pregnancy difference at the p90 of sens_dry"))
  }
  list(entries = e, tables = list(supp_expdist = data.table::rbindlist(rows), supp_precision = data.table::rbindlist(pr)),
       check = data.table::rbindlist(lapply(fits, `[[`, "tab")))
}

# Every run: the refits must reproduce the published primary estimates (relative tol 1e-8).
check_precision <- function(steps, tol = 1e-8) {
  ref <- data.table::rbindlist(lapply(steps, function(s) data.table::rbindlist(list(s$tables$primary_models, s$tables$enso_total), fill = TRUE)), fill = TRUE)
  chk <- steps$precision$check
  for (i in seq_len(nrow(chk))) {
    r <- chk[i]
    p <- ref[model == r$model & term == r$term]
    if (nrow(p) != 1L) stop("check_precision: no unique published row for ", r$model, " / ", r$term, call. = FALSE)
    for (v in c("est", "lo", "hi", "n"))
      if (!isTRUE(all.equal(as.numeric(r[[v]]), as.numeric(p[[v]]), tolerance = tol)))
        stop("check_precision: refit ", r$model, " / ", r$term, " ", v, " differs from the published table", call. = FALSE)
  }
  invisible(TRUE)
}
