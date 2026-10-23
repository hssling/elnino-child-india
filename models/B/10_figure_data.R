# Paper B Phase 3b Task 1: aggregate data behind the figures (papers/B/figures/data/;
# results_keys.md "Amendment 2026-10-23", "Figure data", and the ruling-Q2 note). Aggregates
# only: a small cell (n < 30) keeps its row with its value NA (suppressed = TRUE). No exact counts
# where differencing could recover a suppressed cell: fig2 carries an n band, fig1 weighted birth
# shares. No region_id / birth_ym column: districts are labelled district_gid, months YYYY-MM.

FIG_YEARS <- c(2006L, 2021L)
COHORT_2015 <- c(201506L, 201605L)          # the 2015 El Nino birth cohort, June 2015 - May 2016

# ---------------------------------------------------------------- Figure 1

fig1_timeline <- function(d, enso = data.table::fread(file.path(REPO_ROOT, "data", "derived", "enso_monthly.csv")),
                          years = FIG_YEARS, min_n = MIN_CELL_N) {
  em <- data.table::as.data.table(enso)[year >= years[1] & year <= years[2]]
  out <- em[, .(month = sprintf("%04d-%02d", year, month), oni_anom, roni, enso_phase = phase)]
  if (nrow(out) != 12L * (years[2] - years[1] + 1L) || anyDuplicated(out$month))
    stop("fig1_timeline: enso_monthly.csv does not hold one row per month for ", years[1], "-", years[2], call. = FALSE)
  b <- data.table::data.table(round = as.character(d$round),
                              month = sprintf("%04d-%02d", as.integer(d$birth_ym) %/% 100L, as.integer(d$birth_ym) %% 100L),
                              weight = d$weight)
  tot <- b[, .(w_total = sum(weight)), by = round]
  b <- b[, .(n = .N, w = sum(weight)), by = .(round, month)]
  # ruling Q2: no exact counts (a suppressed month would equal the round total minus the
  # published months); the month's share of the round's weighted births, %, 0.1 decimals
  for (r in sort(unique(b$round))) {
    br <- b[round == r]
    n <- br$n[match(out$month, br$month)]; n[is.na(n)] <- 0L
    w <- br$w[match(out$month, br$month)]
    small <- n < min_n
    share <- round(100 * w / tot[round == r, w_total], 1)
    data.table::set(out, j = paste0("birth_share_pct_", r), value = ifelse(small, NA_real_, share))
    data.table::set(out, j = paste0("suppressed_", r), value = small)
  }
  out[]
}

# ---------------------------------------------------------------- Figure 2

# Long table district x panel x round over the FULL grid (every district in every panel);
# value = weighted mean T2_dry or weighted stunting share; n < min_n (including 0) -> NA;
# n_band (band_n) instead of the exact n.
district_map_data <- function(d, min_n = MIN_CELL_N) {
  x <- data.table::data.table(district_gid = as.character(d$region_id), state = as.character(d$state_gadm),
                              round = as.character(d$round), bym = as.integer(d$birth_ym),
                              t2 = as.numeric(d$T2_dry), w = as.numeric(d$weight),
                              s_haz = d$s_haz2459 %in% TRUE, stunt = as.numeric(d$stunted))
  states <- x[, .(state = names(which.max(table(state)))), by = district_gid]
  x[, cohort := bym >= COHORT_2015[1] & bym <= COHORT_2015[2]]
  wm <- function(v, w) {
    ok <- is.finite(v) & is.finite(w)
    list(n_children = sum(ok), value = if (any(ok)) sum(w[ok] * v[ok]) / sum(w[ok]) else NA_real_)
  }
  vals <- data.table::rbindlist(list(
    x[cohort == TRUE, wm(t2, w), by = district_gid][, `:=`(panel = "T2_dry_cohort_2015_16", round = "pooled")],
    x[cohort == FALSE, wm(t2, w), by = district_gid][, `:=`(panel = "T2_dry_other_cohorts", round = "pooled")],
    x[, wm(t2, w), by = .(district_gid, round)][, panel := "T2_dry"],
    x[s_haz == TRUE, wm(stunt, w), by = .(district_gid, round)][, panel := "stunting_24_59"]), use.names = TRUE)
  rounds <- sort(unique(x$round))
  grid <- rbind(data.table::CJ(district_gid = states$district_gid, panel = c("T2_dry_cohort_2015_16", "T2_dry_other_cohorts"), round = "pooled"),
                data.table::CJ(district_gid = states$district_gid, panel = c("T2_dry", "stunting_24_59"), round = rounds))
  out <- vals[grid, on = c("district_gid", "panel", "round")]
  out[is.na(n_children), n_children := 0L]
  out <- states[out, on = "district_gid"]
  # ruling Q2: an n band, not the exact n (cohort + other = NFHS-4 + NFHS-5 per district, and the
  # by-round panel sums to the state x round counts, so exact n would expose suppressed cells)
  out[, n_band := band_n(n_children, min_n)]
  out <- suppress_small(out, "n_children", "value", min_n)
  out[, n_children := NULL]
  data.table::setcolorder(out, c("district_gid", "state", "panel", "round", "n_band", "value", "suppressed"))
  data.table::setorder(out, panel, round, district_gid)
  out[]
}

# ---------------------------------------------------------------- Figures 3, 4, S1

FIG3_PANELS <- data.table::data.table(
  model = c("H4-bw", "H4-lbw", "H5-haz2459", "H5-haz0023", "H5-stunt2459", "H5-stunt0023", "H6-diarr", "H6-anaem",
            "ENSOxS-bw", "ENSOxS-haz2459", "ENSOxS-haz0023"),
  panel = c("bw", "lbw", "haz", "haz", "stunting", "stunting", "h6", "h6", "ensoxs", "ensoxs", "ensoxs"),
  age_group = c("birth", "birth", "24-59", "0-23", "24-59", "0-23", "recent", "recent", "birth", "24-59", "0-23"))

fig3_coefs <- function(primary) {
  p <- data.table::as.data.table(primary)
  miss <- setdiff(unique(p$model), FIG3_PANELS$model)
  if (length(miss)) stop("fig3_coefs: model(s) without a panel: ", paste(miss, collapse = ", "), call. = FALSE)
  out <- FIG3_PANELS[p, on = "model"]
  out[, .(panel, age_group, model, outcome, family, term, est, lo, hi, p, n, unit, vcov_fixed)]
}

# primary term of each rule-check model code
FIG4_PRIMARY <- data.table::data.table(model_code = c("h4bw", "h4bw", "h5haz", "xsbw", "xshaz"),
                                       term = c("T2_dry", "T3_dry", "dry", "inter", "inter"),
                                       primary_model = c("H4-bw", "H4-bw", "H5-haz2459", "ENSOxS-bw", "ENSOxS-haz2459"))

fig4_sensitivity <- function(sens, verdicts, primary) {
  st <- data.table::as.data.table(sens)[model_code %in% FIG4_PRIMARY$model_code & spec != "placebo"]
  st <- FIG4_PRIMARY[st, on = c("model_code", "term"), nomatch = NULL]       # rule terms only
  pr <- data.table::as.data.table(primary)[, .(primary_model = model, term, primary_est = est, primary_lo = lo, primary_hi = hi)]
  st <- pr[st, on = c("primary_model", "term")]
  v <- data.table::as.data.table(verdicts)[, .(spec, model_code, rule_met, same_as_primary = same)]
  st <- v[st, on = c("spec", "model_code")]
  if (data.table::uniqueN(st[, .(spec, model_code)]) != nrow(v)) stop("fig4_sensitivity: checks do not match sensitivity_verdicts", call. = FALSE)
  if (anyNA(st$primary_est)) stop("fig4_sensitivity: a check has no primary row", call. = FALSE)
  out <- st[, .(spec, model_code, model, term, est, lo, hi, p, n, unit, rule_met, same_as_primary,
                primary_est, primary_lo, primary_hi)]
  data.table::setorder(out, model_code, spec, term)
  out[]
}

figS1_placebo <- function(sens, primary) {
  pr <- data.table::as.data.table(primary)[model == "H4-bw"]
  pl <- data.table::as.data.table(sens)[spec == "placebo"]
  if (nrow(pr) != 6L || nrow(pl) != 4L) stop("figS1_placebo: expected 6 primary and 4 placebo rows", call. = FALSE)
  rbind(pr[, .(source = "primary", term, est, lo, hi, p, n, unit)],
        pl[, .(source = "placebo", term, est, lo, hi, p, n, unit)])
}

# ---------------------------------------------------------------- step

stack_tables <- function(steps, name) data.table::rbindlist(lapply(steps, function(s) s$tables[[name]]), fill = TRUE)

run_figure_data <- function(steps) {
  d <- load_b()
  primary <- stack_tables(steps, "primary_models")
  sens <- stack_tables(steps, "sensitivity")
  ver <- stack_tables(steps, "sensitivity_verdicts")
  list(entries = list(),
       fig_tables = list(fig1_timeline = fig1_timeline(d), fig2_district_map = district_map_data(d),
                         fig3_coefs = fig3_coefs(primary), fig4_sensitivity = fig4_sensitivity(sens, ver, primary),
                         figS1_placebo = figS1_placebo(sens, primary)))
}
