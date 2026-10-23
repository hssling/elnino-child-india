# Paper B Task 7 step 1: samples, descriptive summaries and exposure missingness
# (papers/B/results_keys.md section 1). Aggregates only; cells with n < 30 suppressed.

run_descriptive <- function() {
  S <- "models/B/01_descriptive.R"
  d <- load_b()
  e <- list()
  e[["B.n.children"]] <- ent_int(nrow(d), S, "Rows (children) in analysis_child.parquet")
  for (r in c("nfhs4", "nfhs5"))
    e[[paste0("B.n.children.", r)]] <- ent_int(sum(d$round == r), S, paste("Children,", r))
  e[["B.n.districts"]] <- ent_int(data.table::uniqueN(d$region_id), S, "Distinct point-in-polygon GADM districts (region_id)")

  bw <- d[s_bw == TRUE]
  e[["B.desc.n_bw"]] <- ent_int(nrow(bw), S, "BW sample: singleton births with bw_g (before complete-case)")
  e[["B.desc.n_twins_excluded"]] <- ent_int(sum(d$twin & !is.na(d$bw_g)), S, "Twins with a birth weight, excluded")
  e[["B.desc.n_haz2459"]] <- ent_int(sum(d$s_haz2459), S, "HAZ sample 24-59 months")
  e[["B.desc.n_haz0023"]] <- ent_int(sum(d$s_haz0023), S, "HAZ sample 0-23 months")
  e[["B.desc.n_anaemia"]] <- ent_int(sum(d$s_anaem), S, "Anaemia sample (6-59 months, tested)")
  e[["B.desc.n_diarrhoea"]] <- ent_int(sum(d$s_diarr), S, "Diarrhoea sample (living children)")
  # ruling Q9 (M-16): singleton births without a recorded birth weight; if any part is 1-29, the
  # total and both rounds are suppressed together (no differencing)
  nsb <- d[twin %in% FALSE & is.na(bw_g), .N, by = round]
  nsb <- c(nfhs4 = sum(nsb[round == "nfhs4", N]), nfhs5 = sum(nsb[round == "nfhs5", N]))
  nsb <- c(all = sum(nsb), nsb)
  hide <- any(nsb > 0L & nsb < MIN_CELL_N)
  for (r in names(nsb)) {
    k <- if (r == "all") "B.desc.n_singleton_no_bw" else paste0("B.desc.n_singleton_no_bw.", r)
    desc <- paste("Singleton births without a recorded birth weight (bw_g missing),", if (r == "all") "both rounds" else r)
    e[[k]] <- if (hide) .ent(NA, "suppressed (n < 30)", S, desc) else ent_int(nsb[[r]], S, desc)
  }

  e[["B.desc.bw_mean"]] <- ent_mean_sd(bw$bw_g, bw$weight, 0, S, "Weighted mean birth weight (g), BW sample")
  e[["B.desc.lbw_pct"]] <- ent_pct(wpct(bw$lbw, bw$weight), S, "Weighted share LBW (<2500 g), BW sample")
  e[["B.desc.bw_card_pct"]] <- ent_pct(mean(bw$bw_card), S, "Unweighted share of BW sample weights from a card")
  e[["B.desc.bw_heap500_pct"]] <- ent_pct(mean(bw$bw_g %% 500 == 0), S, "Unweighted share of BW-sample weights at multiples of 500 g")
  h <- d[s_haz2459 == TRUE]
  e[["B.desc.haz2459_mean"]] <- ent_mean_sd(h$haz, h$weight, 2, S, "Weighted mean HAZ, 24-59 months")
  e[["B.desc.stunt2459_pct"]] <- ent_pct(wpct(h$stunted, h$weight), S, "Weighted share stunted, 24-59 months")
  a <- d[s_anaem == TRUE]
  e[["B.desc.anaemia_pct"]] <- ent_pct(wpct(a$anaemia_any, a$weight), S, "Weighted share any anaemia (6-59 months)")
  dd <- d[s_diarr == TRUE]
  e[["B.desc.diarrhoea_pct"]] <- ent_pct(wpct(dd$diarrhoea, dd$weight), S, "Weighted share diarrhoea in last 2 weeks")

  bwe <- bw[!is.na(preg_elnino_dev)]
  e[["B.desc.elnino_dev_pct"]] <- ent_pct(wpct(bwe$preg_elnino_dev, bwe$weight), S,
                                          "Weighted share of BW sample with preg_elnino_dev TRUE")
  he <- h[!is.na(preg_elnino_dev)]
  for (g in c("elnino", "other")) {
    v <- g == "elnino"
    sb <- bwe[preg_elnino_dev == v]; sh <- he[preg_elnino_dev == v]
    e[[paste0("B.desc.bw_mean.", g)]] <- ent_mean_sd(sb$bw_g, sb$weight, 0, S, paste("Weighted mean birth weight, preg_elnino_dev", v))
    e[[paste0("B.desc.lbw_pct.", g)]] <- ent_pct(wpct(sb$lbw, sb$weight), S, paste("Weighted share LBW, preg_elnino_dev", v))
    e[[paste0("B.desc.haz2459_mean.", g)]] <- ent_mean_sd(sh$haz, sh$weight, 2, S, paste("Weighted mean HAZ 24-59, preg_elnino_dev", v))
  }

  sens <- data.table::fread(file.path(REPO_ROOT, "data", "derived", "nfhs", "district_enso_sensitivity.csv"))
  sd_ <- -sens$sens_beta[is.finite(sens$sens_beta)]
  e[["B.desc.n_sens_districts"]] <- ent_int(length(sd_), S, "GADM districts with a 1981-2009 ENSO-rainfall sensitivity")
  e[["B.desc.sens_dry_median"]] <- .ent(stats::median(sd_), fmt_num(stats::median(sd_), 1), S,
                                        "Median sens_dry = -sens_beta (% JJAS rainfall loss per +1 deg C ONI) over districts")
  e[["B.desc.sens_dry_pct_pos"]] <- ent_pct(mean(sd_ > 0), S, "Share of districts with sens_dry > 0 (El Nino dries them)")

  for (r in c("nfhs4", "nfhs5"))
    e[[paste0("B.miss.T2_spi3_pct.", r)]] <- ent_pct(mean(is.na(d$T2_spi3[d$round == r])), S,
                                                     paste("Share of all children with missing T2_spi3,", r))
  tn4 <- d$round == "nfhs4" & d$state_gadm == "Tamil Nadu"
  e[["B.miss.tn4_T2_spi3_pct"]] <- ent_pct(mean(is.na(d$T2_spi3[tn4])), S,
                                           sprintf("Share missing T2_spi3, Tamil Nadu NFHS-4 (2011-12 mask); n=%d", sum(tn4)))

  # ---- tables
  # count columns are named n_* so check_aggregate_table enforces n >= 30 on each (P14, M9)
  flow <- d[, .(n_children = .N, n_bw_sample = sum(s_bw), n_twins_with_bw_excluded = sum(twin & !is.na(bw_g)),
                n_haz_24_59 = sum(s_haz2459), n_haz_0_23 = sum(s_haz0023), n_hb_tested = sum(s_anaem),
                n_diarrhoea_sample = sum(s_diarr)), by = round]
  summ <- function(x) x[, .(
    n_bw_sample = sum(s_bw),
    bw_mean_grams = sum((weight * bw_g)[s_bw]) / sum(weight[s_bw]),
    lbw_share = sum((weight * lbw)[s_bw]) / sum(weight[s_bw]),
    n_haz_24_59 = sum(s_haz2459),
    haz_mean_24_59 = sum((weight * haz)[s_haz2459]) / sum(weight[s_haz2459]),
    stunting_share_24_59 = sum((weight * stunted)[s_haz2459]) / sum(weight[s_haz2459]),
    n_hb_tested_sample = sum(s_anaem),
    share_any_anaemia = sum((weight * anaemia_any)[s_anaem]) / sum(weight[s_anaem]),
    n_diarrhoea_sample = sum(s_diarr),
    diarrhoea_share = sum((weight * diarrhoea)[s_diarr]) / sum(weight[s_diarr]))]
  by_round <- d[, summ(.SD), by = round]
  by_enso <- d[!is.na(preg_elnino_dev), summ(.SD), by = .(preg_elnino_dev)]
  # every summary must rest on >= 30 children
  for (tb in list(by_round, by_enso))
    if (any(unlist(tb[, grep("^n_", names(tb)), with = FALSE]) < MIN_CELL_N))
      stop("run_descriptive: a descriptive cell has n < 30", call. = FALSE)
  miss_cols <- c("T1_spi3", "T2_spi3", "T3_spi3", "T2_tmax_anom")
  missf <- function(x) {
    out <- x[, c(list(n_children = .N), lapply(.SD, function(v) mean(is.na(v)))), .SDcols = miss_cols]
    data.table::setnames(out, miss_cols, paste0("missing_share_", miss_cols))
    out
  }
  ms <- d[, missf(.SD), by = .(state = state_gadm, round)]
  data.table::setorder(ms, state, round)
  ms <- suppress_small(ms, "n_children", paste0("missing_share_", miss_cols))
  d[, birth_year := birth_ym %/% 100L]
  tn <- rbind(d[state_gadm == "Tamil Nadu", missf(.SD), by = .(round, birth_year)][, area := "Tamil Nadu"],
              d[, missf(.SD), by = .(round, birth_year)][, area := "India"])
  data.table::setcolorder(tn, c("area", "round", "birth_year"))
  data.table::setorder(tn, area, round, birth_year)
  # n band instead of exact n (ruling Q2): a suppressed TN birth-year cell must not be recoverable
  # from the state x round counts of missingness_state_round
  tn[, n_band := band_n(n_children)]
  tn <- suppress_small(tn, "n_children", paste0("missing_share_", miss_cols))
  tn[, n_children := NULL]
  # m-1 (Q2 re-review): unrounded fractions next to an n band could be inverted into counts
  for (c in paste0("missing_share_", miss_cols)) data.table::set(tn, j = c, value = round(tn[[c]], 3))
  data.table::setcolorder(tn, c("area", "round", "birth_year", "n_band"))
  d[, birth_year := NULL]

  list(entries = e, tables = list(sample_flow = flow, descriptive_by_round = by_round,
                                  descriptive_by_enso = by_enso, missingness_state_round = ms,
                                  missingness_tn_birthyear = tn))
}
