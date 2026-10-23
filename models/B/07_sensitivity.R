# Paper B Task 7 step 7: pre-specified sensitivity analyses (results_keys.md "Sensitivity
# analyses" and section 7). Each analysis changes only what its contract row states.

H4_DRY <- c("T1_dry", "T2_dry", "T3_dry")
H4_TMAX <- c("T1_tmax_anom", "T2_tmax_anom", "T3_tmax_anom")

# Rows complete on the model variables (for the mother-FE restriction).
complete_on <- function(x, vars) x[stats::complete.cases(x[, unique(vars), with = FALSE])]

# H4-bw variant -> rows T2_dry / T3_dry (canonical names) plus any extra terms.
sens_h4 <- function(x, id, dry = H4_DRY, tmax = H4_TMAX, fe = c("region_id", "birth_ym"),
                    weights = ~weight, cluster = B_CLUSTER, extra = character(), mother = FALSE,
                    lhs = "bw_g", family = "gaussian") {
  rhs <- c(dry, tmax, extra, COV)
  if (mother) {
    x <- complete_on(x, c(lhs, rhs, fe, "weight"))
    x <- x[multi_birth_mothers(mother_id)]
  }
  t <- run_model(id, lhs, rhs, fe, x, c(dry[2:3], extra), family = family, weights = weights, cluster = cluster)
  t[term == dry[2], term := "T2_dry"]
  t[term == dry[3], term := "T3_dry"]
  t[]
}

sens_h5 <- function(x, id, dry = "dry1000_2459", fe = c("region_id", "birth_ym", "age_fe"),
                    weights = ~weight, cluster = B_CLUSTER, mother = FALSE, lhs = "haz", family = "gaussian") {
  rhs <- c(dry, COV)
  if (mother) {
    x <- complete_on(x, c(lhs, rhs, fe, "weight"))
    x <- x[multi_birth_mothers(mother_id)]
  }
  t <- run_model(id, lhs, rhs, fe, x, dry, family = family, weights = weights, cluster = cluster)
  t[term == dry, term := "dry"]
  t[]
}

run_sensitivity <- function(primary_rules) {
  S <- "models/B/07_sensitivity.R"
  d <- load_b()
  bw <- d[s_bw == TRUE]
  hz <- d[s_haz2459 == TRUE]
  tabs <- list()
  add <- function(spec, model, t) {
    .s <- spec; .m <- model                      # t has a `model` column: avoid data.table scoping
    t[, `:=`(spec = .s, model_code = .m)]
    tabs[[length(tabs) + 1L]] <<- t
  }
  # row filters as functions of the frame (nonmover / precovid / nogaz contract rows)
  rowfilters <- list(nonmover = function(x) x$nonmover, precovid = function(x) x$v008 < 1450L,
                     nogaz = function(x) x$source != "GAZ")
  xw_cl <- ~region_id_xw + birth_ym

  # ---- h4bw
  add("card", "h4bw", sens_h4(bw[bw_card == TRUE], "S-card-h4bw"))
  add("lastbirth", "h4bw", sens_h4(bw[bidx == 1L], "S-lastbirth-h4bw"))
  add("unweighted", "h4bw", sens_h4(bw, "S-unweighted-h4bw", weights = NULL))
  for (s in names(rowfilters)) add(s, "h4bw", sens_h4(bw[rowfilters[[s]](bw)], paste0("S-", s, "-h4bw")))
  add("motherfe", "h4bw", sens_h4(bw, "S-motherfe-h4bw", fe = c("mother_id", "birth_ym"), mother = TRUE))
  for (s in c("sm1", "sp1"))
    add(c(sm1 = "shiftm1", sp1 = "shiftp1")[[s]], "h4bw",
        sens_h4(bw, paste0("S-shift", s, "-h4bw"), dry = paste0(H4_DRY, "_", s), tmax = paste0(H4_TMAX, "_", s)))
  add("district", "h4bw", sens_h4(bw, "S-district-h4bw", dry = paste0("dist_", H4_DRY), tmax = paste0("dist_", H4_TMAX),
                                  fe = c("region_id_xw", "birth_ym"), cluster = xw_cl))
  add("buf10", "h4bw", sens_h4(bw, "S-buf10-h4bw", dry = paste0(H4_DRY, "_b10r")))
  add("noflag", "h4bw", sens_h4(bw[!any_flag(bw, c("T1", "T2", "T3"))], "S-noflag-h4bw"))
  add("placebo", "h4bw", sens_h4(bw, "S-placebo-h4bw", extra = c("placebo_dry", "placebo_tmax_anom")))

  # ---- h5haz (24-59)
  add("unweighted", "h5haz", sens_h5(hz, "S-unweighted-h5haz", weights = NULL))
  for (s in names(rowfilters)) add(s, "h5haz", sens_h5(hz[rowfilters[[s]](hz)], paste0("S-", s, "-h5haz")))
  add("motherfe", "h5haz", sens_h5(hz, "S-motherfe-h5haz", fe = c("mother_id", "birth_ym", "age_fe"), mother = TRUE))
  add("district", "h5haz", sens_h5(hz, "S-district-h5haz", dry = "dist_dry1000_2459",
                                   fe = c("region_id_xw", "birth_ym", "age_fe"), cluster = xw_cl))
  add("buf10", "h5haz", sens_h5(hz, "S-buf10-h5haz", dry = "b10r_dry1000_2459"))
  add("noflag", "h5haz", sens_h5(hz[!any_flag(hz, WINDOWS_1000)], "S-noflag-h5haz"))

  # ---- ENSO x sensitivity
  for (out in c("bw", "haz2459")) {
    mc <- if (out == "bw") "xsbw" else "xshaz"
    add("unweighted", mc, fit_ensoxs(d, out, paste0("S-unweighted-", mc), weights = NULL))
    for (s in names(rowfilters)) add(s, mc, fit_ensoxs(d[rowfilters[[s]](d)], out, paste0("S-", s, "-", mc)))
    add("district", mc, fit_ensoxs(d, out, paste0("S-district-", mc), inter = "elnino_x_sens_dry_xw",
                                   fe_region = "region_id_xw", cluster = xw_cl))
    add("cluster1", mc, fit_ensoxs(d, out, paste0("S-cluster1-", mc), cluster = ~region_id))
  }

  # ---- non-verdict analyses
  add("ensoplacebo", "bw", fit_enso_total(d, "bw", "oni", "S-ensoplacebo-bw", extra_terms = c("post_oni_pos", "post_oni_neg")))
  add("hbplaus", "anaem", run_model("S-hbplaus-anaem", "anaemia_any", c("recent12_dry", "recent12_tmax_anom", COV),
                                    c("region_id", "interview_ym", "age_fe"), d[s_anaem == TRUE & hb_plausible == 1],
                                    c("recent12_dry", "recent12_tmax_anom"), family = "logit"))
  add("lpm", "lbw", sens_h4(bw, "S-lpm-lbw", lhs = "lbw"))
  add("lpm", "stunt2459", sens_h5(hz, "S-lpm-stunt2459", lhs = "stunted"))

  st <- data.table::rbindlist(tabs, fill = TRUE)
  data.table::setcolorder(st, c("spec", "model_code", "model", "outcome", "term"))

  # ---- rule checks against the primary
  checks <- st[model_code %in% c("h4bw", "h5haz", "xsbw", "xshaz") & spec != "placebo"]
  rule_terms <- list(h4bw = c("T2_dry", "T3_dry"), h5haz = "dry", xsbw = "inter", xshaz = "inter")
  ver <- checks[, .(rule_met = any(rule_dry_negative(est[term %in% rule_terms[[model_code[1]]]],
                                                     hi[term %in% rule_terms[[model_code[1]]]])),
                    n = n[1]), by = .(spec, model_code)]
  ver[, primary_rule_met := vapply(model_code, function(m) isTRUE(primary_rules[[m]]), logical(1))]
  ver[, same := rule_met == primary_rule_met]

  # ---- entries
  e <- list()
  dig <- c(h4bw = 1, h5haz = 3, xsbw = 2, xshaz = 4)
  for (i in seq_len(nrow(checks))) {
    r <- checks[i]
    e[[sprintf("B.sens.%s.%s.%s", r$spec, r$model_code, r$term)]] <-
      ent_row(r, dig[[r$model_code]], S, sprintf("Sensitivity %s, model %s (%s), term %s; fe %s; cluster %s; weighted %s",
                                                 r$spec, r$model_code, r$model, r$term, r$fe, r$cluster, r$weighted))
  }
  for (i in seq_len(nrow(ver)))
    e[[sprintf("B.sens.%s.%s.n", ver$spec[i], ver$model_code[i])]] <- ent_int(ver$n[i], S, sprintf("n, sensitivity %s / %s", ver$spec[i], ver$model_code[i]))
  pl <- st[spec == "placebo"]
  e[["B.sens.placebo.h4bw.placebo_dry"]] <- ent_row(row_of(pl, "placebo_dry"), 1, S, "Placebo: grams per +1 SPI unit drier in the pseudo-pregnancy window (birth +1..+9), H4-bw with T1-T3 terms")
  e[["B.sens.placebo.h4bw.placebo_tmax"]] <- ent_row(row_of(pl, "placebo_tmax_anom"), 1, S, "Placebo: grams per +1 deg C, pseudo-pregnancy window")
  e[["B.sens.placebo.h4bw.n"]] <- ent_int(pl$n[1], S, "Placebo H4-bw n")
  ep <- st[spec == "ensoplacebo"]
  e[["B.sens.ensoplacebo.bw.post_elnino"]] <- ent_row(row_of(ep, "post_oni_pos"), 1, S, "ENSO placebo: grams per +1 deg C of max(post-birth ONI +1..+12, 0) in ENSO-total bw (ONI)")
  e[["B.sens.ensoplacebo.bw.post_lanina"]] <- ent_flip(row_of(ep, "post_oni_neg"), 1, S, "ENSO placebo: grams per 1 deg C cooler of min(post-birth ONI, 0)")
  e[["B.sens.ensoplacebo.bw.n"]] <- ent_int(ep$n[1], S, "ENSO placebo n")
  hb <- st[spec == "hbplaus"]
  e[["B.sens.hbplaus.anaem.recent12_dry"]] <- ent_row(row_of(hb, "recent12_dry"), 2, S, "H6-anaem restricted to hb_plausible: OR per +1 SPI unit drier")
  e[["B.sens.hbplaus.anaem.recent12_tmax"]] <- ent_row(row_of(hb, "recent12_tmax_anom"), 2, S, "H6-anaem restricted to hb_plausible: OR per +1 deg C")
  e[["B.sens.hbplaus.anaem.n"]] <- ent_int(hb$n[1], S, "H6-anaem hb_plausible n")
  lb <- st[spec == "lpm" & model_code == "lbw"]
  for (k in c("T2_dry", "T3_dry"))
    e[[paste0("B.sens.lpm.lbw.", k)]] <- ent_row(row_of(lb, k), 2, S, paste("LBW linear probability model, percentage points per +1 SPI unit drier,", k), scale = 100)
  e[["B.sens.lpm.lbw.n"]] <- ent_int(lb$n[1], S, "LBW LPM n")
  ls <- st[spec == "lpm" & model_code == "stunt2459"]
  e[["B.sens.lpm.stunt2459.dry"]] <- ent_row(row_of(ls, "dry"), 2, S, "Stunting 24-59 LPM, percentage points per +1 SPI unit drier (1,000-day)", scale = 100)
  e[["B.sens.lpm.stunt2459.n"]] <- ent_int(ls$n[1], S, "Stunting LPM n")

  e[["B.sens.n_checks"]] <- ent_int(nrow(ver), S, "Sensitivity (spec, model) verdict-rule checks")
  e[["B.sens.n_same"]] <- ent_int(sum(ver$same), S, "Checks whose rule outcome equals the primary's")
  diff <- ver[same == FALSE, paste0(spec, "/", model_code)]
  e[["B.sens.differ"]] <- ent_text(if (length(diff)) paste(diff, collapse = ", ") else "none", S,
                                   "Sensitivity checks whose rule outcome differs from the primary")
  flag <- function(r) if (is.finite(r$lo) && is.finite(r$hi) && (r$lo > 0 || r$hi < 0)) "CI excludes 0" else "CI includes 0"
  e[["B.sens.placebo_flag"]] <- ent_text(flag(row_of(pl, "placebo_dry")), S, "Placebo pseudo-pregnancy dry term: does its 95% CI exclude 0?")
  e[["B.sens.ensoplacebo_flag"]] <- ent_text(flag(row_of(ep, "post_oni_pos")), S, "Post-birth El Nino term: does its 95% CI exclude 0?")

  list(entries = e, tables = list(sensitivity = st, sensitivity_verdicts = ver))
}
