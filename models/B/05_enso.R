# Paper B Task 7 step 5: ENSO x district sensitivity (verdict) and total ENSO associations
# (no verdict); results_keys.md section 5.

ENSO_OUTCOMES <- list(bw = list(lhs = "bw_g", sample = "s_bw", age = FALSE),
                      haz2459 = list(lhs = "haz", sample = "s_haz2459", age = TRUE),
                      haz0023 = list(lhs = "haz", sample = "s_haz0023", age = TRUE))

# ENSO x sensitivity fit for one outcome; `inter` / `fe_region` / `cluster` / `weights` and a
# row filter let 07_sensitivity reuse it.
fit_ensoxs <- function(d, out, id, inter = "elnino_x_sens_dry", fe_region = "region_id",
                       cluster = B_CLUSTER, weights = ~weight) {
  o <- ENSO_OUTCOMES[[out]]
  x <- d[get(o$sample) == TRUE]
  fe <- c(fe_region, "birth_ym", if (o$age) "age_fe")
  t <- run_model(id, o$lhs, c(inter, COV), fe, x, inter, weights = weights, cluster = cluster)
  t[term == inter, term := "inter"]
  t[]
}

fit_enso_total <- function(d, out, ix, id, extra = character(), extra_terms = character()) {
  o <- ENSO_OUTCOMES[[out]]
  x <- d[get(o$sample) == TRUE]
  fe <- c("region_id", "birth_month_of_year", "round", if (o$age) "age_fe")
  terms <- c(paste0("preg_", ix, c("_pos", "_neg")), extra_terms)
  run_model(id, o$lhs, c(terms, extra, COV, CMC_NS), fe, x, terms)
}

run_enso <- function() {
  S <- "models/B/05_enso.R"
  d <- load_b()
  e <- list(); tabs <- list(); rules <- list()
  digits <- c(bw = 2, haz2459 = 4, haz0023 = 4)
  for (out in names(ENSO_OUTCOMES)) {
    t <- fit_ensoxs(d, out, paste0("ENSOxS-", out))
    tabs[[length(tabs) + 1]] <- t
    r <- row_of(t, "inter")
    unit <- if (out == "bw") "grams" else "HAZ"
    e[[paste0("B.ensoxs.", out)]] <- ent_row(r, digits[[out]], S, sprintf(
      "ENSOxS %s: %s per unit of preg_elnino_dev x sens_dry (sens_dry = -sens_beta, %% JJAS rain loss per deg C ONI); FE region_id+birth_ym%s; cluster region_id+birth_ym",
      out, unit, if (out == "bw") "" else "+age_fe"))
    e[[paste0("B.ensoxs.", out, ".p")]] <- ent_p(r$p, S, paste("ENSOxS p,", out))
    e[[paste0("B.ensoxs.", out, ".n")]] <- ent_int(r$n, S, paste("ENSOxS n,", out))
    rules[[out]] <- rule_dry_negative(r$est, r$hi)
  }
  e[["B.ensoxs.verdict"]] <- ent_text(verdict_text(rules$bw || rules$haz2459), S,
    "ENSO x sensitivity rule: interaction negative with 95% CI upper bound < 0 for birth weight or HAZ 24-59")

  idx_lab <- c(oni = "ONI (oni_anom)", roni = "RONI", nino34 = "Nino-3.4 anomaly (monthly)")
  for (ix in names(idx_lab)) for (out in names(ENSO_OUTCOMES)) {
    t <- fit_enso_total(d, out, ix, sprintf("ENSO-total-%s-%s", ix, out))
    tabs[[length(tabs) + 1]] <- t
    dg <- if (out == "bw") 1 else 3
    base <- sprintf("B.enso.%s.%s", ix, out)
    desc <- sprintf("ENSO-total %s, %s: FE region_id+birth_month_of_year+round%s, ns(birth_cmc,4), weighted, cluster region_id+birth_ym",
                    idx_lab[[ix]], out, if (out == "bw") "" else "+age_fe")
    rp <- row_of(t, paste0("preg_", ix, "_pos"))
    e[[paste0(base, ".elnino")]] <- ent_row(rp, dg, S, paste(desc, "; per +1 deg C of the El Nino term max(x,0)"))
    e[[paste0(base, ".elnino.p")]] <- ent_p(rp$p, S, paste(desc, "; El Nino term p"))
    rn <- row_of(t, paste0("preg_", ix, "_neg"))
    e[[paste0(base, ".lanina")]] <- ent_flip(rn, dg, S, paste(desc, "; per 1 deg C cooler of the La Nina term min(x,0) (-beta, bounds swapped)"))
    e[[paste0(base, ".lanina.p")]] <- ent_p(rn$p, S, paste(desc, "; La Nina term p"))
    e[[paste0(base, ".n")]] <- ent_int(rp$n, S, paste(desc, "; n"))
  }
  tt <- data.table::rbindlist(tabs)
  list(entries = e, tables = list(primary_models = tt[startsWith(model, "ENSOxS")],
                                  enso_total = tt[startsWith(model, "ENSO-total")]),
       rule = list(xsbw = rules$bw, xshaz = rules$haz2459))
}
