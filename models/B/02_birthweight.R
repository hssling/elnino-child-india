# Paper B Task 7 step 2: H4, birth weight (results_keys.md section 2). Primary local-climate
# model: weighted, district + birth year-month FE, two-way clustered (region_id, birth_ym).

H4_RHS <- c("T1_dry", "T2_dry", "T3_dry", "T1_tmax_anom", "T2_tmax_anom", "T3_tmax_anom")
H4_TERMS <- H4_RHS

run_birthweight <- function() {
  S <- "models/B/02_birthweight.R"
  d <- load_b()
  bw <- d[s_bw == TRUE]
  fe <- c("region_id", "birth_ym")
  t_bw <- run_model("H4-bw", "bw_g", c(H4_RHS, COV), fe, bw, H4_TERMS)
  t_lbw <- run_model("H4-lbw", "lbw", c(H4_RHS, COV), fe, bw, H4_TERMS, family = "logit")

  e <- list()
  for (k in c("T1", "T2", "T3")) {
    r <- row_of(t_bw, paste0(k, "_dry"))
    e[[sprintf("B.h4.bw.%s_dry", k)]] <- ent_row(r, 1, S, sprintf("H4-bw: grams per +1 SPI unit drier (-%s_spi3), weighted, FE region_id+birth_ym, cluster region_id+birth_ym", k))
    e[[sprintf("B.h4.bw.%s_dry.p", k)]] <- ent_p(r$p, S, sprintf("H4-bw p, %s_dry", k))
    e[[sprintf("B.h4.bw.%s_spi3", k)]] <- ent_flip(r, 1, S, sprintf("H4-bw: grams per +1 SPI-3 (wetter), %s; same fit, sign flipped", k))
    rt <- row_of(t_bw, paste0(k, "_tmax_anom"))
    e[[sprintf("B.h4.bw.%s_tmax", k)]] <- ent_row(rt, 1, S, sprintf("H4-bw: grams per +1 deg C Tmax anomaly, %s", k))
    e[[sprintf("B.h4.bw.%s_tmax.p", k)]] <- ent_p(rt$p, S, sprintf("H4-bw p, %s_tmax_anom", k))
    r <- row_of(t_lbw, paste0(k, "_dry"))
    e[[sprintf("B.h4.lbw.%s_dry", k)]] <- ent_row(r, 2, S, sprintf("H4-lbw logit: OR of LBW per +1 SPI unit drier, %s", k))
    e[[sprintf("B.h4.lbw.%s_dry.p", k)]] <- ent_p(r$p, S, sprintf("H4-lbw p, %s_dry", k))
    rt <- row_of(t_lbw, paste0(k, "_tmax_anom"))
    e[[sprintf("B.h4.lbw.%s_tmax", k)]] <- ent_row(rt, 2, S, sprintf("H4-lbw logit: OR per +1 deg C, %s", k))
    e[[sprintf("B.h4.lbw.%s_tmax.p", k)]] <- ent_p(rt$p, S, sprintf("H4-lbw p, %s_tmax_anom", k))
  }
  e[["B.h4.bw.n"]] <- ent_int(t_bw$n[1], S, "H4-bw observations used")
  e[["B.h4.lbw.n"]] <- ent_int(t_lbw$n[1], S, "H4-lbw observations used")
  ok2 <- rule_dry_negative(row_of(t_bw, "T2_dry")$est, row_of(t_bw, "T2_dry")$hi)
  ok3 <- rule_dry_negative(row_of(t_bw, "T3_dry")$est, row_of(t_bw, "T3_dry")$hi)
  e[["B.h4.verdict"]] <- ent_text(verdict_text(ok2 || ok3), S,
                                  "H4 rule: T2_dry or T3_dry coefficient negative with 95% CI upper bound < 0 (H4-bw)")
  e[["B.h4.verdict_terms"]] <- ent_text(verdict_terms(ok2, ok3), S, "Terms meeting the H4 rule")
  list(entries = e, tables = list(primary_models = rbind(t_bw, t_lbw)),
       rule = list(h4bw = ok2 || ok3))
}
