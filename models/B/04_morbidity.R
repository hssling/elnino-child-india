# Paper B Task 7 step 4: H6, recent morbidity, descriptive only (results_keys.md section 4).
# FE region_id + interview_ym + age_fe; logit; weighted; cluster region_id + birth_ym.

run_morbidity <- function() {
  S <- "models/B/04_morbidity.R"
  d <- load_b()
  fe <- c("region_id", "interview_ym", "age_fe")
  td <- run_model("H6-diarr", "diarrhoea", c("recent1_dry", "recent1_tmax_anom", COV), fe,
                  d[s_diarr == TRUE], c("recent1_dry", "recent1_tmax_anom"), family = "logit")
  ta <- run_model("H6-anaem", "anaemia_any", c("recent12_dry", "recent12_tmax_anom", COV), fe,
                  d[s_anaem == TRUE], c("recent12_dry", "recent12_tmax_anom"), family = "logit")
  e <- list()
  for (spec in list(list("diarr", td, "recent1", "month before interview"),
                    list("anaem", ta, "recent12", "12 months before interview"))) {
    k <- spec[[1]]; t <- spec[[2]]; w <- spec[[3]]
    r <- row_of(t, paste0(w, "_dry"))
    e[[sprintf("B.h6.%s.%s_dry", k, w)]] <- ent_row(r, 2, S, sprintf("H6-%s logit (descriptive): OR per +1 SPI unit drier, %s", k, spec[[4]]))
    e[[sprintf("B.h6.%s.%s_dry.p", k, w)]] <- ent_p(r$p, S, sprintf("H6-%s p, %s_dry", k, w))
    r <- row_of(t, paste0(w, "_tmax_anom"))
    e[[sprintf("B.h6.%s.%s_tmax", k, w)]] <- ent_row(r, 2, S, sprintf("H6-%s logit (descriptive): OR per +1 deg C Tmax anomaly, %s", k, spec[[4]]))
    e[[sprintf("B.h6.%s.%s_tmax.p", k, w)]] <- ent_p(r$p, S, sprintf("H6-%s p, %s_tmax_anom", k, w))
    e[[sprintf("B.h6.%s.n", k)]] <- ent_int(r$n, S, sprintf("H6-%s n", k))
  }
  list(entries = e, tables = list(primary_models = rbind(td, ta)))
}
