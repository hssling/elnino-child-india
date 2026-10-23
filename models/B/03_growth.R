# Paper B Task 7 step 3: H5, HAZ and stunting by age group (results_keys.md section 3).
# Exposure: 1,000-day dryness (month-weighted, ruling P13 Q3); FE region_id + birth_ym + age_fe.

run_growth <- function() {
  S <- "models/B/03_growth.R"
  d <- load_b()
  fe <- c("region_id", "birth_ym", "age_fe")
  e <- list(); tabs <- list()
  for (g in c("2459", "0023")) {
    x <- d[get(paste0("s_haz", g)) == TRUE]
    x[, dry := get(paste0("dry1000_", g))]
    th <- run_model(paste0("H5-haz", g), "haz", c("dry", COV), fe, x, "dry")
    ts <- run_model(paste0("H5-stunt", g), "stunted", c("dry", COV), fe, x, "dry", family = "logit")
    tabs[[length(tabs) + 1]] <- th; tabs[[length(tabs) + 1]] <- ts
    lab <- if (g == "2459") "24-59 months, 1,000-day (3,3,3,12,12)/33" else "0-23 months, (3,3,3,12)/21 where complete"
    r <- row_of(th, "dry")
    e[[sprintf("B.h5.haz%s.dry", g)]] <- ent_row(r, 3, S, paste("H5-haz: HAZ per +1 SPI unit drier,", lab))
    e[[sprintf("B.h5.haz%s.dry.p", g)]] <- ent_p(r$p, S, paste("H5-haz p,", lab))
    e[[sprintf("B.h5.haz%s.spi3", g)]] <- ent_flip(r, 3, S, paste("H5-haz: HAZ per +1 SPI-3 (wetter),", lab))
    e[[sprintf("B.h5.haz%s.n", g)]] <- ent_int(r$n, S, paste("H5-haz n,", lab))
    rs <- row_of(ts, "dry")
    e[[sprintf("B.h5.stunt%s.dry", g)]] <- ent_row(rs, 2, S, paste("H5-stunt logit: OR per +1 SPI unit drier,", lab))
    e[[sprintf("B.h5.stunt%s.dry.p", g)]] <- ent_p(rs$p, S, paste("H5-stunt p,", lab))
    e[[sprintf("B.h5.stunt%s.n", g)]] <- ent_int(rs$n, S, paste("H5-stunt n,", lab))
    if (g == "2459") ok <- rule_dry_negative(r$est, r$hi)
  }
  e[["B.h5.verdict"]] <- ent_text(verdict_text(ok), S,
                                  "H5 rule: dry1000 coefficient for HAZ 24-59 negative with 95% CI upper bound < 0")
  list(entries = e, tables = list(primary_models = data.table::rbindlist(tabs)), rule = list(h5haz = ok))
}
