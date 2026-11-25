# Paper B additional (post hoc) analyses A1-A3, as fixed in papers/B/envres/analysis_plan_addendum.md
# (committed e60b3e9 before any of these models was fitted). They do not change any pre-specified
# verdict. Writes papers/B/envres/results_additional.json (keys B.x.*) and the aggregate table
# papers/B/envres/additional_models.csv. Aggregates only.
#
# Run from <repo>/models:  Rscript B/13_additional.R
# Needs data/external/nfhs_subset/maternal_height.parquet (python -I models/B/py/maternal_height.py).

source("common/R/load.R")
source("B/R/fe.R")
source("B/R/prep.R")
source("B/R/io.R")
source("B/R/display.R")
# step files as in run_all.R (definitions only); 12_precision.R gives wquantile, exp_stats, rescale_ci
for (f in c("01_descriptive.R", "02_birthweight.R", "03_growth.R", "04_morbidity.R", "05_enso.R",
            "06_mediation.R", "07_sensitivity.R", "08_results.R", "09_supplementary.R", "10_figure_data.R",
            "11_covariates.R", "12_precision.R"))
  source(file.path("B", f))

S <- "models/B/13_additional.R"
OUT_DIR <- file.path(REPO_ROOT, "papers", "B", "envres")
B_MHEIGHT <- file.path(B_SUBSET, "maternal_height.parquet")

H4_DRY <- c("T1_dry", "T2_dry", "T3_dry")
H4_TMAX <- c("T1_tmax_anom", "T2_tmax_anom", "T3_tmax_anom")
PLACEBO <- c("placebo_dry", "placebo_tmax_anom")
FE0 <- c("region_id", "birth_ym")
SEAS <- list(base = character(), state_cm = "state_gadm^birth_month_of_year",
             district_cm = "region_id^birth_month_of_year")

# mother's height for the rows of d (row-aligned like load_hw57: round, v001, bidx, b3 checked)
load_mheight <- function(d, path = B_MHEIGHT) {
  if (!file.exists(path)) stop(path, " not found: python -I models/B/py/maternal_height.py", call. = FALSE)
  k <- data.table::as.data.table(arrow::read_parquet(path))
  same <- function(a, b) identical(as.integer(a), as.integer(b))
  ok <- nrow(k) == nrow(d) && identical(as.character(k$round), as.character(d$round)) &&
    same(k$v001, d$v001) && same(k$bidx, d$bidx) && same(k$b3, d$birth_cmc)
  if (!ok) stop("maternal_height.parquet is not row-aligned with analysis_child.parquet", call. = FALSE)
  as.numeric(k$mheight_cm)
}

main <- function() {
  set.seed(SEED_B)
  d <- load_b()
  d[, mheight_cm := load_mheight(d)]
  d[, state_gadm := factor(state_gadm)]
  bw <- d[s_bw == TRUE]
  hz <- d[s_haz2459 == TRUE]
  tabs <- list()
  keep <- function(t, analysis, spec) {
    .a <- analysis; .s <- spec
    t[, `:=`(analysis = .a, spec = .s)]
    tabs[[length(tabs) + 1L]] <<- t
    t
  }

  # ---- A1: seasonality and the heat placebo
  for (s in names(SEAS)) {
    fe <- c(FE0, SEAS[[s]])
    keep(run_model(paste0("X-placebo-", s), "bw_g", c(H4_DRY, H4_TMAX, PLACEBO, COV), fe, bw,
                   c(H4_DRY, H4_TMAX, PLACEBO)), "A1_placebo", s)
    keep(run_model(paste0("X-h4bw-", s), "bw_g", c(H4_DRY, H4_TMAX, COV), fe, bw, c(H4_DRY, H4_TMAX)), "A1_h4bw", s)
    keep(run_model(paste0("X-h5haz-", s), "haz", c("dry1000_2459", COV), c(fe, "age_fe"), hz, "dry1000_2459"),
         "A1_h5haz", s)
  }

  # ---- A2: mother's height as a negative-control outcome
  bwh <- bw[!is.na(mheight_cm)]
  for (s in c("base", "state_cm")) {
    fe <- c(FE0, SEAS[[s]])
    keep(run_model(paste0("X-nc-height-", s), "mheight_cm", c(H4_DRY, H4_TMAX, COV), fe, bwh, c(H4_DRY, H4_TMAX)),
         "A2_height", s)
    keep(run_model(paste0("X-nc-height-placebo-", s), "mheight_cm", c(H4_DRY, H4_TMAX, PLACEBO, COV), fe, bwh,
                   c(H4_DRY, H4_TMAX, PLACEBO)), "A2_height_placebo", s)
    keep(run_model(paste0("X-bw-on-height-", s), "bw_g", c(H4_DRY, H4_TMAX, "mheight_cm", COV), fe, bwh, "mheight_cm"),
         "A2_gradient", s)
  }
  tab <- data.table::rbindlist(tabs, fill = TRUE)

  # exposure IQRs (weighted) in the birth-weight sample with height, for the A2 conversion
  iqr <- vapply(c(H4_DRY, H4_TMAX, PLACEBO), function(v) exp_stats(bwh[[v]], bwh$weight)[["iqr"]], numeric(1))

  # ---- entries
  e <- list()
  get_row <- function(an, sp, term) row_of(tab[analysis == an & spec == sp], term)
  lab <- c(base = "district and birth year-month FE", state_cm = "plus state x birth calendar-month FE",
           district_cm = "plus district x birth calendar-month FE")
  for (s in names(SEAS)) {
    for (term in c(H4_DRY, H4_TMAX, PLACEBO)) {
      r <- get_row("A1_placebo", s, term)
      dg <- 1
      e[[sprintf("B.x.a1.placebo.%s.%s", s, term)]] <- ent_row(r, dg, S, paste("birth weight (g per unit), placebo model,", lab[[s]]))
    }
    e[[sprintf("B.x.a1.placebo.%s.n", s)]] <- ent_int(get_row("A1_placebo", s, "placebo_tmax_anom")$n, S, "n, placebo model")
    for (term in c("T2_dry", "T3_dry", H4_TMAX))
      e[[sprintf("B.x.a1.h4bw.%s.%s", s, term)]] <- ent_row(get_row("A1_h4bw", s, term), 1, S,
                                                           paste("birth weight (g per unit), primary H4 model,", lab[[s]]))
    e[[sprintf("B.x.a1.h4bw.%s.n", s)]] <- ent_int(get_row("A1_h4bw", s, "T2_dry")$n, S, "n, H4 model")
    e[[sprintf("B.x.a1.h5haz.%s.dry1000", s)]] <- ent_row(get_row("A1_h5haz", s, "dry1000_2459"), 3, S,
                                                         paste("HAZ 24-59 m per unit of 1,000-day dryness,", lab[[s]]))
  }
  # decision rule A1 (model a = state_cm)
  base_pl <- get_row("A1_placebo", "base", "placebo_tmax_anom")
  a_pl <- get_row("A1_placebo", "state_cm", "placebo_tmax_anom")
  explained <- a_pl$lo <= 0 && a_pl$hi >= 0 && abs(a_pl$est) < 7.65
  e[["B.x.a1.rule.explained"]] <- ent_text(if (explained) "explained by seasonality" else "not explained by seasonality", S,
    sprintf("A1 rule: state x calendar-month placebo heat CI includes 0 and |est| < 7.65 g; base est %.2f, model (a) est %.2f [%.2f, %.2f]",
            base_pl$est, a_pl$est, a_pl$lo, a_pl$hi))

  # A2: height per IQR, implied birth-weight bias per IQR, flags
  flags <- character()
  for (s in c("base", "state_cm")) {
    g <- get_row("A2_gradient", s, "mheight_cm")
    e[[sprintf("B.x.a2.gradient.%s", s)]] <- ent_row(g, 1, S, paste("birth weight (g) per cm of mother's height,", lab[[s]]))
    for (an in c("A2_height", "A2_height_placebo")) {
      terms <- if (an == "A2_height") c(H4_DRY, H4_TMAX) else PLACEBO
      for (term in terms) {
        r <- get_row(an, s, term)
        q <- iqr[[term]]
        rc <- rescale_ci(r$est, r$lo, r$hi, q)
        key <- sprintf("B.x.a2.height.%s.%s", s, term)
        e[[key]] <- ent_ci(rc$est, rc$lo, rc$hi, 3, S, sprintf("mother's height (cm) per IQR (%.4g) of %s, %s", q, term, lab[[s]]))
        bias <- rc$est * g$est
        e[[paste0(key, ".bwbias")]] <- .ent(bias, fmt_num(bias, 1), S, "implied birth-weight bias, g per IQR (height est x gradient)")
        if ((rc$lo > 0 || rc$hi < 0) && abs(bias) > 2) flags <- c(flags, sprintf("%s:%s", s, term))
      }
    }
    e[[sprintf("B.x.a2.n.%s", s)]] <- ent_int(get_row("A2_height", s, "T2_dry")$n, S, "n, height model")
  }
  e[["B.x.a2.rule.flags"]] <- ent_text(if (length(flags)) paste(flags, collapse = "; ") else "none", S,
                                      "A2 rule: terms whose height CI excludes 0 and implied bias > 2 g per IQR")

  # ---- A3: smallest effects of interest, from the primary rescaled estimates (results.json)
  res <- jsonlite::fromJSON(file.path(REPO_ROOT, "papers", "B", "results.json"), simplifyVector = FALSE)
  lo_of <- function(k) as.numeric(sub(".*lo=([-0-9.eE]+).*", "\\1", res[[k]]$description))
  bench <- list(
    list(k = "B.precision.h4bw.T2_dry.per_iqr", b = 15, id = "t2_bw"),
    list(k = "B.precision.h4bw.T3_dry.per_iqr", b = 15, id = "t3_bw"),
    list(k = "B.precision.h5haz2459.dry1000.per_iqr", b = 0.05, id = "dry1000_haz"),
    list(k = "B.precision.ensoxs.bw.at_p90", b = 15, id = "ensoxs_bw_p90"),
    list(k = "B.precision.ensoxs.haz2459.at_p90", b = 0.05, id = "ensoxs_haz_p90"))
  for (x in bench) {
    lo <- lo_of(x$k)
    e[[paste0("B.x.a3.", x$id)]] <- ent_text(if (lo > -x$b) "ruled out" else "not ruled out", S,
      sprintf("A3: harmful effect >= %s ruled out if lower 95%% limit (%.4g, from %s) > -%s", x$b, lo, x$k, x$b))
  }

  dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
  out_tab <- tab[, .(analysis, spec, model, outcome, term, est, lo, hi, n, fe, vcov_fixed)]
  check_aggregate_table("additional_models", out_tab)
  data.table::fwrite(out_tab, file.path(OUT_DIR, "additional_models.csv"))
  write_results(e, file.path(OUT_DIR, "results_additional.json"))
  log_msg("13_additional: wrote", length(e), "entries")
  invisible(e)
}

if (sys.nframe() == 0L) main()
