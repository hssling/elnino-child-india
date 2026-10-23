# Paper B Task 7 step 6: ENSO -> T2/T3 SPI-3 -> birth weight mediation, SECONDARY
# (results_keys.md section 6; rulings P10, R-1, P13 Q4). ENSO-year (Jun-May) block bootstrap,
# B = 1000, 6 PSOCK workers, no birth_ym FE; time-clustered Wald check of the a path.
#
# Standalone (background) use from <repo>/models:
#   Rscript B/06_mediation.R            # full run, cached to B/_cache/mediation.rds
#   Rscript B/06_mediation.R --resume   # reuse finished calls (B/_cache/med_*.rds)
#   Rscript B/06_mediation.R --time 12  # timing only: B = 12 draws of the first call, nothing cached

MED_FE <- c("region_id", "birth_month_of_year", "round")
MED_CALLS <- list(list(spec = "oni", exposure = "preg_oni", extra = character()),
                  list(spec = "onipos", exposure = "preg_oni_pos", extra = "preg_oni_neg"))
MED_B <- 1000L
MED_CORES <- 6L

fmt_pct_ci <- function(est, lo, hi) {
  neg <- round(100 * lo) < 0 | round(100 * hi) < 0
  paste0(fmt_num(100 * est, 0), "% (95% CI ", fmt_num(100 * lo, 0), if (neg) " to " else EN_DASH,
         fmt_num(100 * hi, 0), "%)")
}

med_data <- function() {
  d <- load_b()
  d[s_bw == TRUE, c("bw_g", "preg_oni", "preg_oni_pos", "preg_oni_neg", "T2_spi3", "T3_spi3", COV, CMC_NS,
                    MED_FE, "enso_year", "weight"), with = FALSE]
}

# resume = TRUE reuses a call's cached bootstrap (B/_cache/med_<spec>_<m>_B<B>.rds) if present;
# every finished call is cached, so an interrupted background run can be resumed.
run_mediation <- function(B = MED_B, cores = MED_CORES, calls = MED_CALLS, mediators = c("T2", "T3"),
                          resume = FALSE) {
  S <- "models/B/06_mediation.R"
  bw <- med_data()
  e <- list(); rows <- list(); timing <- list()
  for (cl in calls) for (m in mediators) {
    med <- paste0(m, "_spi3")
    covars <- c(cl$extra, COV, CMC_NS)
    log_msg("mediation", cl$spec, m, "B =", B, "cores =", cores)
    cf <- file.path(B_CACHE, sprintf("med_%s_%s_B%d.rds", cl$spec, m, B))
    fp <- b_fingerprint()
    if (resume && file.exists(cf)) {
      cached <- readRDS(cf); r <- cached$r; secs <- cached$secs
      # P14 (M5): the cached bootstrap must come from the current analysis parquet. Caches written
      # before the fingerprint existed (no data_md5) are accepted only if the parquet is older.
      same <- if (!is.null(cached$data_md5)) identical(cached$data_md5, fp$md5) else fp$mtime < file.mtime(cf)
      if (!same) stop("mediation cache ", basename(cf), " was computed on a different analysis_child.parquet; ",
                      "delete it (or run Rscript B/run_all.R without --use-cache)", call. = FALSE)
      log_msg("  cached call reused:", basename(cf),
              if (is.null(cached$data_md5)) "(pre-fingerprint cache; parquet older than cache)" else "(data md5 matches)")
    } else {
      t0 <- Sys.time()
      r <- mediate_fe(bw, outcome = "bw_g", exposure = cl$exposure, mediator = med, fe = MED_FE,
                      covars = covars, B = B, cores = cores, seed = SEED_B, boot_unit = "enso_year")
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      r$boot <- NULL                                  # aggregates only in the cache
      saveRDS(list(r = r, secs = secs, data_md5 = fp$md5), cf)
      log_msg("  done in", round(secs / 60, 1), "min;", round(secs / B, 2), "s per draw (wall, all workers)")
    }
    timing[[paste(cl$spec, m)]] <- secs
    # time-clustered Wald check of the a path, same complete rows
    vars <- unique(c("bw_g", cl$exposure, med, covars, MED_FE, "enso_year", "weight"))
    cc <- bw[stats::complete.cases(bw[, vars, with = FALSE])]
    ta <- run_model(sprintf("MED-a-%s-%s", cl$spec, m), med, c(cl$exposure, covars), MED_FE, cc,
                    cl$exposure, cluster = ~enso_year)
    ra <- row_of(ta, cl$exposure)
    n_blocks <- data.table::uniqueN(cc$enso_year)
    rows[[length(rows) + 1]] <- data.table::data.table(
      spec = cl$spec, exposure = cl$exposure, mediator = med,
      quantity = c("a", "b", "c_total", "c_direct", "indirect", "prop"),
      est = c(r$a, r$b, r$c_total, r$c_direct, r$indirect, r$prop),
      lo = vapply(c("a", "b", "c_total", "c_direct", "indirect", "prop"), function(q) unname(r$ci[[q]][1]), 0),
      hi = vapply(c("a", "b", "c_total", "c_direct", "indirect", "prop"), function(q) unname(r$ci[[q]][2]), 0),
      n = r$n, B = r$B, B_ok = r$B_ok, flag = r$flag, prop_flag = r$prop_flag, enso_year_blocks = n_blocks,
      a_wald_est = ra$est, a_wald_lo = ra$lo, a_wald_hi = ra$hi, a_wald_p = ra$p, a_wald_df = ra$df,
      a_wald_vcov_fixed = ra$vcov_fixed, minutes = secs / 60)
    base <- sprintf("B.med.%s.%s", cl$spec, m)
    desc <- sprintf("Mediation (secondary) %s -> %s -> bw_g; FE %s; covars COV + ns(birth_cmc,4)%s; weighted; ENSO-year block bootstrap B=%d seed=%d (%d blocks)",
                    cl$exposure, med, paste(MED_FE, collapse = "+"),
                    if (length(cl$extra)) paste0(" + ", cl$extra) else "", B, SEED_B, n_blocks)
    boot_ok <- r$flag == "ok"
    bci <- function(q, digits, what) {
      if (!boot_ok) return(.ent(r[[q]], "not estimated (too many failed bootstrap fits)", S, paste(desc, what)))
      ent_ci(r[[q]], r$ci[[q]][1], r$ci[[q]][2], digits, S, paste(desc, what))
    }
    e[[paste0(base, ".a")]] <- bci("a", 3, "; a path (SPI-3 per +1 deg C), bootstrap percentile CI")
    e[[paste0(base, ".b")]] <- bci("b", 1, "; b path (grams per +1 SPI-3), bootstrap CI")
    e[[paste0(base, ".c_total")]] <- bci("c_total", 1, "; total effect (grams per +1 deg C), bootstrap CI")
    e[[paste0(base, ".c_direct")]] <- bci("c_direct", 1, "; direct effect, bootstrap CI")
    e[[paste0(base, ".indirect")]] <- bci("indirect", 1, "; indirect effect a*b (grams per +1 deg C), bootstrap CI")
    e[[paste0(base, ".prop")]] <- if (boot_ok && r$prop_flag == "ok")
      .ent(r$prop, fmt_pct_ci(r$prop, r$ci$prop[1], r$ci$prop[2]), S,
           sprintf("%s; proportion mediated; lo=%s hi=%s", desc, signif(r$ci$prop[1], 6), signif(r$ci$prop[2], 6)))
    else .ent(r$prop, sprintf("not estimated (%s)", if (boot_ok) r$prop_flag else "too many failed bootstrap fits"),
              S, paste(desc, "; proportion mediated"))
    e[[paste0(base, ".a_wald")]] <- ent_row(ra, 3, S, sprintf("a path %s -> %s, Wald CI clustered by enso_year (%d clusters, t df=%s)",
                                                              cl$exposure, med, n_blocks, ra$df))
    e[[paste0(base, ".a_wald.p")]] <- ent_p(ra$p, S, "a path, enso_year-clustered Wald p")
    e[[paste0(base, ".flag")]] <- ent_text(r$flag, S, "mediate_fe flag")
    e[[paste0(base, ".prop_flag")]] <- ent_text(r$prop_flag, S, "mediate_fe prop_flag")
    e[[paste0(base, ".n")]] <- ent_int(r$n, S, "Mediation rows used (complete cases)")
    e[[paste0(base, ".B_ok")]] <- ent_int(r$B_ok, S, sprintf("Successful bootstrap resamples of %d", B))
    if (cl$spec == "oni" && m == "T2")
      e[["B.med.n_blocks"]] <- ent_int(n_blocks, S, "ENSO-year (Jun-May of birth month - 5) blocks, primary T2 mediation sample")
  }
  e[["B.med.B"]] <- ent_int(B, S, "Bootstrap resamples requested per mediation call")
  list(entries = e, tables = list(mediation = data.table::rbindlist(rows)), timing = timing)
}

if (sys.nframe() == 0L) {
  source("common/R/load.R"); source("B/R/fe.R"); source("B/R/prep.R"); source("B/R/io.R")
  set.seed(SEED_B)
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 2L && args[1] == "--time") {
    nb <- as.integer(args[2])
    t0 <- Sys.time()
    r <- run_mediation(B = nb, calls = MED_CALLS[1], mediators = "T2")
    secs <- r$timing[[1]]
    log_msg(sprintf("TIMING: B=%d, cores=%d: %.1f s wall => %.2f s per draw; projected 4 calls x %d draws = %.1f h",
                    nb, MED_CORES, secs, secs / nb, MED_B, 4 * MED_B * secs / nb / 3600))
  } else {
    res <- cached_step("mediation", function() run_mediation(resume = "--resume" %in% args), use_cache = FALSE)
  }
}
