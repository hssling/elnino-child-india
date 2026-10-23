# Paper B Phase 3b Task 1: pre-specified supplementary analyses deferred by spec §14 (c, d, f);
# results_keys.md "Amendment 2026-10-23: Phase 3b supplementary keys", section 8. No verdicts.
#   (c) ENSO-trim: total-ENSO model with trimester ONI terms (bw, HAZ 24-59)
#   (d) hw57 anaemia cross-check (descriptive; hw57 from the row-aligned kr_harmonised.parquet)
#   (f) WAZ / WHZ, underweight / wasted (H5 family) and moderate or severe anaemia (H6 family)

TRIM_ONI <- c("T1_oni", "T2_oni", "T3_oni")
B_KR <- file.path(B_SUBSET, "kr_harmonised.parquet")

SUPP_GROWTH <- list(waz = list(lhs = "waz", family = "gaussian", digits = 3),
                    whz = list(lhs = "whz", family = "gaussian", digits = 3),
                    underweight = list(lhs = "underweight", family = "logit", digits = 2),
                    wasted = list(lhs = "wasted", family = "logit", digits = 2))
SUPP_AGE <- list(`2459` = c(24L, 59L), `0023` = c(0L, 23L))

# ENSO-trim: ENSO-total with T1_oni + T2_oni + T3_oni in place of preg_oni_pos / preg_oni_neg;
# same sample, FE, spline, weights and cluster (ENSO_OUTCOMES from 05_enso.R).
fit_trim_oni <- function(d, out, id, covars = c(COV, CMC_NS)) {
  o <- ENSO_OUTCOMES[[out]]
  if (is.null(o)) stop("fit_trim_oni: unknown outcome '", out, "'", call. = FALSE)
  x <- d[get(o$sample) == TRUE]
  fe <- c("region_id", "birth_month_of_year", "round", if (o$age) "age_fe")
  run_model(id, o$lhs, c(TRIM_ONI, covars), fe, x, TRIM_ONI)
}

# H5 family with the outcome replaced; sample = outcome non-missing and age_fe in the group.
fit_supp_growth <- function(d, outcome, g) {
  o <- SUPP_GROWTH[[outcome]]; age <- SUPP_AGE[[g]]
  x <- d[!is.na(get(o$lhs)) & age_fe >= age[1] & age_fe <= age[2]]
  x[, dry := get(paste0("dry1000_", g))]
  run_model(sprintf("S-%s%s", outcome, g), o$lhs, c("dry", COV), c("region_id", "birth_ym", "age_fe"),
            x, "dry", family = o$family)
}

# ---------------------------------------------------------------- hw57 cross-check

# DHS anaemia level: 1 severe, 2 moderate, 3 mild -> TRUE; 4 not anaemic -> FALSE; else NA.
hw57_anaemic <- function(h) {
  h <- as.numeric(h)
  out <- rep(NA, length(h))
  out[h %in% c(1, 2, 3)] <- TRUE
  out[h %in% 4] <- FALSE
  out
}

# Unweighted Cohen's kappa of two binary classifications: (po - pe) / (1 - pe); 1 if pe = po = 1.
cohen_kappa <- function(x, y) {
  x <- as.logical(x); y <- as.logical(y)
  if (length(x) != length(y)) stop("cohen_kappa: x and y differ in length", call. = FALSE)
  if (anyNA(x) || anyNA(y)) stop("cohen_kappa: missing values", call. = FALSE)
  po <- mean(x == y)
  px <- mean(x); py <- mean(y)
  pe <- px * py + (1 - px) * (1 - py)
  if (isTRUE(all.equal(pe, 1))) return(if (po == 1) 1 else NA_real_)
  (po - pe) / (1 - pe)
}

# hw57 for the rows of d, read from kr_harmonised.parquet (same rows, same order: Task 5
# excluded none); stops unless round, v001, bidx and b3 == birth_cmc agree on every row.
load_hw57 <- function(d, path = B_KR) {
  if (!file.exists(path)) stop("load_hw57: ", path, " not found (python run_all.py --skip-download)", call. = FALSE)
  k <- data.table::as.data.table(arrow::read_parquet(path, col_select = c("round", "caseid", "v001", "bidx", "b3", "hw57")))
  # (round, caseid, bidx) identifies a child; caseid is compared through mother_id, the dense
  # 0-based rank of (round, caseid) assigned by Task 5 (pandas groupby sort = C-locale order)
  if (anyDuplicated(k[, .(round, caseid, bidx)]))
    stop("load_hw57: (round, caseid, bidx) is not unique in kr_harmonised.parquet", call. = FALSE)
  k_mother <- data.table::frankv(k, c("round", "caseid"), ties.method = "dense") - 1L
  same <- function(a, b) identical(as.integer(a), as.integer(b))
  ok <- nrow(k) == nrow(d) && identical(as.character(k$round), as.character(d$round)) &&
    same(k$v001, d$v001) && same(k$bidx, d$bidx) && same(k$b3, d$birth_cmc) && same(k_mother, d$mother_id)
  if (!ok) stop("load_hw57: kr_harmonised.parquet is not row-aligned with analysis_child.parquet", call. = FALSE)
  as.numeric(k$hw57)
}

# Per-round aggregates on the anaemia sample with a valid hw57 (1-4).
hw57_table <- function(d, hw57) {
  out <- list()
  for (r in sort(unique(as.character(d$round)))) {
    sel <- as.character(d$round) == r & d$s_anaem %in% TRUE
    h <- hw57_anaemic(hw57[sel])
    keep <- !is.na(h)
    a <- as.logical(d$anaemia_any[sel][keep]); b <- h[keep]; w <- d$weight[sel][keep]
    s56 <- sum(w * a) / sum(w); s57 <- sum(w * b) / sum(w)
    out[[r]] <- data.table::data.table(round = r, n_children = sum(keep), hw56_share = s56, hw57_share = s57,
                                       diff_pp = 100 * (s57 - s56), agreement = mean(a == b),
                                       kappa = cohen_kappa(a, b))
  }
  tab <- data.table::rbindlist(out)
  if (any(tab$n_children < MIN_CELL_N)) stop("hw57_table: a round has n < 30", call. = FALSE)
  tab
}

# ---------------------------------------------------------------- step

run_supplementary <- function() {
  S <- "models/B/09_supplementary.R"
  d <- load_b()
  e <- list(); tabs <- list()

  # (c) ENSO-trim
  for (out in c("bw", "haz2459")) {
    t <- fit_trim_oni(d, out, paste0("ENSO-trim-", out))
    tabs[[length(tabs) + 1L]] <- t
    dg <- if (out == "bw") 1 else 3
    unit <- if (out == "bw") "grams" else "HAZ"
    for (k in c("T1", "T2", "T3")) {
      r <- row_of(t, paste0(k, "_oni"))
      desc <- sprintf("ENSO-trim %s: %s per +1 deg C of %s_oni (trimester mean oni_anom), with T1-T3_oni + COV + ns(birth_cmc,4); FE region_id+birth_month_of_year+round%s; weighted; cluster region_id+birth_ym; no verdict",
                      out, unit, k, if (out == "bw") "" else "+age_fe")
      e[[sprintf("B.supp.trim_oni.%s.%s", out, k)]] <- ent_row(r, dg, S, desc)
      e[[sprintf("B.supp.trim_oni.%s.%s.p", out, k)]] <- ent_p(r$p, S, paste(desc, "; p"))
    }
    e[[sprintf("B.supp.trim_oni.%s.n", out)]] <- ent_int(t$n[1], S, paste("ENSO-trim n,", out))
  }

  # (f) secondary growth outcomes (H5 family)
  for (o in names(SUPP_GROWTH)) for (g in names(SUPP_AGE)) {
    t <- fit_supp_growth(d, o, g)
    tabs[[length(tabs) + 1L]] <- t
    r <- row_of(t, "dry")
    lab <- if (g == "2459") "24-59 months, 1,000-day (3,3,3,12,12)/33" else "0-23 months, (3,3,3,12)/21 where complete"
    what <- if (SUPP_GROWTH[[o]]$family == "logit") sprintf("logit OR of %s", o) else sprintf("%s (z)", toupper(o))
    desc <- sprintf("S-%s: %s per +1 SPI unit drier, %s; FE region_id+birth_ym+age_fe; weighted; cluster region_id+birth_ym; no verdict", o, what, lab)
    base <- sprintf("B.supp.%s%s.dry1000", o, g)
    e[[base]] <- ent_row(r, SUPP_GROWTH[[o]]$digits, S, desc)
    e[[paste0(base, ".p")]] <- ent_p(r$p, S, paste(desc, "; p"))
    e[[paste0(base, ".n")]] <- ent_int(r$n, S, paste(desc, "; n"))
  }

  # (f) moderate or severe anaemia (H6 family)
  ta <- run_model("S-anaemia-modsev", "anaemia_modsev", c("recent12_dry", "recent12_tmax_anom", COV),
                  c("region_id", "interview_ym", "age_fe"), d[s_anaem == TRUE & !is.na(anaemia_modsev)],
                  c("recent12_dry", "recent12_tmax_anom"), family = "logit")
  tabs[[length(tabs) + 1L]] <- ta
  for (k in c("dry", "tmax")) {
    r <- row_of(ta, if (k == "dry") "recent12_dry" else "recent12_tmax_anom")
    desc <- sprintf("S-anaemia-modsev logit: OR of moderate or severe anaemia (hw56 raw < 100) per +1 %s over the 12 months before interview; FE region_id+interview_ym+age_fe; weighted; cluster region_id+birth_ym; no verdict",
                    if (k == "dry") "SPI unit drier" else "deg C Tmax anomaly")
    e[[sprintf("B.supp.anaemia_modsev.recent12_%s", k)]] <- ent_row(r, 2, S, desc)
    e[[sprintf("B.supp.anaemia_modsev.recent12_%s.p", k)]] <- ent_p(r$p, S, paste(desc, "; p"))
  }
  e[["B.supp.anaemia_modsev.n"]] <- ent_int(ta$n[1], S, "S-anaemia-modsev n")

  # (d) hw57 cross-check
  hw <- hw57_table(d, load_hw57(d))
  for (i in seq_len(nrow(hw))) {
    r <- hw[i]; b <- sprintf("B.supp.hw57.%s", r$round)
    lab <- sprintf("%s anaemia sample (6-59 months, tested) with hw57 in 1-4", r$round)
    e[[paste0(b, ".hw56_pct")]] <- ent_pct(r$hw56_share, S, paste("Weighted share any anaemia by hw56 (raw < 110),", lab))
    e[[paste0(b, ".hw57_pct")]] <- ent_pct(r$hw57_share, S, paste("Weighted share any anaemia by hw57 (levels 1-3), same rows,", lab))
    e[[paste0(b, ".diff_pp")]] <- .ent(r$diff_pp, fmt_num(r$diff_pp, 2), S, paste("hw57 share minus hw56 share, percentage points,", lab))
    e[[paste0(b, ".agree_pct")]] <- ent_pct(r$agreement, S, paste("Unweighted share classified the same by hw56 and hw57,", lab))
    k <- r$kappa
    e[[paste0(b, ".kappa")]] <- if (is.finite(k)) .ent(k, fmt_num(k, 3), S, paste("Unweighted Cohen's kappa, hw56 vs hw57 any anaemia,", lab))
                                else .ent(NA, "not estimable (pe = 1, po < 1)", S, paste("Unweighted Cohen's kappa,", lab))
    e[[paste0(b, ".n")]] <- ent_int(r$n_children, S, paste("Children in the hw57 cross-check,", lab))
  }

  list(entries = e, tables = list(supp_models = data.table::rbindlist(tabs), supp_hw57 = hw))
}
