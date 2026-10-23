# Paper B Task 7 helpers: exposure construction, covariate bands, verdict rules, contract key
# expansion and small-cell suppression (papers/B/results_keys.md "Definitions").
# Source after common/R/load.R (working directory <repo>/models). Pure functions only;
# the data preparation for the real dataset is prep_b() at the end.

# ---------------------------------------------------------------- time

# YYYYMM integer shifted by k months.
ym_add <- function(ym, k) {
  ym <- as.integer(ym); k <- as.integer(k)
  m0 <- (ym %/% 100L) * 12L + (ym %% 100L - 1L) + k
  as.integer((m0 %/% 12L) * 100L + m0 %% 12L + 1L)
}

# June-May ENSO year of the T2 midpoint month (birth month - 5), labelled by its June's year.
enso_year_of <- function(birth_ym) {
  mid <- ym_add(birth_ym, -5L)
  y <- mid %/% 100L; m <- mid %% 100L
  as.integer(ifelse(m >= 6L, y, y - 1L))
}

# ---------------------------------------------------------------- exposures

# Month-weighted mean of window means; NaN (NA) if any component is missing.
mw_mean <- function(parts, w) {
  if (length(parts) != length(w)) stop("mw_mean: one weight per component (weights)", call. = FALSE)
  num <- Reduce(`+`, Map(function(x, wi) wi * as.numeric(x), parts, w))
  num / sum(w)
}

# 1,000-day SPI-3 (ruling P13 Q3): 24-59 -> (3,3,3,12,12)/33 over T1,T2,T3,inf0_11,inf12_23;
# 0-23 -> (3,3,3,12)/21 over T1-T3 and inf0_11 (inf12_23 ignored). NaN if a used part is NaN.
spi3_1000d <- function(T1, T2, T3, inf0_11, inf12_23, group = c("2459", "0023")) {
  group <- match.arg(group)
  if (group == "2459") mw_mean(list(T1, T2, T3, inf0_11, inf12_23), c(3, 3, 3, 12, 12))
  else mw_mean(list(T1, T2, T3, inf0_11), c(3, 3, 3, 12))
}

# ---------------------------------------------------------------- covariates

band_bord <- function(bord) {
  b <- as.integer(bord)
  factor(ifelse(b >= 4L, "4+", as.character(b)), levels = c("1", "2", "3", "4+"))
}

band_v012 <- function(v012) {
  a <- as.integer(v012)
  if (any(!is.na(a) & a < 15L)) stop("band_v012: v012 below 15", call. = FALSE)
  cut(a, c(14, 19, 24, 29, 34, Inf), labels = c("15-19", "20-24", "25-29", "30-34", "35+"))
}

caste_f <- function(s116) {
  lab <- c(`1` = "SC", `2` = "ST", `3` = "OBC", `4` = "none")
  x <- lab[as.character(s116)]
  x[is.na(x)] <- "missing"
  factor(unname(x), levels = c("SC", "ST", "OBC", "none", "missing"))
}

religion_f <- function(v130) {
  lab <- c(`1` = "Hindu", `2` = "Muslim", `3` = "Christian", `4` = "Sikh", `5` = "Buddhist")
  x <- lab[as.character(v130)]
  x[is.na(x)] <- "other"
  factor(unname(x), levels = c("Hindu", "Muslim", "Christian", "Sikh", "Buddhist", "other"))
}

# Spec §6: v104 == 95 (always) or years in place >= child's age in years + 1; visitors (96)
# and codes > 96 excluded.
nonmover <- function(v104, age_m) {
  v <- as.integer(v104)
  need <- floor(as.numeric(age_m) / 12) + 1
  out <- v == 95L | (v < 95L & v >= need)
  out[is.na(out)] <- FALSE
  out
}

multi_birth_mothers <- function(mother_id) {
  n <- stats::ave(rep(1L, length(mother_id)), mother_id, FUN = length)
  n >= 2L
}

# ---------------------------------------------------------------- rules and display

# Dry-orientation decision rule (P13): negative estimate with the 95% CI upper bound < 0.
rule_dry_negative <- function(est, hi) {
  r <- est < 0 & hi < 0
  r[is.na(r)] <- FALSE
  r
}

verdict_text <- function(ok) ifelse(ok, "consistent with", "not consistent with")

verdict_terms <- function(t2, t3) {
  ifelse(t2 & t3, "T2 and T3", ifelse(t2, "T2", ifelse(t3, "T3", "neither")))
}

# Coefficient on the opposite orientation (dry <-> SPI, or La Nina per 1 degree cooler).
flip_ci <- function(est, lo, hi) list(est = -est, lo = -hi, hi = -lo)

# ---------------------------------------------------------------- contract keys

# Shell-style brace expansion ("{a,b}", empty alternatives allowed); no nesting.
expand_keys <- function(pattern) {
  out <- character()
  todo <- pattern
  while (length(todo)) {
    p <- todo[1]; todo <- todo[-1]
    m <- regexpr("\\{[^{}]*\\}", p)
    if (m < 0) { out <- c(out, p); next }
    inner <- substr(p, m + 1L, m + attr(m, "match.length") - 2L)
    alts <- strsplit(inner, ",", fixed = TRUE)[[1]]
    if (endsWith(inner, ",")) alts <- c(alts, "")
    if (!length(alts)) alts <- ""
    pre <- substr(p, 1L, m - 1L); post <- substr(p, m + attr(m, "match.length"), nchar(p))
    todo <- c(paste0(pre, alts, post), todo)
  }
  unique(out)
}

# Keys of a results_keys.md: backticked "B." tokens in the FIRST cell of markdown table rows.
contract_keys <- function(path) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  rows <- grep("^\\|\\s*`B\\.", lines, value = TRUE)
  first <- vapply(strsplit(rows, "|", fixed = TRUE), function(x) x[2], character(1))
  toks <- unlist(regmatches(first, gregexpr("`B\\.[^`]+`", first)))
  toks <- gsub("`", "", toks)
  unique(unlist(lapply(toks, expand_keys)))
}

# ---------------------------------------------------------------- small cells

MIN_CELL_N <- 30L

# Published size band of a count (fix round 1, ruling Q2): exact n is not published where a
# suppressed cell could be recovered by differencing. NA counts as 0 ("<30").
N_BANDS <- c("<30", "30\u201399", "100\u2013499", "500+")
band_n <- function(n, min_n = MIN_CELL_N) {
  n <- as.numeric(n); n[is.na(n)] <- 0
  out <- ifelse(n < min_n, N_BANDS[1], ifelse(n < 100, N_BANDS[2], ifelse(n < 500, N_BANDS[3], N_BANDS[4])))
  enc2utf8(out)
}

suppress_small <- function(dt, n_col, value_cols, min_n = MIN_CELL_N) {
  dt <- data.table::copy(data.table::as.data.table(dt))
  small <- !is.na(dt[[n_col]]) & dt[[n_col]] < min_n
  for (v in value_cols) data.table::set(dt, i = which(small), j = v, value = NA)
  data.table::set(dt, i = which(small), j = n_col, value = NA)
  dt[, suppressed := small]
  dt[]
}

# ---------------------------------------------------------------- ruling P14 polish

# M4: evaluate expr; fixest's "VCOV matrix is not positive definite and was 'fixed'" warning is
# muffled and recorded (vcov_fixed = TRUE); every other warning passes through.
with_vcov_fix_flag <- function(expr) {
  fixed <- FALSE
  value <- withCallingHandlers(expr, warning = function(w) {
    if (grepl("not positive definite", conditionMessage(w), fixed = TRUE)) {
      fixed <<- TRUE
      invokeRestart("muffleWarning")
    }
  })
  list(value = value, vcov_fixed = fixed)
}

# M8: scale and per-unit of a coefficient as tabulated by fe_table (logit rows are odds ratios).
unit_of <- function(lhs, family, term) {
  scale <- if (identical(family, "logit")) "odds ratio"
           else switch(lhs, bw_g = "grams", haz = "HAZ (z)", waz = "WAZ (z)", whz = "WHZ (z)",
                       lbw = , stunted = "probability (x100 = percentage points)", lhs)
  # the ENSO x sensitivity interaction is tested FIRST: its name ends in "_dry" (ruling Q2b); the
  # text matches the B.ensoxs.* key descriptions
  per <- if (term %in% c("inter", "elnino_x_sens_dry", "elnino_x_sens_dry_xw"))
           "per unit of preg_elnino_dev x sens_dry (sens_dry = -sens_beta, % JJAS rain loss per deg C ONI)"
         else if (grepl("(^|_)dry($|_|1000)", term)) "per +1 SPI unit drier"
         else if (grepl("tmax", term)) "per +1 deg C"
         else if (grepl("_pos$", term)) "per +1 deg C of max(x, 0) (El Nino term)"
         else if (grepl("_neg$", term)) "per +1 deg C of min(x, 0) (raw coefficient; La Nina keys report -beta per 1 deg C cooler)"
         else if (grepl("_spi3$", term)) "per +1 SPI-3"
         else if (grepl("^T[123]_oni$", term)) "per +1 deg C of the trimester mean ONI"
         else paste("per unit of", term)
  paste(scale, per)
}

# M5: conception_shift.parquet metadata (b_source_md5 / b_source_size, written by
# models/B/py/conception_shift.py) against the analysis parquet's md5 and size.
fingerprint_matches <- function(meta, fp) {
  md5 <- meta[["b_source_md5"]]; size <- meta[["b_source_size"]]
  !is.null(md5) && !is.null(size) && identical(as.character(md5), as.character(fp$md5)) &&
    identical(as.character(size), format(fp$size, scientific = FALSE))
}

# M9: names never allowed in a committed table, aligned with tests/test_no_controlled_data.py
# (P12b): FORBIDDEN_COLUMNS | COORDINATE_COLUMNS | CHILD_LEVEL_COLUMNS, the "anaemia" prefix, and
# birth_ym / region_id (a combination there; banned singly here).
B_FORBIDDEN_COLUMNS <- c(
  "v001", "v002", "v003", "caseid", "hw1", "b3", "dhsclust", "dhsid", "dhscc",
  "lat", "lon", "latnum", "longnum", "latitude", "longitude", "lng", "long",
  "bw_g", "haz", "waz", "whz", "lbw", "stunted", "wasted", "underweight", "hb_adj", "diarrhoea",
  "mother_id", "birth_cmc", "interview_cmc", "hw70", "hw71", "hw72", "hw56", "m19",
  "birth_ym", "region_id")

# Every numeric count-like column (n, N, n_*, count, count_*) must be >= MIN_CELL_N or NA
# (suppressed).
check_aggregate_table <- function(name, dt) {
  nm <- tolower(names(dt))
  bad <- c(intersect(nm, B_FORBIDDEN_COLUMNS), grep("^anaemia", nm, value = TRUE))
  if (length(bad)) stop("08_results: table ", name, " has forbidden column(s): ", paste(bad, collapse = ", "), call. = FALSE)
  counts <- names(dt)[grepl("^(n|count)(_|$)", names(dt), ignore.case = TRUE) &
                        vapply(dt, is.numeric, logical(1))]
  for (c in counts) {
    v <- dt[[c]]
    if (any(!is.na(v) & v < MIN_CELL_N))
      stop("08_results: table ", name, " column ", c, " has a cell with n < ", MIN_CELL_N, call. = FALSE)
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------- dataset

SEED_B <- 20261022L
COV <- c("sex", "bord_band", "v012_band", "edu", "wealth", "residence", "caste", "religion")
CMC_NS <- paste0("cmc_ns", 1:4)
WINDOWS_1000 <- c("T1", "T2", "T3", "inf0_11", "inf12_23")

# All Task 7 derived columns on the analysis frame (data.table, modified and returned).
prep_b <- function(d) {
  d <- data.table::as.data.table(d)
  d[, sex := factor(b4)]
  d[, bord_band := band_bord(bord)]
  d[, v012_band := band_v012(v012)]
  d[, edu := factor(v106)]
  d[, wealth := factor(v190)]
  d[, residence := factor(v025)]
  d[, caste := caste_f(s116)]
  d[, religion := religion_f(v130)]
  d[, region_id := factor(region_id)]
  d[, region_id_xw := factor(region_id_xw)]
  d[, round := factor(round)]
  d[, birth_ym := as.integer(birth_ym)]
  d[, interview_ym := as.integer(interview_ym)]
  d[, birth_month_of_year := as.integer(birth_month_of_year)]
  d[, age_fe := as.integer(ifelse(is.na(hw1), age_m, hw1))]
  d[, enso_year := enso_year_of(birth_ym)]
  d[, twin := as.logical(twin)]
  d[, urban := as.logical(urban)]
  d[, preg_elnino_dev := as.logical(preg_elnino_dev)]
  # dry orientation (P13)
  for (w in c("T1", "T2", "T3", "placebo", "recent1", "recent12")) {
    d[, (paste0(w, "_dry")) := -as.numeric(get(paste0(w, "_spi3")))]
    d[, (paste0("dist_", w, "_dry")) := -as.numeric(get(paste0("dist_", w, "_spi3")))]
  }
  # buffer (rural uses buf10, urban plain)
  for (w in WINDOWS_1000) {
    d[, (paste0(w, "_spi3_b10r")) := ifelse(urban %in% TRUE, as.numeric(get(paste0(w, "_spi3"))),
                                            as.numeric(get(paste0(w, "_spi3_buf10"))))]
  }
  for (w in c("T1", "T2", "T3")) d[, (paste0(w, "_dry_b10r")) := -get(paste0(w, "_spi3_b10r"))]
  mk1000 <- function(prefix, suffix, group) {
    g <- function(w) as.numeric(d[[paste0(prefix, w, suffix)]])
    -spi3_1000d(g("T1"), g("T2"), g("T3"), g("inf0_11"), g("inf12_23"), group)
  }
  d[, dry1000_2459 := mk1000("", "_spi3", "2459")]
  d[, dry1000_0023 := mk1000("", "_spi3", "0023")]
  d[, dist_dry1000_2459 := mk1000("dist_", "_spi3", "2459")]
  d[, b10r_dry1000_2459 := mk1000("", "_spi3_b10r", "2459")]
  d[, sens_dry := -sens_beta]
  d[, sens_dry_xw := -sens_beta_xw]
  d[, elnino_x_sens_dry := as.numeric(preg_elnino_dev) * sens_dry]
  d[, elnino_x_sens_dry_xw := as.numeric(preg_elnino_dev) * sens_dry_xw]
  for (ix in c("oni", "roni", "nino34")) {
    x <- as.numeric(d[[paste0("preg_", ix)]])
    d[, (paste0("preg_", ix, "_pos")) := pmax(x, 0)]
    d[, (paste0("preg_", ix, "_neg")) := pmin(x, 0)]
  }
  ns <- splines::ns(d$birth_cmc, df = 4)
  for (k in 1:4) d[, (CMC_NS[k]) := as.numeric(ns[, k])]
  # samples
  d[, s_bw := !twin & !is.na(bw_g)]
  d[, s_haz2459 := !is.na(haz) & age_fe >= 24L & age_fe <= 59L]
  d[, s_haz0023 := !is.na(haz) & age_fe >= 0L & age_fe <= 23L]
  d[, s_anaem := !is.na(anaemia_any)]
  d[, s_diarr := !is.na(diarrhoea)]
  d[, nonmover := nonmover(v104, age_m)]
  d[]
}

# Flags for the noflag sensitivity: TRUE if any listed window has a fallback-rain or MIS flag.
any_flag <- function(d, windows) {
  out <- rep(FALSE, nrow(d))
  for (w in windows) for (f in c("_fallback_rain_any", "_mis_district_any")) {
    v <- d[[paste0(w, f)]]
    out <- out | (!is.na(v) & as.logical(v))
  }
  out
}
