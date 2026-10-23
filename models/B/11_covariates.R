# Paper B ruling Q3: covariate distribution and missingness table (results_keys.md "Amendment
# 2026-10-23 (ruling Q3)", section 10). Descriptive, no verdict. Samples: all children (cov) and
# the BW analytic sample (covbw), by round. Levels as modelled (prep_b / COV); "missing" = raw NA
# or a code outside the levels (caste: code 8 or NA, the modelled "missing" category).
# Disclosure (ruling Q2 rule): cells of 1-29 suppressed; complementary suppression within each
# sample x round x variable group (its cells sum to a published total); a covbw cell whose
# difference from the cov cell is 1-29 is suppressed; true zeros published (zero = TRUE, n NA).

COV_VARS <- c("sex", "birth_order", "mother_age", "education", "wealth", "residence", "caste", "religion")
COV_MAPS <- list(
  sex = c(`1` = "male", `2` = "female"),
  birth_order = c(`1` = "1", `2` = "2", `3` = "3", `4+` = "4plus"),
  mother_age = c(`15-19` = "15_19", `20-24` = "20_24", `25-29` = "25_29", `30-34` = "30_34", `35+` = "35plus"),
  education = c(`0` = "none", `1` = "primary", `2` = "secondary", `3` = "higher"),
  wealth = c(`1` = "poorest", `2` = "poorer", `3` = "middle", `4` = "richer", `5` = "richest"),
  residence = c(`1` = "urban", `2` = "rural"),
  caste = c(`1` = "sc", `2` = "st", `3` = "obc", `4` = "none"),
  religion = c(`1` = "hindu", `2` = "muslim", `3` = "christian", `4` = "sikh", `5` = "buddhist",
               `6` = "other", `7` = "other", `8` = "other", `9` = "other", `96` = "other"))
COV_SOURCE <- c(sex = "b4", birth_order = "bord_band", mother_age = "v012_band", education = "v106",
                wealth = "v190", residence = "v025", caste = "s116", religion = "v130")

# Level of each row for one variable; "missing" for NA or a code outside the levels.
cov_level <- function(d, var) {
  map <- COV_MAPS[[var]]
  if (is.null(map)) stop("cov_level: unknown variable '", var, "'", call. = FALSE)
  x <- d[[COV_SOURCE[[var]]]]
  key <- if (is.factor(x)) as.character(x) else format(x, scientific = FALSE, trim = TRUE)
  key[is.na(x)] <- NA_character_
  out <- unname(map[key])
  out[is.na(out)] <- "missing"
  out
}

# Unweighted n and weighted sum per sample x round x variable x level (full level grid).
cov_cells <- function(d, samples) {
  out <- list()
  for (s in names(samples)) {
    x <- d[samples[[s]]]
    rounds <- sort(unique(as.character(x$round)))
    for (v in COV_VARS) {
      lv <- cov_level(x, v)
      agg <- data.table::data.table(round = as.character(x$round), level = lv, w = x$weight)[
        , .(n = .N, w = sum(w)), by = .(round, level)]
      grid <- data.table::CJ(round = rounds, level = c(unique(unname(COV_MAPS[[v]])), "missing"), sorted = FALSE)
      g <- agg[grid, on = c("round", "level")]
      g[is.na(n), `:=`(n = 0L, w = 0)]
      out[[length(out) + 1L]] <- g[, `:=`(sample = s, variable = v)]
    }
  }
  res <- data.table::rbindlist(out)
  data.table::setcolorder(res, c("sample", "round", "variable", "level", "n", "w"))
  res[]
}

# Complementary rule for one group: a primary suppression (a 1-29 cell, or a covbw cell suppressed
# by the nested rule) is protected by at least one NON-primary suppressed cell with n >= min_n;
# if there is none, the smallest unsuppressed cell with n >= min_n is suppressed too.
.cov_complement <- function(n, s, primary, min_n) {
  if (any(primary) && !any(s & !primary & n >= min_n)) {
    cand <- which(!s & n >= min_n)
    if (length(cand)) s[cand[which.min(n[cand])]] <- TRUE
    else if (sum(s & n > 0L) < 2L)               # every non-zero cell is suppressed: need >= 2 of them
      stop("suppress_cov: a group cannot be protected (one non-zero cell)", call. = FALSE)
  }
  s
}

suppress_cov <- function(cells, min_n = MIN_CELL_N) {
  x <- data.table::copy(data.table::as.data.table(cells))
  x[, zero := n == 0L]
  x[, primary := n > 0L & n < min_n]
  x[, suppressed := primary]
  repeat {
    before <- sum(x$suppressed)
    x[, suppressed := .cov_complement(n, suppressed, primary, min_n), by = .(sample, round, variable)]
    if (all(c("cov", "covbw") %in% x$sample)) {
      m <- merge(x[sample == "cov", .(round, variable, level, n_cov = n, s_cov = suppressed)],
                 x[sample == "covbw", .(round, variable, level, n_bw = n, s_bw = suppressed)],
                 by = c("round", "variable", "level"))
      bad <- m[!s_cov & !s_bw & n_cov - n_bw > 0L & n_cov - n_bw < min_n, paste(round, variable, level)]
      if (length(bad)) x[sample == "covbw" & paste(round, variable, level) %in% bad, `:=`(suppressed = TRUE, primary = TRUE)]
    }
    if (sum(x$suppressed) == before) break
  }
  x[, primary := NULL]
  x[]
}

covariate_table <- function(d, S = "models/B/11_covariates.R", min_n = MIN_CELL_N) {
  d <- data.table::as.data.table(d)
  samples <- list(cov = rep(TRUE, nrow(d)), covbw = d$s_bw %in% TRUE)
  x <- suppress_cov(cov_cells(d, samples), min_n)
  x[, `:=`(N = sum(n), W = sum(w)), by = .(sample, round, variable)]
  x[, `:=`(weighted_pct = 100 * w / W, unweighted_pct = 100 * n / N)]
  e <- list()
  sup <- function(desc) .ent(NA, "suppressed (n < 30)", S, desc)
  lab <- c(cov = "all children in analysis_child.parquet", covbw = "BW analytic sample (singleton births with bw_g)")
  for (i in seq_len(nrow(x))) {
    r <- x[i]
    k <- sprintf("B.supp.%s.%s.%s.%s", r$sample, r$variable, r$level, r$round)
    desc <- sprintf("%s = %s, %s, %s", r$variable, r$level, lab[[r$sample]], r$round)
    e[[k]] <- if (r$suppressed) sup(desc)
      else if (r$level == "missing")
        .ent(as.integer(r$n), sprintf("%s (%s)", fmt_int(r$n), fmt_pct(r$n / r$N, 1)), S,
             paste("Unweighted n missing (unweighted % of the sample x round);", desc))
      else ent_pct(r$w / r$W, S, paste("Weighted share of the sample x round;", desc, "; n =", r$n))
  }
  # survey round: children per round in the pooled sample
  rr <- list()
  for (s in names(samples)) {
    nr <- table(as.character(d$round[samples[[s]]]))
    for (r in names(nr)) {
      e[[sprintf("B.supp.%s.round.%s.all", s, r)]] <- .ent(as.integer(nr[[r]]),
        sprintf("%s (%s)", fmt_int(nr[[r]]), fmt_pct(nr[[r]] / sum(nr), 1)), S,
        sprintf("Children in %s, %s (unweighted %% of the pooled sample)", r, lab[[s]]))
      rr[[length(rr) + 1L]] <- data.table::data.table(sample = s, round = "all", variable = "round", level = r,
        n = as.integer(nr[[r]]), zero = FALSE, weighted_pct = NA_real_, unweighted_pct = 100 * nr[[r]] / sum(nr),
        suppressed = FALSE)
    }
  }
  tab <- x[, .(sample, round, variable, level, n, zero, weighted_pct, unweighted_pct, suppressed)]
  tab[suppressed == TRUE, `:=`(n = NA_integer_, weighted_pct = NA_real_, unweighted_pct = NA_real_)]
  tab[zero == TRUE, n := NA_integer_]          # a true zero: published as zero = TRUE (n column: >= 30 or NA)
  tab <- rbind(tab, data.table::rbindlist(rr))
  list(entries = e, table = tab[])
}

run_covariates <- function() {
  r <- covariate_table(load_b())
  list(entries = r$entries, tables = list(supp_covariates = r$table))
}
