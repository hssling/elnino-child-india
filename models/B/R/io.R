# Paper B Task 7 shared plumbing: paths, cache, data loading, model wrappers and results
# entries. Source after common/R/load.R, B/R/fe.R and B/R/prep.R (working directory
# <repo>/models).

B_DIR <- file.path(MODELS_ROOT, "B")
B_CACHE <- file.path(B_DIR, "_cache")
B_PAPER <- file.path(REPO_ROOT, "papers", "B")
B_OUT <- file.path(B_PAPER, "output")
B_SUBSET <- file.path(REPO_ROOT, "data", "external", "nfhs_subset")
for (p in c(B_CACHE, B_OUT, file.path(B_DIR, "logs"))) dir.create(p, recursive = TRUE, showWarnings = FALSE)

B_CLUSTER <- ~region_id + birth_ym

log_msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
  flush.console()
}

# ---------------------------------------------------------------- cache

cache_path <- function(step) file.path(B_CACHE, paste0(step, ".rds"))

# Run fun() unless use_cache and the step's cache exists; the result is cached either way.
cached_step <- function(step, fun, use_cache = FALSE) {
  f <- cache_path(step)
  if (use_cache && file.exists(f)) {
    log_msg("step", step, ": cached")
    return(readRDS(f))
  }
  t0 <- Sys.time()
  data_md5 <- b_fingerprint()$md5                 # ruling P14 (M5/M6): data the step was computed on
  res <- fun()
  res$runtime_min <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  res$data_md5 <- data_md5
  saveRDS(res, f)
  log_msg("step", step, ": done in", round(res$runtime_min, 1), "min")
  res
}

# ---------------------------------------------------------------- data

.B_DATA <- new.env()
B_ANALYSIS <- file.path(B_SUBSET, "analysis_child.parquet")
B_SHIFT <- file.path(B_SUBSET, "conception_shift.parquet")
B_SHIFT_CMD <- "cd <repo> && python models/B/py/conception_shift.py   (or: python run_all.py --skip-download)"

# Fingerprint of the analysis parquet: md5 (base R), size in bytes, parquet row count, mtime.
b_fingerprint <- function(path = B_ANALYSIS) {
  if (!file.exists(path)) stop("b_fingerprint: ", path, " not found (build it: python run_all.py --skip-download)", call. = FALSE)
  list(md5 = unname(tools::md5sum(path)), size = file.size(path),
       rows = arrow::ParquetFileReader$create(path)$num_rows, mtime = file.mtime(path))
}

# Metadata of conception_shift.parquet (b_source_sha256 / md5 / size / rows).
b_shift_meta <- function(path = B_SHIFT) {
  if (!file.exists(path)) stop("conception_shift.parquet not found; build it: ", B_SHIFT_CMD, call. = FALSE)
  as.list(arrow::read_parquet(path, as_data_frame = FALSE)$metadata)
}

# Stop unless conception_shift.parquet was built from the current analysis parquet (M5).
check_shift_fingerprint <- function(fp = b_fingerprint()) {
  meta <- b_shift_meta()
  if (!fingerprint_matches(meta, fp))
    stop("conception_shift.parquet was built from a different analysis_child.parquet (source md5 ",
         meta[["b_source_md5"]], " vs current ", fp$md5, "); rebuild it: ", B_SHIFT_CMD, call. = FALSE)
  invisible(meta)
}

# Analysis frame with the Task 7 derived columns (prep_b) and the conception-shift / post-birth
# ONI columns; cached in memory for the session.
load_b <- function() {
  if (!is.null(.B_DATA$d)) return(.B_DATA$d)
  check_shift_fingerprint()
  d <- data.table::as.data.table(arrow::read_parquet(B_ANALYSIS))
  x <- data.table::as.data.table(arrow::read_parquet(B_SHIFT))
  if (nrow(x) != nrow(d) || !identical(as.integer(x$birth_cmc), as.integer(d$birth_cmc)) ||
      !identical(as.integer(x$row), seq_len(nrow(d)) - 1L))
    stop("load_b: conception_shift.parquet is not aligned with analysis_child.parquet; rebuild it: ",
         B_SHIFT_CMD, call. = FALSE)
  x[, c("row", "birth_cmc") := NULL]
  for (c in names(x)) data.table::set(d, j = c, value = as.numeric(x[[c]]))
  d <- prep_b(d)
  for (s in c("sm1", "sp1")) for (w in c("T1", "T2", "T3"))
    d[, (paste0(w, "_dry_", s)) := -get(paste0(w, "_spi3_", s))]
  d[, post_oni_pos := pmax(post_oni, 0)]
  d[, post_oni_neg := pmin(post_oni, 0)]
  .B_DATA$d <- d
  d
}

# ---------------------------------------------------------------- models

fml <- function(lhs, rhs, fe) stats::as.formula(
  paste(lhs, "~", paste(rhs, collapse = " + "), "|", paste(fe, collapse = " + ")), env = globalenv())

# Fit through fit_fe and tabulate through fe_table with the SAME cluster (R-2).
run_model <- function(id, lhs, rhs, fe, data, terms, family = "gaussian", weights = ~weight,
                      cluster = B_CLUSTER) {
  t0 <- Sys.time()
  # M4: fixest's non-PD vcov fix (eigenvalue correction) is recorded per fit, not just warned
  r <- with_vcov_fix_flag({
    fit <- fit_fe(fml(lhs, rhs, fe), data, weights = weights, cluster = cluster, family = family)
    fe_table(fit, terms, cluster = cluster)
  })
  tab <- r$value
  .fam <- family; .lhs <- lhs                     # data.table scoping: tab has no such columns, but be explicit
  tab[, `:=`(model = id, outcome = .lhs, family = .fam, weighted = !is.null(weights),
             cluster = paste(all.vars(cluster), collapse = "+"), fe = paste(fe, collapse = "+"),
             vcov_fixed = r$vcov_fixed)]
  tab[, unit := vapply(term, function(t) unit_of(.lhs, .fam, t), character(1))]
  data.table::setcolorder(tab, c("model", "outcome", "term"))
  log_msg("  model", id, "n =", tab$n[1], sprintf("(%.1f s)", as.numeric(difftime(Sys.time(), t0, units = "secs"))),
          if (r$vcov_fixed) "[vcov not positive definite: fixest eigenvalue fix applied]" else "")
  tab[]
}

# Row of a table for one term (est, lo, hi, p, n).
row_of <- function(tab, term) {
  .want <- term                                   # data.table evaluates i in the columns' scope
  r <- tab[tab[["term"]] %in% .want, ]
  if (nrow(r) != 1L) stop("row_of: term '", .want, "' not found exactly once", call. = FALSE)
  r
}

# ---------------------------------------------------------------- results entries

.ent <- function(value, display, script, description) res_entry(value, display, script, description)

ent_ci <- function(est, lo, hi, digits, script, description) {
  if (!is.finite(est) || !is.finite(lo) || !is.finite(hi))
    return(.ent(NA, "not estimable (non-finite estimate)", script, description))
  .ent(est, fmt_ci(est, lo, hi, digits = digits), script,
       sprintf("%s; lo=%s hi=%s", description, signif(lo, 6), signif(hi, 6)))
}

ent_row <- function(r, digits, script, description, scale = 1) {
  ent_ci(scale * r$est, scale * r$lo, scale * r$hi, digits, script, description)
}

ent_flip <- function(r, digits, script, description) {
  f <- flip_ci(r$est, r$lo, r$hi)
  ent_ci(f$est, f$lo, f$hi, digits, script, description)
}

ent_p <- function(p, script, description) {
  if (!is.finite(p)) return(.ent(NA, "not estimable", script, description))
  .ent(p, fmt_p(p), script, description)
}

ent_int <- function(n, script, description) .ent(as.integer(n), fmt_int(n), script, description)

ent_pct <- function(x, script, description, digits = 1) {
  if (!is.finite(x)) return(.ent(NA, "not estimable", script, description))
  .ent(x, fmt_pct(x, digits), script, description)
}

ent_text <- function(x, script, description) .ent(x, x, script, description)

# Weighted mean and SD display, e.g. "2,812 (SD 580)".
wmean_sd <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  x <- x[ok]; w <- w[ok]
  m <- sum(w * x) / sum(w)
  s <- sqrt(sum(w * (x - m)^2) / sum(w))
  c(mean = m, sd = s)
}

ent_mean_sd <- function(x, w, digits, script, description) {
  ms <- wmean_sd(x, w)
  big <- if (digits == 0) "," else ""
  .ent(unname(ms["mean"]), sprintf("%s (SD %s)", fmt_num(ms["mean"], digits, big.mark = big),
                                   fmt_num(ms["sd"], digits, big.mark = big)), script,
       sprintf("%s; weighted SD=%s; n=%d", description, signif(ms["sd"], 6), sum(is.finite(x))))
}

wpct <- function(x, w) {
  ok <- !is.na(x) & is.finite(w)
  sum(w[ok] * as.numeric(x[ok])) / sum(w[ok])
}
