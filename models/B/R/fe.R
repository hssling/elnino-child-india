# Paper B fixest layer: survey-weighted fixed-effects models with two-way clustered SEs, and
# a fixest mediation with a cluster (unit) bootstrap.
#
# Source after common/R/load.R (working directory <repo>/models):
#   source("common/R/load.R"); source("B/R/fe.R")
# It reuses .with_seed, .pct, .na_ci and .prop_ci from common/R.
#
# fit_fe(formula, data, weights = ~weight, cluster = ~region_id + birth_ym, family)
#   formula is a fixest formula, e.g. bw ~ spi3_t2 + covars | region_id + birth_ym.
#   gaussian -> fixest::feols; logit -> fixest::feglm(binomial("logit")). The weights column
#   (survey weight, already rescaled by the caller) must be finite and > 0; weights = NULL fits
#   unweighted. The clustered vcov is requested at estimation and verified afterwards by a
#   semantic check, vcov(fit) == vcov(fit, cluster = cluster) (version-independent; a forgotten
#   argument would otherwise silently give iid or one-way SEs); a weighted call must also have
#   stored its weights. fixest's default fixef.rm = "perfect_fit" applies: singleton FE levels
#   (and, for logit, FE levels with an all-0 or all-1 outcome) are dropped and n is post-drop.
#
# fe_table(fit, terms, cluster = ~region_id + birth_ym)
#   -> data.table(term, est, lo, hi, p, n, vcov, df). Estimates, 95% CIs and p from the fit's
#   own vcov via fixest's coeftable/confint. Unless cluster = NULL (an explicit opt-out for a
#   deliberately different specification), the fit's vcov must equal the vcov clustered by
#   `cluster` (same semantic check), else an error. vcov is fixest's label (e.g.
#   "Clustered (region_id & birth_ym)"); df is the reference distribution's degrees of freedom.
#   Reference distribution (fixest defaults, kept deliberately): gaussian (feols) uses
#   t with df = min(G) - 1 over the cluster dimensions (ssc G.df = "min", t.df = "min");
#   logit (feglm) uses the normal (z), reported as df = Inf. At NFHS sizes min(G) >= 36, so
#   the difference is negligible. For a logit fit est/lo/hi are odds ratios (exp of the
#   coefficient and of the CI bounds).
#
# mediate_fe(data, outcome, exposure, mediator, fe, covars, B, cores, seed, weights, fitter,
#            boot_unit = "region_id")
#   Linear product-of-coefficients mediation with fixed effects `fe` (character vector of
#   column names) and the same survey weights in all three fits:
#     a-path: mediator ~ exposure + covars | fe
#     b-path: outcome  ~ exposure + mediator + covars | fe   (b; c_direct = exposure coef)
#     total:  outcome  ~ exposure + covars | fe              (c_total)
#   indirect = a * b, prop = indirect / c_total (for OLS on the same rows and weights
#   c_total = c_direct + indirect exactly).
#   Scope: linear only - a binary outcome or mediator is fitted as a linear probability model
#   and the indirect effect is then on the risk-difference scale; one mediator per call (T2
#   and T3 SPI-3 need separate calls); `covars` must be existing column names, so spline or
#   other transformed terms (e.g. ns(birth_cmc, 3)) must be precomputed as columns.
#   All three models use the same rows: complete on the union of their variables plus
#   boot_unit and the weight (dropped counts per first missing variable in `dropped`). The
#   fits use fixef.rm = "none" (fe_ols), so singleton FE levels are kept (they do not change
#   OLS coefficients) and every model has exactly the complete-row n.
#   Bootstrap unit (boot_unit, any column): whole units are resampled with replacement; draw k
#   of unit u becomes the distinct unit "u_k", and every FE in `fe` that is equal to or nested
#   in the unit (each FE level within one unit) is relabelled the same way ("<level>_k"), so
#   each copy has its own FE levels; FEs that cross units are left as they are.
#   CHOICE OF UNIT (binding for Task 7):
#     - local-climate exposures (vary by district and month): boot_unit = "region_id".
#     - ENSO / any time-varying exposure common to all districts (e.g. ONI by birth month):
#       the district bootstrap is anti-conservative (it treats the months as fixed and the
#       a-path SE almost vanishes). Task 7 MUST use a time-block unit - the ENSO year
#       (June-May) - as boot_unit, and must NOT include birth_ym in `fe` (ENSO is collinear
#       with the birth year-month FE), e.g. fe = c("region_id", "birth_month") or region_id
#       plus a trend column in covars.
#   Resample i uses seed + i - 1 (applied locally; the caller's RNG stream is restored).
#   Flags and proportion rules are those of common/R/mediation.R: a resample fails when a fit
#   errors or leaves a needed coefficient unidentified; if more than 10% of the B resamples
#   fail, every CI is NA, flag = "too_many_failed_fits", prop_flag = "not_assessed"; otherwise
#   CIs are 2.5/97.5% percentiles over successful resamples and prop_flag is set by .prop_ci
#   ("ok", "total_effect_ci_includes_0", "opposite_signs", "nonfinite_resample_proportion").
#   cores > 1: PSOCK cluster (also on Windows); each worker gets the master's library paths,
#   setwd(MODELS_ROOT), sources common/R/load.R and B/R/fe.R, single-threaded fixest and
#   data.table; only the columns the models use are shipped; results are identical to the
#   serial run. The fitter is shipped by value, so a closure is refused with cores > 1.
#   The bootstrap fits use fe_ols (iid vcov: only the coefficients are used).

if (!exists(".prop_ci", mode = "function") || !exists(".with_seed", mode = "function"))
  stop("B/R/fe.R: source common/R/load.R first", call. = FALSE)

.fe_weight_var <- function(weights) {
  if (is.null(weights)) return(NULL)
  if (!inherits(weights, "formula") || length(all.vars(weights)) != 1L)
    stop("fit_fe: weights must be a one-sided formula naming one column (e.g. ~weight) or NULL", call. = FALSE)
  all.vars(weights)
}

.fe_check_weights <- function(weights, d, who = "fit_fe") {
  wv <- .fe_weight_var(weights)
  if (is.null(wv)) return(invisible(NULL))
  if (!wv %in% names(d)) stop(who, ": weight column '", wv, "' not found", call. = FALSE)
  w <- d[[wv]]
  if (!is.numeric(w) || any(!is.finite(w)) || any(w <= 0))
    stop(who, ": weight column '", wv, "' must be finite and > 0 (no NA)", call. = FALSE)
  invisible(wv)
}

# fixest's label for the fit's vcov, e.g. "Clustered (region_id & birth_ym)" or "IID".
.fe_vcov_label <- function(fit) {
  type <- tryCatch(attr(stats::vcov(fit), "vcov_type"), error = function(e) NULL)
  if (is.null(type)) type <- attr(fit$cov.scaled, "vcov_type")
  if (is.null(type)) "" else paste(type, collapse = " ")
}

# Semantic, version-independent check: the fit's vcov equals the vcov clustered by `cluster`.
.fe_check_vcov <- function(fit, cluster, who) {
  if (!inherits(cluster, "formula"))
    stop(who, ": cluster must be a one-sided formula, e.g. ~region_id + birth_ym", call. = FALSE)
  same <- tryCatch(isTRUE(all.equal(stats::vcov(fit), stats::vcov(fit, cluster = cluster),
                                    check.attributes = FALSE)),
                   error = function(e) FALSE)
  if (!same)
    stop(who, ": the fit's vcov ('", .fe_vcov_label(fit), "') is not the vcov clustered by ",
         paste(all.vars(cluster), collapse = " + "), call. = FALSE)
  invisible(fit)
}

# fit_fe post-check: clustered vcov, and a weighted call must have stored its weights.
.fe_check_fit <- function(fit, cluster, weights) {
  .fe_check_vcov(fit, cluster, "fit_fe")
  if (!is.null(weights) && is.null(fit$weights))
    stop("fit_fe: the weights were not applied to the fit", call. = FALSE)
  invisible(fit)
}

fit_fe <- function(formula, data, weights = ~weight, cluster = ~region_id + birth_ym,
                   family = c("gaussian", "logit")) {
  family <- match.arg(family)
  if (!inherits(cluster, "formula"))
    stop("fit_fe: cluster must be a one-sided formula, e.g. ~region_id + birth_ym", call. = FALSE)
  .fe_check_weights(weights, data)
  fit <- if (family == "gaussian") {
    if (is.null(weights)) fixest::feols(formula, data = data, cluster = cluster)
    else fixest::feols(formula, data = data, weights = weights, cluster = cluster)
  } else {
    if (is.null(weights)) fixest::feglm(formula, data = data, family = stats::binomial("logit"), cluster = cluster)
    else fixest::feglm(formula, data = data, family = stats::binomial("logit"), weights = weights, cluster = cluster)
  }
  .fe_check_fit(fit, cluster, weights)
}

.fe_is_logit <- function(fit) {
  !is.null(fit$family) && is.list(fit$family) && identical(fit$family$link, "logit")
}

fe_table <- function(fit, terms, cluster = ~region_id + birth_ym) {
  if (!is.null(cluster)) .fe_check_vcov(fit, cluster, "fe_table")
  ct <- fixest::coeftable(fit)
  absent <- setdiff(terms, rownames(ct))
  if (length(absent)) stop("fe_table: term(s) not in the fit: ", paste(absent, collapse = ", "), call. = FALSE)
  ci <- stats::confint(fit, level = 0.95)
  est <- unname(ct[terms, 1]); lo <- unname(ci[terms, 1]); hi <- unname(ci[terms, 2])
  if (.fe_is_logit(fit)) { est <- exp(est); lo <- exp(lo); hi <- exp(hi) }
  df <- if (identical(colnames(ct)[4], "Pr(>|t|)")) as.numeric(fixest::degrees_freedom(fit, "t")) else Inf
  data.table::data.table(term = terms, est = est, lo = lo, hi = hi, p = unname(ct[terms, 4]),
                         n = as.integer(stats::nobs(fit)), vcov = .fe_vcov_label(fit), df = df)
}

# ---- mediation ----

# Default mediation fitter: weighted feols, iid vcov (only coefficients are used), no FE-based
# row removal (fixef.rm = "none": singletons kept, n = complete rows).
fe_ols <- function(formula, data, weights) {
  if (is.null(weights)) fixest::feols(formula, data = data, vcov = "iid", fixef.rm = "none", notes = FALSE)
  else fixest::feols(formula, data = data, weights = weights, vcov = "iid", fixef.rm = "none", notes = FALSE)
}

.fe_formula <- function(lhs, rhs, fe) {
  f <- paste(lhs, "~", paste(rhs, collapse = " + "))
  if (length(fe)) f <- paste(f, "|", paste(fe, collapse = " + "))
  stats::as.formula(f, env = globalenv())
}

# Rows complete on vars; dropped = count of rows by the first missing variable (in vars order).
.fe_complete_rows <- function(d, vars) {
  absent <- setdiff(vars, names(d))
  if (length(absent)) stop("mediate_fe: missing column(s): ", paste(absent, collapse = ", "), call. = FALSE)
  miss <- rep(NA_character_, nrow(d))
  for (v in rev(vars)) {
    x <- d[[v]]
    bad <- if (is.numeric(x)) !is.finite(x) else is.na(x)
    miss[bad] <- v
  }
  .fe_keep <- is.na(miss)
  dropped <- if (any(!.fe_keep)) data.table::data.table(variable = miss[!.fe_keep])[, .(n = .N), by = variable]
             else data.table::data.table(variable = character(), n = integer())
  list(data = d[.fe_keep], dropped = dropped)
}

# FEs in `fe` equal to or nested in `unit` (every FE level within a single unit).
.fe_nested_fe <- function(d, unit, fe) {
  u <- as.character(d[[unit]])
  Filter(function(f) f == unit || all(tapply(u, as.character(d[[f]]), function(v) length(unique(v))) == 1L), fe)
}

# Resampled copy of d: units of `unit` drawn with replacement; draw k of unit u relabelled
# "u_k", and each column in `relabel` (the unit and FEs nested in it) gets "<level>_k".
.fe_bootstrap_units <- function(d, unit = "region_id", seed = NULL, relabel = unit) {
  d <- data.table::as.data.table(d)
  units <- unique(as.character(d[[unit]]))
  pick <- .with_seed(seed, sample(units, length(units), replace = TRUE))
  rows <- split(seq_len(nrow(d)), as.character(d[[unit]]))[pick]
  # dot-prefixed local: data.table evaluates i in the scope of d's columns first
  .fe_ix <- unlist(rows, use.names = FALSE)
  draw <- rep(seq_along(pick), lengths(rows))
  out <- d[.fe_ix]
  for (col in unique(c(unit, relabel))) {
    lab <- paste0(as.character(out[[col]]), "_", draw)
    data.table::set(out, j = col, value = factor(lab, levels = unique(lab)))
  }
  out[]
}

.fe_fit_or_stop <- function(formula, d, weights, fitter, what) {
  m <- tryCatch(fitter(formula, d, weights), error = function(e) e)
  if (inherits(m, "error")) stop("mediate_fe: ", what, " fit failed: ", conditionMessage(m), call. = FALSE)
  m
}

.fe_coef_or_stop <- function(model, term, what) {
  v <- stats::coef(model)[term]
  if (is.na(v)) stop("mediate_fe: ", what, " coefficient '", term, "' is unidentified (NA)", call. = FALSE)
  unname(v)
}

.fe_med_fits <- function(d, outcome, exposure, mediator, covars, fe, weights, fitter) {
  fa <- .fe_fit_or_stop(.fe_formula(mediator, c(exposure, covars), fe), d, weights, fitter, "a-path")
  fb <- .fe_fit_or_stop(.fe_formula(outcome, c(exposure, mediator, covars), fe), d, weights, fitter, "b-path")
  ft <- .fe_fit_or_stop(.fe_formula(outcome, c(exposure, covars), fe), d, weights, fitter, "total-effect")
  list(a = .fe_coef_or_stop(fa, exposure, "a-path"),
       b = .fe_coef_or_stop(fb, mediator, "b-path"),
       c_direct = .fe_coef_or_stop(fb, exposure, "b-path"),
       c_total = .fe_coef_or_stop(ft, exposure, "total-effect"),
       n_by_model = c(a = stats::nobs(fa), b = stats::nobs(fb), total = stats::nobs(ft)))
}

# One bootstrap resample (top level so parLapply ships it to PSOCK workers).
.fe_boot_one <- function(i, d, seed, outcome, exposure, mediator, covars, fe, weights, fitter,
                         boot_unit, relabel) {
  di <- .fe_bootstrap_units(d, unit = boot_unit, seed = seed + i - 1L, relabel = relabel)
  r <- tryCatch(.fe_med_fits(di, outcome, exposure, mediator, covars, fe, weights, fitter),
                error = function(e) e)
  if (inherits(r, "error"))
    return(list(i = as.integer(i), a = NA_real_, b = NA_real_, c_direct = NA_real_, c_total = NA_real_,
                error = conditionMessage(r)))
  list(i = as.integer(i), a = r$a, b = r$b, c_direct = r$c_direct, c_total = r$c_total, error = NA_character_)
}

.fe_worker_init <- function(root, libs) {
  .libPaths(libs)
  setwd(root)
  source("common/R/load.R")
  source("B/R/fe.R")
  fixest::setFixest_nthreads(1)
  data.table::setDTthreads(1)
  invisible(TRUE)
}

mediate_fe <- function(data, outcome, exposure, mediator, fe, covars = character(), B = 1000,
                       cores = 6, seed = 1, weights = ~weight, fitter = fe_ols,
                       boot_unit = "region_id") {
  if (length(B) != 1L || is.na(B) || B < 1) stop("mediate_fe: B must be a positive integer", call. = FALSE)
  B <- as.integer(B)
  if (length(cores) != 1L || is.na(cores) || cores < 1) stop("mediate_fe: cores must be a positive integer", call. = FALSE)
  cores <- as.integer(cores)
  if (cores > 1L) {
    fenv <- environment(fitter)
    if (!is.null(fenv) && !identical(fenv, globalenv()) && !isNamespace(fenv))
      stop("mediate_fe: with cores > 1 the fitter must be a top-level function (global environment or a ",
           "package), not a closure: its state would not be shared with the workers. Use cores = 1.", call. = FALSE)
  }
  if (!is.character(boot_unit) || length(boot_unit) != 1L || is.na(boot_unit))
    stop("mediate_fe: boot_unit must be one column name", call. = FALSE)
  d0 <- data.table::as.data.table(data)
  if (!boot_unit %in% names(d0))
    stop("mediate_fe: boot_unit column '", boot_unit, "' not found (bootstrap unit)", call. = FALSE)
  wv <- .fe_weight_var(weights)
  vars <- unique(c(outcome, exposure, mediator, covars, fe, boot_unit, wv))
  cr <- .fe_complete_rows(d0, vars)
  d <- droplevels(cr$data[, vars, with = FALSE])        # only the used columns (also to workers)
  if (nrow(d) == 0L) stop("mediate_fe: no complete rows remain", call. = FALSE)
  .fe_check_weights(weights, d, "mediate_fe")
  relabel <- .fe_nested_fe(d, boot_unit, fe)

  est <- .fe_med_fits(d, outcome, exposure, mediator, covars, fe, weights, fitter)
  if (any(est$n_by_model != nrow(d)))
    stop("mediate_fe: a model dropped rows: n = ",
         paste(names(est$n_by_model), est$n_by_model, collapse = ", "), " against ", nrow(d),
         " complete rows (a custom fitter must keep every row, e.g. fixef.rm = 'none')", call. = FALSE)
  indirect <- est$a * est$b

  args <- list(d = d, seed = seed, outcome = outcome, exposure = exposure, mediator = mediator,
               covars = covars, fe = fe, weights = weights, fitter = fitter,
               boot_unit = boot_unit, relabel = relabel)
  if (cores > 1L) {
    cl <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, .fe_worker_init, MODELS_ROOT, .libPaths())
    res <- do.call(parallel::parLapply, c(list(cl = cl, X = seq_len(B), fun = .fe_boot_one), args))
  } else {
    res <- do.call(lapply, c(list(X = seq_len(B), FUN = .fe_boot_one), args))
  }
  boot <- data.table::rbindlist(res)
  boot[, indirect := a * b]
  boot[, prop := indirect / c_total]
  data.table::setcolorder(boot, c("i", "a", "b", "c_direct", "c_total", "indirect", "prop", "error"))
  ok <- is.na(boot$error)
  B_ok <- sum(ok); B_failed <- B - B_ok
  prop_nonfinite <- ok & !is.finite(boot$prop)

  if (B_failed > 0.1 * B) {
    flag <- "too_many_failed_fits"
    prop_flag <- "not_assessed"
    ci <- list(indirect = .na_ci(), prop = .na_ci(), a = .na_ci(), b = .na_ci(),
               c_direct = .na_ci(), c_total = .na_ci())
  } else {
    flag <- "ok"
    bo <- boot[ok]
    ci <- list(indirect = .pct(bo$indirect), prop = .na_ci(), a = .pct(bo$a), b = .pct(bo$b),
               c_direct = .pct(bo$c_direct), c_total = .pct(bo$c_total))
    pc <- .prop_ci(indirect, est$c_total, bo, ci$c_total)
    prop_flag <- pc$prop_flag
    ci$prop <- pc$ci
  }
  list(a = est$a, b = est$b, c_total = est$c_total, c_direct = est$c_direct,
       indirect = indirect, prop = indirect / est$c_total, ci = ci, flag = flag, prop_flag = prop_flag,
       n = nrow(d), n_by_model = est$n_by_model, dropped = cr$dropped,
       B = B, B_ok = as.integer(B_ok), B_failed = as.integer(B_failed),
       B_prop_nonfinite = as.integer(sum(prop_nonfinite)),
       boot = boot[], fail_messages = unique(boot$error[!ok]), cores = cores,
       boot_unit = boot_unit, relabelled = relabel)
}
