# Paper B display layer (ruling Q9; results_keys.md "Amendment 2026-10-23 (ruling Q9)").
# Applied by write_b_results at assembly: it rewrites `display` strings only (p values to 1-2
# significant digits, machine codes to prose) and never a `value`, so the cached steps (including
# the mediation bootstrap) are reused unchanged. models/common/R/results.R is not touched
# (Papers A and C use it).

# "<0.001" below 0.001; ">0.99" when p rounds to 1 at 2 significant digits; otherwise 2
# significant digits (0.27, 0.049, 0.0021, 0.50). NA stays NA.
fmt_p_sig <- function(p) {
  if (any(!is.na(p) & (p < 0 | p > 1))) stop("fmt_p_sig: p values must lie in [0, 1]", call. = FALSE)
  s <- signif(p, 2)
  out <- trimws(formatC(s, digits = 2, format = "fg", flag = "#"))
  out[!is.na(p) & s >= 1] <- ">0.99"
  out[!is.na(p) & p < 0.001] <- "<0.001"
  out[is.na(p)] <- NA_character_
  out
}

VERDICT_TERMS_PROSE <- c(
  neither = "neither the second- nor the third-trimester term met the rule",
  T2 = "the second-trimester term met the rule",
  T3 = "the third-trimester term met the rule",
  `T2 and T3` = "both the second- and the third-trimester terms met the rule")

MED_FLAG_PROSE <- c(ok = "bootstrap completed (at most 10% of resamples failed)",
                    too_many_failed_fits = "bootstrap not usable (more than 10% of resamples failed)")

PROP_FLAG_PROSE <- c(ok = "estimated",
                     total_effect_ci_includes_0 = "not estimated (total-association CI includes 0)",
                     opposite_signs = "not estimated (indirect and total associations of opposite sign)",
                     nonfinite_resample_proportion = "not estimated (non-finite proportion in a bootstrap resample)",
                     not_assessed = "not assessed (too many failed bootstrap fits)")

SENS_SPEC_PROSE <- c(card = "card-recorded birth weight", lastbirth = "last births", unweighted = "unweighted",
                     nonmover = "non-movers", precovid = "pre-COVID interviews", motherfe = "mother fixed effects",
                     shiftm1 = "conception shifted 1 month earlier", shiftp1 = "conception shifted 1 month later",
                     district = "district linkage", buf10 = "10 km rural buffer", noflag = "flagged exposures excluded",
                     nogaz = "gazetteer clusters excluded", cluster1 = "district-only clustering")

SENS_MODEL_PROSE <- c(h4bw = "trimester dryness, birth weight", h5haz = "1,000-day dryness, HAZ 24\u201359 months",
                      xsbw = "El Ni\u00f1o \u00d7 sensitivity, birth weight",
                      xshaz = "El Ni\u00f1o \u00d7 sensitivity, HAZ 24\u201359 months")

# Look a code up in a prose table; an unknown code is an error, never passed through silently.
prose_of <- function(table, code, what) {
  if (length(code) != 1L || is.na(code) || !code %in% names(table))
    stop("display: no prose for ", what, " code '", paste(code, collapse = " "), "'", call. = FALSE)
  enc2utf8(unname(table[[code]]))
}

verdict_terms_prose <- function(code) prose_of(VERDICT_TERMS_PROSE, code, "B.h4.verdict_terms")

# "nonmover/xsbw, card/h4bw" -> "non-movers, El Nino x sensitivity, birth weight; card-recorded ..."
differ_prose <- function(code) {
  if (identical(code, "none")) return("none")
  pairs <- trimws(strsplit(code, ",", fixed = TRUE)[[1]])
  enc2utf8(paste(vapply(pairs, function(p) {
    sm <- strsplit(p, "/", fixed = TRUE)[[1]]
    if (length(sm) != 2L) stop("display: B.sens.differ pair '", p, "' is not spec/model", call. = FALSE)
    paste0(prose_of(SENS_SPEC_PROSE, sm[1], "sensitivity spec"), ", ", prose_of(SENS_MODEL_PROSE, sm[2], "sensitivity model"))
  }, character(1)), collapse = "; "))
}

apply_b_display <- function(entries) {
  for (k in names(entries)) {
    e <- entries[[k]]
    v <- e$value
    if (endsWith(k, ".p")) {
      if (is.numeric(v) && length(v) == 1L && is.finite(v)) e$display <- fmt_p_sig(v)
    } else if (grepl("^B\\.med\\..+\\.flag$", k)) {
      e$display <- prose_of(MED_FLAG_PROSE, v, k)
    } else if (grepl("^B\\.med\\..+\\.prop_flag$", k)) {
      e$display <- prose_of(PROP_FLAG_PROSE, v, k)
    } else if (grepl("^B\\.med\\..+\\.prop$", k) && grepl("^not estimated \\((.+)\\)$", e$display)) {
      code <- sub("^not estimated \\((.+)\\)$", "\\1", e$display)
      if (code %in% names(PROP_FLAG_PROSE)) e$display <- prose_of(PROP_FLAG_PROSE, code, k)
    } else if (k == "B.sens.differ") {
      e$display <- differ_prose(v)
    } else if (k == "B.h4.verdict_terms") {
      e$display <- verdict_terms_prose(v)
    }
    entries[[k]] <- e
  }
  entries
}
