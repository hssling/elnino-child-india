# Paper B Task 7 step 8: papers/B/results.json (key set must equal papers/B/results_keys.md)
# and the aggregate CSVs in papers/B/output/. The forbidden-name list and the count-column check
# (check_aggregate_table) live in B/R/prep.R (ruling P14, M9).

# ---- generated text keys (results_keys.md section 9, Phase 3b amendment): fixed templates
# filled from the verdict keys and the total-ENSO HAZ estimate and CI; no causal wording.

VERDICTS <- c("consistent with", "not consistent with")
.TXT_LEAD <- "In Indian children surveyed in NFHS-4 and NFHS-5, "

.check_verdict <- function(v, what) {
  if (length(v) != 1L || is.na(v) || !v %in% VERDICTS)
    stop("text keys: ", what, " is not a verdict string ('", paste(v, collapse = " "), "')", call. = FALSE)
  invisible(v)
}

# Sign clause of the total-ENSO HAZ estimate (ruling Q2, M-1): "opposite in sign" / "in the
# direction" only when the 95% CI excludes 0; neutral wording otherwise.
sign_clause <- function(x, lo = NA_real_, hi = NA_real_) {
  ok <- function(v) !is.null(v) && length(v) == 1L && !is.na(v) && is.finite(v)
  if (!ok(x)) return("not estimable")
  if (x == 0) return("zero")
  if (!ok(lo) || !ok(hi)) return("not estimable")
  if (x > 0) { if (lo > 0) "positive (95% CI excluding 0), opposite in sign to the hypothesis"
               else "positive, with a 95% CI including 0" }
  else { if (hi < 0) "negative (95% CI excluding 0), in the direction of the hypothesis"
         else "negative, with a 95% CI including 0" }
}

# lo / hi of an ent_ci entry, from its description ("...; lo=<6 sig. digits> hi=<...>").
ci_from_description <- function(desc) {
  m <- regmatches(desc, regexec("lo=([^ ;]+) hi=([^ ;]+)", desc))[[1]]
  if (length(m) != 3L) return(c(lo = NA_real_, hi = NA_real_))
  c(lo = suppressWarnings(as.numeric(m[2])), hi = suppressWarnings(as.numeric(m[3])))
}

H4_ORDINAL <- c(T2 = "second trimester", T3 = "third trimester", `T2 and T3` = "second and third trimesters")

conclusion_sentence <- function(h4, h4_terms, h5, xs, enso_haz, enso_lo = NA_real_, enso_hi = NA_real_) {
  .check_verdict(h4, "B.h4.verdict"); .check_verdict(h5, "B.h5.verdict"); .check_verdict(xs, "B.ensoxs.verdict")
  ok <- c(h4 = h4 == VERDICTS[1], h5 = h5 == VERDICTS[1], xs = xs == VERDICTS[1])
  if (ok[["h4"]] && !h4_terms %in% names(H4_ORDINAL))
    stop("text keys: B.h4.verdict is 'consistent with' but B.h4.verdict_terms is '", h4_terms, "'", call. = FALSE)
  tail <- paste0("the total El Ni\u00f1o association with height-for-age at 24\u201359 months was ",
                 sign_clause(enso_haz, enso_lo, enso_hi), ".")
  if (!any(ok))
    return(paste0(.TXT_LEAD, "drier conditions in pregnancy and early childhood and El Ni\u00f1o exposure in pregnancy ",
                  "did not show the pre-specified associations with lower birth weight or lower height-for-age; ", tail))
  parts <- c(h4 = sprintf("drier conditions in the %s of pregnancy were associated with lower birth weight",
                          if (ok[["h4"]]) H4_ORDINAL[[h4_terms]] else ""),
             h5 = "drier conditions over the first 1,000 days were associated with lower height-for-age at 24\u201359 months",
             xs = "El Ni\u00f1o exposure in pregnancy was associated with lower birth weight or height-for-age where El Ni\u00f1o dries the district")
  labels <- c(h4 = "the birth-weight hypothesis (H4)", h5 = "the height-for-age hypothesis (H5)",
              xs = "the El Ni\u00f1o \u00d7 district-sensitivity hypothesis")
  s <- paste0(.TXT_LEAD, paste(parts[ok], collapse = "; "), ", consistent with the pre-specified ",
              if (sum(ok) == 1L) "hypothesis" else "hypotheses", "; ")
  if (any(!ok)) s <- paste0(s, paste(labels[!ok], collapse = " and "), " did not meet the pre-specified rule; ")
  paste0(s, tail)
}

verdict_summary <- function(h4, h4_terms, h5, xs) {
  .check_verdict(h4, "B.h4.verdict"); .check_verdict(h5, "B.h5.verdict"); .check_verdict(xs, "B.ensoxs.verdict")
  sprintf(paste0("By the pre-specified rules, the birth-weight results were %s H4 (%s), ",
                 "the height-for-age results at 24\u201359 months were %s H5, and the results %s the pre-specified ",
                 "El Ni\u00f1o \u00d7 district-sensitivity hypothesis."), h4, verdict_terms_prose(h4_terms), h5,
          if (xs == VERDICTS[1]) "were consistent with" else "did not support")
}

build_text_keys <- function(entries) {
  S <- "models/B/08_results.R"
  need <- c("B.h4.verdict", "B.h4.verdict_terms", "B.h5.verdict", "B.ensoxs.verdict", "B.enso.oni.haz2459.elnino")
  absent <- setdiff(need, names(entries))
  if (length(absent)) stop("build_text_keys: missing key(s) ", paste(absent, collapse = ", "), call. = FALSE)
  v <- function(k) entries[[k]]$value
  ci <- ci_from_description(entries[["B.enso.oni.haz2459.elnino"]]$description)
  cs <- conclusion_sentence(v("B.h4.verdict"), v("B.h4.verdict_terms"), v("B.h5.verdict"), v("B.ensoxs.verdict"),
                            v("B.enso.oni.haz2459.elnino"), ci[["lo"]], ci[["hi"]])
  vs <- verdict_summary(v("B.h4.verdict"), v("B.h4.verdict_terms"), v("B.h5.verdict"), v("B.ensoxs.verdict"))
  src <- "generated by template from B.h4.verdict, B.h4.verdict_terms, B.h5.verdict, B.ensoxs.verdict"
  list(B.abstract.conclusion_sentence = res_entry(cs, cs, S, paste("Abstract conclusion sentence;", src,
                                                                   "and the estimate and 95% CI of B.enso.oni.haz2459.elnino")),
       B.discussion.verdict_summary = res_entry(vs, vs, S, paste("Discussion verdict summary;", src)))
}

write_b_results <- function(steps, keys_md = file.path(B_PAPER, "results_keys.md"),
                            json = file.path(B_PAPER, "results.json"), extra_tables = list(),
                            fig_dir = file.path(B_PAPER, "figures", "data")) {
  entries <- list()
  for (s in steps) {
    dup <- intersect(names(entries), names(s$entries))
    if (length(dup)) stop("08_results: duplicate key(s): ", paste(dup, collapse = ", "), call. = FALSE)
    entries <- c(entries, s$entries)
  }
  txt <- build_text_keys(entries)                 # Phase 3b section 9: from the verdict keys
  dup <- intersect(names(entries), names(txt))
  if (length(dup)) stop("08_results: duplicate key(s): ", paste(dup, collapse = ", "), call. = FALSE)
  entries <- c(entries, txt)
  entries <- apply_b_display(entries)             # ruling Q9: display strings only (B/R/display.R)
  want <- contract_keys(keys_md)
  missing <- setdiff(want, names(entries)); extra <- setdiff(names(entries), want)
  if (length(missing) || length(extra))
    stop("08_results: key set differs from the contract. Missing: ", paste(head(missing, 40), collapse = ", "),
         if (length(missing) > 40) " ..." else "", ". Extra: ", paste(head(extra, 40), collapse = ", "), call. = FALSE)
  write_results(entries, json)
  log_msg("results.json:", length(entries), "keys (contract:", length(want), ")")

  # aggregate tables: same-named tables from several steps are stacked
  tabs <- list()
  for (s in steps) for (nm in names(s$tables)) tabs[[nm]] <- data.table::rbindlist(list(tabs[[nm]], s$tables[[nm]]), fill = TRUE)
  for (nm in names(extra_tables)) tabs[[nm]] <- data.table::as.data.table(extra_tables[[nm]])
  for (nm in names(tabs)) {
    dt <- tabs[[nm]]
    check_aggregate_table(nm, dt)
    data.table::fwrite(dt, file.path(B_OUT, paste0(nm, ".csv")))
  }
  log_msg("output CSVs:", paste(sort(names(tabs)), collapse = ", "))

  # figure data (Phase 3b): same aggregate checks, written to papers/B/figures/data/
  figs <- list()
  for (s in steps) for (nm in names(s$fig_tables)) {
    if (!is.null(figs[[nm]])) stop("08_results: duplicate figure table ", nm, call. = FALSE)
    figs[[nm]] <- data.table::as.data.table(s$fig_tables[[nm]])
  }
  if (length(figs)) dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(figs)) {
    check_aggregate_table(nm, figs[[nm]])
    data.table::fwrite(figs[[nm]], file.path(fig_dir, paste0(nm, ".csv")))
  }
  if (length(figs)) log_msg("figure data CSVs:", paste(sort(names(figs)), collapse = ", "))
  invisible(list(n_keys = length(entries), tables = names(tabs), figures = names(figs)))
}
