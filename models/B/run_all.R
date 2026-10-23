# Paper B: full Task 7 analysis, one command from <repo>/models:
#   Rscript B/run_all.R              # fresh run (clears models/B/_cache first; mediation ~3-4 h)
#   Rscript B/run_all.R --use-cache  # reuse cached step results (e.g. a finished background
#                                    # mediation run: Rscript B/06_mediation.R)
# Needs data/external/nfhs_subset/analysis_child.parquet (Task 5) and conception_shift.parquet;
# both are built by `python run_all.py --skip-download` (step nfhs; conception_shift.parquet is
# rebuilt whenever missing or stale, ruling P14 M7), or the latter alone by
# `python models/B/py/conception_shift.py`. load_b() and the provenance check stop with that
# command if conception_shift.parquet is missing or was built from other analysis bytes.
# Writes papers/B/results.json (key set = papers/B/results_keys.md), the aggregate tables in
# papers/B/output/, the figure data in papers/B/figures/data/ and papers/B/output/provenance.csv
# (analysis parquet sha256/md5/rows, build commit, code commit and dirty flag; every cached step
# must carry the current data md5).
# Seed: SEED_B (B/R/prep.R) for the global RNG and the mediation bootstrap.
# Env B_THREADS: fixest threads for the main models (default 0 = fixest default).

source("common/R/load.R")
source("B/R/fe.R")
source("B/R/prep.R")
source("B/R/io.R")
source("B/R/display.R")
for (f in c("01_descriptive.R", "02_birthweight.R", "03_growth.R", "04_morbidity.R", "05_enso.R",
            "06_mediation.R", "07_sensitivity.R", "08_results.R", "09_supplementary.R", "10_figure_data.R",
            "11_covariates.R", "12_precision.R"))
  source(file.path("B", f))

args <- commandArgs(trailingOnly = TRUE)
use_cache <- "--use-cache" %in% args
if (!use_cache) {
  old <- list.files(B_CACHE, full.names = TRUE)
  if (length(old)) unlink(old, recursive = TRUE)
  log_msg("run_all: cache cleared (", length(old), " files)")
}
th <- as.integer(Sys.getenv("B_THREADS", "0"))
if (th > 0L) fixest::setFixest_nthreads(th)
set.seed(SEED_B)
t_all <- Sys.time()

steps <- list()
steps$descriptive <- cached_step("descriptive", run_descriptive, use_cache)
steps$birthweight <- cached_step("birthweight", run_birthweight, use_cache)
steps$growth <- cached_step("growth", run_growth, use_cache)
steps$morbidity <- cached_step("morbidity", run_morbidity, use_cache)
steps$enso <- cached_step("enso", run_enso, use_cache)
steps$mediation <- cached_step("mediation", function() run_mediation(resume = use_cache), use_cache)
primary_rules <- c(steps$birthweight$rule, steps$growth$rule, steps$enso$rule)
steps$sensitivity <- cached_step("sensitivity", function() run_sensitivity(primary_rules), use_cache)
# Phase 3b (results_keys.md amendment 2026-10-23): deferred supplementary models, then the figure
# aggregates built from the step tables above (papers/B/figures/data/). figure_data is NOT cached
# (ruling Q2, M-4): it is cheap and must always reflect the current step tables.
steps$supplementary <- cached_step("supplementary", run_supplementary, use_cache)
steps$covariates <- cached_step("covariates", run_covariates, use_cache)     # ruling Q3 (Table S2)
steps$precision <- cached_step("precision", run_precision, use_cache)        # ruling Q5
check_precision(steps)                            # every run: refits reproduce the published estimates
t_fig <- Sys.time()
steps$figure_data <- run_figure_data(steps)
steps$figure_data$runtime_min <- as.numeric(difftime(Sys.time(), t_fig, units = "mins"))
steps$figure_data$data_md5 <- b_fingerprint()$md5
unlink(cache_path("figure_data"))                 # a cache from before this rule is never read

# ---- provenance (ruling P14, M5/M6): every step must come from the current analysis parquet,
# and conception_shift.parquet must have been built from it.
fp <- b_fingerprint()
meta <- check_shift_fingerprint(fp)
step_md5 <- vapply(steps, function(s) if (is.null(s$data_md5)) NA_character_ else s$data_md5, "")
if (anyNA(step_md5) || any(step_md5 != fp$md5))
  stop("run_all: cached step(s) ", paste(names(step_md5)[is.na(step_md5) | step_md5 != fp$md5], collapse = ", "),
       " were computed on a different (or unrecorded) analysis_child.parquet; run without --use-cache", call. = FALSE)
git1 <- function(...) tryCatch(system2("git", c("-C", REPO_ROOT, ...), stdout = TRUE, stderr = FALSE),
                               error = function(e) NA_character_)
# Code commit (ruling P15, I-1): the last commit touching the code and contract that produce the
# tables, so a re-run on a clean tree reproduces provenance.csv byte for byte (a later commit of
# the regenerated outputs alone does not change it). "dirty" = uncommitted changes in those paths.
CODE_PATHS <- c("models/B", "models/common/R", "papers/B/results_keys.md")
code_commit <- git1("log", "-1", "--format=%h", "--", CODE_PATHS)[1]
code_dirty <- length(git1("status", "--porcelain", "--", CODE_PATHS)) > 0L
log_msg("provenance: analysis_child.parquet mtime", format(fp$mtime, "%Y-%m-%d %H:%M:%S %z"),
        "(log only); code commit", code_commit, if (code_dirty) "(DIRTY tree)" else "(clean)")
prov <- data.table::data.table(
  item = c("analysis_parquet", "analysis_sha256", "analysis_md5", "analysis_bytes", "analysis_rows",
           "analysis_build_commit", "conception_shift_source_match", "steps_data_match",
           "results_code_commit", "results_code_dirty", "mediation_B", "seed"),
  value = c("data/external/nfhs_subset/analysis_child.parquet (git-ignored, controlled)", meta[["b_source_sha256"]],
            fp$md5, format(fp$size, scientific = FALSE), as.character(fp$rows),
            git1("log", "-1", "--format=%h", "--", "pipeline")[1],
            "TRUE", "TRUE", code_commit, as.character(code_dirty), as.character(MED_B), as.character(SEED_B)),
  note = c("Task 5 analysis dataset",
           "sha256 recorded by models/B/py/conception_shift.py from the same bytes (md5 and size verified here)",
           "tools::md5sum at results assembly", "file size", "parquet row count",
           "last commit touching pipeline/ (the code that builds the parquet; a rebuild from it was byte-identical, 2026-10-23)",
           "conception_shift.parquet b_source_md5/size equal the analysis parquet's",
           paste("every step cache records the same data md5:", paste(names(step_md5), collapse = "/")),
           paste("last commit touching", paste(CODE_PATHS, collapse = ", ")),
           "TRUE if those paths had uncommitted changes when the tables were assembled",
           "bootstrap resamples per mediation call", "SEED_B"))
write_b_results(steps, extra_tables = list(provenance = prov))
rt <- vapply(steps, function(s) if (is.null(s$runtime_min)) NA_real_ else s$runtime_min, 0)
log_msg("run_all: step runtimes (min):", paste(names(rt), round(rt, 1), sep = "=", collapse = ", "))
log_msg("run_all: total", round(as.numeric(difftime(Sys.time(), t_all, units = "mins")), 1), "min")
