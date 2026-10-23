# El Niño, local climate and child health in India (NFHS-4 and NFHS-5): code and aggregate data

Code and aggregate results for a study of in-utero and early-life climate exposures and child health in India. The study links children under five in the National Family Health Survey (NFHS-4, 2015–16, and NFHS-5, 2019–21) to the rainfall deficit (SPI-3 dryness) and maximum-temperature (Tmax) anomaly at their survey cluster in each trimester of pregnancy and in early childhood, and to El Niño–Southern Oscillation (ENSO) conditions during pregnancy. Outcomes are birth weight and low birth weight (H4), height-for-age and stunting (H5), and, descriptively, recent diarrhoea and anaemia (H6). Models are survey-weighted fixed-effects regressions (R, fixest) with district and birth-month fixed effects.

Repository: https://github.com/hssling/elnino-child-india

**Paper (manuscript in preparation):**
Siddalingaiah H S. Rainfall deficit, heat and El Niño in pregnancy and early childhood and the birth weight and growth of children under five in India: an observational study of two national household surveys linked to gridded climate data. Manuscript in preparation.

Please cite the paper once it is published, and this repository through `CITATION.cff`.

Companion repositories (same climate and ENSO pipeline): https://github.com/hssling/elnino-malaria-india (https://doi.org/10.5281/zenodo.23424043) and https://github.com/hssling/elnino-heat-india (https://doi.org/10.5281/zenodo.23556278).

## NFHS data: controlled, not redistributed

The analysis uses the NFHS-4 and NFHS-5 children's recode (**KR**) files and the cluster GPS (**GE**) files:

| Survey | KR (children under 5) | GE (cluster GPS) |
|---|---|---|
| NFHS-4 (India DHS 2015–16) | `IAKR74DT.zip` | `IAGE71FL.zip` |
| NFHS-5 (India DHS 2019–21) | `IAKR7EDT.zip` | `IAGE7AFL.zip` |

These files come from the DHS Program (https://dhsprogram.com/) on registration and are **not redistributed** here. To reproduce the analysis, register a project with the DHS Program, request the India KR datasets and, separately, the GPS (GE) datasets (GPS access needs its own approval), and download the four zip files listed above. The GE coordinates are displaced by the DHS Program (up to 2 km for urban clusters and up to 5 km for rural clusters, 10 km for 1% of them) and must not be published or mapped at point level.

This repository holds **no record-level or cluster-level NFHS data**: no child, mother or household rows, no cluster identifiers (`v001`), household or line numbers (`v002`, `v003`), `caseid`, `mother_id` or coordinates, and no parquet, cache or log files. Every NFHS-derived file is an aggregate (state, district or national counts, shares and model coefficients); table cells based on fewer than 30 children are suppressed (left empty and flagged). Before this repository was assembled, every file was checked by the same controlled-data guard that the analysis repository runs on every commit (forbidden NFHS column names, coordinates, child-level tables of more than 20 rows, DHS recode file names, unscanned binary formats).

**Cite the surveys.** International Institute for Population Sciences (IIPS) and ICF. 2017. *National Family Health Survey (NFHS-4), 2015-16: India.* Mumbai: IIPS. International Institute for Population Sciences (IIPS) and ICF. 2021. *National Family Health Survey (NFHS-5), 2019-21: India, Volume I.* Mumbai: IIPS. The microdata and GPS datasets are produced by IIPS and ICF and distributed by the DHS Program (ICF, Rockville, Maryland, USA).

**Ethics.** The study is a secondary analysis of anonymised NFHS data obtained from the DHS Program with permission. The NFHS surveys were approved by the institutional review boards of the International Institute for Population Sciences, Mumbai, and ICF, and informed consent was obtained from all participants by the survey teams.

## Pre-specified analysis

The analysis was pre-specified. The design specification is dated 2026-10-22, and the results contract, which fixes the samples, the models, the decision rules for the three tested hypotheses, the sensitivity analyses and the primary result keys, was committed before any model was fitted. Four dated amendments added keys later (listed below): the supplementary analyses that the design had pre-specified but deferred, committed before their models were run; and, after the primary results were known, covariate and exposure-distribution tables, rescaled precision summaries of the pre-specified estimates, and display formats. The amendments fit no new primary models and carry no verdicts. The three documents are included unchanged except as noted below:

- `docs/superpowers/specs/2026-10-22-paper-b-nfhs-child-health-design.md`: design specification (hypotheses H4–H6, exposures, outcomes, models, decision rules). §14 records the amendments made during execution, each with its reason.
- `docs/superpowers/plans/2026-10-22-paper-b-nfhs-analysis.md`: data and analysis plan.
- `papers/B/results_keys.md`: the results contract. Every number reported in the paper is a key of `papers/B/results.json` listed in this file; keys added later are under the four dated amendments listed below.

Timeline, from the history of the author's working repository (first commit of each file; times are IST, UTC+05:30):

| Step | File | First commit (IST) |
|---|---|---|
| Design specification | `docs/superpowers/specs/2026-10-22-paper-b-nfhs-child-health-design.md` | 2026-10-22 14:29 (`62e26e8`) |
| Data and analysis plan | `docs/superpowers/plans/2026-10-22-paper-b-nfhs-analysis.md` | 2026-10-22 14:44 (`cbd8aea`) |
| Results contract (`papers/B/results_keys.md`) | `papers/B/results_keys.md` | 2026-10-23 01:24 (`f6bedff`) |
| First model results (`papers/B/results.json`) | `papers/B/results.json` | 2026-10-23 05:04 (`3685cc0`) |
| Contract amendment for the Phase 3b supplementary keys | `papers/B/results_keys.md` | 2026-10-23 09:53 (`1ae9c79`) |
| Supplementary models (`models/B/09_supplementary.R`) | `models/B/09_supplementary.R` | 2026-10-23 10:11 (`8f6785e`) |

Contract amendments of 2026-10-23 (each a section of `papers/B/results_keys.md`; commits of the working repository, IST):

| Amendment | Commit | What it added | Relative to the results |
|---|---|---|---|
| Phase 3b supplementary keys | 09:53 (`1ae9c79`) | Trimester-specific ENSO terms, anaemia-classification cross-check, secondary outcomes (pre-specified in the design, deferred), figure data, text keys built from the verdicts | Before the supplementary models were run; after the primary results |
| Covariate distribution and missingness table | 11:12 (`e9d1b22`) | Descriptive covariate table | After the primary results; descriptive, no verdicts |
| Exposure distributions and rescaled precision | 11:34 (`d30c9c0`) | Exposure percentiles; pre-specified estimates rescaled per IQR or SD of exposure and at percentiles of dry-sensitivity | After the primary results; restatements, no new models or verdicts |
| Reader-facing displays | 11:51 (`8ffdde8`) | Display formats (2-significant-digit p values, prose codes), one sample-flow count | After the primary results; display only |

The copies here differ from the committed originals only in three ways: local drive paths are replaced by repository-relative ones (for example, the folder holding the NFHS zip files becomes `data/external/nfhs/`); e-mail addresses are removed; and three journal-selection notes in the design specification (its target-journal lines) are replaced by the marked text "[journal-selection note redacted from the public copy; not part of the analysis protocol]". These notes concerned only where the paper would be submitted. The hypotheses, exposures, outcomes, models, decision rules, sensitivity analyses and amendments are verbatim. Internal decision identifiers (for example "ruling P13") refer to the author's analysis log; the decisions that affect the analysis are described in the documents themselves.

## Contents

```
./                              11 files
data/derived/                   7 files
                                # ENSO and district climate tables: inputs of the exposure linkage
data/derived/nfhs/              child_counts_by_state_round.csv, district_enso_sensitivity.csv
                                # NFHS aggregates: state-round counts, district ENSO sensitivity
docs/superpowers/plans/         2026-10-22-paper-b-nfhs-analysis.md
                                # data and analysis plan, 2026-10-22
docs/superpowers/specs/         2026-10-22-paper-b-nfhs-child-health-design.md
                                # design specification, 2026-10-22
models/                         .Rprofile, DESCRIPTION, renv.lock, setup.R
                                # renv lock file and setup, R 4.5.3
models/B/                       14 files
                                # R analysis of this study: run_all.R and steps 01-10
models/B/R/                     display.R, fe.R, io.R, prep.R
                                # helpers: fixed-effects fits, data preparation, input/output
models/B/py/                    conception_shift.py
                                # conception-shift exposure windows, Python
models/common/R/                8 files
                                # shared R functions: results formatting, mediation and others
models/common/tests/            run.R
                                # testthat suite
models/common/tests/testthat/   14 files
models/renv/                    .gitignore, activate.R, settings.json
papers/B/                       results.json, results_keys.md
                                # every reported number, and the results contract
papers/B/figures/               make_figures.py
papers/B/figures/data/          fig1_timeline.csv, fig2_district_map.csv, fig3_coefs.csv, fig4_sensitivity.csv, figS1_placebo.csv
                                # figure data
papers/B/output/                16 files
                                # aggregate tables written by models/B
pipeline/                       13 files
                                # shared climate and ENSO modules imported by pipeline/nfhs
pipeline/nfhs/                  7 files
                                # NFHS extraction, outcomes, crosswalk, cluster and district exposures, dataset
pipeline/outcomes/              __init__.py, states.py
                                # state frame, shared module
pipeline/reference/             nfhs_district_crosswalk.csv, nfhs_district_review.csv, nfhs_report_estimates.csv, qc_reviewed.csv, state_frame.csv
                                # QC masks, state frame, NFHS district crosswalk, NFHS report estimates
tests/                          __init__.py
tests/nfhs/                     8 files
                                # pytest suite for pipeline/nfhs, synthetic data only
```

No manuscript text and no figure images are included.

## Aggregate data

All tables are UTF-8 CSV; missing is an empty cell (NaN), never 0.

- `papers/B/results.json`: every number reported in the paper, as `{"<key>": {"value", "display", "script", "description"}}` (690 keys). `display` is the formatted text used in the paper; `description` records the model, the sample, the interval bounds and the seed where relevant. Verdicts (`consistent with` / `not consistent with`) are applied only by the decision rules in `results_keys.md`.
- `papers/B/output/`: the aggregate tables written by the R analysis: sample flow, descriptive tables by survey round and by El Niño exposure in pregnancy, missingness of the exposures by state and round, the primary, ENSO, sensitivity and supplementary model coefficients (estimate, 95% CI, p value, n, vcov, fixed effects, unit), the mediation results and `provenance.csv` (checksums of the analysis dataset from which the results were computed, the code commit and the seed).
- `papers/B/figures/data/`: the data behind each figure. District map values are suppressed (`suppressed` = True, value empty) where fewer than 30 children contribute.
- `data/derived/nfhs/child_counts_by_state_round.csv`: by survey round and NFHS state, the number of children and, per indicator (lbw, stunting, underweight, wasting, anaemia_any, anaemia_modsev, diarrhoea), the unweighted denominator `n` and the weighted prevalence `prev_pct`, blank where `n` < 30.
- `data/derived/nfhs/district_enso_sensitivity.csv`: per GADM district, the OLS slope of June–September rainfall percentage departure on June–September ONI, 1981–2009 (`sens_beta`, % per ONI unit; `sens_se`; `sens_n` valid years). Climate only; NaN with fewer than 20 valid years and for the island districts.
- `data/derived/`: the climate and ENSO inputs of the exposure linkage, as in the companion repositories: monthly and annual ENSO indices (`enso_monthly.csv`, `enso_annual.csv`, `enso_episodes.csv`), the GADM district-month climate table from India Meteorological Department (IMD) gridded data (`climate_district_month.csv.gz`), the district-year table (`climate_district_year.csv`) and the QC-mask tables. `pipeline/reference/` holds the reviewed QC masks, the state frame, the NFHS district crosswalk and its review, and the NFHS report estimates used to check the extraction.

## Reproduction

Requirements: Python 3.10 or later (tested with 3.14), R 4.5.3, and the four DHS zip files.

1. **Place the DHS files** (never commit them; `data/external/` is git-ignored):

   ```
   data/external/nfhs/IAKR74DT.zip   data/external/nfhs/IAGE71FL.zip
   data/external/nfhs/IAKR7EDT.zip   data/external/nfhs/IAGE7AFL.zip
   ```

2. **GADM boundaries.** Download the GADM 4.1 GeoPackage for India, `gadm41_IND.gpkg` (https://gadm.org/download_country.html; direct link https://geodata.ucdavis.edu/gadm/gadm4.1/gpkg/gadm41_IND.gpkg), and save it as `data/external/gadm41_IND.gpkg`. It is used for the point-in-polygon district assignment of the survey clusters, for the district weights and for the maps. It is not redistributed here (GADM licence: free for academic and other non-commercial use).

   For Figure 2 only, India's external outline also includes the areas claimed by India that GADM 4.1 keeps in its Pakistan and China files (disputed-area units Z06, Z02, Z03 and Z08). Download the two GADM 4.1 level-0 GeoJSON files and save them at the paths `papers/B/figures/make_figures.py` expects (`CLAIMED_FILES`):

   ```
   https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_PAK_0.json -> data/external/gadm41_extra/pak0/gadm41_PAK_0.json
   https://geodata.ucdavis.edu/gadm/gadm4.1/json/gadm41_CHN_0.json -> data/external/gadm41_extra/chn0/gadm41_CHN_0.json
   ```

   sha256 of the files used (downloaded 2026-10-23): `gadm41_PAK_0.json` b7e852a489cd2cb55f782bed8656372aec463ae68f34791159e8b693e694f49b (158,483 bytes); `gadm41_CHN_0.json` 3f82d0d1c021acfb83d457e7e79547f90620953dcb4b94a059269f106e7660ac (1,246,198 bytes). The outline so drawn is not authenticated by the Survey of India.

3. **Climate grids for the cluster-level exposures.** The cluster (GPS) exposures are computed from the IMD 0.25° daily rainfall and 1° daily temperature grids and the district grid weights, which are too large to distribute. Build them with the climate pipeline of a companion repository (for example `python run_all.py` in https://github.com/hssling/elnino-heat-india) and copy its `data/derived/grids/`, `data/derived/_parts/` and `data/raw/imd/` folders into this repository. The district-level tables in `data/derived/` are included and need no rebuild.

4. **Analysis dataset** (Python, from the repository root):

   ```
   python -m pip install -r requirements.txt
   python -m pipeline.nfhs.assemble              # KR extraction, outcomes, cluster-GPS and district exposures, gates
   python models/B/py/conception_shift.py        # conception-shift windows used by the sensitivity analyses
   ```

   The Stata files are streamed from the zips (selected columns only) and never extracted. All record-level and cluster-level outputs are written to the git-ignored `data/external/nfhs_subset/` (under 150 MB); only the aggregates in `data/derived/nfhs/` are written elsewhere. Code comments and messages that mention `run_all.py` refer to the pipeline orchestrator of the analysis repository (the `nfhs` step runs the same two commands).

5. **Models** (R):

   ```
   cd models
   Rscript setup.R          # restores the package library from renv.lock (renv)
   Rscript B/run_all.R      # all steps; clears models/B/_cache first
   ```

   `B/run_all.R` writes `papers/B/results.json` (it stops unless the key set equals `results_keys.md`), `papers/B/output/` and `papers/B/figures/data/`. The mediation bootstrap (B = 1,000, resampling ENSO years) takes about 3–4 hours on 6 workers; `Rscript B/run_all.R --use-cache` reuses finished steps. The seed is fixed (`SEED_B` in `models/B/R/prep.R`). `papers/B/output/provenance.csv` records the checksums of the analysis dataset; in the author's environment a rebuild from the same DHS files gave the same checksums (they also depend on the package versions, see `requirements.txt` and `models/renv.lock`).

**Tests.** `python -m pytest` runs the pipeline tests on synthetic data (tests that need the DHS files are marked `nfhs` and deselected by default; `-m nfhs` runs them). `cd models && Rscript common/tests/run.R` runs the R testthat suite.

## Software

- R 4.5.3; packages pinned in `models/renv.lock` (renv 1.3.1; fixest 0.14.2, arrow 25.0.1, data.table 1.18.6.1, jsonlite 2.0.0, testthat 3.3.2).
- Python 3.14 (3.10 or later required) with pandas 2.3, numpy 2.3, pyarrow 24, scipy 1.16, geopandas 1.1, shapely 2.1, pyogrio 0.12, pyproj, requests 2.32 and pytest.

## Licences

- Code: MIT (`LICENSE`).
- Aggregate data and figure data (`papers/B/results.json`, `papers/B/output/`, `papers/B/figures/data/`, `data/derived/`, `pipeline/reference/`) and any figures drawn from them: CC BY 4.0 (`DATA_LICENSE.md`), subject to the third-party terms listed there. The NFHS microdata and GPS data are not covered by these licences and are not included.
- GADM-derived identifiers, names and boundaries (the GADM district and state codes such as `region_id` and `district_gid`, the district and state names taken from GADM, and the district and national outlines drawn in Figure 2, including the disputed-area units from the GADM Pakistan and China files) remain subject to the GADM licence (free for academic and other non-commercial use; redistribution or commercial use needs GADM's permission, https://gadm.org/license.html), not to CC BY 4.0.

## Contact

Dr Siddalingaiah H S, Department of Community Medicine, Shridevi Institute of Medical Sciences and Research Hospital, Tumakuru, Karnataka, India. ORCID: https://orcid.org/0000-0002-4771-8285.
