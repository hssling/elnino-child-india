# Paper B — NFHS Child Health × Local Climate and ENSO: Data and Analysis Plan

> **For agentic workers:** REQUIRED SUB-SKILL: use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task by task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** turn NFHS-4 and NFHS-5 child (KR) records into a child-level analysis dataset linked to district-month climate and ENSO exposures. Then fit the pre-specified Paper B models and write `papers/B/results.json`.

**Architecture:**
- **Python, `pipeline/nfhs/`:**
  - streams only the needed columns from the zipped Stata files;
  - cleans outcomes with fixed unit rules;
  - maps NFHS district codes to GADM districts;
  - builds exposure windows from the Phase 1 district-month climate and ENSO tables.
- **Storage:** the child-level dataset is controlled data. It stays git-ignored under `data/external/nfhs_subset/`, and only aggregate tables are ever committed.
- **R, `models/B/`, using `fixest`:** survey-weighted fixed-effects models with two-way clustered standard errors, mediation and sensitivity analyses, written through `models/common/R/results.R`.
- **Later plan:** manuscript drafting, journal selection and the build come in a separate Phase 3b plan, as for Papers A and C.

**Tech Stack:** Python 3 (pandas, numpy, pytest); R 4.5.3 with renv (`models/`), adding `fixest` from CRAN and reusing `models/common`.

**Spec:** `docs/superpowers/specs/2026-10-22-paper-b-nfhs-child-health-design.md` (approved 2026-10-22; §12–§13 are binding, as amended by spec §14).

## Amendment 2026-10-23 (rulings P4–P15)

The plan was executed with these controller rulings, which supersede the original text where they conflict. The stale lines below have been corrected in place. Spec §14 has the full deviation list.

- **P4: cluster-GPS linkage is primary.**
  - The author confirmed DHS GPS dataset access.
  - Task 4b was added. GE files are opened read-only, streamed from the zip, and used only in Task 4b.
  - Each cluster point is linked to its 0.25° rain cell and 1° temperature cell, and assigned to a GADM district by point in polygon.
  - The district linkage of Tasks 3–4 (`dist_*` exposures, `region_id_xw`) is kept as a pre-specified sensitivity analysis.
  - Coordinates and cluster-level exposures stay git-ignored.
- **P8:**
  - the primary district FE and clusters use the point-in-polygon `region_id`;
  - the Punjab NFHS-4 codes 39 (Sangrur) and 53 (SBS Nagar) are swapped on the DHS side, and an explicit, tested override applies to the district-linkage sensitivity only.
- **P9:**
  - the pre-fit QC mask of a rain cell is the union over the dominant district and the own districts of every cluster using that cell;
  - the flags `qc_cell_masked` and `mis_district_any` are carried.
- **P10: ENSO mediation uses an ENSO-year block bootstrap.**
  - The bootstrap is `mediate_fe(..., boot_unit = "enso_year")`, where the ENSO year is the June–May year of birth month − 5.
  - There is no `birth_ym` FE; the FE are region_id + birth month of year + survey, with ns(birth_cmc, 4).
  - The number of blocks and an a-path Wald CI clustered by ENSO year are reported.
  - Local-climate mediation, had it been run, would use region_id.
- **P11:**
  - mother's age `v012` was added to the extract and enters the covariates as 5-year bands;
  - the Task 7 sensitivity inputs are carried in the analysis dataset.
- **P12 / P12b:** the controlled-data guard fails closed and checks derived identifiers, coordinates and binary formats.
- **P13: dry orientation and verdict rules.**
  - Every verdict model uses dry-oriented terms: `T*_dry = −T*_spi3`, `dry1000 = −spi3_1000d` and `sens_dry = −sens_beta`.
  - H4, H5 and ENSO × sensitivity are "consistent with" when the dry-term coefficient is negative and its 95% CI upper bound is below 0.
  - HAZ verdicts use the 24–59-month group.
  - The 1,000-day exposure is month-weighted: (3,3,3,12,12)/33 for 24–59 months; (3,3,3,12)/21 for 0–23 months, where complete.
  - Mediation is `preg_oni` → T2 or T3 SPI-3 → birth weight, with a `preg_oni_pos` sensitivity version.
  - The conception shift is built by `models/B/py/conception_shift.py`.
- **P14:**
  - non-PD vcov fixes are flagged per fit;
  - `conception_shift.parquet` carries the analysis parquet's fingerprint;
  - `run_all.py` rebuilds that file when it is missing or stale;
  - `papers/B/output/provenance.csv` is added.
- **P15:** the plan, spec and README were brought into line with these rulings. Spec §14 records the deviations and deferrals; deferrals go to Phase 3b.

## Global Constraints

- **NFHS sources.** Zips in `data/external/nfhs/IAKR74DT.zip` and `data/external/nfhs/IAKR7EDT.zip`.
  - These folders are READ-ONLY. Never extract whole files.
  - D: has about 3.8 GB free. The total of `data/external/nfhs_subset/` must stay under 150 MB.
- **GPS (GE) files** (corrected 2026-10-23, P4). The author confirmed DHS GPS dataset access.
  - GE files are opened read-only, streamed from the zip, in Task 4b only.
  - Cluster-GPS linkage is primary; district-level linkage is a sensitivity analysis.
  - Coordinates and cluster-level tables never enter git.
- **Controlled data.** No NFHS record-level or cluster-level data, and no cluster IDs (`v001`), in any committed or public file.
  - Committed outputs are aggregates only: state or district counts and means, and model coefficients.
  - Tests use synthetic records only, built in the test.
- **Units, all verified in the spec §12:**
  - `m19` is in **grams**; 9996 and 9998 → missing; values > 6000 are excluded.
  - `hw56` is **g/dL × 10**: any anaemia = raw < 110; moderate or severe = raw < 100; codes 994–999 → missing.
  - `hw70/71/72` are z × 100; ≥ 9996 → missing; WHO flags are excluded.
  - `h11`: 2 = yes, 0 = no, 8 → missing.
- **Missing is not zero.** QC-masked climate is NaN. Use `pipeline.climate_annual.apply_qc_to_month(df, "district")`.
- **Primary model:**
  - local-climate exposures, with district fixed effects plus birth year-month fixed effects;
  - survey-weighted (`v005`/1e6, rescaled to each round's sum);
  - standard errors two-way clustered by district and birth year-month;
  - estimated with `fixest`.
- **Decision rules (spec §13), fixed before results:**
  - **H4:** consistent with if the T2 **or** T3 SPI-3 coefficient for birth weight is positive (wetter means heavier; equivalently, drought lowers weight) with a 95% CI excluding 0. Report the drought direction in the text.
  - **H5:** the same rule for HAZ against 1,000-day SPI-3 in the 24–59-month group.
  - **H6:** descriptive only.
  - **ENSO × sensitivity:** consistent with if the interaction lowers birth weight or HAZ with a 95% CI excluding 0.
  - **As implemented (P13):**
    - one dry-positive orientation for every verdict: a negative coefficient on `T2_dry`/`T3_dry`, `dry1000` (HAZ 24–59) or `preg_elnino_dev × sens_dry`, with 95% CI upper bound < 0;
    - `sens_dry = −sens_beta`;
    - HAZ verdicts use 24–59 months.
- **Commits:** reuse existing modules; stage by explicit path; trailer exactly `Co-Authored-By: Claude Opus 5.5 <[e-mail removed]>`. Never touch `papers/0`, `papers/A`, `papers/C`, `release/`.

## Review Focus

1. **Stata value labels versus codes.** `read_stata` with `convert_categoricals=True` would turn codes into label strings and break the numeric rules. Read with `convert_categoricals=False`, and fetch labels separately only for districts. Test in Task 1.
2. **A child whose birth month falls before the climate series, or whose window runs past interview or after 2025.** The window exposure is NaN, never 0 or a partial mean. Test in Task 4.
3. **District names that differ between NFHS rounds and GADM** (spelling, split districts, new districts). These must be matched or listed in the review file, never silently dropped. Test in Task 3.
4. **The SE clustering and weights in fixest.** Two-way clustering and weights must actually be applied; a forgotten argument silently gives iid SEs. Test in Task 6.
5. **Leakage of controlled data into git:**
   - a committed file containing `v001`, `caseid`, or more than 20 rows of child-level data fails a check;
   - the analysis parquet is git-ignored.

   Test in Task 5.

---

### Task 1: Stream and harmonise KR columns from both rounds

**Files:**
- Create: `pipeline/nfhs/__init__.py`, `pipeline/nfhs/extract.py`, `tests/nfhs/test_extract.py`
- Modify: `.gitignore` (add `data/external/nfhs_subset/`), `pyproject.toml` (marker `nfhs`: needs local controlled NFHS data; deselected by default)

**Interfaces:**
- Produces:
  - `KR_SOURCES = {"nfhs4": (zip_path, "IAKR74FL.DTA"), "nfhs5": (zip_path, "IAKR7EFL.DTA")}`
  - `HARMONISED`: a dict mapping each output name to the source variable per round, covering `v001, v005, v008, v016, v024, district (sdistri|sdist), v025, v104, v106, v130, v190, s116, b0, b1, b2, b3, b4, b5, bord, bidx, hw1, hw53, hw56, hw57, hw70, hw71, hw72, h11, m19, m19a, sphase` (NFHS-5 only; NaN for NFHS-4) and `caseid`.
  - `read_kr(round_id: str, columns: list[str], zip_path=None) -> pd.DataFrame`. It streams the zip member and keeps only the listed columns, with `convert_categoricals=False`.
  - `district_labels(round_id: str, zip_path=None) -> dict[int, str]`, the value labels of the district variable.
  - `extract_all(out_dir=EXTERNAL/"nfhs_subset") -> Path`. It writes `kr_harmonised.parquet` with a `round` column ("nfhs4"/"nfhs5") and returns the path.

- [ ] **Step 1: Write failing tests**
  - Build a synthetic DTA in a tmp zip (`pandas.DataFrame.to_stata` with value labels, zipped) containing all the needed columns.
  - `test_read_kr_returns_numeric_codes`: `h11` stays numeric 2/0/8, not label strings.
  - `test_read_kr_selects_columns_only`: the returned columns equal the requested set.
  - `test_district_labels_roundtrip`: returns `{1: "Nicobars", 2: "North & Middle Andaman"}` for the synthetic labels.
  - `test_extract_all_harmonises_district_name_and_round`: `sdistri`/`sdist` both become `district`; `sphase` is NaN for nfhs4.
- [ ] **Step 2:** Run `python -m pytest tests/nfhs/test_extract.py -q`. Expected: ImportError.
- [ ] **Step 3: Implement.** Stream with `zipfile.ZipFile(...).open(member)` passed to `pd.read_stata(..., columns=..., convert_categoricals=False)`. If pandas cannot seek the stream, fall back to `pd.read_stata(..., iterator=True, chunksize=50000)` while filtering columns. Never extract to disk.
- [ ] **Step 4:** Run the tests (PASS). Then the real-data smoke test, marked `nfhs`: `python -m pytest -m nfhs tests/nfhs/test_extract.py -q`. It checks NFHS-4 rows = 259,627 and NFHS-5 rows = 232,920, and that the parquet is under 150 MB.
- [ ] **Step 5: Commit** the code and tests. The parquet is git-ignored and is never added.

### Task 2: Outcome cleaning and the report-estimate check

**Files:**
- Create: `pipeline/nfhs/outcomes.py`, `pipeline/reference/nfhs_report_estimates.csv`, `tests/nfhs/test_outcomes.py`

**Interfaces:**
- Consumes: `kr_harmonised.parquet` (Task 1).
- Produces:
  - `clean_outcomes(kr: pd.DataFrame) -> pd.DataFrame`, which adds the following columns:
    - `bw_g`, plus `bw_card` (bool, `m19a` = card);
    - `lbw`;
    - `haz`, `waz`, `whz` (z, float);
    - `stunted`, `underweight`, `wasted`;
    - `hb_adj` (g/dL, float);
    - `anaemia_any`, `anaemia_modsev`;
    - `diarrhoea`;
    - `age_m`, `alive`, `twin`;
    - `birth_cmc` (= `b3`), `interview_cmc` (= `v008`);
    - `weight` (= `v005`/1e6, rescaled so each round's weights sum to that round's row count).
  - `national_estimates(df) -> pd.DataFrame`: weighted national prevalence per round, computed with the same denominators the NFHS reports use.
  - `check_against_report(est, report_csv) -> list[str]`: returns the failures beyond ±1.0 percentage point.

- [ ] **Step 1: Write failing tests**, using synthetic rows:
  - `m19` 3000 → 3000 g; 9996 → NaN; 9998 → NaN; 7000 → NaN (> 6000).
  - `hw56` 105 → `hb_adj` 10.5, `anaemia_any` = 1 and `anaemia_modsev` = 0; 95 → modsev 1; 997 → NaN.
  - `hw70` −250 → `haz` −2.5, `stunted` = 1; 9998 → NaN.
  - `h11` 2 → 1, 0 → 0, 8 → NaN.
  - Twins (`b0` = 1) are flagged.
  - Weights rescale to each round's row count.
  - `check_against_report` flags a 1.5-point gap and passes a 0.5-point gap.
- [ ] **Step 2:** Run (FAIL). Then **Step 3: implement.** WHO flags: z-score outside [−6, 6] for HAZ, [−6, 5] for WAZ, [−5, 5] for WHZ → NaN.
- [ ] **Step 4: Report estimates.**
  - Fill `nfhs_report_estimates.csv` with the columns `round, indicator, denominator_note, value_pct, source, page`. The values come from the NFHS-4 national report (`data/external/nfhs/FR339.pdf`, read-only) and the NFHS-5 national report (public PDF from the DHS Program or IIPS; download it to `data/raw/nfhs/` with a MANIFEST row).
  - Indicators: stunting, underweight and wasting (< 5); anaemia any (6–59 months); diarrhoea in the last 2 weeks; LBW among children weighed.
  - Use the report tables whose denominators match KR, not the PR-based factsheet values. Record any KR-vs-report denominator difference in `denominator_note`.
- [ ] **Step 5:** The real-data test, marked `nfhs`: every indicator is within ±1.0 point of the report. A failure stops the plan for review; do not change the tolerance.
- [ ] **Step 6: Commit** the code, tests and the estimates CSV, which contains published aggregates only.

### Task 3: NFHS district → GADM district crosswalk

**Files:**
- Create: `pipeline/nfhs/districts.py`, `pipeline/reference/nfhs_district_crosswalk.csv`, `pipeline/reference/nfhs_district_review.csv`, `tests/nfhs/test_districts.py`

**Interfaces:**
- Consumes:
  - `district_labels(round_id)` (Task 1);
  - the GADM district names and `region_id`/`state_id` from `data/derived/climate_district_year.csv`;
  - `pipeline.outcomes.states` for state normalisation.
- Produces:
  - `build_crosswalk() -> pd.DataFrame` with columns `round, nfhs_code, nfhs_name, state_code, region_id, match_type` (exact | normalised | fuzzy | manual | split | unmatched).
  - `crosswalk(round_id) -> dict[int, str]`, mapping each code to `region_id`.

**Rules:**
- Match only within the same state.
- Normalise case, punctuation, "&"/"and", and the old/new spellings already in `pipeline/reference/state_frame.csv` aliases.
- Fuzzy matches with ratio ≥ 0.9 are accepted but logged.
- New districts carved after 2011 map to their parent GADM district, with `match_type = split`.
- Anything else goes to `nfhs_district_review.csv`. The task is done when there are 0 unmatched rows that have children. Manual rows give a reason.

- [ ] **Step 1: Failing tests:**
  - an exact match;
  - "Bangalore" → "Bengaluru Urban" through an alias;
  - a split district → parent;
  - a cross-state name collision is NOT matched ("Aurangabad" in Bihar vs Maharashtra);
  - an unmatched name goes to review.
- [ ] **Step 2:** FAIL. **Step 3:** implement. **Step 4:** Build both crosswalks on the real labels (marked `nfhs`). Gate: every district code that occurs in KR is mapped.
- [ ] **Step 5: Commit.** The crosswalk CSVs hold only district names and codes, no records.

### Task 4: District climate and ENSO exposure windows

**Files:**
- Create: `pipeline/nfhs/exposure.py`, `tests/nfhs/test_exposure.py`

**Interfaces:**
- Consumes:
  - `data/derived/climate_district_month.csv.gz`, masked with `apply_qc_to_month(df, "district")`;
  - `pipeline.climate_annual.spi(monthly, scale, base)`;
  - `data/derived/enso_monthly.csv` (ONI, RONI, Niño-3.4 by central month);
  - `enso_annual.csv` (`el_nino_developing`).
- Produces:
  - `district_monthly_metrics() -> pd.DataFrame`, keyed `region_id, year, month`, with:
    - `spi3`: per district, on the masked monthly rain, base 1981–2010;
    - `tmax_anom`: against the district's 1981–2010 calendar-month mean;
    - `heatwave_days`, `tmax_ge40_days`.
  - `window_exposure(children: pd.DataFrame, metrics, windows: dict[str, tuple[int, int]]) -> pd.DataFrame`.
    - Windows are month offsets relative to the birth month: `T1 = (-9, -7)`, `T2 = (-6, -4)`, `T3 = (-3, -1)`, `inf0_11 = (0, 11)`, `inf12_23 = (12, 23)`, `placebo = (1, 9)`, plus interview-relative `recent1 = (-1, -1)`, `recent3 = (-3, -1)` and `recent12 = (-12, -1)`.
    - Output columns: `<window>_spi3` (the mean of the monthly SPI-3 over the window), `<window>_tmax_anom`, `<window>_hw_days` (sum).
    - **A window is NaN if any month is masked or missing, or if it extends past the interview month.**
  - `enso_exposure(children, enso_monthly, enso_annual) -> pd.DataFrame` with:
    - `preg_oni` (mean ONI over T1–T3), `preg_oni_pos` / `preg_oni_neg` (pmax/pmin);
    - `preg_roni`, `preg_nino34`;
    - `preg_elnino_dev` (any pregnancy month in JJAS of an El Niño-developing year);
    - trimester ONI means.
  - `district_enso_sensitivity(metrics, enso_annual) -> pd.DataFrame` (`region_id, sens_beta, sens_se`). It is the per-district OLS slope of the JJAS rainfall percent departure on `oni_jjas`, over 1981–2009 only, so it uses no data from the birth window.

- [ ] **Step 1: Failing tests:**
  - a birth in Mar 2015 gives T3 = Dec 2014–Feb 2015;
  - one masked month inside T2 makes `T2_spi3` NaN;
  - an `inf12_23` window past the interview is NaN;
  - `preg_elnino_dev` is True for a pregnancy overlapping Jul 2015 and False for one entirely in 2017;
  - the sensitivity on a synthetic district with rain = −10 × ONI recovers −10 ± 1e-9;
  - the sensitivity uses no year ≥ 2010.
- [ ] **Step 2:** FAIL. **Step 3:** implement, with vectorised month-offset joins. **Step 4:** PASS. Then a real-data check (`nfhs`): T2/T3 exposure is non-missing for ≥ 90% of children outside A&N and Lakshadweep.
- [ ] **Step 5: Commit.**

### Task 4b: Cluster-GPS exposure linkage (primary; added 2026-10-23, ruling P4)

The author confirmed DHS GPS approval, so cluster-GPS linkage is now the **primary** exposure linkage. District-level exposures from Tasks 3–4 stay as a sensitivity analysis.

**Files:**
- Create: `pipeline/nfhs/gps.py`, `tests/nfhs/test_gps.py`.

**Interfaces:**
- Consumes:
  - the GE zips (shapefiles, read-only), read with geopandas or pyogrio straight from the zip without extracting it;
  - the Phase 1 0.25° rain grids (`data/derived/grids/rain_monthly_<year>.npz`), with the grid definition taken from the Phase 1 build code;
  - the raw daily 1° temperature `.grd` files and the Phase 1 heatwave logic;
  - the GADM districts (the `gadm41_IND.gpkg` used in Phase 1);
  - the QC masks;
  - `pipeline.climate_annual.spi`.
- Produces:
  - `cluster_table() -> pd.DataFrame`, with columns `round, v001, lat, lon, urban, source (GPS|GAZ|MIS), region_id`. It is saved git-ignored as `data/external/nfhs_subset/clusters.parquet`.
    - `region_id` is assigned by point-in-polygon against GADM.
    - MIS clusters at (0,0) get `region_id` from the GE district name or the KR district crosswalk, with no point.
  - `cluster_cell_monthly() -> pd.DataFrame`, keyed `round, v001, year, month`, with:
    - `spi3`: per cell, on the masked monthly rain, base 1981–2010, using the cell that contains the point;
    - `tmax_anom`: against the nearest 1° cell's 1981–2010 calendar-month normal;
    - `heatwave_days` and `tmax_ge40_days`: cell-level, computed with the Phase 1 heatwave code.

    Masking rules:
    - a cluster in a district with a reviewed QC mask for that month → NaN;
    - A&N and Lakshadweep rain → NaN;
    - MIS clusters fall back to district-month values, flagged.
  - The Task 4 `window_exposure` logic is reused with the cluster key `["round", "v001"]`. Generalise its join key rather than duplicating the code.
  - **Sensitivity:** a 10 km buffer mean (area-weighted over the 0.25° cells) for rural clusters.

**Steps:**
- [ ] **Step 1: Write failing tests on synthetic inputs:**
  - a point at a known lon/lat maps to the expected 0.25° and 1° cells;
  - a point on a cell edge is handled deterministically;
  - a (0,0) MIS cluster gets no point and takes its district's values, with `source=MIS`;
  - point-in-polygon assigns a synthetic square district;
  - a masked district-month gives a NaN cluster exposure;
  - the 10 km buffer mean equals the area-weighted mean on a synthetic grid.
- [ ] **Step 2: Run the tests.** Expected: FAIL.
- [ ] **Step 3: Implement.** Never write coordinates to any committed file.
- [ ] **Step 4: Real-data checks** (`nfhs` marker; report aggregates only):
  - counts per round of GPS, GAZ and MIS clusters;
  - share of clusters whose point-in-polygon GADM district equals the Task 3 crosswalk district (expected ≥ 95%, because displacement stays within the district);
  - share of children with non-missing T2/T3 SPI-3 (≥ 90% outside A&N and Lakshadweep).
- [ ] **Step 5: Commit** code and tests only. Nothing at cluster level is committed.

The Task 5 guard (`test_no_controlled_data.py`) must also fail on any tracked file that has `lat`/`lon` columns together with `v001`.

### Task 5: Analysis dataset and the controlled-data guard

**Files:**
- Create: `pipeline/nfhs/assemble.py`, `tests/nfhs/test_assemble.py`, `tests/test_no_controlled_data.py`
- Modify: `run_all.py` (step `nfhs` that runs only if the NFHS zips exist; not a download step), `data/README.md` (Paper B data section: sources, rules, controlled-data note)

**Interfaces:**
- Consumes: Tasks 1–4.
- Produces:
  - `build(out_dir=EXTERNAL/"nfhs_subset") -> dict[str, Path]`. It writes:
    - `analysis_child.parquet` (git-ignored), with the `region_id` crosswalk, outcomes, exposures, covariates, `birth_ym` (int YYYYMM) and the weights;
    - committed aggregate `data/derived/nfhs/child_counts_by_state_round.csv` (n and weighted prevalences by state and round);
    - committed `data/derived/nfhs/district_enso_sensitivity.csv`, which is climate-only.
  - `tests/test_no_controlled_data.py` scans every git-tracked file and fails if any of these is present: a column named `v001`, `caseid`, `hw1` or `b3`, or a CSV with more than 20 rows that has a `bw_g` or `haz` column.

- [ ] **Steps:**
  1. Write the failing tests, covering synthetic assembly, the guard on a planted file, and `birth_ym` construction from CMC (CMC 1385 → 201505).
  2. Implement.
  3. Real build (`nfhs`). Gates: under 150 MB; row counts equal Task 1 minus documented exclusions.
  4. Run `python run_all.py --skip-download`; all gates must PASS.
  5. Commit.

### Task 6: R fixest layer — weighted two-way-clustered FE models and fixest mediation

**Files:**
- Create: `models/B/R/fe.R`, `models/common/tests/testthat/test-fe.R`. These are shared helpers, sourced by `models/B`; keep them in `models/B/R` and source them explicitly.
- Modify: `models/DESCRIPTION` (add `fixest`) and `models/renv.lock` (`renv::install("fixest")` then `renv::snapshot()`).

**Interfaces:**
- Produces:
  - `fit_fe(formula, data, weights = ~weight, cluster = ~region_id + birth_ym, family = c("gaussian", "logit")) -> fixest`.
  - `fe_table(fit, terms) -> data.table(term, est, lo, hi, p, n)`. For logit, `est`/`lo`/`hi` are odds ratios.
  - `mediate_fe(data, outcome, exposure, mediator, fe, covars, B = 1000, cores = 6, seed = 1) -> list(a, b, c_total, indirect, prop, ci, flag, prop_flag)`. It uses the same flags and prop rules as `models/common/R/mediation.R` and a cluster bootstrap by `region_id`, with a PSOCK cluster seeded `seed + i - 1`.
  - **Corrected 2026-10-23 (P10):** `mediate_fe` takes `boot_unit` (default `region_id`). ENSO and other time-varying exposures must use a time-block unit, `enso_year`, and no `birth_ym` FE. `fe_table(fit, terms, cluster)` checks that the fit's vcov is the vcov clustered by `cluster`.

- [ ] **Step 1: Failing tests on simulated data:**
  - Linear: a 40-district × 60-month panel with a district effect, a month effect and a known β = −25. `fit_fe` recovers β within 3 SE.
  - The two-way-clustered SE differs from the iid SE, and equals `fixest::summary(fit, cluster = ~region_id + birth_ym)`.
  - Supplying weights changes the estimate when the weights correlate with the treatment effect.
  - Logit `fe_table` returns ORs.
  - `mediate_fe` recovers a × b within tolerance on simulated data with B = 40.
- [ ] **Step 2:** FAIL. **Step 3:** implement. **Step 4:** PASS. Run the full R suite with `SKIP_SLOW_R` allowed.
- [ ] **Step 5: Commit** the R code, tests, DESCRIPTION and renv.lock.

### Task 7: Paper B analyses → `papers/B/results.json`

**Files:**
- Create:
  - `models/B/run_all.R`, `models/B/01_descriptive.R`, `models/B/02_birthweight.R`, `models/B/03_growth.R`, `models/B/04_morbidity.R`, `models/B/05_enso.R`, `models/B/06_mediation.R`, `models/B/07_sensitivity.R`, `models/B/08_results.R`;
  - `papers/B/results_keys.md`, `papers/B/results.json`;
  - `papers/B/output/*.csv`, which are aggregate tables only.

**Interfaces:**
- Consumes: `analysis_child.parquet` (read with `arrow`; add to renv if needed), plus the Task 6 helpers and `models/common/R/results.R`.

**Pre-specified models.** Every model uses weights and two-way clustering unless stated otherwise. **Corrected 2026-10-23 (P10, P13):** `papers/B/results_keys.md`, the contract committed before any model ran, is the binding specification. Items 1, 2, 4 and 6 below are as implemented.

1. **H4, birth weight.**
   - Formula: `bw_g ~ T1_dry + T2_dry + T3_dry + T1_tmax_anom + T2_tmax_anom + T3_tmax_anom + sex + bord + mother covariates | region_id + birth_ym`, with `T*_dry = −T*_spi3` (P13).
   - Also LBW with logit.
2. **H5, HAZ and stunting.** Run by age group (0–23 months; 24–59 months).
   - The exposure is the 1,000-day SPI-3, as the month-weighted mean (P13):
     - 24–59 months: (3·T1 + 3·T2 + 3·T3 + 12·inf0_11 + 12·inf12_23)/33;
     - 0–23 months: (3·T1 + 3·T2 + 3·T3 + 12·inf0_11)/21, where complete.

     It enters dry-oriented as `dry1000`.
   - Fixed effects: `region_id + birth_ym`, plus age in months.
3. **H6, recent windows.**
   - Diarrhoea on `recent1`; anaemia on `recent12`.
   - Fixed effects: `region_id + interview year-month`.
   - Descriptive only.
4. **ENSO × sensitivity.**
   - Formula: `outcome ~ preg_elnino_dev:sens_dry + covariates | region_id + birth_ym`, with `sens_dry = −sens_beta` (P13). The main effects are absorbed: `preg_elnino_dev` by `birth_ym` and `sens_dry` by `region_id`.
   - Outcomes: birth weight and HAZ. The verdict uses birth weight, or HAZ in the 24–59-month group only.
5. **Total ENSO association.**
   - Formula: `outcome ~ preg_oni_pos + preg_oni_neg + covariates + spline(birth_cmc) | region_id + birth_month_of_year + survey`.
   - RONI and Niño-3.4 as alternatives.
6. **Mediation (secondary).** ENSO → T2/T3 SPI-3 → birth weight, via `mediate_fe`, with B = 1000 and cores = 6.
   - **Corrected (P10, P13):**
     - exposure `preg_oni`, with a `preg_oni_pos` + `preg_oni_neg` sensitivity version;
     - one call per mediator;
     - `boot_unit = "enso_year"`;
     - FE region_id + birth month of year + survey, with ns(birth_cmc, 4) columns and no `birth_ym` FE;
     - a time-clustered Wald CI for the a path.
7. **Sensitivity analyses:**
   - card-only birth weight; last births only;
   - unweighted;
   - non-movers;
   - excluding `v008` ≥ 1450;
   - mother fixed effects (mothers with ≥ 2 births);
   - placebo window;
   - conception shift ±1 month.

**Results output:**
- `08_results.R` writes `papers/B/results.json` covering descriptive Ns, every primary estimate (`fmt_ci`), the H4 and H5 verdicts by rule, the ENSO-interaction verdict, mediation, and a sensitivity summary.
- `results_keys.md` documents every key: meaning, model and display.

- [ ] **Steps:**
  1. Write `results_keys.md` first, as the contract.
  2. Implement the scripts.
  3. Run `Rscript B/run_all.R`; the mediation step runs in the background with ≤ 6 workers.
  4. Check that the key set equals the contract.
  5. Commit the scripts, contract, results and aggregate outputs. The guard in Task 5 must still pass.

---

## Self-review notes

- **Spec coverage:**
  - Data and units (§3, §5, §12) → Tasks 1–2.
  - District linkage (§13.1) → Task 3.
  - Windows, SPI-3, ENSO and sensitivity (§4, §12, §13.2) → Task 4.
  - Controlled-data rules (§8) → Task 5.
  - Weighting and clustering (§13.4) and fixest (§6) → Task 6.
  - H4–H6, age groups (§13.3), the ENSO pathway, mediation and sensitivity analyses (§6, §13) → Task 7.
- **Deferred to the Phase 3b plan:** the manuscript, figures, journal choice, and the public code and aggregate repository.
  - **Corrected 2026-10-23:** the cluster-GPS upgrade was not deferred. It became primary under P4 and was built in Task 4b.
  - **Also deferred (P15, spec §14):**
    - trimester-specific ONI terms, as a pre-specified supplementary model reported whatever the result;
    - the `hw57` anaemia cross-check, descriptive;
    - WAZ/WHZ, underweight, wasting and moderate/severe anaemia, as supplementary tables using the same model families.
  - **Deviations recorded in spec §14:**
    - rainfall % departure is not a window exposure;
    - the location rule is the point cell plus a 10 km rural buffer as a sensitivity analysis;
    - interview-month FE:
      - none in the birth-weight models, by design;
      - near-absorbed by birth year-month plus `hw1` age FE in the HAZ, stunting and anaemia models (exact in NFHS-4; in NFHS-5 about 45% of children differ by one month);
      - H6 uses interview year-month FE.
