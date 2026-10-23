# Paper B — El Niño, local climate and child health in India (NFHS-4 and NFHS-5): Design

**Date:** 2026-10-22
**Status:** approved 2026-10-22, and executed in Phase 3 (data and analyses, branch phase3-paperB) on 2026-10-23. §14 records the amendments made during execution (rulings P4–P15).
**Parent specs:**
- `docs/superpowers/specs/2026-10-14-el-nino-health-program-design.md` (program; Paper B);
- `docs/superpowers/specs/2026-10-15-phase2-papers-a-c-design.md` §6 (roadmap: Paper B adds birth weight with trimester-specific prenatal exposure).

**Built on:**
- the Phase 1 climate layer (IMD grids and district-month tables, ENSO tables);
- the Phase 2b R modelling layer (`models/common`);
- the house rules for results contracts and submission builds used for Papers A and C.

---

## 1. Purpose and why

Papers A and C examined El Niño through state-year counts. They found real but modest associations, and found that ENSO adds little predictive skill beyond local weather. Paper B moves to the **individual child**. It asks whether El Niño-related climate conditions experienced **in the womb and in early life** relate to:
- birth weight;
- later growth (stunting, wasting, underweight);
- anaemia;
- diarrhoea.

Exposure is local where possible. It is linked to each child's survey cluster and birth month, not to a state-year.

The policy hook is the first 1,000 days. If El Niño-related drought or heat during pregnancy lowers birth weight or growth, nutrition programmes have a reason to target pregnant women when an El Niño develops. The 2026–27 event is the timely case.

## 2. Research questions and pre-specified hypotheses

- **H4 (birth weight):** El Niño-developing conditions during pregnancy are associated with lower birth weight and higher odds of low birth weight (LBW, < 2.5 kg). The association runs mainly through local rainfall deficit and heat in the second and third trimesters.
- **H5 (growth):** Children whose first 1,000 days (conception to age 2) overlapped El Niño-related drought have lower height-for-age z-scores and higher odds of stunting.
- **H6 (recent morbidity, secondary):**
  - diarrhoea in the two weeks before interview is associated with local rainfall and heat in the month before interview;
  - anaemia (6–59 months) is associated with exposure over the preceding 12 months.

  These are secondary and descriptive, because seasonal confounding is strong.

**Decision rules.** Fix them before results, in the analysis plan, in the same style as Papers A and C: estimate, CI, and the direction and significance needed for "consistent with".

## 3. Data

**Survey data.** NFHS-4 (2015–16) and NFHS-5 (2019–21) DHS recode files are held locally in `data/external/nfhs/`.
- **Controlled access:** the DHS data-use agreement forbids redistribution, so they never enter git or any public repository.
- **KR (children under 5).** NFHS-4 `IAKR74DT.zip`, NFHS-5 `IAKR7EDT.zip`. Variables already confirmed present:
  - identifiers: `v001` cluster, `v005` weight, `v024` state, `sdistri`/`sdist` district;
  - dates: `v008` and `v016` interview date, `b1`/`b2`/`b3` birth month, year and CMC; `b17` day of birth (NFHS-5 only); `hw1` age;
  - child: `b4` sex;
  - outcomes: `m19`/`m19a` birth weight and its source (card or recall); `hw70`/`hw71`/`hw72` HAZ/WAZ/WHZ; `hw53`/`hw56`/`hw57` haemoglobin and anaemia; `h11` diarrhoea;
  - covariates: `v106` education, `v190` wealth, `v025` residence, `m14` antenatal visits, `m15` place of delivery.
- **GE (cluster GPS).** NFHS-4 `IAGE71FL.zip`, NFHS-5 `IAGE7AFL.zip`. Clusters are displaced by up to 2 km (urban) and 5 km, or 10 km for 1%, (rural).
- **BR (births).** NFHS-5 only (`IABR7EDT.zip`); optional for the sibling sensitivity analysis. NFHS-4 has no BR file in the folder, but KR holds multiple under-5 siblings per mother.

**Exposures** (all already built in this repo):
- **ENSO:** `data/derived/enso_monthly.csv` (monthly ONI, RONI, Niño-3.4 and DMI) and `enso_episodes.csv`.
- **Local climate:**
  - monthly 0.25° rain grids, `data/derived/grids/rain_monthly_<year>.npz` (R19-masked; git-ignored);
  - temperature exists only as raw daily 1° IMD `.grd` files in `data/raw/imd/`. Cluster-level Tmax and heatwave days need cell-level processing that reuses the Phase 1 heatwave code;
  - `data/derived/climate_district_month.csv.gz`, giving rain, rainy days, Tmax, Tmin, days ≥ 40 °C and heatwave days per GADM district-month, 1981–2025. QC masks are **not** stored in it; apply them with `pipeline.climate_annual.apply_qc_to_month(df, "district")`.

**Coverage.** Births in the five years before each survey run from about Feb 2010 to 2021. ENSO events that reach the in-utero and early-life windows, from `data/derived/enso_episodes.csv`:
- El Niño 2009–10 (strong; 2009 is a developing year and reaches NFHS-4 births from Feb 2010), 2014–16 (very strong), 2018–19 and 2019–20;
- La Niña 2010–11 (strong), 2011–12, 2017–18 and 2020–21. The 2021–23 episode starts after the last NFHS-5 interview in May 2021.

By `enso_annual.csv`, `el_nino_developing` is True only for 2009 and 2015 in this window. An El Niño-developing indicator is therefore essentially a single birth-cohort contrast. This is the binding limitation. The design leans on **local climate anomalies**, which vary across space and carry far more information, and on the program spec's **ENSO × district rainfall-sensitivity** design. ENSO terms are reported as total associations, with honest uncertainty.

## 4. Exposure construction

- **Location: cluster GPS → IMD grid.**
  - Average the grid cells within a buffer matching the displacement: 2 km urban, 5 km rural, and 10 km as a sensitivity analysis.
  - Clusters at (0, 0) or missing coordinates are mapped to their district through the district-month table and flagged.
  - Fallback / sensitivity: use the district-month table for all children. This needs an NFHS-district ↔ GADM-district crosswalk, because NFHS-5 has about 707 districts against GADM's 666.
- **Timing:**
  - Conception month is approximately birth month − 9. NFHS has no gestational age, so this is stated as a limitation, with a sensitivity analysis on a −8 or −10 shift.
  - Trimester windows: T1 = months −9 to −7, T2 = −6 to −4, T3 = −3 to −1 before the birth month.
  - Infancy windows: 0–11 and 12–23 months after birth, truncated at interview.
  - Recent windows, for diarrhoea and anaemia: the month of interview and the preceding 1, 3 and 12 months.
- **Local metrics per window:**
  - rainfall anomaly as a percentage departure from the cluster's 1981–2010 monthly normal, plus SPI-3 for the window;
  - Tmax anomaly against the same baseline;
  - heatwave days and days ≥ 40 °C.

  QC-masked cells or months are set to NaN, never 0, following the R19 and Task 0 rules.
- **ENSO metrics per window:**
  - mean ONI, with RONI and Niño-3.4 as alternatives (the Amendment A1 comparison carried over);
  - El Niño-developing exposure, meaning any trimester in a JJAS of an El Niño-developing year;
  - El Niño and La Niña entered separately as `pmax`/`pmin` terms.

## 5. Outcomes and samples

| Outcome | Sample | Definition |
|---|---|---|
| Birth weight (g), LBW | All births in the KR 5-year window with a recorded weight. KR holds all births in the window, including dead children. Sensitivity: last births only; card only. | `m19` is already in **grams**; DHS floors it at 500 g. Codes 9996 (not weighed) and 9998 (don't know) → missing. Exclude values > 6,000 g as implausible. Primary: card **and** recall (`m19a`); sensitivity: card only. Check heaping at round values. Twins (`b0` > 0) excluded. |
| HAZ, WAZ, WHZ; stunting, underweight, wasting | Children 0–59 months with a valid measurement | `hw70/71/72` ÷ 100. Exclude WHO flags (≥ 9996, \|z\| outside WHO limits). Binary outcomes use z < −2. |
| Anaemia (any; moderate or severe) | 6–59 months, tested | `hw56` altitude-adjusted Hb is stored as **g/dL × 10**: any = raw < 110; moderate or severe = raw < 100. Exclude codes 994–999, including NFHS-5's 997 ("inconsistent"). `hw57` as a cross-check. |
| Diarrhoea, last 2 weeks | Living children under 5 | `h11` = 2 (yes, last two weeks); 0 = no; 8 (don't know) → missing. Neither round uses code 1. |

**Missing is not zero.** Records with missing outcomes are excluded and counted, never filled.

## 6. Statistical models

**Primary:** individual-level fixed-effects regression with the R package **fixest**.
- Linear models for birth weight and z-scores; logit, or a linear probability model as sensitivity, for binary outcomes.
- **District fixed effects** absorb stable geography.
- **Birth-month-of-year fixed effects** absorb seasonality. **Survey fixed effects** and **interview-month fixed effects** are also included.
- Smooth terms: child age in months (spline) and a linear or spline trend in birth date.
- Child sex and birth order; mother's age, education, wealth quintile, residence, caste and religion.
- **Standard errors clustered by district.**
- Survey weights: primary analyses are unweighted with design variables included; weighted analyses are a sensitivity. This is a choice, documented in Methods.

**Model families** (estimand clarity, as in Amendment A1.4):
1. **Total ENSO association:** outcome ~ trimester ENSO terms + fixed effects + covariates, with no local climate.
2. **Local climate:** outcome ~ trimester local rainfall and heat anomalies + fixed effects + covariates. This is the main, better-identified exposure.
3. **Mediation:** ENSO → trimester local anomaly → outcome, using product-of-coefficients with a cluster bootstrap by district. Reuse the logic of `models/common/R/mediation.R`, adapted to the fixest fits; the adaptation is the builder's call.

**H6** uses recent-window exposures, with interview-month and district fixed effects. These results are reported as secondary.

**Sensitivity analyses:**
- card-only birth weight;
- the conception-month shift;
- a 10 km buffer;
- district-table exposure instead of GPS;
- non-movers only: `v104` == 95 ("always"), or `v104` ≥ the child's age in years + 1 to cover pregnancy. Visitors (96) excluded;
- mother fixed effects with siblings in KR, for birth weight and HAZ;
- excluding post-COVID-onset NFHS-5 interviews (`v008` ≥ 1450, October 2020 onward; `sphase` == 2 also contains pre-COVID Jan–Mar 2020 interviews);
- RONI and Niño-3.4 instead of ONI;
- placebo windows: ENSO in the 12 months after birth for birth weight; a future-pregnancy window.

**No out-of-sample forecasting** in Paper B. The question is association with a well-identified exposure, not prediction. Papers A and C already cover the forecasting question.

## 7. Outputs and manuscript

- **Results contract:** `papers/B/results.json` with `results_keys.md`, the same schema as Papers A and C. Display formats come from `models/common/R/results.R` (U+2212 minus, " to ", `fmt_ci`/`fmt_p`).
- **Figures:**
  - maps of exposure and outcome, at state or district level only;
  - trimester-specific coefficient plots;
  - an exposure timeline of ENSO against birth cohorts.
- **Tables:**
  - descriptive by ENSO exposure;
  - main models;
  - mediation;
  - sensitivity analyses.
- **Reporting guideline:** STROBE, with the RECORD extension for routinely collected data. The DHS survey design is described.
- **Ethics wording:** use the program spec's prescribed NFHS text ("Not applicable … anonymised … approved by the IIPS and ICF institutional review boards … informed consent obtained by the survey teams"), as in the author's prior NFHS papers.
- **Data availability:**
  - NFHS microdata are available from the DHS Program on registration; the paper does not redistribute them.
  - Public **code** plus **derived exposure tables** go into a public repository: cluster-free aggregates, with no individual or cluster records.
- [journal-selection note redacted from the public copy; not part of the analysis protocol]
- **Related-work disclosure:** the author's NFHS anaemia papers are in `<author's manuscripts folder>`. The differentiation (exposure focus on ENSO and in-utero climate) must be stated, and overlapping analyses avoided.

## 8. Constraints

- **Disk:** D: has about **3.8 GB free**. Stream from the zip files and extract only the needed columns into a git-ignored subset under `data/external/nfhs_subset/`, targeting under 150 MB. Never unzip whole IR or HR files.
- **Read-only sources:** the NFHS folders are read-only; work from copies or by streaming.
- **Privacy:**
  - no cluster coordinates, cluster IDs or individual records in any committed or public file;
  - cluster-level exposure tables stay git-ignored;
  - only state, district or aggregate summaries leave the machine.
- **Shared code:** reuse `models/common` where it fits. New shared R helpers get tests. Add fixest to renv, from CRAN.
- **House rules:** missing is not zero; seeds fixed; one command per stage (`Rscript B/run_all.R`); commits staged by explicit path, with the standard trailer; papers/0, papers/A and papers/C are never touched.

## 9. Work breakdown (the plan will detail it)

1. **NFHS extraction and harmonisation.** Python, streaming KR and GE. Covers variable mapping across the two rounds, outcome cleaning rules, and checks against the published NFHS factsheet national estimates within tolerance.
2. **Exposure linkage.** Cluster GPS buffers to the grids, the district crosswalk fallback, the trimester, infancy and recent windows, and the ENSO metrics.
3. **Models** in R with fixest: main models, mediation and sensitivity analyses.
4. **Results contract and manuscript.** Drafted in parallel with step 3, as for Papers A and C. Then journal research, formatting and the build.
5. **Public code and derived-aggregate repository,** plus a Zenodo DOI.

## 10. Left to the builder

- exact spline degrees of freedom;
- logit versus linear probability as the primary model for each binary outcome;
- buffer implementation (a cell-area-weighted mean is suggested);
- crosswalk construction;
- the factsheet-check tolerance (suggested ±1 percentage point national);
- figure styling within journal rules.

## 11. Risks

| Risk | Mitigation |
|---|---|
| Few ENSO events (2–3 El Niño) | Lead with local climate, report ENSO as total association with an honest CI, and run the RONI/Niño-3.4 comparison. |
| Recall bias and heaping in birth weight | Card-only sensitivity analysis and a heaping check. |
| Selection: only weighed births, only living children | Describe it; sensitivity analysis on all births with birth weight. |
| Gestational age unknown | Approximate conception, run the shift sensitivity, and state it. |
| GPS displacement | Buffers and the district fallback. |
| Disk space | Column-subset streaming, with a size check. |
| Overlap with the author's anaemia papers | Different exposure and question, disclosed in the cover letter. |

## 12. Corrections and the assistant's calls after the builder check (2026-10-22)

**Verified against the files (read-only):**
- KR holds **all births** in the 5-year window (`bidx` ≥ 2 is about 25% of records).
- `v104` exists in both KR files.
- NFHS-4 KR has 259,627 rows, with births from Feb 2010 to Nov 2016.
- NFHS-5 KR has 232,920 rows, with births 2014–2021 and interviews Jun 2019 to May 2021 (none Apr–Sep 2020).
- District codes differ between rounds (640 vs 707 codes; GADM has 666).
- GE files: about 12–14% of cluster locations are gazetteer (GAZ) and about 0.5% are missing (MIS). Displacement is kept within the NFHS district.

**Assistant's calls** (minor gaps; the author may override any of them):
- **Location.**
  - Each cluster point is assigned to a GADM district, giving a consistent 666 districts across rounds. District fixed effects and the district fallback both use this.
  - MIS clusters are placed through the GE district-name crosswalk.
  - GAZ clusters are included, with a sensitivity analysis that excludes them.
  - Exposure is the grid cell containing the point. Area-weighting is used only where a 5 km buffer crosses a cell edge, with a 10 km sensitivity analysis. For 1° temperature, the nearest cell is used.
- **Rain metric.**
  - SPI-3 per trimester is primary.
  - Percent departure is used only when the window's normal is ≥ 50 mm.
  - Masked cells and months are set to NaN. The Tamil Nadu 2011–12 artefact is masked through the reviewed district masks.
  - Andaman & Nicobar and Lakshadweep rain is NaN.
- **ENSO terms.**
  - Primary: a single pregnancy-mean ONI, split into El Niño and La Niña terms.
  - Secondary: trimester-specific ENSO terms, which are collinear with one another.
  - "El Niño-developing exposure" means any overlap of the pregnancy with JJAS of an El Niño-developing year.
- **Time controls.**
  - Calendar month-of-interview (1–12) fixed effects, not year-month, which avoids exact collinearity with age and birth date.
  - Age in months as fixed effects or a spline.
  - A birth-date trend.
- **Covariates.**
  - Included: child sex, birth order (`bord`), mother's age, education, wealth, residence, caste (`s116`, with a "none/missing" level) and religion (`v130`).
  - Excluded from total-effect models: `m14`/`m15`, which are recorded for the last birth only and lie on the causal path.
- **H6 diarrhoea:** dead children are excluded; recent-window exposures are used.
- **Placebo:** a pseudo-pregnancy in months +1 to +9 after birth, set to NaN if incomplete at interview.
- **Mother fixed effects:** only for mothers with at least 2 births with the outcome in the window.
- **Factsheet check.**
  - Compare with NFHS report tables built on the same KR denominators; factsheet anthropometry and anaemia figures come from PR.
  - Tolerance is ±1 percentage point, national.
  - Document any KR–PR gap.
- **Weighted analyses:** `v005`/1e6, rescaled to each round's sum of weights.
- **Data reading:** stream the zipped Stata files reading selected columns only, either with `pandas.read_stata` on the zip member stream (to verify) or with a fixed-width reader. Never extract full files, because of disk space.
- **Gestational age:** absent (`b20` is all not-applicable in NFHS-5). Conception is approximated, and Methods says so.
- **Mediation:** the shared `mediate()` relabels units for glmmTMB, so the bootstrap needs a fixest-based fitter wrapper. Budget hours on parallel workers.
- **Deviations from the program spec, stated openly.**
  - Buffer sizes: the program said 5/10 km. This design uses the cell containing the point, with a 5 km edge buffer and a 10 km sensitivity analysis.
  - [journal-selection note redacted from the public copy; not part of the analysis protocol]
  - GC covariates: not used in the primary model.

**Open questions for the author** (see the conversation): DHS registration and GPS approval; the fixed-effects design choice for the local-climate model; the ENSO × district sensitivity design; H5 age groups; weighting; the primary birth-weight sample; the target journal; decision rules.

## 13. Author decisions (2026-10-22) and resulting design

1. **DHS access.** *Superseded by §14 (ruling P4, 2026-10-23): cluster-GPS linkage is primary, and district linkage is a sensitivity analysis.* The author will check whether the DHS project registration covers this study and includes the GPS (GE) datasets.
   - Until that is confirmed, the **primary linkage is district-level**. Each child is linked to the district-month climate of the GADM district matched to their NFHS district code (the crosswalk is built in the plan) and to their birth month.
   - GE files are **not opened** until the author confirms access.
   - Cluster-GPS linkage, described in §4 and §12, becomes an **upgrade step** that runs only after confirmation, as a sensitivity analysis or a replacement for the primary analysis.
2. **Design.**
   - **Local climate is primary.** Trimester rainfall (SPI-3) and heat anomalies at the child's location, with **district fixed effects plus birth year-month fixed effects**. This is the strongest design, because it absorbs every national shock, including ENSO itself, and identifies effects from spatial variation in local anomalies.
   - **El Niño is secondary.** It is studied as (a) **ENSO × district rainfall sensitivity**: does a district's historical ENSO–rainfall sensitivity, estimated 1981–2009 before the birth window, modify the El Niño-year association with outcomes? This is identified from spatial contrasts. And (b) a total association, using pregnancy-mean ONI with district and season fixed effects, reported with honest uncertainty.
   - **Mediation** links (b) to the local-climate pathway, using a shared specification without birth year-month fixed effects. It is therefore secondary.
3. **H5 growth.** Analysed in **two age groups**:
   - **0–23 months:** in-utero exposure plus exposure since birth, truncated at interview;
   - **24–59 months:** the full first 1,000 days.
4. **Weighting.**
   - **Survey-weighted primary.** Weights are `v005`/1e6, rescaled to each round.
   - Standard errors are **two-way clustered by district and birth year-month**, which respects the shared ENSO timing.
   - Unweighted analyses are a sensitivity check.

**Assistant's calls**, still open to override:
- **Birth-weight primary sample:** all births with a recorded weight. Last births and card-only are sensitivity analyses.
- [journal-selection note redacted from the public copy; not part of the analysis protocol]
- **Decision rules,** to be fixed in the plan before any model runs:
  - **H4** is "consistent with" if the T2 **or** T3 rainfall-deficit (SPI-3 < 0) coefficient for birth weight is negative, with a 95% CI excluding 0, in the primary local-climate model.
  - **H5** applies the same rule to HAZ, for the 1,000-day drought exposure in the 24–59-month group.
  - **H6** is descriptive only, with no verdict.
  - The **ENSO × sensitivity** interaction is "consistent with" an El Niño pathway if the interaction is negative for birth weight or HAZ, with a 95% CI excluding 0.

## 14. Amendment 2026-10-23: execution rulings (P4–P15), deviations and deferrals

This section summarises the controller rulings made while the plan was executed. Where they conflict with §1–§13, they supersede those sections. The binding analysis specification is `papers/B/results_keys.md`, committed before any model was run. The data build is `analysis_child.parquet`, built by c6f4210, sha256 09ced2f2…5ab0; see `papers/B/output/provenance.csv`.

**Design as executed.**
- **Linkage (P4; supersedes §13.1).**
  - The author confirmed DHS GPS dataset access.
  - Cluster-GPS linkage is the primary exposure linkage: the cluster point's 0.25° rain cell and 1° temperature cell, with district fixed effects and clustering on the point-in-polygon GADM district.
  - The NFHS-district crosswalk linkage (`dist_*` exposures, `region_id_xw`) is kept as a pre-specified sensitivity analysis.
  - Coordinates and cluster-level tables stay git-ignored.
- **Punjab codes (P8).** NFHS-4 Punjab codes 39 (Sangrur) and 53 (SBS Nagar) are swapped on the DHS side. An explicit, tested override corrects the district-linkage sensitivity. The primary analysis uses point-in-polygon districts.
- **QC mask (P9).** A rain cell's pre-fit QC mask is the union of the masks of its dominant district and of every district whose clusters use the cell.
- **ENSO mediation (P10).** Common-in-time ENSO exposures make a district bootstrap anti-conservative.
  - ENSO mediation therefore uses an ENSO-year (June–May, of birth month − 5) block bootstrap with no birth year-month FE.
  - The FE are district + birth month of year + survey, with a birth-date spline.
  - The number of blocks (12) and an a-path Wald CI clustered by ENSO year are reported.
  - Mediation stays secondary.
- **Mother's age (P11).** Mother's age (`v012`) was added to the extract and enters the model in 5-year bands.
- **Dry orientation and verdict rules (P13; makes §13's rules operational).**
  - Every verdict model uses dry-positive terms: `−SPI-3` per trimester, 1,000-day dryness, and `sens_dry = −sens_beta`.
  - H4, H5 and ENSO × sensitivity are "consistent with" when the dry-term coefficient is negative and the 95% CI upper bound is below 0.
  - HAZ verdicts use the 24–59-month group.
  - The 1,000-day exposure is the month-weighted mean of the window SPI-3 values: (3,3,3,12,12)/33, or (3,3,3,12)/21 for 0–23 months, where complete.
- **Reproducibility (P14).** Non-PD vcov corrections are flagged per fit. Derived files carry data fingerprints, and `provenance.csv` records the data and code commits.

**Deviations from §4–§6 and §12 (to be stated in Methods).**
- **(a) Rainfall % departure.** Standardised SPI-3 is the pre-specified window exposure. The rainfall percent departure is used only for the district ENSO-rainfall sensitivity slope (`sens_beta`, 1981–2009 JJAS). It is not a window exposure.
- **(b) Location buffer.** Exposure is the 0.25° cell containing the displaced point. A 10 km area-weighted buffer for rural clusters is a sensitivity analysis. This replaces the §12 "5 km edge buffer" and the §4 2 km urban / 5 km rural buffers. Rationale: DHS displaces rural clusters by ≤ 5 km, and 1% by ≤ 10 km, so the 10 km buffer bounds the displacement. Methods must say this.
- **(c) Trimester ONI terms.** The approved plan specified a single pregnancy-mean ONI split into positive and negative parts, and that was run, with RONI and Niño-3.4 as alternatives. The secondary trimester-specific ONI terms of §12 are **deferred to Phase 3b** as a pre-specified supplementary model, to be reported whatever the result.
- **(e) Interview-month FE (deviation).** Which models carry an interview-month FE:
  - **HAZ, stunting and anaemia models** (birth year-month FE plus age-in-months FE, `hw1`): these nearly absorb the interview year-month, since interview = birth + age.
    - In NFHS-4 the absorption is exact: `hw1` equals the CMC age for 99.3% of children.
    - In NFHS-5 it is only near-complete: `hw1` is day-based and is one month below the CMC age for about 45% of children.
    - No separate interview-month FE is added.
  - **Birth-weight models** (H4, ENSO × sensitivity, total ENSO): no interview-month FE, by design. Birth weight is fixed at birth, so interview timing is not part of the exposure–outcome path. This is a deviation from §6.
  - **H6:** uses interview year-month FE.

**Deferred to Phase 3b (descriptive or supplementary, same model families).**
- **(d)** The `hw57` anaemia cross-check (descriptive).
- **(f)** WAZ/WHZ, underweight, wasting and moderate/severe anaemia, as supplementary tables. The §2 hypotheses do not require them.
- The manuscript, figures, STROBE/RECORD, journal choice, and the public code and aggregate repository, as already planned.
