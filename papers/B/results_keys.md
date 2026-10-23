# Paper B: results keys (contract)

This contract was written and committed before any Paper B model was run (Task 7; rulings P4, P8, P10, P13). It lists every key that `models/B/08_results.R` writes to `papers/B/results.json`. Entries use the Paper 0 schema `{"<key>": {"value", "display", "script", "description"}}`, written with `res_entry()` / `write_results()` (`models/common/R/results.R`). One command, run from `<repo>/models`, builds the file: `Rscript B/run_all.R`.

**Key notation.** In the key tables, `{a,b,c}` expands shell-style to one key per alternative, and an empty alternative (`{,.p}`) means the key without that suffix. For example, `B.h4.bw.{T2,T3}_dry{,.p}` stands for `B.h4.bw.T2_dry`, `B.h4.bw.T2_dry.p`, `B.h4.bw.T3_dry` and `B.h4.bw.T3_dry.p`. `08_results.R` expands every key in the first column of the tables below. It stops with an error unless the expanded set equals the set of keys written to `results.json`.

## Conventions

| Format | Display example | Rule |
|---|---|---|
| `ci1` | `−12.3 (95% CI −20.1 to −4.5)` | `fmt_ci(est, lo, hi, digits = 1)`. Grams, or grams per °C. |
| `ci2` | `1.04 (95% CI 1.01–1.07)` | `fmt_ci(digits = 2)`. Odds ratios, and the ENSO × sensitivity birth-weight interaction. |
| `ci3` | `−0.021 (95% CI −0.040 to −0.002)` | `fmt_ci(digits = 3)`. HAZ coefficients, and mediation a paths. |
| `ci4` | `−0.0012 (95% CI −0.0020 to −0.0004)` | `fmt_ci(digits = 4)`. The ENSO × sensitivity HAZ interaction. |
| `pct_ci` | `34% (95% CI 12–60%)` | Percent with 0 decimals. Negative bounds use U+2212 and `" to "`. |
| `p` | `0.012`, `<0.001` | `fmt_p`. |
| `int` | `12,345` | `fmt_int`. |
| `pct1` | `18.4%` | `fmt_pct(x, 1)`. `x` is a proportion. |
| `mean_sd` | `2,812 (SD 580)` | Weighted mean and weighted SD. Grams use 0 decimals; HAZ uses 2. |
| `verdict` | `consistent with` / `not consistent with` | Fixed strings. |
| text | as given per key | Fixed vocabulary only. |

All intervals in the key tables are 95% Wald CIs from the fit's own clustered vcov (`fe_table`), except where a key says bootstrap. Estimation on the t scale uses fixest defaults (t with df = min(G) − 1 for feols; z for feglm). Odds ratios are `exp(coef)`, with the CI bounds exponentiated.

`description` carries the model, the sample, `lo`/`hi` to 6 significant digits and the seed where relevant. `value` is the point estimate on the display scale, the count, the proportion or the verdict string. A non-estimable value is `null`, with the display `not estimable (<reason>)`.

## Definitions (fixed before any model was run)

**Data.** `data/external/nfhs_subset/analysis_child.parquet` (Task 5, b73e703; git-ignored; read with `arrow`). It has one row per NFHS-4/5 KR child: 492,547 rows. No row is excluded except by the sample rules below. *Note 2026-10-23 (ruling P15):* the file actually analysed was built by the pipeline at c6f4210 (the Task 5 review fix round; same rows and values as b73e703 plus that round's checks), sha256 `09ced2f28341727bfc7a7fed8476cf5dcfc9e9fc00518cb8d21c8be94b8c5ab0`; see `papers/B/output/provenance.csv`.

**Samples.**

| Sample | Rule |
|---|---|
| BW | `twin == FALSE` and `bw_g` non-missing. All births in the KR window are included, alive or dead. Card and recall both count (`m19a`). |
| HAZ 24–59 | `haz` non-missing and `age_fe` 24–59. |
| HAZ 0–23 | `haz` non-missing and `age_fe` 0–23. |
| Stunting 24–59 / 0–23 | As the HAZ samples, with `stunted` non-missing. |
| Anaemia | `anaemia_any` non-missing (6–59 months, tested; Task 2 rule). |
| Diarrhoea | `diarrhoea` non-missing (living children; Task 2 rule). |

- `age_fe` is `hw1`, the child's age in months as recorded with the measurement. Where `hw1` is missing, `age_fe` is `age_m` (interview CMC − birth CMC).
- Each model also drops rows with a missing model variable (complete cases). fixest's default `fixef.rm` then removes singleton and perfect-fit FE levels. A model's `n` is `nobs(fit)` after both steps.

**Weights.** `weight` is `v005`/1e6, rescaled to each round's sum (Task 5). Every model is weighted unless stated otherwise.

**Covariates (`COV`)**, all entered as factors:
- sex (`b4`);
- birth order band (`bord`: 1, 2, 3, 4+);
- mother's age band (`v012`: 15–19, 20–24, 25–29, 30–34, 35+);
- education (`v106`);
- wealth quintile (`v190`);
- residence (`v025`);
- caste (`s116`: SC, ST, OBC, none of these, and missing for code 8 or NA);
- religion (`v130`: Hindu, Muslim, Christian, Sikh, Buddhist, and other for codes 6, 7, 8, 9 and 96).

`m14` and `m15` are not used (spec §12).

**Dry orientation (ruling P13).** Every verdict model uses one convention: a positive value means drier.
- **Trimester dryness:** `T1_dry = −T1_spi3`, `T2_dry = −T2_spi3`, `T3_dry = −T3_spi3`. The SPI-3 values are the cluster-GPS window means (P4).
- **ENSO sensitivity:** `sens_dry = −sens_beta`, the district's 1981–2009 JJAS rainfall loss in percent per +1 °C ONI. A positive value means El Niño dries the district.
- **1,000-day dryness:** `dry1000 = −spi3_1000d`.
- **Coefficient direction:** a negative coefficient on a dry term means drought lowers the outcome.
- **SPI scale:** keys ending `_spi3` or `.spi3` report the same coefficient on the SPI-3 scale (the sign flipped and the bounds swapped). They are descriptive only.

**1,000-day SPI-3 (`spi3_1000d`, P13).** This is a month-weighted mean of the window means:
- **24–59 months:** `(3·T1 + 3·T2 + 3·T3 + 12·inf0_11 + 12·inf12_23) / 33`, using the `_spi3` columns. It is NaN if any component is NaN.
- **0–23 months:** `(3·T1 + 3·T2 + 3·T3 + 12·inf0_11) / 21`. It is defined only where all four are non-NaN.
  - Windows that run past the interview are NaN (Task 4). This group is therefore in practice the children whose 0–11-month window is complete at interview, aged about 11–23 months.
  - This group has no verdict.

**Recent windows (H6).** Diarrhoea uses `recent1` (the month before the interview). Anaemia uses `recent12` (the 12 months before the interview). Each has `_dry = −_spi3` and `_tmax_anom`.

**ENSO year (`enso_year`, P10 and R-1).** This is the June–May year containing the T2 midpoint month, the birth month − 5. It is labelled by the calendar year of its June. For a midpoint month of June–December it equals that month's year; for January–May it equals the year − 1.

**El Niño-developing exposure.** `preg_elnino_dev` (Task 4): any pregnancy month in June–September of an `el_nino_developing` year. Rows where it is NA are dropped from models that use it.

**Birth-date spline.** `ns(birth_cmc, df = 4)` from `splines`. The basis is computed once on all 492,547 rows and stored as the columns `cmc_ns1`–`cmc_ns4`.

**Seed.** `SEED_B = 20261022`. It is passed to `mediate_fe` and also set globally.

## Pre-specified models

All models are fitted with `fit_fe` and tabulated with `fe_table`, passing the same `cluster` to both (R-2). Unless a row says otherwise, models are weighted and the cluster is `~region_id + birth_ym`. `region_id` is the point-in-polygon GADM district (P8). HAZ, stunting, anaemia and diarrhoea models always add `age_fe` (age in months) as a fixed effect.

| ID | Outcome (family) | Right-hand side | Fixed effects |
|---|---|---|---|
| H4-bw | `bw_g` (gaussian) | `T1_dry + T2_dry + T3_dry + T1_tmax_anom + T2_tmax_anom + T3_tmax_anom + COV` | `region_id + birth_ym` |
| H4-lbw | `lbw` (logit) | as H4-bw | `region_id + birth_ym` |
| H5-haz | `haz` (gaussian), 24–59 and 0–23 separately | `dry1000 + COV` | `region_id + birth_ym + age_fe` |
| H5-stunt | `stunted` (logit), 24–59 and 0–23 | `dry1000 + COV` | `region_id + birth_ym + age_fe` |
| H6-diarr | `diarrhoea` (logit) | `recent1_dry + recent1_tmax_anom + COV` | `region_id + interview_ym + age_fe` |
| H6-anaem | `anaemia_any` (logit) | `recent12_dry + recent12_tmax_anom + COV` | `region_id + interview_ym + age_fe` |
| ENSO×S | `bw_g` (BW sample); `haz` (24–59, 0–23) | `elnino_x_sens_dry + COV`, where `elnino_x_sens_dry = preg_elnino_dev × sens_dry`. Both main effects are absorbed: `sens_dry` by `region_id`, `preg_elnino_dev` by `birth_ym`. | `region_id + birth_ym` (+ `age_fe` for HAZ) |
| ENSO-total | `bw_g`; `haz` (24–59, 0–23) | `<idx>_pos + <idx>_neg + COV + cmc_ns1..4`, with `<idx>` = `preg_oni` (primary), `preg_roni` or `preg_nino34`. `_pos = max(x, 0)`, `_neg = min(x, 0)`. | `region_id + birth_month_of_year + round` (+ `age_fe` for HAZ) |

The H6 models are descriptive only. The ENSO-total models report total associations and have no verdict. For La Niña the display is per 1 °C cooler: `−β` on `_neg`, with the bounds swapped; the p value is unchanged.

**Mediation (SECONDARY; P10, R-1, P13 Q4).** `mediate_fe(data, outcome = "bw_g", exposure, mediator, fe = c("region_id", "birth_month_of_year", "round"), covars = c(COV, cmc_ns1..4), B = 1000, cores = 6, seed = SEED_B, boot_unit = "enso_year")`.
- The data are the BW sample, restricted to the complete rows of the call's variables.
- The mediator is `T2_spi3` in one call and `T3_spi3` in another, on the SPI scale (no verdict).
- **Primary (`oni`):** the exposure is `preg_oni`.
- **Sensitivity (`onipos`):** the exposure is `preg_oni_pos`, with `preg_oni_neg` added to `covars`.
- There is no `birth_ym` FE.
- **Time-clustered Wald check of the a path:** `fit_fe(mediator ~ exposure + COV + cmc_ns1..4 | region_id + birth_month_of_year + round, cluster = ~enso_year)`, tabulated by `fe_table(cluster = ~enso_year)`.
- No other mediation analysis is run.

## Decision rules (spec §13, ruling P13; fixed before results)

- **H4, birth weight:** `consistent with` if, in H4-bw, the coefficient on `T2_dry` **or** on `T3_dry` is negative with a 95% CI excluding 0 (`hi < 0`). Otherwise `not consistent with`.
- **H5, HAZ:** `consistent with` if, in H5-haz for the 24–59-month group, the coefficient on `dry1000` is negative with `hi < 0`. Otherwise `not consistent with`. The 0–23 group is reported without a verdict.
- **ENSO × sensitivity:** `consistent with` an El Niño pathway if the coefficient on `elnino_x_sens_dry` is negative with `hi < 0` for birth weight **or** for HAZ 24–59. Otherwise `not consistent with`. HAZ 0–23 is reported without a verdict.
- **H6:** descriptive only, with no verdict.
- **Mediation and ENSO-total:** no verdict.

## Sensitivity analyses (each pre-specified; each gets rows in `sensitivity.csv` and keys)

Each analysis changes only what its row states; everything else stays as in the primary model.

| ID | Applied to | Change |
|---|---|---|
| `card` | h4bw | `bw_card == TRUE` (birth weight from a card). |
| `lastbirth` | h4bw | `bidx == 1` (the most recent birth). |
| `unweighted` | h4bw, h5haz, xsbw, xshaz | `weights = NULL`. |
| `nonmover` | h4bw, h5haz, xsbw, xshaz | `v104 == 95` (always lived here), or `v104 < 95` and `v104 ≥ floor(age_m / 12) + 1`. Visitors (`v104 == 96`) and codes > 96 are excluded. |
| `precovid` | h4bw, h5haz, xsbw, xshaz | `v008 < 1450` (interviews before October 2020). |
| `motherfe` | h4bw, h5haz | FE `mother_id + birth_ym` (+ `age_fe` for HAZ), replacing `region_id`. Restricted to mothers with ≥ 2 children in that model's complete-case sample. The cluster is unchanged. |
| `shiftm1` | h4bw | Conception at birth − 10: T1 (−10, −8), T2 (−7, −5), T3 (−4, −2). Applied to the SPI-3 and tmax windows. |
| `shiftp1` | h4bw | Conception at birth − 8: T1 (−8, −6), T2 (−5, −3), T3 (−2, 0). |
| `district` | h4bw, h5haz, xsbw, xshaz | District linkage: `dist_*` window exposures; `sens_dry_xw = −sens_beta_xw`; FE and cluster use `region_id_xw` in place of `region_id`. |
| `buf10` | h4bw, h5haz | Rural clusters (`urban == FALSE`) use `<window>_spi3_buf10`; urban clusters keep `<window>_spi3`. tmax is unchanged. |
| `noflag` | h4bw, h5haz | Excludes children with `<window>_fallback_rain_any` or `<window>_mis_district_any` in any window the model uses. |
| `nogaz` | h4bw, h5haz, xsbw, xshaz | Excludes gazetteer clusters (`source == "GAZ"`; spec §12). |
| `cluster1` | xsbw, xshaz | SEs clustered by `region_id` only (`cluster = ~region_id` in `fit_fe` and `fe_table`). |
| `placebo` | h4bw | H4-bw plus `placebo_dry = −placebo_spi3` and `placebo_tmax_anom`, the pseudo-pregnancy window from birth +1 to +9 months (NaN if incomplete at interview). |
| `ensoplacebo` | ENSO-total bw (ONI) | Adds `post_oni_pos` and `post_oni_neg`, where `post_oni` is the mean `oni_anom` over birth +1 to +12 months (by central month), built with `pipeline.nfhs.exposure.enso_exposure`. |
| `hbplaus` | H6-anaem | `hb_plausible == 1` (Hb 4.0–18.0 g/dL; ruling P3). |
| `lpm` | H4-lbw, H5-stunt 24–59 | A linear probability model (feols). Coefficients × 100 are percentage points. |

Model codes: `h4bw` = H4-bw (terms `T2_dry`, `T3_dry`); `h5haz` = H5-haz 24–59 (term `dry`); `xsbw` and `xshaz` = ENSO×S for birth weight and HAZ 24–59 (term `inter`).

The shifted windows and `post_oni` come from `models/B/py/conception_shift.py` (P13 Q5). It imports `pipeline.nfhs.gps` and `pipeline.nfhs.exposure` read-only and writes the git-ignored file `data/external/nfhs_subset/conception_shift.parquet`. That file is aligned to the analysis parquet's row order, and the alignment is checked by `birth_cmc`.

**Sensitivity verdict check.** For every (spec, model) pair with model code h4bw, h5haz, xsbw or xshaz, `07_sensitivity.R` applies that model's rule:
- h4bw: `T2_dry` or `T3_dry` negative with `hi < 0`;
- h5haz: `dry` negative with `hi < 0`;
- xsbw and xshaz: `inter` negative with `hi < 0`.

It compares the result with the same rule applied to the primary fit of that model. A pair whose rule outcome differs from the primary is a sensitivity failure (listed in `B.sens.differ`).

## Aggregate outputs (`papers/B/output/`, committed)

These are aggregate tables only. No cell is based on n < 30; any such cell is suppressed (values NA, `suppressed = TRUE`). There are no `v001`, `caseid`, `mother_id` or coordinate columns, and no child-level rows.

| File | Content |
|---|---|
| `sample_flow.csv` | Counts per round: all children, BW sample (and twins excluded), HAZ 24–59 and 0–23, anaemia, diarrhoea. |
| `descriptive_by_round.csv` | Weighted outcome summaries by round. |
| `descriptive_by_enso.csv` | Weighted outcome summaries by `preg_elnino_dev` (exposed or not). |
| `missingness_state_round.csv` | Percent missing of T1/T2/T3 SPI-3 and T2 tmax anomaly, by GADM state × round. |
| `missingness_tn_birthyear.csv` | The same for Tamil Nadu by birth year × round (the NFHS-4 2011–12 mask), plus national rows. |
| `primary_models.csv` | Every coefficient row of the H4, H5, H6 and ENSO×S models (term, est, lo, hi, p, n, vcov, df). |
| `enso_total.csv` | ENSO-total models for ONI, RONI and Niño-3.4. |
| `mediation.csv` | Mediation estimates, bootstrap CIs, flags and the time-clustered a-path check. |
| `sensitivity.csv` | Every sensitivity row (spec, model, term, est, lo, hi, p, n). |
| `sensitivity_verdicts.csv` | The rule check per (spec, model) against the primary. |
| `provenance.csv` | Added by the amendment below. The analysis parquet's sha256, md5, size, row count and build commit, plus the results code commit. |

**Amendment 2026-10-23 (ruling P14). Provenance and output columns only. No result key and no value changes.**
- `provenance.csv` was added.
- Every coefficient table gained two columns:
  - `vcov_fixed`: TRUE where fixest corrected a non-positive-definite clustered vcov;
  - `unit`: the scale and per-unit of `est`. For example, the LPM rows in `sensitivity.csv` are probabilities; their keys are × 100, in percentage points.
- `sample_flow.csv` count columns were renamed `n_*`. The mediation `n_blocks` column was renamed `enso_year_blocks`.

## 1. Samples and descriptive (`01_descriptive.R`)

| Key | Meaning | Display |
|---|---|---|
| `B.n.children{,.nfhs4,.nfhs5}` | Rows in the analysis dataset: total, then by round | int |
| `B.n.districts` | Distinct `region_id` (GADM districts) among all children | int |
| `B.desc.n_bw` | BW sample size (before complete-case) | int |
| `B.desc.n_twins_excluded` | Twins with a birth weight, excluded from BW | int |
| `B.desc.{n_haz2459,n_haz0023,n_anaemia,n_diarrhoea}` | Sample sizes (before complete-case) | int |
| `B.desc.bw_mean` | Weighted mean (SD) of `bw_g`, BW sample | mean_sd (g, 0 decimals) |
| `B.desc.lbw_pct` | Weighted % LBW, BW sample | pct1 |
| `B.desc.bw_card_pct` | Unweighted % of the BW sample with a card weight | pct1 |
| `B.desc.bw_heap500_pct` | Unweighted % of BW-sample weights that are multiples of 500 g (heaping check) | pct1 |
| `B.desc.haz2459_mean` | Weighted mean (SD) HAZ, 24–59 | mean_sd (2 decimals) |
| `B.desc.stunt2459_pct` | Weighted % stunted, 24–59 | pct1 |
| `B.desc.anaemia_pct` | Weighted % any anaemia | pct1 |
| `B.desc.diarrhoea_pct` | Weighted % diarrhoea in the last 2 weeks | pct1 |
| `B.desc.elnino_dev_pct` | Weighted % of the BW sample with `preg_elnino_dev` TRUE | pct1 |
| `B.desc.bw_mean.{elnino,other}` | Weighted mean (SD) of `bw_g` by `preg_elnino_dev` | mean_sd |
| `B.desc.lbw_pct.{elnino,other}` | Weighted % LBW by `preg_elnino_dev` | pct1 |
| `B.desc.haz2459_mean.{elnino,other}` | Weighted mean (SD) HAZ 24–59 by `preg_elnino_dev` | mean_sd |
| `B.desc.sens_dry_median` | Median `sens_dry` over the districts with a value (% per °C ONI) | 1 decimal, U+2212 |
| `B.desc.sens_dry_pct_pos` | % of those districts with `sens_dry > 0` (El Niño dries them) | pct1 |
| `B.desc.n_sens_districts` | Districts with a non-missing `sens_beta` | int |
| `B.miss.T2_spi3_pct.{nfhs4,nfhs5}` | % of all children with missing `T2_spi3`, by round | pct1 |
| `B.miss.tn4_T2_spi3_pct` | % missing `T2_spi3`, Tamil Nadu NFHS-4 (2011–12 mask) | pct1 |

## 2. H4: birth weight (`02_birthweight.R`)

| Key | Meaning | Display |
|---|---|---|
| `B.h4.bw.{T1,T2,T3}_dry{,.p}` | H4-bw: grams per +1 SPI unit drier in the trimester; p | ci1; p |
| `B.h4.bw.{T1,T2,T3}_spi3` | The same coefficient on the SPI-3 scale (grams per +1 SPI, wetter) | ci1 |
| `B.h4.bw.{T1,T2,T3}_tmax{,.p}` | H4-bw: grams per +1 °C Tmax anomaly; p | ci1; p |
| `B.h4.bw.n` | H4-bw n | int |
| `B.h4.lbw.{T1,T2,T3}_dry{,.p}` | H4-lbw: OR of LBW per +1 SPI unit drier; p | ci2; p |
| `B.h4.lbw.{T1,T2,T3}_tmax{,.p}` | H4-lbw: OR per +1 °C; p | ci2; p |
| `B.h4.lbw.n` | H4-lbw n | int |
| `B.h4.verdict` | H4 verdict by rule | verdict |
| `B.h4.verdict_terms` | Which term(s) met the rule | text: `T2`, `T3`, `T2 and T3`, or `neither` |

## 3. H5: growth (`03_growth.R`)

| Key | Meaning | Display |
|---|---|---|
| `B.h5.{haz2459,haz0023}.dry{,.p}` | H5-haz: HAZ per +1 SPI unit drier over the 1,000-day (or 0–23) exposure; p | ci3; p |
| `B.h5.{haz2459,haz0023}.spi3` | The same on the SPI-3 scale | ci3 |
| `B.h5.{haz2459,haz0023}.n` | n | int |
| `B.h5.{stunt2459,stunt0023}.dry{,.p}` | H5-stunt: OR of stunting per +1 SPI unit drier; p | ci2; p |
| `B.h5.{stunt2459,stunt0023}.n` | n | int |
| `B.h5.verdict` | H5 verdict by rule (24–59 HAZ) | verdict |

## 4. H6: recent morbidity, descriptive (`04_morbidity.R`)

| Key | Meaning | Display |
|---|---|---|
| `B.h6.diarr.{recent1_dry,recent1_tmax}{,.p}` | H6-diarr: OR per +1 SPI unit drier, or per +1 °C, in the month before interview; p | ci2; p |
| `B.h6.diarr.n` | n | int |
| `B.h6.anaem.{recent12_dry,recent12_tmax}{,.p}` | H6-anaem: OR per +1 SPI unit drier, or per +1 °C, over the 12 months before interview; p | ci2; p |
| `B.h6.anaem.n` | n | int |

## 5. ENSO (`05_enso.R`)

| Key | Meaning | Display |
|---|---|---|
| `B.ensoxs.bw{,.p,.n}` | ENSO×S, birth weight: grams per unit of `preg_elnino_dev × sens_dry` (per 1% per °C of El Niño rainfall loss, in El Niño-developing pregnancies); p; n | ci2; p; int |
| `B.ensoxs.{haz2459,haz0023}{,.p,.n}` | ENSO×S, HAZ per unit of the interaction; p; n | ci4; p; int |
| `B.ensoxs.verdict` | ENSO × sensitivity verdict by rule | verdict |
| `B.enso.{oni,roni,nino34}.bw.{elnino,lanina}{,.p}` | ENSO-total: grams per +1 °C of the El Niño term, or per 1 °C cooler of the La Niña term; p | ci1; p |
| `B.enso.{oni,roni,nino34}.{haz2459,haz0023}.{elnino,lanina}{,.p}` | The same for HAZ | ci3; p |
| `B.enso.{oni,roni,nino34}.{bw,haz2459,haz0023}.n` | n | int |

## 6. Mediation, secondary (`06_mediation.R`)

`{oni,onipos}` = primary exposure `preg_oni` / sensitivity `preg_oni_pos` (with `preg_oni_neg` as a covariate). `{T2,T3}` = mediator `T2_spi3` / `T3_spi3`.

| Key | Meaning | Display |
|---|---|---|
| `B.med.{oni,onipos}.{T2,T3}.a` | a path: SPI-3 per +1 °C of the exposure; ENSO-year block-bootstrap percentile CI | ci3 |
| `B.med.{oni,onipos}.{T2,T3}.a_wald{,.p}` | a path, time-clustered Wald CI (cluster `enso_year`); p | ci3; p |
| `B.med.{oni,onipos}.{T2,T3}.b` | b path: grams per +1 SPI-3 (adjusted for the exposure); bootstrap CI | ci1 |
| `B.med.{oni,onipos}.{T2,T3}.{c_total,c_direct,indirect}` | Total, direct and indirect (a × b) effect, grams per +1 °C; bootstrap CI | ci1 |
| `B.med.{oni,onipos}.{T2,T3}.prop` | Proportion mediated; bootstrap CI | pct_ci when `prop_flag == "ok"`, otherwise `not estimated (<prop_flag>)` |
| `B.med.{oni,onipos}.{T2,T3}.{flag,prop_flag}` | mediate_fe flags | text, as returned |
| `B.med.{oni,onipos}.{T2,T3}.{n,B_ok}` | Rows used; successful resamples | int |
| `B.med.n_blocks` | ENSO-year blocks (distinct `enso_year`) in the primary T2 mediation sample | int |
| `B.med.B` | Bootstrap resamples requested | int |

If `flag == "too_many_failed_fits"`, every bootstrap-CI display is `not estimated (too many failed bootstrap fits)`.

## 7. Sensitivity (`07_sensitivity.R`)

| Key | Meaning | Display |
|---|---|---|
| `B.sens.{card,lastbirth,unweighted,nonmover,precovid,motherfe,shiftm1,shiftp1,district,buf10,noflag,nogaz}.h4bw.{T2_dry,T3_dry,n}` | H4-bw under each analysis | ci1; int |
| `B.sens.{unweighted,nonmover,precovid,motherfe,district,buf10,noflag,nogaz}.h5haz.{dry,n}` | H5-haz 24–59 under each analysis | ci3; int |
| `B.sens.{unweighted,nonmover,precovid,district,nogaz,cluster1}.xsbw.{inter,n}` | ENSO×S birth weight under each analysis | ci2; int |
| `B.sens.{unweighted,nonmover,precovid,district,nogaz,cluster1}.xshaz.{inter,n}` | ENSO×S HAZ 24–59 under each analysis | ci4; int |
| `B.sens.placebo.h4bw.{placebo_dry,placebo_tmax,n}` | Placebo pseudo-pregnancy window terms in H4-bw | ci1; int |
| `B.sens.ensoplacebo.bw.{post_elnino,post_lanina,n}` | Post-birth ONI terms (grams per +1 °C; La Niña per 1 °C cooler) in ENSO-total bw | ci1; int |
| `B.sens.hbplaus.anaem.{recent12_dry,recent12_tmax,n}` | H6-anaem restricted to plausible Hb | ci2; int |
| `B.sens.lpm.lbw.{T2_dry,T3_dry,n}` | LBW as an LPM, percentage points per +1 SPI unit drier | ci2; int |
| `B.sens.lpm.stunt2459.{dry,n}` | Stunting 24–59 as an LPM, percentage points | ci2; int |
| `B.sens.n_checks` | (spec, model) rule checks | int |
| `B.sens.n_same` | Checks whose rule outcome equals the primary's | int |
| `B.sens.differ` | The (spec/model) pairs that differ, comma-separated | text, or `none` |
| `B.sens.{placebo_flag,ensoplacebo_flag}` | Whether the placebo `placebo_dry` / `post_elnino` CI excludes 0 | text: `CI includes 0` or `CI excludes 0` |

## Amendment 2026-10-23: Phase 3b supplementary keys

This amendment was committed before any of the models below was run (Phase 3b plan, Task 1). It adds keys and files only. Every key and value above stays as it is. The added items are the pre-specified deferrals of spec §14: (c) trimester ONI terms, (d) the `hw57` cross-check and (f) the secondary outcomes. They are reported whatever the result, and **none of them carries a verdict**. Two text keys, generated by code from the existing verdicts, and the figure data files are added as well.

`Rscript B/run_all.R --use-cache` builds everything: the cached Task 7 steps and the mediation bootstrap are reused, and two new cached steps run, `supplementary` (`09_supplementary.R`) and `figure_data` (`10_figure_data.R`). The conventions, samples, weights, covariates (`COV`), dry orientation, clustering (`~region_id + birth_ym`) and small-cell rule are those defined above.

### Supplementary models (`09_supplementary.R`)

| ID | Outcome (family) | Right-hand side | Fixed effects |
|---|---|---|---|
| ENSO-trim | `bw_g` (gaussian), BW sample; `haz` (gaussian), HAZ 24–59 | `T1_oni + T2_oni + T3_oni + COV + cmc_ns1..4`. These replace `preg_oni_pos + preg_oni_neg` in ENSO-total. `<T>_oni` is the mean `oni_anom` over the trimester's months (by central month; Task 4 `enso_exposure`), linear, not split into El Niño and La Niña parts. | `region_id + birth_month_of_year + round` (+ `age_fe` for HAZ), as ENSO-total |
| S-waz, S-whz | `waz` or `whz` (gaussian), 24–59 and 0–23 separately | `dry1000 + COV`. This is H5-haz with the outcome replaced. The exposure is `dry1000_2459` for 24–59 and `dry1000_0023` for 0–23. | `region_id + birth_ym + age_fe` |
| S-underweight, S-wasted | `underweight` or `wasted` (logit), 24–59 and 0–23 | as S-waz / S-whz | `region_id + birth_ym + age_fe` |
| S-anaemia-modsev | `anaemia_modsev` (logit), anaemia sample | `recent12_dry + recent12_tmax_anom + COV`. This is H6-anaem with the outcome replaced (moderate or severe: `hw56` raw < 100). | `region_id + interview_ym + age_fe` |

Samples for the secondary outcomes follow the HAZ rule: the outcome non-missing and `age_fe` 24–59 or 0–23. `underweight` is defined exactly where `waz` is, and `wasted` exactly where `whz` is (Task 2). Complete-case and `fixef.rm` rules are as above.

**`hw57` cross-check (descriptive).** `hw57` is the DHS anaemia level (1 severe, 2 moderate, 3 mild, 4 not anaemic). It is not in `analysis_child.parquet`, so it is read from the git-ignored `data/external/nfhs_subset/kr_harmonised.parquet`. That file has the same 492,547 rows in the same order (Task 5 excluded none). Its alignment is checked row by row on `round`, `v001`, `bidx` and `b3 == birth_cmc`, and the step stops on any mismatch. Only `hw57` is attached.
- **Sample:** the anaemia sample (`anaemia_any` non-missing) with `hw57` in 1–4, by round.
- **Measures:**
  - any anaemia by `hw56` is `anaemia_any` (raw < 110);
  - any anaemia by `hw57` is `hw57` in 1–3;
  - both are weighted shares on the same rows;
  - the difference is the `hw57` share minus the `hw56` share, in percentage points;
  - agreement is the unweighted share of children classified the same way by both;
  - kappa is the unweighted Cohen's kappa of the two binary classifications, `(po − pe) / (1 − pe)`. When `pe = 1`, it is 1 if `po = 1`, otherwise not estimable.

### 8. Supplementary keys (`09_supplementary.R`; no verdicts)

| Key | Meaning | Display |
|---|---|---|
| `B.supp.trim_oni.bw.{T1,T2,T3}{,.p}` | ENSO-trim, birth weight: grams per +1 °C of the trimester mean ONI (`T1_oni`, `T2_oni`, `T3_oni`); p | ci1; p |
| `B.supp.trim_oni.haz2459.{T1,T2,T3}{,.p}` | ENSO-trim, HAZ 24–59: HAZ per +1 °C of the trimester mean ONI; p | ci3; p |
| `B.supp.trim_oni.{bw,haz2459}.n` | ENSO-trim n | int |
| `B.supp.hw57.{nfhs4,nfhs5}.hw56_pct` | Weighted % any anaemia by `hw56` (raw < 110) on the cross-check sample | pct1 |
| `B.supp.hw57.{nfhs4,nfhs5}.hw57_pct` | Weighted % any anaemia by `hw57` (levels 1–3) on the same rows | pct1 |
| `B.supp.hw57.{nfhs4,nfhs5}.diff_pp` | `hw57` share minus `hw56` share, percentage points | 2 decimals, U+2212 |
| `B.supp.hw57.{nfhs4,nfhs5}.agree_pct` | Unweighted % of children classified the same way by both | pct1 |
| `B.supp.hw57.{nfhs4,nfhs5}.kappa` | Unweighted Cohen's kappa of the two classifications | 3 decimals |
| `B.supp.hw57.{nfhs4,nfhs5}.n` | Children in the cross-check sample | int |
| `B.supp.{waz2459,waz0023,whz2459,whz0023}.dry1000{,.p,.n}` | S-waz / S-whz: z-score per +1 SPI unit drier over the 1,000-day exposure (24–59) or the (3,3,3,12)/21 exposure (0–23); p; n | ci3; p; int |
| `B.supp.{underweight2459,underweight0023,wasted2459,wasted0023}.dry1000{,.p,.n}` | S-underweight / S-wasted logit: OR per +1 SPI unit drier; p; n | ci2; p; int |
| `B.supp.anaemia_modsev.{recent12_dry,recent12_tmax}{,.p}` | S-anaemia-modsev logit: OR of moderate or severe anaemia per +1 SPI unit drier, or per +1 °C, over the 12 months before interview; p. The key says `recent12`, not `dry1000`, because the exposure is the H6 window. | ci2; p |
| `B.supp.anaemia_modsev.n` | S-anaemia-modsev n | int |

### 9. Generated text keys (`08_results.R`; built from the verdict keys, no new model)

These are fixed templates filled by code, so the manuscript states the verdicts and signs exactly as `results.json` holds them. They contain no causal wording ("cause", "effect", "impact", "led to", "reduce"). Their inputs are `B.h4.verdict`, `B.h4.verdict_terms`, `B.h5.verdict`, `B.ensoxs.verdict`, and the sign of the value of `B.enso.oni.haz2459.elnino`.

| Key | Meaning | Display |
|---|---|---|
| `B.abstract.conclusion_sentence` | One-sentence conclusion (template below) | text |
| `B.discussion.verdict_summary` | One-sentence summary of the three verdicts (template below) | text |

**`B.abstract.conclusion_sentence`.**
- **All three verdicts `not consistent with`:** "In Indian children surveyed in NFHS-4 and NFHS-5, drier conditions in pregnancy and early childhood and El Niño exposure in pregnancy did not show the pre-specified associations with lower birth weight or lower height-for-age; the total El Niño association with height-for-age at 24–59 months was {sign clause}."
- **Otherwise:** "In Indian children surveyed in NFHS-4 and NFHS-5, {parts meeting their rule, joined by "; "}, consistent with the pre-specified hypotheses; {labels of the parts not meeting their rule, joined by " and "} did not meet the pre-specified rule; the total El Niño association with height-for-age at 24–59 months was {sign clause}." If every part meets its rule, the "did not meet" clause is omitted.
- **Parts meeting their rule:**
  - H4: "drier conditions in the {B.h4.verdict_terms} trimester of pregnancy were associated with lower birth weight". When `B.h4.verdict_terms` is "T2 and T3", "trimesters" is used.
  - H5: "drier conditions over the first 1,000 days were associated with lower height-for-age at 24–59 months".
  - ENSO × sensitivity: "El Niño exposure in pregnancy was associated with lower birth weight or height-for-age where El Niño dries the district".
- **Labels:** "the birth-weight hypothesis (H4)", "the height-for-age hypothesis (H5)" and "the El Niño × district-sensitivity hypothesis".
- **Sign clause:**
  - value > 0: "positive, opposite in sign to the hypothesis";
  - value < 0: "negative, in the direction of the hypothesis";
  - value 0: "zero";
  - null: "not estimable".

**`B.discussion.verdict_summary`.** "By the pre-specified rules, the birth-weight results were {B.h4.verdict} H4 (trimester terms meeting the rule: {B.h4.verdict_terms}), the height-for-age results at 24–59 months were {B.h5.verdict} H5, and the El Niño × district-sensitivity results were {B.ensoxs.verdict} an El Niño pathway."

### Supplementary aggregate outputs (`papers/B/output/`, committed)

| File | Content |
|---|---|
| `supp_models.csv` | Every coefficient row of ENSO-trim, S-waz, S-whz, S-underweight, S-wasted and S-anaemia-modsev, in the `primary_models.csv` columns |
| `supp_hw57.csv` | Per round: `n_children`, the `hw56` and `hw57` weighted shares, the difference in percentage points, agreement and kappa |

### Figure data (`papers/B/figures/data/`, committed; `10_figure_data.R`)

These are aggregates only, under the same rules as `papers/B/output/`:
- every count column (`n`, `n_*`) is ≥ 30 or NA, and `check_aggregate_table` enforces this;
- a cell with n < 30 keeps its row, with its values and n set to NA and `suppressed = TRUE`. It is never dropped;
- there are no `v001`, `caseid`, `mother_id`, coordinate, `birth_ym` or `region_id` columns, and no child-level rows.

| File | Content |
|---|---|
| `fig1_timeline.csv` | One row per calendar month, January 2006 to December 2021 (192 rows). The columns are: `month` (YYYY-MM), `oni_anom` and `roni` (`data/derived/enso_monthly.csv`, by central month), and `enso_phase` (that file's NOAA episode phase: `el_nino`, `la_nina` or `neutral`). Then, per round (`nfhs4`, `nfhs5`): `n_births_<round>` (children in `analysis_child.parquet` born that month), `weighted_births_<round>` (sum of `weight`) and `suppressed_<round>`. A month with n < 30, including 0, is suppressed. |
| `fig2_district_map.csv` | Long format, one row per GADM district × panel × round. The full grid is used: every `region_id` in the analysis dataset appears in every panel. Columns: `district_gid` (the GADM GID_2 held in `region_id`), `state`, `panel`, `round`, `n_children`, `value` and `suppressed`. **Panels:** <br>• `T2_dry_cohort_2015_16`: weighted mean `T2_dry` of children born June 2015–May 2016 (the 2015 El Niño cohort), with `round` = `pooled`; <br>• `T2_dry_other_cohorts`: the same for all other birth months, `pooled`; <br>• `T2_dry`: weighted mean `T2_dry` of all children with a value, by round; <br>• `stunting_24_59`: weighted share stunted in the HAZ 24–59 sample, by round. <br>`n_children` counts the children behind `value`. When it is < 30 (including 0), `value` and `n_children` are NA and `suppressed` is TRUE. |
| `fig3_coefs.csv` | Every row of `primary_models.csv` (term, est, lo, hi, p, n, unit and the model columns), plus `panel` and `age_group`. `panel` is `bw`, `lbw`, `haz`, `stunting`, `h6` or `ensoxs`. `age_group` is `birth`, `24-59`, `0-23` or `recent`. |
| `fig4_sensitivity.csv` | The 32 (spec, model) rule checks of `sensitivity.csv` (model codes h4bw, h5haz, xsbw and xshaz; placebo excluded), one row per rule term. Each row has `rule_met`, `same_as_primary`, and the primary fit's `primary_est`, `primary_lo` and `primary_hi` for the same model and term. |
| `figS1_placebo.csv` | The H4-bw primary rows (T1–T3 dry and Tmax), with `source = primary`, and the placebo fit's rows (`T2_dry`, `T3_dry`, `placebo_dry`, `placebo_tmax_anom`), with `source = placebo`. Columns: term, est, lo, hi, p, n, unit. |

### Note 2026-10-23 (Task 1 review, fix round 1, ruling Q2): differencing, text templates, display columns

This note was committed before the files were regenerated. It changes no value key: every model, estimate, count and verdict key keeps its value. Two kinds of thing change. The two generated text keys of section 9 get new templates, and the aggregate files and display columns listed below are amended.

**1. Suppressed cells must not be recoverable by differencing (I-1).** Per-cell suppression is not enough when a suppressed cell's n equals a published total minus published cells. The following changes close the routes found.
- **`fig2_district_map.csv`.**
  - The exact `n_children` column is removed. It is replaced by `n_band`: `<30`, `30–99`, `100–499` or `500+`, using U+2013.
  - `value` and `suppressed` stay. A suppressed row has `n_band` `<30` (including 0) and `value` NA.
  - **Routes closed:** the cohort panels against the by-round panels (cohort + other = NFHS-4 + NFHS-5 per district), and the by-round panel against the state × round counts of `missingness_state_round.csv`.
- **`fig1_timeline.csv`.**
  - `n_births_<round>` and `weighted_births_<round>` are removed.
  - They are replaced by `birth_share_pct_<round>`: the month's share of the round's weighted births (sum of `weight` in the month over the round's sum), in %, rounded to 0.1.
  - A month with n < 30, including 0, has NA, and `suppressed_<round>` stays.
  - **Route closed:** a suppressed month's n equalled `sample_flow.n_children` minus the published months.
- **`missingness_tn_birthyear.csv`.**
  - `n_children` is replaced by `n_band`, with the same bands.
  - **Route closed:** the one suppressed Tamil Nadu NFHS-4 birth-year cell equalled the Tamil Nadu NFHS-4 `n_children` of `missingness_state_round.csv` minus the published birth years.
  - The India rows are banded too, for one rule per file.
- **Not changed, because they leave no suppressed cell to recover:**
  - `missingness_state_round.csv` (no state × round cell is suppressed, and its rows sum to `sample_flow`);
  - `sample_flow.csv`, `descriptive_by_round.csv`, `descriptive_by_enso.csv` and `supp_hw57.csv` (no suppressed cells; every n is at least 30);
  - the model tables (model n only).
- `prep.R` holds the band rule (`band_n`). A test checks that a synthetic district or month with one suppressed cell cannot have its n recovered from the published columns.

**2. Text templates (section 9; M-1, M-2, M-3).** These replace the section 9 templates.
- **Sign clause.** It now uses the 95% CI of `B.enso.oni.haz2459.elnino`, read from its description (`lo=`, `hi=`):
  - estimate > 0 and lo > 0: "positive (95% CI excluding 0), opposite in sign to the hypothesis";
  - estimate < 0 and hi < 0: "negative (95% CI excluding 0), in the direction of the hypothesis";
  - estimate > 0 with a CI including 0: "positive, with a 95% CI including 0";
  - estimate < 0 with a CI including 0: "negative, with a 95% CI including 0";
  - estimate 0: "zero";
  - missing estimate or bounds: "not estimable".
- **H4 part of the alternative template.** "drier conditions in the {second | third} trimester of pregnancy were associated with lower birth weight", or "…in the second and third trimesters…". The trimester is written as an ordinal; "T2" and "T3" are not used.
- **Plural.** "consistent with the pre-specified hypothesis" is used when one part meets its rule, and "…hypotheses" when more than one does.
- **`B.discussion.verdict_summary`.** "By the pre-specified rules, the birth-weight results were {B.h4.verdict} H4 (trimester terms meeting the rule: {B.h4.verdict_terms}), the height-for-age results at 24–59 months were {B.h5.verdict} H5, and the results {did not support | supported} the pre-specified El Niño × district-sensitivity hypothesis." The last clause is "did not support" when `B.ensoxs.verdict` is `not consistent with`, and "supported" when it is `consistent with`.

**3. ENSO-trim CIs are approximate (I-2).** The intervals of `B.supp.trim_oni.*` should be read as approximate, for three reasons.
- The fixest eigenvalue fix was applied to the two-way clustered vcov, which was not positive definite (`vcov_fixed = TRUE` in `supp_models.csv`).
- `T1_oni`, `T2_oni` and `T3_oni` are functions of the birth month alone, common to all districts. They are collinear: r ≈ 0.8 between adjacent trimesters.
- The terms rest on about 12 ENSO events, the same caveat as the mediation.

The fix was also applied to S-waz0023, S-whz0023, S-underweight0023, S-wasted2459 and S-wasted0023, so seven fits in all.

**4. Other display and pipeline changes.**
- **Unit labels (M-6).** `unit_of` labels 1,000-day dryness terms (`dry1000_*`, `b10r_dry1000_*`, `dist_dry1000_*`) "per +1 SPI unit drier". This corrects the `unit` display column of the h5haz and stunting-LPM rows in `sensitivity.csv`, and therefore in `fig4_sensitivity.csv`.
- **Figure data not cached (M-4).** `figure_data` is computed on every run and is not cached, so it always reflects the current step tables. It still records the data md5 for `provenance.csv`.
- **`hw57` alignment key (M-5).** The key adds the woman's `caseid`, which is compared through `mother_id`, the dense rank of (`round`, `caseid`) that Task 5 assigns. With `bidx`, it identifies each child uniquely, and the step stops unless (`round`, `caseid`, `bidx`) is unique.

## Amendment 2026-10-23 (ruling Q3): covariate distribution and missingness table

This amendment was committed before the step that produces its keys was run. It adds one descriptive step, `covariates` (`11_covariates.R`), one output file, `papers/B/output/supp_covariates.csv`, and the keys below. Every existing key and value is unchanged. The table is descriptive and carries no verdict. It serves the journal's request for a table of all covariates with their missing counts (supplement Table S2).

**Samples.**
- `cov`: all children in `analysis_child.parquet`, 492,547 rows.
- `covbw`: the BW analytic sample (`s_bw`: singleton births with `bw_g`) before complete-case.

Both are reported by round (`nfhs4`, `nfhs5`).

**Variables and levels, as modelled** (`prep_b`, `COV`).

| `<var>` | Source | `<level>` |
|---|---|---|
| `sex` | `b4` | `male` (1), `female` (2) |
| `birth_order` | `bord` band | `1`, `2`, `3`, `4plus` |
| `mother_age` | `v012` band | `15_19`, `20_24`, `25_29`, `30_34`, `35plus` |
| `education` | `v106` | `none` (0), `primary` (1), `secondary` (2), `higher` (3) |
| `wealth` | `v190` | `poorest`, `poorer`, `middle`, `richer`, `richest` (1–5) |
| `residence` | `v025` | `urban` (1), `rural` (2) |
| `caste` | `s116` | `sc` (1), `st` (2), `obc` (3), `none` (4: none of these) |
| `religion` | `v130` | `hindu` (1), `muslim` (2), `christian` (3), `sikh` (4), `buddhist` (5), `other` (6, 7, 8, 9 and 96) |

`missing` is the children whose raw code is NA or outside the levels.
- For `caste` this is code 8 (don't know) or NA. It is the "missing" category that the models keep as its own level.
- For every other variable, any missing child would be counted here. In this dataset only `caste` has missing values, so the other `missing` cells are true zeros.
- The levels and `missing` partition each sample × round.

**Measures.**
- **Level:** the weighted share of the sample × round, with `weight` and the denominator including `missing`. Across a variable, the level shares and the weighted share of `missing` sum to 100%.
- **`missing`:** the unweighted n and the unweighted %.
- **`round`:** the unweighted n per round in the pooled sample and its share (equal to `sample_flow.csv`).

**Disclosure rule** (ruling Q2, applied to this table).
- A cell with 1–29 children is suppressed.
- A true zero (no child, e.g. no missing value) describes no child and is published as 0. In the CSV it has `n` = NA and `zero = TRUE`, because `check_aggregate_table` admits only n ≥ 30 or NA in count columns.
- **Complementary suppression.** The cells of one sample × round × variable sum to a published total (`sample_flow.csv`, `descriptive_by_round.csv`). So if the suppressed cells of a group include no cell of n ≥ 30, the smallest cell with n ≥ 30 in that group is also suppressed. Differencing then gives only the sum of at least two cells.
- **Nested samples.** `covbw` is a subset of `cov`, so a published pair of cells gives the non-BW count by difference. If that difference is 1–29, the `covbw` cell is suppressed, and the complementary rule is applied again until nothing changes.
- **Weighted shares.** A suppressed cell's weighted share is NA. Published weighted shares give at most a suppressed cell's weighted share, never its n, as in `fig1_timeline.csv`.
- A test checks that a synthetic group with one small cell, and a nested pair with a small difference, cannot have that n recovered from the published columns and totals.

**`supp_covariates.csv`** has the columns `sample`, `round`, `variable`, `level`, `n`, `zero`, `weighted_pct`, `unweighted_pct` and `suppressed`. A suppressed cell has `n`, `weighted_pct` and `unweighted_pct` NA.

### 10. Covariate keys (`11_covariates.R`; descriptive, no verdict)

`{cov,covbw}` is the sample (all children / BW analytic sample).

| Key | Meaning | Display |
|---|---|---|
| `B.supp.{cov,covbw}.sex.{male,female}.{nfhs4,nfhs5}` | Weighted % of the sample × round in the level | pct1, or `suppressed (n < 30)` |
| `B.supp.{cov,covbw}.birth_order.{1,2,3,4plus}.{nfhs4,nfhs5}` | as above | pct1 |
| `B.supp.{cov,covbw}.mother_age.{15_19,20_24,25_29,30_34,35plus}.{nfhs4,nfhs5}` | as above | pct1 |
| `B.supp.{cov,covbw}.education.{none,primary,secondary,higher}.{nfhs4,nfhs5}` | as above | pct1 |
| `B.supp.{cov,covbw}.wealth.{poorest,poorer,middle,richer,richest}.{nfhs4,nfhs5}` | as above | pct1 |
| `B.supp.{cov,covbw}.residence.{urban,rural}.{nfhs4,nfhs5}` | as above | pct1 |
| `B.supp.{cov,covbw}.caste.{sc,st,obc,none}.{nfhs4,nfhs5}` | as above | pct1 |
| `B.supp.{cov,covbw}.religion.{hindu,muslim,christian,sikh,buddhist,other}.{nfhs4,nfhs5}` | as above | pct1 |
| `B.supp.{cov,covbw}.{sex,birth_order,mother_age,education,wealth,residence,caste,religion}.missing.{nfhs4,nfhs5}` | Missing. `value` is the unweighted n missing. | `1,234 (0.5%)` (unweighted %), `0 (0.0%)`, or `suppressed (n < 30)` |
| `B.supp.{cov,covbw}.round.{nfhs4,nfhs5}.all` | Children in the round, in the pooled sample. `value` is n. | `259,627 (52.7%)` (unweighted % of the pooled sample) |

A suppressed key has `value` null.

### Note 2026-10-23 (rulings Q3 and Q2 minors m-1, m-3): disclosure details and one template word

This note was committed before the files were regenerated. It changes no value key, and the current text keys keep their wording.
- **Q3, complementary rule, made precise.** A *primary* suppression is either a 1–29 cell or a `covbw` cell suppressed by the nested rule; the latter can have n ≥ 30, but the group total would otherwise give it back. A group that contains a primary suppression must also contain at least one suppressed cell that is not primary and has n ≥ 30. If it does not, the smallest unsuppressed cell with n ≥ 30 is suppressed. If no such cell is left, at least two non-zero cells must be suppressed, or the step stops. This replaces the sentence "if the suppressed cells of a group include no cell of n ≥ 30…" of the Q3 amendment.
- **m-1.** The `missing_share_*` columns of `missingness_tn_birthyear.csv` are rounded to 3 decimals. Unrounded fractions next to an n band could be inverted into counts.
- **m-3.** In `B.discussion.verdict_summary`, the alternative wording for `B.ensoxs.verdict` = `consistent with` is "…and the results were consistent with the pre-specified El Niño × district-sensitivity hypothesis." The current output ("did not support") is unchanged.

## Amendment 2026-10-23 (ruling Q5): exposure distributions and rescaled precision

This amendment was committed before the step that produces its keys was run. It adds one cached step, `precision` (`12_precision.R`), two output files, and the keys below. Every existing key and value is unchanged. These are descriptive restatements of existing estimates and carry no verdict. The manuscript uses them to say which effect sizes the null results exclude.

**Analytic samples.**
- Each exposure's distribution is taken in the estimation sample of its main model: the rows fixest actually used (`fixest::obs`), after complete cases and singleton or perfect-fit removal.
- The step refits each main model with the same formula, data, weights and cluster as in 02–05, through `fit_fe` and `fe_table`.
- Every run checks that the refit estimate, lo and hi equal the `primary_models.csv` row to a relative tolerance of 1e-8, and the run stops if they do not. So the samples are those behind the published estimates. The n of every sample equals that model's `.n` key.

| `<exposure>` | Variable | Main model (sample) |
|---|---|---|
| `T1_dry`, `T2_dry`, `T3_dry` | `T*_dry` (SPI units, drier positive) | H4-bw |
| `T1_tmax`, `T2_tmax`, `T3_tmax` | `T*_tmax_anom` (°C) | H4-bw |
| `dry1000_2459` | `dry1000_2459` | H5-haz 24–59 |
| `dry1000_0023` | `dry1000_0023` | H5-haz 0–23 |
| `recent1_dry`, `recent1_tmax` | `recent1_dry`, `recent1_tmax_anom` | H6-diarr |
| `recent12_dry`, `recent12_tmax` | `recent12_dry`, `recent12_tmax_anom` | H6-anaem |
| `preg_oni` | `preg_oni` (°C) | ENSO-total ONI, birth weight |
| `sens_dry_bw` | `sens_dry` (% JJAS rain loss per °C ONI) | ENSO×S birth weight |
| `sens_dry_haz2459` | `sens_dry` | ENSO×S HAZ 24–59 |

**Statistics.**
- Weighted percentiles (`weight`): the smallest value whose cumulative normalised weight is ≥ p, i.e. the inverse of the weighted ECDF.
- `sd` is the weighted SD, as in `wmean_sd`.
- `min` and `max` are unweighted. They are exposure grid values, not counts.
- `n` is the sample size, which equals the model's `.n` key.
- Exposures in SPI units or °C are displayed with 2 decimals, and `sens_dry` with 1, using U+2212.

**Rescaled precision.** For a primary null estimate β with 95% CI (lo, hi) and a positive scale s, the rescaled estimate is s·β, with CI (s·lo, s·hi). For a negative s the bounds are swapped. Display: `fmt_ci` with 1 decimal for grams and 3 for HAZ.
- **H4-bw `T2_dry` and `T3_dry`:** grams. s = the IQR (p75 − p25) or SD of that trimester's dryness in the H4-bw sample.
- **H5-haz 24–59 `dry`:** HAZ. s = the IQR or SD of `dry1000_2459` in the H5-haz 24–59 sample.
- **ENSO×S, birth weight and HAZ 24–59.** The coefficient is the El Niño-developing-pregnancy difference per unit of `sens_dry`.
  - `per_iqr` and `per_sd` scale it by the IQR or SD of `sens_dry` in that model's sample. They give the change in the El Niño difference across the IQR (or 1 SD) of district sensitivity.
  - `at_p50` and `at_p90` scale it by the p50 or p90 of `sens_dry` in that sample. They give the implied El Niño-pregnancy difference at a median or a 90th-percentile sensitive child's district.
  - The interaction term itself (`preg_elnino_dev × sens_dry`) is zero outside El Niño pregnancies, so its own IQR is not a useful scale.

**Outputs.**
- `papers/B/output/supp_expdist.csv`: one row per exposure, with `exposure`, `model`, `n`, `p5`, `p10`, `p25`, `p50`, `p75`, `p90`, `p95`, `min`, `max`, `sd` and `iqr`.
- `papers/B/output/supp_precision.csv`: one row per rescaled key, with `model`, `term`, `scale`, `s`, `est`, `lo`, `hi` and `unit`.

### 11. Exposure distributions and precision keys (`12_precision.R`; descriptive, no verdict)

| Key | Meaning | Display |
|---|---|---|
| `B.supp.expdist.{T1_dry,T2_dry,T3_dry,T1_tmax,T2_tmax,T3_tmax,dry1000_2459,dry1000_0023,recent1_dry,recent1_tmax,recent12_dry,recent12_tmax,preg_oni,sens_dry_bw,sens_dry_haz2459}.{p5,p10,p25,p50,p75,p90,p95,min,max,sd}` | Weighted percentile, unweighted min or max, or weighted SD of the exposure in its main model's sample | number: 2 decimals (1 for `sens_dry_*`) |
| `B.supp.expdist.{T1_dry,T2_dry,T3_dry,T1_tmax,T2_tmax,T3_tmax,dry1000_2459,dry1000_0023,recent1_dry,recent1_tmax,recent12_dry,recent12_tmax,preg_oni,sens_dry_bw,sens_dry_haz2459}.n` | Sample size (equals the model's `.n` key) | int |
| `B.precision.h4bw.{T2_dry,T3_dry}.{per_iqr,per_sd}` | H4-bw: grams per IQR or per SD of the trimester dryness; rescaled 95% CI | ci1 |
| `B.precision.h5haz2459.dry1000.{per_iqr,per_sd}` | H5-haz 24–59: HAZ per IQR or per SD of `dry1000_2459`; rescaled 95% CI | ci3 |
| `B.precision.ensoxs.bw.{per_iqr,per_sd,at_p50,at_p90}` | ENSO×S birth weight, grams: the interaction × IQR or SD of `sens_dry`, or the El Niño-pregnancy difference at `sens_dry` p50 or p90; scaled 95% CI | ci1 |
| `B.precision.ensoxs.haz2459.{per_iqr,per_sd,at_p50,at_p90}` | ENSO×S HAZ 24–59: the same, in HAZ | ci3 |

## Amendment 2026-10-23 (ruling Q9): reader-facing displays, p-value format, one sample-flow key

This amendment was committed before the files were regenerated. **No estimate, CI, p value or verdict value changes**, and no key's `value` changes, except that the two generated sentences (section 9) are text keys whose value is their display. Only `display` strings change, plus the three keys added in item 4.

The changes are applied at assembly, by a Paper-B display layer (`models/B/R/display.R`, called by `write_b_results`), so the cached steps, including the mediation bootstrap, are reused unchanged. `models/common/R/results.R` is not changed; Papers A and C use it.

1. **p values.**
   - The display of every `B.*` key ending in `.p` uses `fmt_p_sig(p)`, as the journal asks for 1–2 significant digits:
     - "<0.001" when p < 0.001;
     - ">0.99" when p rounds to 1 at 2 significant digits;
     - otherwise 2 significant digits, e.g. 0.27, 0.049, 0.0021, 0.50.
   - A non-estimable p keeps "not estimable".
   - This replaces the `p` row of the Conventions table (`fmt_p`, 3 decimals) for Paper B displays. The `value` stays the unrounded p.
2. **Text keys in prose.** The `value` keeps the machine code; the `display` becomes prose.

   **`B.h4.verdict_terms`**

   | Code | Display |
   |---|---|
   | `neither` | "neither the second- nor the third-trimester term met the rule" |
   | `T2` | "the second-trimester term met the rule" |
   | `T3` | "the third-trimester term met the rule" |
   | `T2 and T3` | "both the second- and the third-trimester terms met the rule" |

   **`B.med.*.flag`**

   | Code | Display |
   |---|---|
   | `ok` | "bootstrap completed (at most 10% of resamples failed)" |
   | `too_many_failed_fits` | "bootstrap not usable (more than 10% of resamples failed)" |

   **`B.med.*.prop_flag`, and the `not estimated (<prop_flag>)` display of `B.med.*.prop`**

   | Code | Display |
   |---|---|
   | `ok` | "estimated" |
   | `total_effect_ci_includes_0` | "not estimated (total-association CI includes 0)" |
   | `opposite_signs` | "not estimated (indirect and total associations of opposite sign)" |
   | `nonfinite_resample_proportion` | "not estimated (non-finite proportion in a bootstrap resample)" |
   | `not_assessed` | "not assessed (too many failed bootstrap fits)" |

   The `.prop` displays read the same text as their flag, e.g. "not estimated (total-association CI includes 0)".

   **`B.sens.differ`.** Each `spec/model` pair is written as "<analysis>, <model>", and pairs are joined by "; ". `none` is displayed as "none".
   - Analyses:

     | Code | Prose |
     |---|---|
     | `card` | card-recorded birth weight |
     | `lastbirth` | last births |
     | `unweighted` | unweighted |
     | `nonmover` | non-movers |
     | `precovid` | pre-COVID interviews |
     | `motherfe` | mother fixed effects |
     | `shiftm1` | conception shifted 1 month earlier |
     | `shiftp1` | conception shifted 1 month later |
     | `district` | district linkage |
     | `buf10` | 10 km rural buffer |
     | `noflag` | flagged exposures excluded |
     | `nogaz` | gazetteer clusters excluded |
     | `cluster1` | district-only clustering |

   - Models:

     | Code | Prose |
     |---|---|
     | `h4bw` | trimester dryness, birth weight |
     | `h5haz` | 1,000-day dryness, HAZ 24–59 months |
     | `xsbw` | El Niño × sensitivity, birth weight |
     | `xshaz` | El Niño × sensitivity, HAZ 24–59 months |

   - Example: `nonmover/xsbw` → "non-movers, El Niño × sensitivity, birth weight".

   `B.sens.placebo_flag` and `B.sens.ensoplacebo_flag` are already prose and are unchanged.
3. **`B.discussion.verdict_summary`.** The parenthesis "(trimester terms meeting the rule: {B.h4.verdict_terms})" becomes "({display of B.h4.verdict_terms})", e.g. "(neither the second- nor the third-trimester term met the rule)". The rest of the section 9 template is unchanged.
4. **Sample flow (Table S1; manuscript review M-16).** These keys are added in `01_descriptive.R`.

| Key | Meaning | Display |
|---|---|---|
| `B.desc.n_singleton_no_bw{,.nfhs4,.nfhs5}` | Singleton births (`twin == FALSE`) in the analysis dataset without a recorded birth weight (`bw_g` missing: not weighed, don't know, or > 6,000 g); total, then by round | int, or `suppressed (n < 30)` with value null when n is 1–29 |

The two round counts sum to the total. Differencing is safe because no part is suppressed: if either round count were 1–29, both round keys and the total would be suppressed together.


## Note 2026-10-23 (ruling Q12): description of `at_p50` and `at_p90` (description only)

Wording only. No key, value, display, model or verdict changes. In the Q5 amendment above, read the `at_p50` and `at_p90` bullet as follows:

- `at_p50` and `at_p90` scale the interaction by the p50 or p90 of `sens_dry` in that sample. They give the implied El Niño-pregnancy difference at a median or a 90th-percentile sensitive child's district, **relative to a district with `sens_dry` = 0 (El Niño does not alter rainfall)**. The birth year-month fixed effects absorb the El Niño main effect, so the contrast is always with such a district. The percentiles are of children, not of districts.

The `description` strings of these keys in `results.json` are written by `models/B/12_precision.R` (the `add()` calls for `at_p50` and `at_p90`), not by this file. This note does not regenerate them. To carry the reference into `results.json`, edit those two description templates and rerun the results step (`Rscript B/run_all.R --use-cache`). The values, displays and verdicts would not change.
