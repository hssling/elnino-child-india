# Paper B: additional analyses (plan fixed 25 November 2026) [journal name redacted from the public copy]

**Status.** These analyses were planned on 25 November 2026, after the primary results were known and after the first journal declined the manuscript at pre-review screening. They are **post hoc**. They do not change the pre-specified hypotheses, models or verdicts (design specification 2026-10-22; results contract `papers/B/results_keys.md`), and the manuscript will report them as additional analyses, labelled as such. This plan is committed before any of the models below is fitted.

**Data.** The NFHS-4 and NFHS-5 children's recode (the same 492,547 children), the same climate linkage, the same survey weights and the same clustering (`region_id` + `birth_ym`). Mother's height (`v438`, cm × 10) is newly extracted from the raw children's recode files and joined by row alignment with `kr_harmonised.parquet` (the same checks as `load_hw57`). Only aggregate estimates leave the machine.

## A1. Region-specific seasonality and the heat placebo

**Question.** Post-birth heat, a placebo exposure, was associated with birth weight: 15.3 g (95% CI 6.1 to 24.4) per 1 °C. Is this explained by region-specific seasonal patterns of births that the district and birth year-month fixed effects do not absorb?

**Models.**
- The birth-weight model with the placebo windows, as in the sensitivity analysis `S-placebo-h4bw`: T1–T3 dryness, T1–T3 maximum-temperature anomaly, post-birth dryness and post-birth heat, the primary covariates, and fixed effects `region_id` + `birth_ym`.
- (a) Adds **state × birth calendar month** fixed effects (`state_gadm` × `birth_month_of_year`).
- (b) Instead adds **district × birth calendar month** fixed effects (`region_id` × `birth_month_of_year`).
- The same two extensions of the primary H4 model (without the placebo terms) and of the primary H5 model for height-for-age at 24–59 months (1,000-day dryness).

**Decision rule (fixed now).** In model (a), the heat placebo failure is judged **explained by seasonality** if the post-birth heat estimate's 95% CI includes 0 **and** its absolute value is below half the primary placebo estimate (7.65 g per 1 °C).
- If it is explained, the trimester heat estimates from model (a) are reported as the interpretable heat estimates, still observational.
- Otherwise the heat terms remain uninterpretable, and the paper says so.
- Model (b) is a check on (a).

## A2. Negative-control outcome: mother's height

**Rationale.** Adult height is set long before the index pregnancy, so the weather in that pregnancy cannot affect it. An association would reveal confounding, for example by socioeconomic or regional patterns that also track the weather of birth timing.

**Model.**
- Outcome: mother's height in cm. Valid values are 130.0–200.0 cm; other values and DHS special codes are set to missing.
- Regressed on T1–T3 dryness and T1–T3 heat, with the primary covariates and fixed effects, in the birth-weight sample. Also under model (a) of A1.
- A second model regresses it on the post-birth placebo windows.

**Reporting.**
- Estimates per interquartile range (IQR) of each exposure, in cm.
- Converted into an implied birth-weight bias: the height coefficient multiplied by the within-model association of mother's height with birth weight (g per cm, from the primary birth-weight model with mother's height added).

**Rule.** A term **flags confounding** if its 95% CI excludes 0 **and** the implied birth-weight bias exceeds 2 g per IQR. This is the bias size that would matter against the effects in A3.

## A3. Smallest effects of interest

The paper will state which effects the data rule out, not only which estimates include 0. These benchmarks are chosen now, from sources external to these data:
- **Birth weight:** 15 g per IQR of trimester dryness. This is the smallest well-estimated whole-pregnancy environmental association in a very large birth cohort: about 15 g lower term birth weight with heat above the 90th percentile in 29.6 million US births (manuscript reference 12).
- **Height-for-age at 24–59 months:** 0.05 z per IQR of 1,000-day dryness. This is near the 0.034 z per SPEI unit for in-utero rainfall reported in NFHS-4 (manuscript reference 21), rounded up.
- **El Niño × dry-sensitivity at the 90th-percentile district:** the same external benchmarks, 15 g of birth weight and 0.05 z of height-for-age. These benchmarks are not chosen from this paper's own estimates.

**Rule.** An effect of the benchmark size or larger is **ruled out** (in the harmful direction) if the 95% CI's lower limit lies above minus the benchmark. The harmful direction is lower birth weight or lower height-for-age with more dryness. Each estimate is reported as ruled out or not. No estimate is re-estimated for this step; it uses the primary estimates and their rescaling in Table S3.

## Outputs and labelling

- Keys `B.x.*`, written by `models/B/13_additional.R` to `papers/B/envres/results_additional.json`, with a traceability table.
- Aggregates only, cells under the minimum size suppressed.
- The manuscript reports A1–A3 in a separate "Additional analyses (not pre-specified)" subsection, with this plan cited by commit.
