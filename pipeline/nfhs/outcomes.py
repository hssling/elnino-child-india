"""Clean NFHS KR child outcomes and check national estimates against the reports.

Units (verified, spec §12): m19 grams (9996/9998 missing, > 6000 excluded);
hw56 haemoglobin g/dL x 10, altitude-adjusted (994-999 missing); hw70/71/72
z-scores x 100 (>= 9996 missing); h11 2 = yes, 0 = no, 8 = missing.
Missing is NaN, never 0.

Two rules depart from DHS convention on purpose (per the spec): h11 = 8
(don't know) is NaN, whereas DHS counts it as "no" in the denominator; and
m19 > 6000 g is excluded, whereas DHS keeps every m19 < 9996. Checked on the
real KR files, each moves the national estimate by <= 0.03 pp (diarrhoea
9.19/7.34 vs DHS-style 9.17/7.33; LBW 18.22/18.27 vs 18.21/18.26, NFHS-4/5).

hb_plausible (4.0 <= hb_adj <= 18.0 g/dL) is a sensitivity flag only; the
primary anaemia rule uses every hw56 < 994.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from pipeline.paths import EXTERNAL, ROOT

SUBSET_PARQUET = EXTERNAL / "nfhs_subset" / "kr_harmonised.parquet"
REPORT_CSV = ROOT / "pipeline" / "reference" / "nfhs_report_estimates.csv"
TOLERANCE_PP = 1.0

BW_MAX_G = 6000
HB_MISSING_MIN = 994          # hw56 994-999 are missing codes
HB_PLAUSIBLE_G_DL = (4.0, 18.0)  # hb_plausible sensitivity flag, inclusive
Z_MISSING_MIN = 9996          # hw70/71/72 >= 9996 are missing / flagged codes
# WHO flag limits (z): outside -> NaN
WHO_LIMITS = {"haz": (-6.0, 6.0), "waz": (-6.0, 5.0), "whz": (-5.0, 5.0)}
Z_SOURCE = {"haz": "hw70", "waz": "hw71", "whz": "hw72"}
Z_BINARY = {"haz": "stunted", "waz": "underweight", "whz": "wasted"}


def _binary(cond: pd.Series, valid: pd.Series) -> pd.Series:
    """1.0/0.0 where valid, NaN elsewhere."""
    return cond.astype(float).where(valid)


def clean_outcomes(kr: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of `kr` with cleaned outcome, timing and weight columns added."""
    df = kr.copy()

    m19 = pd.to_numeric(df["m19"], errors="coerce").astype(float)
    df["bw_g"] = m19.where(m19 <= BW_MAX_G)  # 9996/9998 and > 6000 g -> NaN
    df["bw_card"] = (df["m19a"] == 1).fillna(False).astype(bool)
    df["lbw"] = _binary(df["bw_g"] < 2500, df["bw_g"].notna())

    for z, src in Z_SOURCE.items():
        raw = pd.to_numeric(df[src], errors="coerce").astype(float)
        val = (raw / 100.0).where(raw < Z_MISSING_MIN)
        lo, hi = WHO_LIMITS[z]
        df[z] = val.where((val >= lo) & (val <= hi))
        df[Z_BINARY[z]] = _binary(df[z] < -2, df[z].notna())

    hw56 = pd.to_numeric(df["hw56"], errors="coerce").astype(float)
    hb_raw = hw56.where(hw56 < HB_MISSING_MIN)
    df["hb_adj"] = hb_raw / 10.0
    df["anaemia_any"] = _binary(hb_raw < 110, hb_raw.notna())
    df["anaemia_modsev"] = _binary(hb_raw < 100, hb_raw.notna())
    # sensitivity flag only (Task 7); the primary rule above keeps every hw56 < 994
    df["hb_plausible"] = _binary(df["hb_adj"].between(HB_PLAUSIBLE_G_DL[0], HB_PLAUSIBLE_G_DL[1]),
                                 df["hb_adj"].notna())

    h11 = pd.to_numeric(df["h11"], errors="coerce").astype(float)
    df["diarrhoea"] = h11.map({0.0: 0.0, 1.0: 1.0, 2.0: 1.0})  # 8 / NaN -> NaN

    df["birth_cmc"] = df["b3"]
    df["interview_cmc"] = df["v008"]
    df["age_m"] = df["v008"] - df["b3"]
    df["alive"] = (df["b5"] == 1).astype(bool)
    df["twin"] = (df["b0"] > 0).astype(bool)

    w = df["v005"].astype(float) / 1e6
    df["weight"] = w * df.groupby("round")["round"].transform("size") / w.groupby(df["round"]).transform("sum")
    return df


def _age_measured(df: pd.DataFrame) -> pd.Series:
    """Age in months at measurement (hw1) where recorded, else from CMC dates."""
    if "hw1" in df.columns:
        return df["hw1"].astype(float).fillna(df["age_m"].astype(float))
    return df["age_m"].astype(float)


# indicator -> (outcome column, denominator mask builder). Denominators follow
# the NFHS national report tables (see nfhs_report_estimates.csv).
_INDICATORS = {
    "stunting": ("stunted", lambda d: _age_measured(d) < 60),
    "underweight": ("underweight", lambda d: _age_measured(d) < 60),
    "wasting": ("wasted", lambda d: _age_measured(d) < 60),
    "anaemia_any": ("anaemia_any", lambda d: _age_measured(d).between(6, 59)),
    "diarrhoea": ("diarrhoea", lambda d: d["alive"] & (d["age_m"] < 60)),
    "lbw": ("lbw", lambda d: d["age_m"] < 60),
}


def national_estimates(df: pd.DataFrame) -> pd.DataFrame:
    """Weighted national prevalence (%) per round and indicator.

    Columns: round, indicator, value_pct, n_unweighted, n_weighted.
    """
    rows = []
    for rnd, d in df.groupby("round", sort=True):
        for ind, (col, mask_fn) in _INDICATORS.items():
            sel = mask_fn(d) & d[col].notna()
            w = d.loc[sel, "weight"]
            y = d.loc[sel, col]
            value = 100.0 * float((w * y).sum() / w.sum()) if len(w) else np.nan
            rows.append({"round": rnd, "indicator": ind, "value_pct": value,
                         "n_unweighted": int(sel.sum()), "n_weighted": float(w.sum())})
    return pd.DataFrame(rows)


def check_against_report(est: pd.DataFrame, report_csv=REPORT_CSV,
                         tol: float = TOLERANCE_PP) -> list[str]:
    """Failures where the estimate differs from the report by more than `tol` points.

    A report row with no matching estimate is also a failure.
    """
    ref = pd.read_csv(Path(report_csv))
    merged = ref.merge(est[["round", "indicator", "value_pct"]], on=["round", "indicator"],
                       how="left", suffixes=("_report", "_est"))
    failures = []
    for r in merged.itertuples(index=False):
        if pd.isna(r.value_pct_est):
            failures.append(f"{r.round} {r.indicator}: no estimate (report {r.value_pct_report:.1f})")
            continue
        diff = r.value_pct_est - r.value_pct_report
        if abs(diff) > tol:
            failures.append(f"{r.round} {r.indicator}: estimate {r.value_pct_est:.2f} vs report "
                            f"{r.value_pct_report:.1f} (diff {diff:+.2f} pp > {tol})")
    return failures
