"""Climate and ENSO exposure windows for NFHS children (Paper B task 4).

Months are DHS century-month codes (CMC = (year - 1900) * 12 + month): a
child's birth month is b3 (`birth_cmc`), the interview month v008
(`interview_cmc`).

Windows are inclusive month offsets (lo, hi) from the birth month, or from the
interview month when given as (lo, hi, "interview"). A window value is NaN when
any of its months is missing or masked for that variable, falls outside the
metrics series (before 1981 or after 2025 for the district table), or lies after
the interview month. The interview month itself counts as inside the observable
period: a window ending at the interview month is kept. There are no partial
means. Missing rain makes only the SPI-3 window NaN; temperature windows are
judged on temperature alone, because QC masks blank rain only.

Island proxies (ruling P7): Andaman & Nicobar and Lakshadweep have no IMD 0.25 deg
rain cells, so their rain is a fallback proxy (Lakshadweep: Kerala coast); their
rain is set NaN before SPI (SPI windows and sens_beta NaN). Where their
temperature is a fallback (A&N: nearest mainland 1 deg cell), tmax_mean,
heatwave_days and tmax_ge40_days are NaN. This is done here, not in
qc_reviewed.csv (shared with Papers A/C). Mainland fallbacks (Saiha/Lawngtlai,
Anjaw/Dibang Valley) are kept and flagged: `<window>_fallback_rain_any`.

`window_exposure` takes any metrics table keyed by `key` + year + month. The
district linkage uses key "region_id"; the cluster-GPS linkage (Task 4b) uses
key ["round", "v001"] with cluster-cell metrics.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from pipeline.climate_annual import BASE, QC_REVIEWED, apply_qc_to_month, mask_months, spi
from pipeline.paths import DERIVED, month_table

ENSO_MONTHLY_CSV = DERIVED / "enso_monthly.csv"
ENSO_ANNUAL_CSV = DERIVED / "enso_annual.csv"

INTERVIEW = "interview"
# Conception is approximated as birth month - 9 (spec §4); shift sensitivity = pass shifted windows.
PREG_WINDOWS = {"T1": (-9, -7), "T2": (-6, -4), "T3": (-3, -1)}
WINDOWS = {
    **PREG_WINDOWS,
    "inf0_11": (0, 11),
    "inf12_23": (12, 23),
    "placebo": (1, 9),
    "recent1": (-1, -1, INTERVIEW),
    "recent3": (-3, -1, INTERVIEW),
    "recent12": (-12, -1, INTERVIEW),
}
# metric column -> (output suffix, aggregation over the window)
WINDOW_VARIABLES = {"spi3": ("spi3", "mean"), "tmax_anom": ("tmax_anom", "mean"),
                    "heatwave_days": ("hw_days", "sum")}
# boolean metric column -> output suffix; True if any window month has it (only if the column exists)
WINDOW_FLAGS = {"fallback_rain": "fallback_rain_any"}

SENS_YEARS = (1981, 2009)     # before the birth window: no outcome-period data
SENS_MIN_YEARS = 20           # valid JJAS years needed for a district slope
JJAS = (6, 7, 8, 9)
EXCLUDED_STATES = ("Andaman and Nicobar", "Lakshadweep")   # island proxy climate (ruling P7)
TEMP_COLUMNS = ("tmax_mean", "tmin_mean", "heatwave_days", "tmax_ge40_days")


def cmc(year, month):
    """DHS century-month code."""
    return (np.asarray(year) - 1900) * 12 + np.asarray(month)


def _as_list(key) -> list[str]:
    return [key] if isinstance(key, str) else list(key)


def _parse_window(name, w):
    if len(w) == 2:
        lo, hi, anchor = int(w[0]), int(w[1]), "birth"
    elif len(w) == 3 and w[2] in ("birth", INTERVIEW):
        lo, hi, anchor = int(w[0]), int(w[1]), w[2]
    else:
        raise ValueError(f"window {name!r}: expected (lo, hi) or (lo, hi, 'birth'|'interview'), got {w!r}")
    if lo > hi:
        raise ValueError(f"window {name!r}: lo {lo} > hi {hi}")
    return lo, hi, anchor


# ---------------------------------------------------------------- district monthly metrics

def _mask_island_proxies(month: pd.DataFrame, excluded_states) -> pd.DataFrame:
    """Rain NaN in `excluded_states`; temperature NaN there where fallback_temp. Adds island_proxy."""
    month = month.copy()
    if len(excluded_states) == 0:
        month["island_proxy"] = False
        return month
    if "state_name" not in month:
        raise KeyError("district month table lacks state_name (needed to mask island proxy climate)")
    island = month["state_name"].isin(excluded_states).to_numpy()
    month["island_proxy"] = island
    month.loc[island, [c for c in ("rain_mm", "rainy_days") if c in month]] = np.nan
    if "fallback_temp" in month:
        ft = island & month["fallback_temp"].astype("boolean").fillna(False).to_numpy(dtype=bool)
        month.loc[ft, [c for c in TEMP_COLUMNS if c in month]] = np.nan
    return month


def district_monthly_metrics(month: pd.DataFrame | None = None, base=BASE, qc_path=QC_REVIEWED,
                             excluded_states=EXCLUDED_STATES) -> pd.DataFrame:
    """District-month exposure metrics keyed region_id, year, month.

    `month` defaults to data/derived/climate_district_month.csv.gz. Rain is masked with
    apply_qc_to_month(month, "district") (qc_path=None: no masks, for tests). Each district's
    series is reindexed to a contiguous month range, so a missing month is NaN and never bridged.
    Columns: spi3 (SPI-3 of the masked rain, gamma fits on `base`), tmax_anom (tmax_mean minus the
    district's `base` calendar-month mean), heatwave_days, tmax_ge40_days, rain_mm (masked),
    qc_masked, fallback_rain, fallback_temp (where present), island_proxy.
    Island proxy climate in `excluded_states` is blanked before SPI and the normals (see module doc)."""
    if month is None:
        month = pd.read_csv(month_table("district", DERIVED), dtype={"region_id": str})
    month = month.copy()
    month["region_id"] = month["region_id"].astype(str)
    if month.duplicated(["region_id", "year", "month"]).any():
        raise ValueError("duplicate district-month rows")
    if qc_path is None:
        month = mask_months(month, None)
    else:
        month = apply_qc_to_month(month, "district", path=qc_path)
    month = _mask_island_proxies(month, excluded_states)

    keep = ["tmax_mean", "heatwave_days", "tmax_ge40_days", "rain_mm", "qc_masked",
            "fallback_rain", "fallback_temp", "island_proxy"]
    keep = [c for c in keep if c in month]
    parts = []
    for rid, g in month.groupby("region_id", sort=True):
        idx = pd.PeriodIndex.from_fields(year=g.year.to_numpy(), month=g.month.to_numpy(), freq="M")
        full = pd.period_range(idx.min(), idx.max(), freq="M")
        g = g.set_index(idx)[keep].reindex(full)
        out = pd.DataFrame({"region_id": rid, "year": full.year, "month": full.month})
        out["spi3"] = spi(g["rain_mm"].astype(float), 3, base).to_numpy()
        for c in keep:
            out[c] = g[c].to_numpy()
        parts.append(out)
    met = pd.concat(parts, ignore_index=True)
    for c in ("qc_masked", "fallback_rain", "fallback_temp", "island_proxy"):
        if c in met:
            met[c] = met[c].astype("boolean").fillna(False).astype(bool)
    in_base = met.year.between(*base)
    normal = met[in_base].groupby(["region_id", "month"])["tmax_mean"].mean().rename("_tmax_normal")
    met = met.join(normal, on=["region_id", "month"])
    met["tmax_anom"] = met["tmax_mean"] - met["_tmax_normal"]
    cols = ["region_id", "year", "month", "spi3", "tmax_anom", "heatwave_days", "tmax_ge40_days",
            "rain_mm", "tmax_mean"] + [c for c in ("qc_masked", "fallback_rain", "fallback_temp", "island_proxy")
                                         if c in met]
    return met[cols].sort_values(["region_id", "year", "month"]).reset_index(drop=True)


# ---------------------------------------------------------------- windows

def _gather(values: np.ndarray, row: np.ndarray, col0: np.ndarray, lo: int, hi: int,
            limit: np.ndarray | None, how: str, fill: float = np.nan) -> np.ndarray:
    """Aggregate values[row, col0 + k] for k in lo..hi; NaN if any month is NaN, out of range or
    (col0 + k) > limit. row < 0 means no series."""
    n_rows, n_cols = values.shape
    stack = np.empty((len(row), hi - lo + 1))
    for j, k in enumerate(range(lo, hi + 1)):
        c = col0 + k
        ok = (row >= 0) & (c >= 0) & (c < n_cols)
        if limit is not None:
            ok &= c <= limit
        v = np.full(len(row), fill)
        v[ok] = values[row[ok], c[ok]]
        stack[:, j] = v
    if how == "any":
        return np.nan_to_num(stack, nan=0.0).max(axis=1) > 0
    return stack.mean(axis=1) if how == "mean" else stack.sum(axis=1)   # NaN propagates


def window_exposure(children: pd.DataFrame, metrics: pd.DataFrame, windows: dict = WINDOWS,
                    key="region_id", variables: dict = WINDOW_VARIABLES, flags: dict = WINDOW_FLAGS,
                    birth_col="birth_cmc", interview_col="interview_cmc") -> pd.DataFrame:
    """Window exposures per child, indexed like `children`, columns <window>_<suffix>.

    `metrics` is keyed by `key` (a column name or list, e.g. ["round", "v001"]) + year + month
    and holds the columns named in `variables` ({metric column: (suffix, "mean"|"sum")}).
    A window is NaN if any month in it is missing/NaN, outside the metrics series, after the
    child's interview month (the interview month itself is inside), or if the child's key is not in
    `metrics`. `flags` ({boolean metric column: suffix}) adds <window>_<suffix>, True if any window
    month has the flag; flags whose column is absent from `metrics` are not produced."""
    keys = _as_list(key)
    for c in keys + ["year", "month"]:
        if c not in metrics:
            raise KeyError(f"metrics lacks key column {c!r}")
    for c in keys + [birth_col, interview_col]:
        if c not in children:
            raise KeyError(f"children lacks column {c!r}")
    parsed = {name: _parse_window(name, w) for name, w in windows.items()}
    if metrics.duplicated(keys + ["year", "month"]).any():
        raise ValueError("metrics has duplicate key-year-month rows")
    for v, (_, how) in variables.items():
        if v not in metrics:
            raise KeyError(f"metrics lacks variable {v!r}")
        if how not in ("mean", "sum"):
            raise ValueError(f"variable {v!r}: aggregation must be mean or sum")

    birth = pd.to_numeric(children[birth_col], errors="coerce").to_numpy(dtype=float)
    interview = pd.to_numeric(children[interview_col], errors="coerce").to_numpy(dtype=float)

    # dense (series x month) arrays, restricted to the month range the children need
    m_cmc = cmc(metrics["year"].to_numpy(dtype=int), metrics["month"].to_numpy(dtype=int))
    need_lo, need_hi = [], []
    for lo, hi, anchor in parsed.values():
        a = birth if anchor == "birth" else interview
        if np.isfinite(a).any():
            need_lo.append(np.nanmin(a) + lo)
            need_hi.append(np.nanmax(a) + hi)
    out = pd.DataFrame(index=children.index)
    first = int(max(m_cmc.min(), min(need_lo))) if need_lo else int(m_cmc.min())
    last = int(min(m_cmc.max(), max(need_hi))) if need_hi else int(m_cmc.max())
    sel = (m_cmc >= first) & (m_cmc <= last)
    mkeys = pd.MultiIndex.from_frame(metrics[keys]) if len(keys) > 1 else pd.Index(metrics[keys[0]])
    series = mkeys[sel].unique()
    srow = series.get_indexer(mkeys[sel])
    scol = (m_cmc[sel] - first).astype(int)
    n_cols = max(last - first + 1, 0)

    ckeys = pd.MultiIndex.from_frame(children[keys]) if len(keys) > 1 else pd.Index(children[keys[0]])
    crow = series.get_indexer(ckeys) if len(series) else np.full(len(children), -1)
    crow = np.where(np.isfinite(birth) & np.isfinite(interview), crow, -1)
    birth_i = np.where(np.isfinite(birth), birth, 0).astype(int)
    interview_i = np.where(np.isfinite(interview), interview, 0).astype(int)
    limit = interview_i - first

    for v, (suffix, how) in variables.items():
        values = np.full((len(series), n_cols), np.nan)
        if n_cols and len(series):
            values[srow, scol] = metrics.loc[sel, v].to_numpy(dtype=float)
        for name, (lo, hi, anchor) in parsed.items():
            col0 = (birth_i if anchor == "birth" else interview_i) - first
            out[f"{name}_{suffix}"] = _gather(values, crow, col0, lo, hi, limit, how)
    for f, suffix in flags.items():
        if f not in metrics:
            continue
        values = np.zeros((len(series), n_cols))
        if n_cols and len(series):
            values[srow, scol] = metrics.loc[sel, f].astype("boolean").fillna(False).to_numpy(dtype=float)
        for name, (lo, hi, anchor) in parsed.items():
            col0 = (birth_i if anchor == "birth" else interview_i) - first
            out[f"{name}_{suffix}"] = _gather(values, crow, col0, lo, hi, limit, "any", fill=0.0)
    return out


# ---------------------------------------------------------------- ENSO

def enso_exposure(children: pd.DataFrame, enso_monthly: pd.DataFrame, enso_annual: pd.DataFrame,
                  windows: dict = PREG_WINDOWS, birth_col="birth_cmc") -> pd.DataFrame:
    """Pregnancy ENSO exposures per child (indexed like `children`).

    preg_oni / preg_roni / preg_nino34: mean over all pregnancy months (union of `windows`, by
    default T1-T3 = birth - 9 .. birth - 1) of oni_anom / roni / nino34_anom (by central month).
    ONI and RONI are 3-month running means: the value at centre month m averages SST over m-1..m+1,
    so T3_oni (centres birth-3..birth-1) draws on SST up to the birth month itself and preg_oni on
    birth-10..birth. preg_nino34 is the monthly (unsmoothed) Nino-3.4 anomaly: the strict-window
    alternative.
    preg_oni_pos = max(preg_oni, 0), preg_oni_neg = min(preg_oni, 0); <T>_oni per trimester;
    preg_elnino_dev: any pregnancy month in Jun-Sep of an el_nino_developing year (nullable
    boolean; NA if a pregnancy year is absent from enso_annual and no month is positive)."""
    parsed = {name: _parse_window(name, w) for name, w in windows.items()}
    if any(a != "birth" for _, _, a in parsed.values()):
        raise ValueError("pregnancy windows must be birth-relative")
    offsets = sorted({k for lo, hi, _ in parsed.values() for k in range(lo, hi + 1)})
    if offsets != list(range(offsets[0], offsets[-1] + 1)):
        raise ValueError("pregnancy windows must be contiguous")
    lo, hi = offsets[0], offsets[-1]

    em = enso_monthly.copy()
    em["series"] = "enso"
    kids = pd.DataFrame({"series": "enso", birth_col: children[birth_col].to_numpy(),
                         "_far": cmc(9999, 12)}, index=children.index)
    names = {"oni_anom": "oni", "roni": "roni", "nino34_anom": "nino34"}
    variables = {c: (s, "mean") for c, s in names.items() if c in em}
    win = {"preg": (lo, hi), **{n: (a, b) for n, (a, b, _) in parsed.items()}}
    w = window_exposure(kids, em, win, key="series", variables=variables,
                        birth_col=birth_col, interview_col="_far")
    out = pd.DataFrame(index=children.index)
    out["preg_oni"] = w["preg_oni"]
    out["preg_oni_pos"] = w["preg_oni"].clip(lower=0)
    out["preg_oni_neg"] = w["preg_oni"].clip(upper=0)
    for c in ("roni", "nino34"):
        if f"preg_{c}" in w:
            out[f"preg_{c}"] = w[f"preg_{c}"]
    for n in parsed:
        out[f"{n}_oni"] = w[f"{n}_oni"]

    dev = dict(zip(enso_annual["year"].astype(int), enso_annual["el_nino_developing"].astype(bool)))
    b = pd.to_numeric(children[birth_col], errors="coerce").to_numpy(dtype=float)
    any_true = np.zeros(len(b), dtype=bool)
    unknown = ~np.isfinite(b)
    for k in range(lo, hi + 1):
        c = np.where(np.isfinite(b), b, 0).astype(int) + k
        yr, mo = (c - 1) // 12 + 1900, (c - 1) % 12 + 1
        d = pd.Series(yr).map(dev).astype("boolean")
        unknown |= d.isna().to_numpy() & np.isfinite(b)
        any_true |= np.isin(mo, JJAS) & d.fillna(False).to_numpy(dtype=bool) & np.isfinite(b)
    res = pd.array(any_true, dtype="boolean")
    res[unknown & ~any_true] = pd.NA
    out["preg_elnino_dev"] = res
    return out


# ---------------------------------------------------------------- ENSO-rainfall sensitivity

def district_enso_sensitivity(metrics: pd.DataFrame, enso_annual: pd.DataFrame,
                              years=SENS_YEARS, min_years=SENS_MIN_YEARS) -> pd.DataFrame:
    """Per-district OLS slope of JJAS rainfall % departure on oni_jjas, `years` only (1981-2009).

    JJAS rain = sum of the four (masked) monthly rain_mm values; a year with any of them missing is
    dropped. The departure is from the district's own mean JJAS rain over the same valid years, so
    no data from 2010 on enters. Returns region_id, sens_beta (% per ONI unit), sens_se (classical
    OLS), sens_n; NaN where fewer than `min_years` valid years remain."""
    m = metrics[metrics.year.between(*years) & metrics.month.isin(JJAS)]
    jj = m.groupby(["region_id", "year"])["rain_mm"].agg(["sum", "count"])
    jj = jj["sum"].where(jj["count"] == len(JJAS)).rename("jjas").reset_index()
    full = metrics[["region_id"]].drop_duplicates()
    oni = enso_annual.set_index("year")["oni_jjas"]
    rows = []
    for rid in full.region_id:
        g = jj[(jj.region_id == rid)].copy()
        g["oni"] = g.year.map(oni)
        g = g.dropna(subset=["jjas", "oni"])
        n = len(g)
        beta = se = np.nan
        if n >= min_years:
            normal = g.jjas.mean()
            if normal > 0:
                y = 100 * (g.jjas.to_numpy() - normal) / normal
                x = g.oni.to_numpy()
                xc = x - x.mean()
                sxx = (xc ** 2).sum()
                if sxx > 0:
                    beta = (xc * (y - y.mean())).sum() / sxx
                    resid = y - y.mean() - beta * xc
                    se = np.sqrt((resid ** 2).sum() / (n - 2) / sxx)
        rows.append({"region_id": rid, "sens_beta": beta, "sens_se": se, "sens_n": n})
    return pd.DataFrame(rows, columns=["region_id", "sens_beta", "sens_se", "sens_n"])


# ---------------------------------------------------------------- real-data coverage (aggregates)

def load_children(kr_path=None) -> pd.DataFrame:
    """Cleaned NFHS children with region_id from the Task 3 crosswalk (controlled data, local only)."""
    from pipeline.nfhs import districts
    from pipeline.nfhs.outcomes import SUBSET_PARQUET, clean_outcomes
    kr = pd.read_parquet(kr_path or SUBSET_PARQUET)
    kids = clean_outcomes(kr)
    kids["region_id"] = None
    for rnd in kids["round"].unique():
        xw = districts.crosswalk(rnd)
        sel = kids["round"] == rnd
        kids.loc[sel, "region_id"] = kids.loc[sel, "district"].map(lambda c: xw.get(int(c)) if pd.notna(c) else None)
    return kids


def coverage_check(kids: pd.DataFrame | None = None, metrics: pd.DataFrame | None = None) -> pd.DataFrame:
    """Share of children with non-missing window exposures (aggregates only), by round x group:
    level "group" (A&N/Lakshadweep vs outside) and level "state" (GADM state of the region_id)."""
    kids = load_children() if kids is None else kids
    metrics = district_monthly_metrics() if metrics is None else metrics
    w = window_exposure(kids, metrics, WINDOWS, flags={})
    month = pd.read_csv(month_table("district", DERIVED), usecols=["region_id", "state_name"],
                        dtype={"region_id": str}).drop_duplicates("region_id")
    state = kids["region_id"].map(dict(zip(month.region_id, month.state_name))).fillna("(unlinked)")
    grp = np.where(state.isin(EXCLUDED_STATES), "A&N/Lakshadweep", "outside A&N/Lakshadweep")
    df = w.notna().assign(round=kids["round"].to_numpy(), linked=kids["region_id"].notna().to_numpy())
    parts = []
    for level, g in (("group", grp), ("state", state.to_numpy())):
        d = df.assign(group=g)
        agg = d.groupby(["round", "group"]).mean()
        agg.insert(0, "n_children", d.groupby(["round", "group"]).size())
        parts.append(agg.reset_index().assign(level=level))
    out = pd.concat(parts, ignore_index=True)
    return out[["level"] + [c for c in out.columns if c != "level"]]
