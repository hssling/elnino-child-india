"""Paper B analysis dataset: one row per NFHS-4/5 KR child (task 5).

`build()` joins, for every child of the Task 1 KR subset:

* outcomes, timing and weights from `outcomes.clean_outcomes` (hb_plausible kept);
* identifiers for the models: `mother_id` (int, dense rank of (round, caseid); caseid is the DHS
  concatenation of v001, v002, v003, so this is the (round, v001, v002, v003) mother), `birth_ym`
  and `interview_ym` (int YYYYMM from the CMC: CMC 1385 -> 201505). caseid itself is dropped;
* regions (ruling P8): `region_id` = the point-in-polygon GADM district of the child's cluster
  (gps.cluster_table; used for district fixed effects and clustering in the primary analysis);
  `region_id_xw` = the Task 3 crosswalk district of the KR district code after XW_CODE_OVERRIDES;
* PRIMARY exposures (ruling P4): cluster-GPS windows from gps.cluster_window_exposure (plain names
  <window>_spi3 / _tmax_anom / _hw_days / _spi3_buf10 and flags <window>_fallback_rain_any,
  _fallback_temp_any, _mis_district_any, _qc_masked_any, _qc_cell_masked_any);
* SENSITIVITY exposures: district-linked windows from exposure.window_exposure on region_id_xw,
  prefixed `dist_` (e.g. dist_T2_spi3);
* ENSO pregnancy exposures (exposure.enso_exposure) and the district ENSO-rainfall sensitivity
  (sens_beta / sens_se on region_id; sens_beta_xw / sens_se_xw on region_id_xw);
* covariates (COVARIATES) and cluster attributes (source, urban, region_method, state_gadm,
  mis_district, island_proxy, boundary_flag).

Island children (Andaman & Nicobar, Lakshadweep) stay in with NaN climate (ruling P7). Rows are
excluded only by EXCLUSION_RULES, each counted per round in the exclusions table.

Outputs: the git-ignored data/external/nfhs_subset/analysis_child.parquet (+ the exclusion counts
next to it) and two committed aggregates in data/derived/nfhs/: child_counts_by_state_round.csv
(n and weighted prevalence by state and round; prevalence suppressed where n < MIN_CELL_N) and
district_enso_sensitivity.csv (climate only). No record-level or cluster-level row, cluster id,
caseid or coordinate is ever written outside the git-ignored directory.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from pipeline.nfhs import districts
from pipeline.nfhs.exposure import EXCLUDED_STATES, WINDOWS, enso_exposure, window_exposure
from pipeline.paths import DERIVED, EXTERNAL

SUBSET_DIR = EXTERNAL / "nfhs_subset"
KR_PARQUET = SUBSET_DIR / "kr_harmonised.parquet"
ANALYSIS_PARQUET = SUBSET_DIR / "analysis_child.parquet"
EXCLUSIONS_CSV = SUBSET_DIR / "analysis_child_exclusions.csv"
DERIVED_NFHS = DERIVED / "nfhs"
COUNTS_CSV = DERIVED_NFHS / "child_counts_by_state_round.csv"
SENS_CSV = DERIVED_NFHS / "district_enso_sensitivity.csv"
ENSO_MONTHLY_CSV = DERIVED / "enso_monthly.csv"
ENSO_ANNUAL_CSV = DERIVED / "enso_annual.csv"

SIZE_LIMIT_MB = 150.0
PIP_MIN_SHARE = 0.95
T2_SPI_MIN_COVERAGE = 0.90
MIN_CELL_N = 30

# Ruling P8: in NFHS-4 the KR/GE district codes 39 (labelled Sangrur) and 53 (labelled Shahid
# Bhagat Singh Nagar) are swapped on the DHS side: the GPS points of all ~71 of their clusters fall
# in the other district (task-4b report), and NFHS-5 labels the same two codes the other way round.
# Applied only when looking up the Task 3 crosswalk (region_id_xw and the dist_* exposures); the
# crosswalk CSVs are not edited and the point-in-polygon region_id does not need it.
XW_CODE_OVERRIDES: dict[tuple[str, int], int] = {("nfhs4", 39): 53, ("nfhs4", 53): 39}

# Spec §3 / §12 covariates present in the Task 1 subset (b4 sex, bord birth order, v012 mother's
# age at interview (added by ruling P11), v106 education, v190 wealth, v025 residence, s116 caste,
# v130 religion; v104 years at residence for the non-movers sensitivity).
COVARIATES = ["b4", "bord", "v012", "v106", "v190", "v025", "s116", "v130", "v104"]
# Spec covariates NOT in the Task 1 subset (reported, not extracted: §12 excludes them from the
# total-effect models).
MISSING_COVARIATES = {
    "m14": "antenatal visits (spec §3; excluded from total-effect models by §12) - not extracted",
    "m15": "place of delivery (spec §3; excluded from total-effect models by §12) - not extracted",
}

EXCLUSION_RULES = {
    "missing_birth_or_interview_cmc": "b3 or v008 missing: no birth month / interview month",
    "invalid_weight": "v005 missing or <= 0",
    "cluster_not_in_ge": "(round, v001) absent from the GE cluster table: no region or exposure",
}

PRIMARY_EXTRA_FLAGS = {"qc_masked": "qc_masked_any", "qc_cell_masked": "qc_cell_masked_any"}

OUTCOMES = ["bw_g", "bw_card", "m19a", "lbw", "haz", "waz", "whz", "stunted", "underweight", "wasted",
            "hb_adj", "anaemia_any", "anaemia_modsev", "hb_plausible", "diarrhoea"]
# Task 7 sensitivity inputs (ruling P11): m19a (card-only birth weight; bw_card = m19a == 1), bidx
# (bidx 1 = last birth), v104 (non-movers), v008 (= interview_cmc; post-COVID cut), mother_id,
# hw1 (age in months at measurement), birth_cmc (= b3), birth_month_of_year (from the CMC).
TIMING = ["birth_cmc", "interview_cmc", "v008", "birth_ym", "interview_ym", "birth_month_of_year",
          "b1", "b2", "v016", "age_m",
          "hw1", "sphase", "alive", "twin", "b0", "b5"]

COUNT_INDICATORS = ["lbw", "stunting", "underweight", "wasting", "anaemia_any", "anaemia_modsev",
                    "diarrhoea"]


# ---------------------------------------------------------------- small helpers

def cmc_to_ym(c):
    """DHS century-month code -> int YYYYMM (CMC 1385 -> 201505). Series in -> Int32 Series."""
    if np.isscalar(c):
        c = int(c)
        return ((c - 1) // 12 + 1900) * 100 + (c - 1) % 12 + 1
    s = pd.to_numeric(pd.Series(c), errors="coerce")
    ok = s.notna()
    out = pd.Series(pd.NA, index=s.index, dtype="Int32")
    v = s[ok].astype(np.int64)
    out[ok] = (((v - 1) // 12 + 1900) * 100 + (v - 1) % 12 + 1).astype(np.int32).to_numpy()
    return out


def xw_district_code(round_id: str, code: int) -> int:
    """KR district code to look up in the Task 3 crosswalk (XW_CODE_OVERRIDES applied)."""
    return XW_CODE_OVERRIDES.get((round_id, int(code)), int(code))


def _lookup_xw(rounds, codes, table: dict) -> list:
    out = []
    for r, c in zip(rounds, codes):
        out.append(table.get(r, {}).get(xw_district_code(r, c)) if pd.notna(c) else None)
    return out


def _check_rounds(kids: pd.DataFrame, table: dict, what: str) -> None:
    unknown = sorted(set(kids["round"].astype(str)) - set(table))
    if unknown:
        raise ValueError(f"no {what} table for round(s) {unknown}")


def _check_complete(kids: pd.DataFrame, values, what: str) -> None:
    """Raise if any child lacks `what` (lists the (round, KR district code) pairs, no records)."""
    miss = pd.Series(values, index=kids.index).isna()
    if miss.any():
        pairs = sorted({(r, int(c)) if pd.notna(c) else (r, None)
                        for r, c in zip(kids.loc[miss, "round"], kids.loc[miss, "district"])},
                       key=str)
        raise ValueError(f"{int(miss.sum())} children have no {what}; (round, district code): {pairs[:10]}")


def _exclude(kr: pd.DataFrame, clusters: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Apply EXCLUSION_RULES in order; return (kept rows, counts per round x reason)."""
    keep = pd.Series(True, index=kr.index)
    v005 = pd.to_numeric(kr["v005"], errors="coerce")
    known = pd.MultiIndex.from_arrays([clusters["round"].astype(str), clusters["v001"].astype(np.int64)])
    kid_key = pd.MultiIndex.from_arrays([kr["round"].astype(str),
                                         pd.to_numeric(kr["v001"], errors="coerce").fillna(-1).astype(np.int64)])
    masks = {
        "missing_birth_or_interview_cmc": kr["b3"].isna() | kr["v008"].isna(),
        "invalid_weight": v005.isna() | (v005 <= 0),
        "cluster_not_in_ge": pd.Series(~kid_key.isin(known), index=kr.index),
    }
    rows = []
    for reason in EXCLUSION_RULES:
        hit = masks[reason] & keep
        for rnd in sorted(kr["round"].unique()):
            rows.append({"round": rnd, "reason": reason,
                         "n_excluded": int((hit & (kr["round"] == rnd)).sum()),
                         "rule": EXCLUSION_RULES[reason]})
        keep &= ~hit
    return kr[keep].copy(), pd.DataFrame(rows, columns=["round", "reason", "n_excluded", "rule"])


# ---------------------------------------------------------------- assembly (pure core)

def assemble(kr: pd.DataFrame, clusters: pd.DataFrame, cluster_monthly: pd.DataFrame,
             district_metrics: pd.DataFrame, enso_monthly: pd.DataFrame, enso_annual: pd.DataFrame,
             sens: pd.DataFrame, xw: dict, xw_flags: dict | None = None
             ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(analysis frame, exclusion counts). `xw` / `xw_flags` are {round: {KR district code:
    region_id / boundary_flag}} from the Task 3 crosswalk (overrides applied here)."""
    from pipeline.nfhs.gps import CLUSTER_KEY, CLUSTER_WINDOW_FLAGS, cluster_window_exposure
    from pipeline.nfhs.outcomes import clean_outcomes

    kept, exclusions = _exclude(kr, clusters)
    kids = clean_outcomes(kept.reset_index(drop=True))     # weights rescaled on the kept rows
    kids["v001"] = kids["v001"].astype(np.int64)
    kids["birth_ym"] = cmc_to_ym(kids["birth_cmc"])
    kids["interview_ym"] = cmc_to_ym(kids["interview_cmc"])
    kids["birth_month_of_year"] = (kids["birth_ym"] % 100).astype("Int8")
    kids["mother_id"] = kids.groupby(["round", "caseid"], sort=True).ngroup().astype(np.int32)

    # cluster attributes (primary region = point in polygon)
    cl = clusters.copy()
    cl["v001"] = cl["v001"].astype(np.int64)
    cl["mis_district"] = (cl["source"].astype(str) == "MIS") | cl["lat"].isna() | cl["lon"].isna()
    cl["island_proxy"] = cl["state_name"].isin(EXCLUDED_STATES)
    cl = cl.rename(columns={"state_name": "state_gadm"})
    attrs = ["region_id", "region_method", "source", "urban", "state_gadm", "mis_district", "island_proxy"]
    kids = kids.merge(cl[CLUSTER_KEY + attrs], on=CLUSTER_KEY, how="left", validate="many_to_one")

    # crosswalk region (sensitivity) with the P8 override
    _check_rounds(kids, xw, "crosswalk")
    kids["region_id_xw"] = _lookup_xw(kids["round"], kids["district"], xw)
    _check_complete(kids, kids["region_id_xw"], "crosswalk district (region_id_xw)")
    if xw_flags is None:
        kids["boundary_flag"] = ""
    else:
        _check_rounds(kids, xw_flags, "boundary_flag")
        flag = pd.Series(_lookup_xw(kids["round"], kids["district"], xw_flags), index=kids.index)
        _check_complete(kids, flag, "boundary_flag")
        kids["boundary_flag"] = flag.astype(str)

    # exposures
    cm = cluster_monthly.copy()
    cm["v001"] = cm["v001"].astype(np.int64)
    prim = cluster_window_exposure(kids, cm, WINDOWS, flags={**CLUSTER_WINDOW_FLAGS, **PRIMARY_EXTRA_FLAGS})
    dm = district_metrics.copy()
    dm["region_id"] = dm["region_id"].astype(str)
    dist_flags = {"fallback_rain": "fallback_rain_any"}
    dist = window_exposure(kids.assign(region_id=kids["region_id_xw"]), dm, WINDOWS, key="region_id",
                           flags=dist_flags).add_prefix("dist_")
    enso = enso_exposure(kids, enso_monthly, enso_annual)

    s = sens.set_index(sens["region_id"].astype(str))
    if s.index.duplicated().any():
        raise ValueError("duplicate region_id rows in the district ENSO sensitivity table")
    island = kids["island_proxy"].astype(bool)
    for suffix, col in (("", "region_id"), ("_xw", "region_id_xw")):
        absent = sorted(set(kids[col].astype(str)) - set(s.index))
        if absent:
            raise ValueError(f"{len(absent)} {col} value(s) missing from the district ENSO sensitivity "
                             f"table, e.g. {absent[:5]}")
        for v in ("sens_beta", "sens_se"):
            kids[v + suffix] = kids[col].astype(str).map(s[v])
        bad = kids["sens_beta" + suffix].isna() & ~island
        if bad.any():
            regions = sorted(kids.loc[bad, col].astype(str).unique())
            raise ValueError(f"sens_beta{suffix} NaN outside the island districts (ruling P7) for "
                             f"{col} {regions[:5]}")

    for frame in (prim, dist, enso):
        for c in frame.columns:
            if frame[c].dtype == float:
                frame[c] = frame[c].astype(np.float32)

    ids = ["round", "v001", "bidx", "mother_id", "v024", "district"]
    regions = ["region_id", "region_method", "region_id_xw", "boundary_flag", "source", "urban",
               "state_gadm", "mis_district", "island_proxy"]
    base = kids[ids + regions + TIMING + OUTCOMES + COVARIATES + ["v005", "weight",
                "sens_beta", "sens_se", "sens_beta_xw", "sens_se_xw"]]
    out = pd.concat([base, prim, dist, enso], axis=1)
    for c in ("mis_district", "island_proxy", "urban"):
        out[c] = out[c].astype("boolean")
    return out, exclusions


# ---------------------------------------------------------------- committed aggregates

def child_counts(df: pd.DataFrame, state_labels: dict) -> pd.DataFrame:
    """Unweighted n and weighted prevalence (%) by round x state (v024) x indicator.

    Denominators follow outcomes.national_estimates (same KR rules). prev_pct is NaN where the
    indicator's n < MIN_CELL_N. Aggregates only."""
    from pipeline.nfhs.outcomes import _INDICATORS
    ind = dict(_INDICATORS)
    ind["anaemia_modsev"] = ("anaemia_modsev", ind["anaemia_any"][1])
    rows = []
    for (rnd, st), d in df.groupby(["round", "v024"], sort=True):
        for name in COUNT_INDICATORS:
            col, mask_fn = ind[name]
            sel = mask_fn(d) & d[col].notna()
            n = int(sel.sum())
            w, y = d.loc[sel, "weight"], d.loc[sel, col]
            prev = 100.0 * float((w * y).sum() / w.sum()) if n >= MIN_CELL_N and w.sum() > 0 else np.nan
            rows.append({"round": rnd, "state_code": int(st),
                         "state_name": state_labels.get(rnd, {}).get(int(st)),
                         "n_children": len(d), "indicator": name, "n": n, "prev_pct": prev})
    out = pd.DataFrame(rows, columns=["round", "state_code", "state_name", "n_children", "indicator",
                                      "n", "prev_pct"])
    out["prev_pct"] = out["prev_pct"].round(2)
    return out


# ---------------------------------------------------------------- gates

def gates(df: pd.DataFrame, n_task1: int, exclusions: pd.DataFrame, path) -> list:
    """Task 5 gates (GateResult): rows, point-in-polygon share, T2 SPI coverage, file size."""
    from pipeline.validate import GateResult
    res = []
    n_exc = int(exclusions["n_excluded"].sum())
    ok = len(df) == n_task1 - n_exc
    res.append(GateResult("nfhs_rows", ok, f"analysis rows {len(df):,} = Task 1 {n_task1:,} - exclusions "
                                           f"{n_exc:,}" + ("" if ok else " (MISMATCH)")))
    pip = df.groupby("round")["region_method"].apply(lambda s: float((s == "pip").mean()))
    res.append(GateResult("nfhs_pip_share", bool((pip >= PIP_MIN_SHARE).all()),
                          "point-in-polygon region_id share " + ", ".join(f"{r} {v:.4f}" for r, v in pip.items())
                          + f" (>= {PIP_MIN_SHARE})"))
    main = df[~df["state_gadm"].isin(EXCLUDED_STATES)]
    cov = main.groupby("round")["T2_spi3"].apply(lambda s: float(s.notna().mean()))
    res.append(GateResult("nfhs_t2_spi_coverage", bool((cov >= T2_SPI_MIN_COVERAGE).all()),
                          "T2 cluster SPI-3 non-missing outside A&N/Lakshadweep "
                          + ", ".join(f"{r} {v:.4f}" for r, v in cov.items()) + f" (>= {T2_SPI_MIN_COVERAGE})"))
    path = Path(path)
    mb = path.stat().st_size / 2**20
    total = sum(p.stat().st_size for p in path.parent.iterdir() if p.is_file()) / 2**20
    res.append(GateResult("nfhs_size", mb < SIZE_LIMIT_MB and total < SIZE_LIMIT_MB,
                          f"{path.name} {mb:.1f} MB; {path.parent.name}/ total {total:.1f} MB (< {SIZE_LIMIT_MB:.0f})"))
    return res


# ---------------------------------------------------------------- real build

def missing_sources() -> list[str]:
    """What the build needs and is absent: both KR zips (labels; the KR subset is extracted from
    them when missing) and either both GE zips or the two cluster caches. [] = buildable."""
    from pipeline.nfhs import extract, gps
    out = [f"KR zip {Path(z).name} not found" for z, _ in extract.KR_SOURCES.values() if not Path(z).exists()]
    if not all(Path(z).exists() for z, _ in gps.GE_SOURCES.values()):
        out += [f"GE zips absent and no cache {Path(p).name}"
                for p in (gps.CLUSTERS_PARQUET, gps.CLUSTER_MONTHLY_PARQUET) if not Path(p).exists()]
    return out


def sources_available() -> bool:
    return not missing_sources()


def ensure_kr(kr_path=KR_PARQUET) -> Path:
    """The Task 1 KR subset; extracted from the KR zips (local read, not a download) if missing."""
    from pipeline.nfhs import extract
    kr_path = Path(kr_path)
    if not kr_path.exists():
        kr_path.parent.mkdir(parents=True, exist_ok=True)
        extract.extract_all(kr_path.parent)
    return kr_path


def _sensitivity_table(sens: pd.DataFrame) -> pd.DataFrame:
    from pipeline.nfhs.gps import _gadm_table
    g = _gadm_table()
    out = sens.merge(g[["region_id", "region_name", "state_name"]], on="region_id", how="left")
    return out[["region_id", "region_name", "state_name", "sens_beta", "sens_se", "sens_n"]]


def build(out_dir=SUBSET_DIR, derived_dir=DERIVED_NFHS, kr_path=KR_PARQUET) -> dict[str, Path]:
    """Real build (controlled data; local only). Returns the written paths. Tests pass tmp dirs so
    the committed aggregates are never rewritten as a side effect."""
    from pipeline.nfhs import extract, gps
    from pipeline.nfhs.exposure import district_enso_sensitivity, district_monthly_metrics

    out_dir, derived_dir = Path(out_dir), Path(derived_dir)
    kr = pd.read_parquet(ensure_kr(kr_path))
    clusters = gps.load_or_build_clusters()
    monthly = gps.load_or_build_monthly()          # built with both rounds (union QC mask, P9)
    dm = district_monthly_metrics()
    em, ea = pd.read_csv(ENSO_MONTHLY_CSV), pd.read_csv(ENSO_ANNUAL_CSV)
    sens = district_enso_sensitivity(dm, ea)
    rounds = sorted(kr["round"].unique())
    xw = {r: districts.crosswalk(r) for r in rounds}
    flags = {r: districts.boundary_flags(r) for r in rounds}

    df, exc = assemble(kr, clusters, monthly, dm, em, ea, sens, xw, flags)
    paths = {"analysis": out_dir / ANALYSIS_PARQUET.name, "exclusions": out_dir / EXCLUSIONS_CSV.name,
             "counts": derived_dir / COUNTS_CSV.name, "sensitivity": derived_dir / SENS_CSV.name}
    out_dir.mkdir(parents=True, exist_ok=True)
    derived_dir.mkdir(parents=True, exist_ok=True)
    df.to_parquet(paths["analysis"], index=False, compression="zstd")
    exc.to_csv(paths["exclusions"], index=False)
    labels = {r: extract.state_labels(r) for r in rounds}
    child_counts(df, labels).to_csv(paths["counts"], index=False)
    _sensitivity_table(sens).to_csv(paths["sensitivity"], index=False)
    return paths


def run(out_dir=SUBSET_DIR, derived_dir=DERIVED_NFHS, kr_path=KR_PARQUET) -> list:
    """build() then gates() on the written file (run_all step `nfhs`)."""
    paths = build(out_dir, derived_dir, kr_path)
    n_task1 = len(pd.read_parquet(kr_path, columns=["round"]))
    df = pd.read_parquet(paths["analysis"], columns=["round", "region_method", "state_gadm", "T2_spi3"])
    exc = pd.read_csv(paths["exclusions"])
    return gates(df, n_task1, exc, paths["analysis"])


if __name__ == "__main__":
    for g in run():
        print(("PASS" if g.passed else "FAIL"), g.name, g.detail)
