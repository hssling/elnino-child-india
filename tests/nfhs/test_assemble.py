"""Tests for the Paper B analysis dataset (task 5). Synthetic records only, built here."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pipeline.nfhs import assemble
from pipeline.nfhs.exposure import cmc
from pipeline.nfhs.extract import HARMONISED

ROUNDS = {"nfhs4": 1, "nfhs5": 2}
YEARS = (2012, 2021)


# ---------------------------------------------------------------- birth_ym

def test_birth_ym_from_cmc():
    assert assemble.cmc_to_ym(1385) == 201505
    out = assemble.cmc_to_ym(pd.Series([1381, 1392, 1393, np.nan]))
    assert out.tolist()[:3] == [201501, 201512, 201601]
    assert pd.isna(out.iloc[3])
    assert str(out.dtype) == "Int32"


# ---------------------------------------------------------------- P8 Punjab override

def test_punjab_nfhs4_codes_39_53_are_swapped_for_the_crosswalk_only():
    assert assemble.xw_district_code("nfhs4", 39) == 53
    assert assemble.xw_district_code("nfhs4", 53) == 39
    assert assemble.xw_district_code("nfhs4", 40) == 40
    assert assemble.xw_district_code("nfhs5", 39) == 39            # NFHS-5 labels are right
    assert assemble.xw_district_code("nfhs5", 53) == 53
    assert assemble.XW_CODE_OVERRIDES == {("nfhs4", 39): 53, ("nfhs4", 53): 39}


def test_punjab_override_matches_the_committed_crosswalk_names():
    """NFHS-4 labels say 39 = Sangrur, 53 = SBS Nagar; NFHS-5 labels the same codes the other way
    round, agreeing with the NFHS-4 GPS points (task-4b report) -> the NFHS-4 labels are swapped."""
    xw = pd.read_csv(assemble.districts.CROSSWALK_CSV)
    name = {(r, c): n for r, c, n in zip(xw["round"], xw["nfhs_code"], xw["nfhs_name"])}
    assert name[("nfhs4", 39)] == "sangrur" and name[("nfhs4", 53)] == "shahid bhagat singh nagar"
    assert name[("nfhs5", 39)] == "shahid bhagat singh nagar" and name[("nfhs5", 53)] == "sangrur"


# ---------------------------------------------------------------- synthetic assembly

def _kid(rnd, v001, district, caseid, bidx, b3, v008, **kw):
    row = {c: np.nan for c in HARMONISED}
    row.update(round=rnd, v001=v001, district=district, caseid=caseid, bidx=bidx, b3=b3, v008=v008,
               v005=1_000_000, v024=3, v025=2, v104=95, v106=1, v130=1, v190=3, s116=1, b0=0, b4=1,
               b5=1, bord=1, hw1=v008 - b3, m19=3000, m19a=1, hw70=-150, hw71=-100, hw72=-50,
               hw56=105, h11=0, sphase=np.nan, v012=27)
    row.update(kw)
    return row


def _inputs():
    b = cmc(2015, 5)            # 1385
    kr = pd.DataFrame([
        _kid("nfhs4", 1, 39, "  1  1  2", 1, b, b + 6, v005=2_000_000),   # Punjab code 39 (swapped)
        _kid("nfhs4", 1, 39, "  1  1  2", 2, b - 30, b + 6),             # sibling, same mother
        _kid("nfhs4", 2, 53, "  2  5  1", 1, b, b + 6, m19=9996),        # Punjab code 53 (swapped)
        _kid("nfhs4", 3, 70, "  3  1  1", 1, b, b + 6),                  # island cluster
        _kid("nfhs4", 99, 70, " 99  1  1", 1, b, b + 6),                 # cluster absent from GE
        _kid("nfhs5", 1, 39, "  1  1  2", 1, b + 48, b + 60),            # same v001 / caseid, other round
        _kid("nfhs5", 4, 53, "  4  2  2", 1, b + 48, b + 60),            # MIS cluster
    ])
    clusters = pd.DataFrame({
        "round": ["nfhs4", "nfhs4", "nfhs4", "nfhs5", "nfhs5"],
        "v001": [1, 2, 3, 1, 4],
        "lat": [30.9, 30.2, 11.6, 30.9, np.nan], "lon": [76.1, 75.8, 92.7, 76.1, np.nan],
        "urban": [False, True, False, False, True],
        "source": ["GPS", "GPS", "GPS", "GPS", "MIS"],
        "region_id": ["R_SBS", "R_SANG", "R_AN", "R_SBS", "R_SANG"],
        "district": [39, 53, 70, 39, 53],
        "region_method": ["pip", "pip", "pip", "pip", "district"],
        "region_xw": ["R_SANG", "R_SBS", "R_AN", "R_SBS", "R_SANG"],
        "boundary_flag": ["", "", "", "", ""],
        "state_name": ["Punjab", "Punjab", "Andaman and Nicobar", "Punjab", "Punjab"],
    })
    rows = []
    for r, v, s in zip(clusters["round"], clusters["v001"], clusters["state_name"]):
        val = np.nan if s == "Andaman and Nicobar" else 10 * ROUNDS[r] + v      # encodes the cluster
        for y in range(*YEARS):
            for m in range(1, 13):
                rows.append({"round": r, "v001": v, "year": y, "month": m, "spi3": val,
                             "spi3_buf10": val + 0.5 if np.isfinite(val) else np.nan,
                             "tmax_anom": 0.25, "heatwave_days": 1.0, "tmax_ge40_days": 0.0,
                             "qc_masked": bool(v == 2 and y == 2015 and m == 1), "qc_cell_masked": False,
                             "fallback_rain": False, "fallback_temp": False,
                             "island_proxy": s == "Andaman and Nicobar", "mis_district": v == 4})
    cluster_monthly = pd.DataFrame(rows)
    dvals = {"R_SBS": -1.0, "R_SANG": -2.0, "R_AN": np.nan}
    dm = pd.DataFrame([{"region_id": rid, "year": y, "month": m, "spi3": v, "tmax_anom": 0.75,
                        "heatwave_days": 2.0, "fallback_rain": False}
                       for rid, v in dvals.items() for y in range(*YEARS) for m in range(1, 13)])
    em = pd.DataFrame([{"year": y, "month": m, "oni_anom": 1.0, "roni": 0.5, "nino34_anom": 1.5}
                       for y in range(*YEARS) for m in range(1, 13)])
    ea = pd.DataFrame({"year": list(range(*YEARS)), "el_nino_developing": [y == 2014 for y in range(*YEARS)],
                       "oni_jjas": 0.0})
    sens = pd.DataFrame({"region_id": ["R_SBS", "R_SANG", "R_AN"], "sens_beta": [-5.0, -9.0, np.nan],
                         "sens_se": [1.0, 2.0, np.nan], "sens_n": [29, 29, 0]})
    xw = {"nfhs4": {39: "R_SANG", 53: "R_SBS", 70: "R_AN"}, "nfhs5": {39: "R_SBS", 53: "R_SANG"}}
    flags = {"nfhs4": {39: "", 53: "", 70: ""}, "nfhs5": {39: "", 53: "split"}}
    return dict(kr=kr, clusters=clusters, cluster_monthly=cluster_monthly, district_metrics=dm,
                enso_monthly=em, enso_annual=ea, sens=sens, xw=xw, xw_flags=flags)


@pytest.fixture(scope="module")
def built():
    return assemble.assemble(**_inputs())


def _row(df, rnd, v001, bidx=1):
    sel = df[(df["round"] == rnd) & (df["v001"] == v001) & (df["bidx"] == bidx)]
    assert len(sel) == 1
    return sel.iloc[0]


def test_rows_and_documented_exclusions(built):
    df, exc = built
    assert len(df) == 6                                   # 7 synthetic children - 1 exclusion
    e = exc.set_index(["round", "reason"])["n_excluded"]
    assert e[("nfhs4", "cluster_not_in_ge")] == 1
    assert exc["n_excluded"].sum() == 7 - len(df)
    assert (df["v001"] != 99).all()


def test_island_children_are_kept_with_nan_climate(built):
    df, _ = built
    r = _row(df, "nfhs4", 3)
    assert r["island_proxy"] and np.isnan(r["T2_spi3"]) and np.isnan(r["dist_T2_spi3"])


def test_identifiers_birth_ym_and_weights(built):
    df, _ = built
    assert "caseid" not in df.columns
    assert _row(df, "nfhs4", 1)["birth_ym"] == 201505
    assert df["birth_ym"].dtype == "Int32"
    m = df.set_index(["round", "v001", "bidx"])["mother_id"]
    assert m[("nfhs4", 1, 1)] == m[("nfhs4", 1, 2)]                 # siblings share a mother
    assert m[("nfhs4", 1, 1)] != m[("nfhs5", 1, 1)]                 # same caseid text, other round
    assert df["mother_id"].nunique() == 5 and str(df["mother_id"].dtype) == "int32"
    # v005/1e6 rescaled so each round's weights sum to its row count (after exclusions)
    assert df.groupby("round")["weight"].sum().to_dict() == pytest.approx({"nfhs4": 4.0, "nfhs5": 2.0})
    w4 = df[df["round"] == "nfhs4"].set_index(["v001", "bidx"])["weight"]
    assert w4[(1, 1)] == pytest.approx(2 * w4[(2, 1)])


def test_primary_exposure_is_cluster_gps_and_district_is_dist_prefixed(built):
    df, _ = built
    r = _row(df, "nfhs4", 1)
    assert r["T2_spi3"] == pytest.approx(11.0)                      # cluster (nfhs4, 1)
    assert r["T2_spi3_buf10"] == pytest.approx(11.5)
    assert r["T2_tmax_anom"] == pytest.approx(0.25) and r["T2_hw_days"] == pytest.approx(3.0)
    assert r["dist_T2_spi3"] == pytest.approx(-1.0)                 # R_SBS via the P8 override
    assert r["dist_T2_tmax_anom"] == pytest.approx(0.75)
    assert _row(df, "nfhs5", 1)["T2_spi3"] == pytest.approx(21.0)   # same v001, other round
    for c in ("T2_fallback_rain_any", "T2_fallback_temp_any", "T2_mis_district_any",
              "T2_qc_masked_any", "T2_qc_cell_masked_any", "dist_T2_fallback_rain_any"):
        assert c in df.columns
    assert bool(_row(df, "nfhs5", 4)["T2_mis_district_any"])
    r53 = _row(df, "nfhs4", 2)                                      # born May 2015, Jan 2015 masked
    assert bool(r53["T2_qc_masked_any"]) and bool(r53["recent12_qc_masked_any"])   # Nov-Jan; Nov-Oct
    assert not bool(r53["T3_qc_masked_any"])                         # Feb-Apr


def test_regions_with_punjab_override(built):
    df, _ = built
    r39, r53 = _row(df, "nfhs4", 1), _row(df, "nfhs4", 2)
    assert r39["region_id"] == "R_SBS" and r39["region_id_xw"] == "R_SBS"   # crosswalk says R_SANG
    assert r53["region_id"] == "R_SANG" and r53["region_id_xw"] == "R_SANG"
    assert _row(df, "nfhs5", 1)["region_id_xw"] == "R_SBS"               # NFHS-5: no override
    assert _row(df, "nfhs5", 4)["boundary_flag"] == "split"
    assert _row(df, "nfhs5", 4)["source"] == "MIS" and bool(_row(df, "nfhs5", 4)["mis_district"])


def test_sensitivity_merged_on_both_regions(built):
    df, _ = built
    r = _row(df, "nfhs4", 2)
    assert r["sens_beta"] == pytest.approx(-9.0) and r["sens_se"] == pytest.approx(2.0)
    assert r["sens_beta_xw"] == pytest.approx(-9.0)
    assert np.isnan(_row(df, "nfhs4", 3)["sens_beta"])


def test_enso_outcomes_and_covariates_are_carried(built):
    df, _ = built
    r = _row(df, "nfhs4", 1)
    assert r["preg_oni"] == pytest.approx(1.0) and r["T2_oni"] == pytest.approx(1.0)
    assert bool(r["preg_elnino_dev"])                     # pregnancy Aug 2014 - Apr 2015 hits JJAS 2014
    assert r["bw_g"] == 3000 and np.isnan(_row(df, "nfhs4", 2)["bw_g"])
    assert r["haz"] == pytest.approx(-1.5) and r["anaemia_any"] == 1.0
    for c in assemble.COVARIATES + ["hb_plausible", "interview_ym", "age_m"]:
        assert c in df.columns, c


def test_task7_sensitivity_variables_are_carried(built):
    """Ruling P11: mother's age and every variable the Task 7 sensitivity analyses use."""
    df, _ = built
    for c in ("v012", "m19a", "bw_card", "bidx", "v104", "v008", "interview_cmc", "mother_id",
              "hw1", "birth_cmc", "birth_month_of_year"):
        assert c in df.columns, c
    assert "v012" in assemble.COVARIATES and "v012" not in assemble.MISSING_COVARIATES
    r = _row(df, "nfhs4", 1)                              # born CMC 1385 = May 2015
    assert r["birth_month_of_year"] == 5 and r["v012"] == 27 and r["m19a"] == 1
    assert (df["v008"] == df["interview_cmc"]).all()
    assert str(df["birth_month_of_year"].dtype) == "Int8"


# ---------------------------------------------------------------- committed aggregates

def test_child_counts_suppress_small_cells():
    n_a, n_b = 40, 5
    df = pd.DataFrame({
        "round": "nfhs4", "v024": [1] * n_a + [2] * n_b,
        "weight": 1.0, "lbw": [1.0] * 10 + [0.0] * (n_a - 10) + [1.0] * n_b,
        "bw_g": 2000.0, "age_m": 12.0, "hw1": 12.0, "alive": True,
    })
    for c in ("stunted", "underweight", "wasted", "anaemia_any", "anaemia_modsev", "diarrhoea"):
        df[c] = 0.0
    out = assemble.child_counts(df, {"nfhs4": {1: "state a", 2: "state b"}})
    lbw = out[out["indicator"] == "lbw"].set_index("state_name")
    assert lbw.loc["state a", "n"] == n_a and lbw.loc["state a", "prev_pct"] == pytest.approx(25.0)
    assert lbw.loc["state b", "n"] == n_b and np.isnan(lbw.loc["state b", "prev_pct"])
    assert set(out.columns) == {"round", "state_code", "state_name", "n_children", "indicator", "n", "prev_pct"}
    assert not ({"v001", "caseid", "bw_g", "haz"} & set(out.columns))


def test_gates_on_synthetic_build(built, tmp_path):
    df, exc = built
    p = tmp_path / "a.parquet"
    df.to_parquet(p, index=False)
    res = {g.name: g for g in assemble.gates(df, n_task1=7, exclusions=exc, path=p)}
    assert res["nfhs_rows"].passed
    assert res["nfhs_size"].passed
    assert res["nfhs_pip_share"].passed is False          # 1 MIS of 6 children -> 83% < 95%
    bad = {g.name: g for g in assemble.gates(df, n_task1=8, exclusions=exc, path=p)}
    assert not bad["nfhs_rows"].passed


# ---------------------------------------------------------------- no silent NaN (ruling P12, M2)

def test_unknown_round_raises():
    inp = _inputs()
    del inp["xw"]["nfhs5"]
    with pytest.raises(ValueError, match="crosswalk"):
        assemble.assemble(**inp)


def test_missing_crosswalk_code_raises():
    inp = _inputs()
    del inp["xw"]["nfhs5"][53]
    with pytest.raises(ValueError, match="crosswalk district"):
        assemble.assemble(**inp)


def test_missing_boundary_flag_entry_raises():
    inp = _inputs()
    del inp["xw_flags"]["nfhs5"][53]
    with pytest.raises(ValueError, match="boundary_flag"):
        assemble.assemble(**inp)


def test_region_missing_from_sensitivity_raises():
    inp = _inputs()
    inp["sens"] = inp["sens"][inp["sens"]["region_id"] != "R_SANG"]
    with pytest.raises(ValueError, match="sensitivity"):
        assemble.assemble(**inp)


def test_mainland_nan_sensitivity_raises_but_island_nan_is_allowed():
    inp = _inputs()
    inp["sens"].loc[inp["sens"]["region_id"] == "R_SBS", "sens_beta"] = np.nan
    with pytest.raises(ValueError, match="outside the island"):
        assemble.assemble(**inp)
    df, _ = assemble.assemble(**_inputs())                  # R_AN (island) NaN is fine
    assert df.loc[df["island_proxy"].astype(bool), "sens_beta"].isna().all()


def test_duplicate_cluster_or_sensitivity_rows_raise():
    inp = _inputs()
    inp["clusters"] = pd.concat([inp["clusters"], inp["clusters"].iloc[[0]]], ignore_index=True)
    with pytest.raises(Exception):
        assemble.assemble(**inp)
    inp = _inputs()
    inp["sens"] = pd.concat([inp["sens"], inp["sens"].iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate"):
        assemble.assemble(**inp)


# ---------------------------------------------------------------- sources (ruling P12, M1)

def test_missing_sources_reports_kr_and_ge(monkeypatch, tmp_path):
    from pipeline.nfhs import extract, gps
    monkeypatch.setattr(extract, "KR_SOURCES", {"nfhs4": (tmp_path / "kr4.zip", "m"), "nfhs5": (tmp_path / "kr5.zip", "m")})
    monkeypatch.setattr(gps, "GE_SOURCES", {"nfhs4": (tmp_path / "ge4.zip", "m"), "nfhs5": (tmp_path / "ge5.zip", "m")})
    monkeypatch.setattr(gps, "CLUSTERS_PARQUET", tmp_path / "clusters.parquet")
    monkeypatch.setattr(gps, "CLUSTER_MONTHLY_PARQUET", tmp_path / "cm.parquet")
    miss = assemble.missing_sources()
    assert len(miss) == 4 and not assemble.sources_available()
    for f in ("kr4.zip", "kr5.zip", "ge4.zip", "ge5.zip"):
        (tmp_path / f).write_bytes(b"x")
    assert assemble.missing_sources() == [] and assemble.sources_available()
    (tmp_path / "ge5.zip").unlink()                         # GE absent: the caches must exist
    assert len(assemble.missing_sources()) == 2
    (tmp_path / "clusters.parquet").write_bytes(b"x")
    (tmp_path / "cm.parquet").write_bytes(b"x")
    assert assemble.missing_sources() == []


def test_kr_subset_is_extracted_when_missing(monkeypatch, tmp_path):
    from pipeline.nfhs import extract
    calls = []

    def fake_extract(out_dir):
        calls.append(Path(out_dir))
        (Path(out_dir) / "kr_harmonised.parquet").write_bytes(b"x")
        return Path(out_dir) / "kr_harmonised.parquet"
    monkeypatch.setattr(extract, "extract_all", fake_extract)
    kr = tmp_path / "sub" / "kr_harmonised.parquet"
    assert assemble.ensure_kr(kr) == kr and calls == [tmp_path / "sub"]
    assert assemble.ensure_kr(kr) == kr and len(calls) == 1      # present: not re-extracted


# ---------------------------------------------------------------- real data (aggregates only)

@pytest.mark.nfhs
def test_real_analysis_dataset_gates(tmp_path):
    """Builds into tmp_path (never rewrites the committed CSVs) and compares with them (M4)."""
    paths = assemble.build(out_dir=tmp_path / "subset", derived_dir=tmp_path / "derived")
    for key, committed in (("counts", assemble.COUNTS_CSV), ("sensitivity", assemble.SENS_CSV)):
        pd.testing.assert_frame_equal(pd.read_csv(paths[key]), pd.read_csv(committed))
    df = pd.read_parquet(paths["analysis"])
    n_task1 = len(pd.read_parquet(assemble.KR_PARQUET, columns=["round"]))
    exc = pd.read_csv(paths["exclusions"])
    res = assemble.gates(df, n_task1, exc, paths["analysis"])
    for g in res:
        print(g.name, g.passed, g.detail)
    assert all(g.passed for g in res)
    assert "caseid" not in df.columns and not ({"lat", "lon"} & set(df.columns))
    p = df[(df["round"] == "nfhs4") & df["district"].isin([39, 53])]
    assert (p["region_id"] == p["region_id_xw"]).mean() > 0.9        # P8 override agrees with GPS
