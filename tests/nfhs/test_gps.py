"""Tests for the NFHS cluster-GPS exposure linkage (Paper B task 4b).

Synthetic grids, districts and clusters only (no NFHS records, no real
coordinates). Real-data checks are marked `nfhs` and print aggregates only.
"""
import subprocess
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import box

from pipeline import imd_grid
from pipeline.imd_grid import Grid
from pipeline.nfhs import gps

ROOT = Path(__file__).resolve().parents[2]


# ---------------------------------------------------------------- cells

def test_point_maps_to_expected_rain_and_temp_cells():
    ilat, ilon = gps.cell_index([80.30], [20.10], imd_grid.RAIN)
    # 0.25 deg centres: lon 66.5 + 0.25 j, lat 6.5 + 0.25 i; (80.30, 20.10) is in the cell centred (80.25, 20.0)
    assert (ilat[0], ilon[0]) == (54, 55)
    assert imd_grid.RAIN.lats[ilat[0]] == 20.0 and imd_grid.RAIN.lons[ilon[0]] == 80.25
    ilat, ilon = gps.cell_index([80.30], [20.10], imd_grid.TEMP)
    # 1 deg centres: lon 67.5 + j, lat 7.5 + i; nearest centre (80.5, 20.5)
    assert (ilat[0], ilon[0]) == (13, 13)
    assert imd_grid.TEMP.lats[ilat[0]] == 20.5 and imd_grid.TEMP.lons[ilon[0]] == 80.5


def test_point_on_cell_edge_goes_north_east_deterministically():
    # rain edges at 66.375 + 0.25 k; 80.125 / 20.125 are edges -> the cell to the east / north
    ilat, ilon = gps.cell_index([80.125, 80.125, 80.125 - 1e-9], [20.125, 20.125, 20.125 - 1e-9], imd_grid.RAIN)
    assert imd_grid.RAIN.lons[ilon[0]] == 80.25 and imd_grid.RAIN.lats[ilat[0]] == 20.25
    assert (ilat[0], ilon[0]) == (ilat[1], ilon[1])          # same answer every time
    assert imd_grid.RAIN.lons[ilon[2]] == 80.0 and imd_grid.RAIN.lats[ilat[2]] == 20.0
    # 1 deg edges at integer degrees
    ilat, ilon = gps.cell_index([80.0], [20.0], imd_grid.TEMP)
    assert imd_grid.TEMP.lons[ilon[0]] == 80.5 and imd_grid.TEMP.lats[ilat[0]] == 20.5


def test_point_outside_grid_or_missing_has_no_cell():
    ilat, ilon = gps.cell_index([10.0, np.nan], [10.0, np.nan], imd_grid.RAIN)
    assert list(ilat) == [-1, -1] and list(ilon) == [-1, -1]


def test_nearest_valid_cell_fallback_within_limit_only():
    g = Grid(80.0, 20.0, 0.25, 4, 4)
    valid = np.ones((4, 4), bool)
    valid[0, 0] = False                                      # point's own cell is sea
    ilat, ilon = gps.cell_index([80.02, 81.5], [20.02, 20.0], g)
    li, lj, fb = gps.nearest_valid_cell([80.02, 81.5], [20.02, 20.0], ilat, ilon, valid, g, max_km=40)
    assert fb[0] and valid[li[0], lj[0]] and (li[0], lj[0]) in {(0, 1), (1, 0)}
    assert li[1] == -1 and not fb[1]                          # off the grid and > 40 km from it
    valid[:] = False
    li, lj, fb = gps.nearest_valid_cell([80.02], [20.02], ilat[:1], ilon[:1], valid, g, max_km=40)
    assert li[0] == -1 and not fb[0]


# ---------------------------------------------------------------- GE normalisation and regions

def _squares():
    """Two synthetic districts side by side: A = [80, 81] x [20, 21], B = [81, 82] x [20, 21]."""
    return gpd.GeoDataFrame({"region_id": ["A", "B"], "state_name": ["S", "S"]},
                            geometry=[box(80, 20, 81, 21), box(81, 20, 82, 21)], crs="EPSG:4326")


def _ge(rows):
    df = pd.DataFrame(rows, columns=["DHSCLUST", "LATNUM", "LONGNUM", "URBAN_RURA", "SOURCE", "DHSREGCO"])
    return df


def test_normalise_ge_mis_zero_point_has_no_coordinates():
    ge = _ge([[1, 20.5, 80.5, "R", "GPS", 7.0], [2, 0.0, 0.0, "U", "MIS", 8.0], [3, 20.5, 81.5, "U", "GAZ", 8.0]])
    out = gps.normalise_ge(ge, "nfhs5")
    assert list(out.columns[:7]) == ["round", "v001", "lat", "lon", "urban", "source", "district"]
    mis = out[out.v001 == 2].iloc[0]
    assert mis.source == "MIS" and np.isnan(mis.lat) and np.isnan(mis.lon)
    assert list(out.urban) == [False, True, True] and list(out.district) == [7, 8, 8]
    assert (out["round"] == "nfhs5").all()


def test_point_in_polygon_assigns_synthetic_square_district():
    ge = gps.normalise_ge(_ge([[1, 20.5, 80.5, "R", "GPS", 7.0], [2, 20.5, 81.5, "U", "GAZ", 8.0],
                               [3, 0.0, 0.0, "R", "MIS", 8.0], [4, 20.5, 82.02, "R", "GPS", 8.0],
                               [5, 25.0, 85.0, "R", "GPS", 7.0]]), "nfhs4")
    xw = {"nfhs4": {7: "A", 8: "B"}}
    out = gps.assign_regions(ge, _squares(), xw, max_km=10)
    got = dict(zip(out.v001, zip(out.region_id, out.region_method)))
    assert got[1] == ("A", "pip")
    assert got[2] == ("B", "pip")
    assert got[3] == ("B", "district")           # MIS: from the district code via the crosswalk
    assert got[4] == ("B", "nearest")            # ~2 km east of B, inside the tolerance
    assert got[5] == ("A", "district_fallback")  # far outside every polygon: crosswalk district
    assert list(out.region_xw) == ["A", "B", "B", "B", "A"]


def test_point_on_shared_border_is_deterministic():
    ge = gps.normalise_ge(_ge([[1, 20.5, 81.0, "R", "GPS", 7.0]]), "nfhs4")
    a = gps.assign_regions(ge, _squares(), {"nfhs4": {7: "A"}}, max_km=10)
    b = gps.assign_regions(ge, _squares().iloc[::-1], {"nfhs4": {7: "A"}}, max_km=10)
    assert a.region_id.iloc[0] == b.region_id.iloc[0] == "A"   # lowest region_id wins


# ---------------------------------------------------------------- monthly series

G = Grid(80.125, 20.125, 0.25, 4, 4)       # centres 80.125..80.875 / 20.125..20.875; covers square A
TG = Grid(80.5, 20.5, 1.0, 2, 2)           # 1 deg centres (80.5, 20.5), (81.5, 20.5), ...
YEARS = (1981, 2012)


def _months():
    return pd.period_range(f"{YEARS[0]}-01", f"{YEARS[1]}-12", freq="M")


def _rain(seed=0):
    rng = np.random.default_rng(seed)
    t = len(_months())
    return rng.gamma(2.0, 50.0, size=(t, G.nlat, G.nlon))


def _temp():
    t = len(_months())
    mon = np.tile(np.arange(1, 13), t // 12)
    yr = np.repeat(np.arange(YEARS[0], YEARS[1] + 1), 12)
    tmax = 30.0 + mon[:, None, None] + 0.1 * (yr[:, None, None] - 1981) + np.zeros((t, TG.nlat, TG.nlon))
    tmax[:, 1, 1] = np.nan                                    # one sea cell
    hw = np.where(np.isnan(tmax), np.nan, 2.0)
    ge40 = np.where(np.isnan(tmax), np.nan, 3.0)
    return {"tmax_mean": tmax, "heatwave_days": hw, "tmax_ge40_days": ge40}


def _district_metrics(regions=("A", "B")):
    rows = []
    for r in regions:
        for y in range(YEARS[0], YEARS[1] + 1):
            for m in range(1, 13):
                rows.append({"region_id": r, "year": y, "month": m, "spi3": 0.5, "tmax_anom": 1.5,
                             "heatwave_days": 4.0, "tmax_ge40_days": 5.0, "fallback_rain": False,
                             "fallback_temp": False, "qc_masked": False, "island_proxy": False})
    return pd.DataFrame(rows)


def _clusters():
    return pd.DataFrame({
        "round": ["nfhs4"] * 4,
        "v001": [1, 2, 3, 4],
        "lat": [20.20, 20.60, np.nan, 20.60],
        "lon": [80.20, 80.60, np.nan, 80.60],
        "urban": [False, True, False, False],
        "source": ["GPS", "GAZ", "MIS", "GPS"],
        "region_id": ["A", "A", "B", "I"],      # I: an island district
    })


def _build(clusters=None, qc=None, buffer_w=None, rain=None, island=("I",), dominant="A"):
    return gps.build_cluster_monthly(
        clusters if clusters is not None else _clusters(),
        rain=_rain() if rain is None else rain, temp=_temp(), months=_months(),
        district_metrics=_district_metrics(("A", "B", "I")),
        qc_months=qc if qc is not None else pd.DataFrame(columns=["region_id", "year", "month"]),
        cell_region=np.full((G.nlat, G.nlon), dominant, dtype=object),
        island_regions=set(island), years_by_round={"nfhs4": (2009, 2012)},
        buffer_w=buffer_w, rain_grid=G, temp_grid=TG)


def test_cluster_monthly_uses_containing_cell_spi_and_temperature():
    from pipeline.climate_annual import spi
    out = _build()
    assert list(out.columns[:4]) == ["round", "v001", "year", "month"]
    assert out.year.min() == 2009 and out.year.max() == 2012 and len(out) == 4 * 48
    c1 = out[out.v001 == 1].reset_index(drop=True)
    il, jl = gps.cell_index([80.20], [20.20], G)
    want = spi(pd.Series(_rain()[:, il[0], jl[0]], index=_months()), 3).to_numpy()[-48:]
    np.testing.assert_allclose(c1.spi3.to_numpy(), want, rtol=1e-5)
    # tmax anomaly: value minus the cell's 1981-2010 calendar-month mean (trend 0.1/yr -> mean year 1995.5)
    np.testing.assert_allclose(c1.tmax_anom.to_numpy(), 0.1 * (c1.year.to_numpy() - 1995.5), atol=1e-4)
    assert (c1.heatwave_days == 2.0).all() and (c1.tmax_ge40_days == 3.0).all()
    assert not c1.mis_district.any() and not c1.fallback_rain.any() and not c1.fallback_temp.any()


def test_mis_cluster_takes_district_values_and_is_flagged():
    out = _build()
    mis = out[out.v001 == 3]
    assert mis.mis_district.all()
    assert (mis.spi3 == 0.5).all() and (mis.tmax_anom == 1.5).all() and (mis.heatwave_days == 4.0).all()
    assert (mis.spi3_buf10 == 0.5).all()


def test_masked_district_month_gives_nan_cluster_spi():
    qc = pd.DataFrame({"region_id": ["A"], "year": [2010], "month": [7]})
    out = _build(qc=qc)
    c1 = out[out.v001 == 1].set_index(["year", "month"])
    assert c1.loc[(2010, 7), "qc_masked"]
    # SPI-3 at Jul, Aug, Sep 2010 accumulates the masked July
    for m in (7, 8, 9):
        assert np.isnan(c1.loc[(2010, m), "spi3"]) and np.isnan(c1.loc[(2010, m), "spi3_buf10"])
    assert np.isfinite(c1.loc[(2010, 6), "spi3"]) and np.isfinite(c1.loc[(2010, 10), "spi3"])
    assert np.isfinite(c1.loc[(2010, 7), "tmax_anom"])        # masks blank rain only
    assert out[(out.v001 == 3)].spi3.notna().all()            # district B unaffected


def test_masked_cell_months_stay_out_of_spi_fit():
    from pipeline.climate_annual import spi
    qc = pd.DataFrame({"region_id": ["A"] * 12, "year": [1990] * 12, "month": list(range(1, 13))})
    rain = _rain()
    rain[108:120] = 5000.0                                     # 1990: an artefact that must not enter the fit
    out = _build(qc=qc, rain=rain)
    il, jl = gps.cell_index([80.20], [20.20], G)
    clean = pd.Series(rain[:, il[0], jl[0]], index=_months())
    clean.iloc[108:120] = np.nan
    want = spi(clean, 3).to_numpy()[-48:]
    np.testing.assert_allclose(out[out.v001 == 1].spi3.to_numpy(), want, rtol=1e-5)


def test_own_district_artefact_excluded_from_fit_when_neighbour_dominates_cell():
    from pipeline.climate_annual import spi
    # the cell is dominated by unmasked district B; cluster 1's own district A has an artefact year
    qc = pd.DataFrame({"region_id": ["A"] * 12, "year": [1990] * 12, "month": list(range(1, 13))})
    rain = _rain()
    il, jl = gps.cell_index([80.20], [20.20], G)
    rain[108:120, il[0], jl[0]] = 5000.0
    out = _build(qc=qc, rain=rain, dominant="B")
    clean = pd.Series(rain[:, il[0], jl[0]], index=_months())
    clean.iloc[108:120] = np.nan
    want = spi(clean, 3).to_numpy()[-48:]
    np.testing.assert_allclose(out[out.v001 == 1].spi3.to_numpy(), want, rtol=1e-5)


def test_cell_dominated_by_masked_neighbour_gives_nan_with_cell_flag():
    cl = _clusters()
    cl.loc[0, "region_id"] = "B"                                 # cluster 1 in unmasked district B
    qc = pd.DataFrame({"region_id": ["A"], "year": [2010], "month": [7]})
    out = _build(clusters=cl, qc=qc, dominant="A")
    c1 = out[out.v001 == 1].set_index(["year", "month"])
    for m in (7, 8, 9):
        assert np.isnan(c1.loc[(2010, m), "spi3"])
        assert c1.loc[(2010, m), "qc_cell_masked"] and not c1.loc[(2010, m), "qc_masked"]
    assert not c1.loc[(2010, 6), "qc_cell_masked"] and not c1.loc[(2010, 10), "qc_cell_masked"]
    assert not out[out.v001 == 3].qc_cell_masked.any()          # MIS: district values, no cell


def test_zero_point_with_gps_source_is_treated_as_no_point():
    ge = gps.normalise_ge(_ge([[1, 0.0, 0.0, "R", "GPS", 7.0]]), "nfhs4")
    assert ge.source.iloc[0] == "GPS" and np.isnan(ge.lat.iloc[0]) and np.isnan(ge.lon.iloc[0])
    reg = gps.assign_regions(ge, _squares(), {"nfhs4": {7: "A"}}, max_km=10)
    assert (reg.region_id.iloc[0], reg.region_method.iloc[0]) == ("A", "district")
    cl = _clusters()
    cl.loc[0, ["lat", "lon"]] = np.nan                           # GPS source, no point
    out = _build(clusters=cl)
    c1 = out[out.v001 == 1]
    assert c1.mis_district.all() and (c1.spi3 == 0.5).all() and (c1.tmax_anom == 1.5).all()


def test_island_cluster_rain_nan_and_fallback_temperature_nan():
    cl = _clusters()
    cl.loc[3, ["lat", "lon"]] = [21.60, 81.60]                 # in the sea 1 deg cell -> nearest-cell fallback
    out = _build(clusters=cl)
    isl = out[out.v001 == 4]
    assert isl.spi3.isna().all() and isl.spi3_buf10.isna().all()
    assert isl.fallback_temp.all() and isl.tmax_anom.isna().all() and isl.heatwave_days.isna().all()
    assert isl.island_proxy.all()
    # the same point on the mainland keeps the fallback temperature, flagged
    out2 = _build(clusters=cl, island=())
    m = out2[out2.v001 == 4]
    assert m.fallback_temp.all() and m.tmax_anom.notna().all()


def test_rain_sea_cell_uses_nearest_valid_cell_flagged():
    rain = _rain()
    il, jl = gps.cell_index([80.20], [20.20], G)
    rain[:, il[0], jl[0]] = np.nan                            # point's own 0.25 deg cell has no data
    out = _build(rain=rain)
    c1 = out[out.v001 == 1]
    assert c1.fallback_rain.all() and c1.spi3.notna().all()


# ---------------------------------------------------------------- buffer

def test_buffer_inside_one_cell_has_weight_one():
    cl = pd.DataFrame({"lat": [20.375], "lon": [80.375]})      # a cell centre; 10 km < half cell (~13 km)
    w = gps.buffer_weights(cl, G, km=10)
    assert len(w) == 1 and w.weight.iloc[0] == pytest.approx(1.0)
    il, jl = gps.cell_index([80.375], [20.375], G)
    assert (w.ilat.iloc[0], w.ilon.iloc[0]) == (il[0], jl[0])


def test_buffer_mean_equals_area_weighted_mean_on_synthetic_grid():
    # point at the corner of four cells: each holds ~1/4 of the circle
    cl = pd.DataFrame({"lat": [20.5], "lon": [80.5]})
    w = gps.buffer_weights(cl, G, km=10)
    assert len(w) == 4
    np.testing.assert_allclose(w.weight.to_numpy(), 0.25, atol=0.005)
    assert w.weight.sum() == pytest.approx(1.0)
    vals = np.arange(G.nlat * G.nlon, dtype=float).reshape(1, G.nlat, G.nlon)
    got = gps.buffer_mean(vals, w, n_points=1, valid=np.ones((G.nlat, G.nlon), bool))
    want = sum(r.weight * vals[0, r.ilat, r.ilon] for r in w.itertuples())
    assert got[0, 0] == pytest.approx(want)
    assert got[0, 0] == pytest.approx(vals[0, 1:3, 1:3].mean(), abs=0.05)


def test_buffer_mean_renormalises_over_land_and_is_nan_if_a_land_cell_is_missing():
    w = pd.DataFrame({"row": [0, 0, 0], "ilat": [0, 0, 1], "ilon": [0, 1, 0], "weight": [0.5, 0.3, 0.2]})
    valid = np.ones((2, 2), bool)
    valid[1, 0] = False                                         # sea cell: dropped, weights renormalised
    vals = np.array([[[1.0, 2.0], [100.0, 0.0]], [[1.0, np.nan], [100.0, 0.0]]])
    got = gps.buffer_mean(vals, w, n_points=1, valid=valid)
    assert got[0, 0] == pytest.approx((0.5 * 1 + 0.3 * 2) / 0.8)
    assert np.isnan(got[1, 0])                                  # a land cell missing: no partial mean


def test_buffer_used_for_rural_clusters_point_value_for_urban():
    cl = _clusters()
    w = gps.buffer_weights(cl[cl.source != "MIS"].reset_index(drop=True), G, km=10)
    w = w.assign(row=w.row.map(dict(enumerate(cl.index[cl.source != "MIS"]))))
    out = _build(clusters=cl, buffer_w=w)
    urb = out[out.v001 == 2]
    np.testing.assert_allclose(urb.spi3_buf10.to_numpy(), urb.spi3.to_numpy())
    rur = out[out.v001 == 1]
    assert rur.spi3_buf10.notna().any()
    assert not np.allclose(rur.spi3_buf10.to_numpy(), rur.spi3.to_numpy(), equal_nan=True)


# ---------------------------------------------------------------- window reuse and agreement

def test_cluster_windows_reuse_window_exposure_with_cluster_key():
    out = _build()
    from pipeline.nfhs.exposure import cmc
    kids = pd.DataFrame({"round": ["nfhs4", "nfhs4", "nfhs5"], "v001": [1, 3, 1],
                         "birth_cmc": [cmc(2011, 6)] * 3, "interview_cmc": [cmc(2012, 12)] * 3})
    w = gps.cluster_window_exposure(kids, out)
    c1 = out[out.v001 == 1].set_index(["year", "month"])
    assert w.loc[0, "T3_spi3"] == pytest.approx(c1.loc[[(2011, 3), (2011, 4), (2011, 5)], "spi3"].mean())
    assert w.loc[1, "T3_spi3"] == pytest.approx(0.5)
    assert np.isnan(w.loc[2, "T3_spi3"])                        # same v001, other round: not linked
    assert "T3_spi3_buf10" in w
    assert w.loc[1, "T3_mis_district_any"] and not w.loc[0, "T3_mis_district_any"]


def test_cluster_windows_refuse_windows_before_output_pad():
    from pipeline.nfhs.exposure import cmc
    kids = pd.DataFrame({"round": ["nfhs4"], "v001": [1], "birth_cmc": [cmc(2011, 6)],
                         "interview_cmc": [cmc(2012, 12)]})
    with pytest.raises(ValueError, match="padded"):
        gps.cluster_window_exposure(kids, _build(), {"pre": (-gps.OUTPUT_PAD_MONTHS - 1, -1)})
    gps.cluster_window_exposure(kids, _build(), {"T1s": (-10, -8), "r": (-24, -1, "interview")})


def test_agreement_scores_parent_and_carved_children():
    cl = pd.DataFrame({"round": ["nfhs4"] * 5, "v001": [1, 2, 3, 4, 5],
                       "source": ["GPS", "GPS", "GPS", "GPS", "MIS"],
                       "region_id": ["P", "C", "X", "R", "P"], "region_xw": ["P", "P", "P", "R2", "P"],
                       "boundary_flag": ["remnant", "remnant", "remnant", "", "remnant"]})
    a = gps.agreement(cl, accepted={("nfhs4", "P"): {"P", "C"}})
    row = a.set_index(["round", "group"])
    assert row.loc[("nfhs4", "all"), "n"] == 4                  # MIS excluded (no point)
    assert row.loc[("nfhs4", "all"), "agree"] == pytest.approx(2 / 4)
    assert row.loc[("nfhs4", "remnant"), "agree"] == pytest.approx(2 / 3)
    assert row.loc[("nfhs4", "non_remnant"), "agree"] == pytest.approx(0.0)


# ---------------------------------------------------------------- controlled-data guards

def test_cluster_outputs_are_git_ignored():
    for p in (gps.CLUSTERS_PARQUET, gps.CLUSTER_MONTHLY_PARQUET):
        r = subprocess.run(["git", "check-ignore", "-q", str(p)], cwd=ROOT)
        assert r.returncode == 0, f"{p} is not git-ignored"


def test_no_tracked_file_holds_coordinates_with_cluster_ids():
    files = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True).stdout.split()
    bad = []
    for f in files:
        if not f.lower().endswith((".csv", ".tsv", ".txt")):
            continue
        try:
            head = (ROOT / f).open(encoding="utf-8", errors="ignore").readline().lower()
        except OSError:
            continue
        cols = {c.strip().strip('"') for c in head.replace("\t", ",").split(",")}
        if "v001" in cols and ({"lat", "lon"} & cols or {"latnum", "longnum"} & cols):
            bad.append(f)
    assert not bad, bad


# ---------------------------------------------------------------- real data (aggregates only)

@pytest.mark.nfhs
def test_real_cluster_linkage_aggregates():
    cl = gps.load_or_build_clusters()
    counts = cl.groupby(["round", "source"]).size().unstack(fill_value=0)
    print("\nclusters by source\n" + counts.to_string())
    assert counts.loc["nfhs4"].to_dict() == {"GAZ": 3384, "GPS": 25011, "MIS": 131}
    assert counts.loc["nfhs5"].to_dict() == {"GAZ": 4335, "GPS": 25744, "MIS": 118}
    print(cl.groupby(["round", "region_method"]).size().unstack(fill_value=0).to_string())
    assert cl.region_id.notna().all()
    kr = pd.read_parquet(gps.KR_PARQUET, columns=["round", "v001"]).drop_duplicates()
    missing = kr.merge(cl[["round", "v001"]], how="left", indicator=True)["_merge"].eq("left_only").sum()
    assert missing == 0, f"{missing} KR clusters not in the cluster table"

    a = gps.agreement(cl)
    print(a.to_string(index=False))
    allrow = a[a.group == "all"].set_index("round")
    nonrem = a[a.group == "non_remnant"].set_index("round")
    assert (nonrem.agree >= 0.95).all()
    assert (allrow.agree >= 0.95).all()

    cov = gps.coverage_check(clusters=cl)
    print(cov.to_string(index=False))
    main = cov[cov.group == "outside A&N/Lakshadweep"]
    assert (main["T2_spi3"] >= 0.90).all() and (main["T3_spi3"] >= 0.90).all()
    isl = cov[cov.group == "A&N/Lakshadweep"]
    assert (isl["T2_spi3"] == 0).all()
    size = gps.CLUSTER_MONTHLY_PARQUET.stat().st_size + gps.CLUSTERS_PARQUET.stat().st_size
    print(f"cluster files {size / 1e6:.1f} MB")
    assert size < 100e6
