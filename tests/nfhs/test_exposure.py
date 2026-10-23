"""Tests for NFHS child climate / ENSO exposure windows (Paper B task 4).

Synthetic metrics and children only (no NFHS records). Real-data checks are
marked `nfhs` and report aggregates only.
"""
import numpy as np
import pandas as pd
import pytest

from pipeline.nfhs import exposure
from pipeline.nfhs.exposure import (
    WINDOWS,
    cmc,
    district_enso_sensitivity,
    district_monthly_metrics,
    enso_exposure,
    window_exposure,
)


def _metrics(regions=("R1",), y0=2005, y1=2022, key="region_id"):
    """Metrics whose values encode the calendar month: spi3 = cmc, tmax_anom = yyyymm, hw = 1."""
    rows = []
    for r in regions:
        for y in range(y0, y1 + 1):
            for m in range(1, 13):
                rows.append({key: r, "year": y, "month": m, "spi3": float(cmc(y, m)),
                             "tmax_anom": float(y * 100 + m), "heatwave_days": 1.0})
    return pd.DataFrame(rows)


def _child(birth, interview, region="R1"):
    """birth / interview as (year, month)."""
    return {"region_id": region, "birth_cmc": cmc(*birth), "interview_cmc": cmc(*interview)}


def test_cmc_matches_dhs_definition():
    assert cmc(1900, 1) == 1
    assert cmc(2015, 3) == 115 * 12 + 3


def test_default_windows_match_brief():
    assert WINDOWS["T1"][:2] == (-9, -7) and WINDOWS["T2"][:2] == (-6, -4) and WINDOWS["T3"][:2] == (-3, -1)
    assert WINDOWS["inf0_11"][:2] == (0, 11) and WINDOWS["inf12_23"][:2] == (12, 23)
    assert WINDOWS["placebo"][:2] == (1, 9)
    assert WINDOWS["recent1"] == (-1, -1, "interview") and WINDOWS["recent3"] == (-3, -1, "interview")
    assert WINDOWS["recent12"] == (-12, -1, "interview")


def test_t3_for_march_2015_birth_is_dec2014_to_feb2015():
    kids = pd.DataFrame([_child((2015, 3), (2019, 6))])
    out = window_exposure(kids, _metrics(), {"T3": (-3, -1)})
    months = [cmc(2014, 12), cmc(2015, 1), cmc(2015, 2)]
    assert out.loc[0, "T3_spi3"] == pytest.approx(np.mean(months))
    assert out.loc[0, "T3_tmax_anom"] == pytest.approx(np.mean([201412, 201501, 201502]))
    assert out.loc[0, "T3_hw_days"] == pytest.approx(3.0)       # sum, not mean


def test_one_masked_month_in_t2_makes_only_t2_spi3_nan():
    met = _metrics()
    # birth Mar 2015 -> T2 = Sep, Oct, Nov 2014; mask Oct 2014 rain-derived SPI
    met.loc[(met.year == 2014) & (met.month == 10), "spi3"] = np.nan
    kids = pd.DataFrame([_child((2015, 3), (2019, 6))])
    out = window_exposure(kids, met, WINDOWS)
    assert np.isnan(out.loc[0, "T2_spi3"])
    assert np.isfinite(out.loc[0, "T1_spi3"]) and np.isfinite(out.loc[0, "T3_spi3"])
    assert np.isfinite(out.loc[0, "T2_tmax_anom"])            # temperature not masked by a rain mask


def test_missing_month_row_makes_window_nan():
    met = _metrics()
    met = met[~((met.year == 2014) & (met.month == 10))]        # month absent, not NaN
    out = window_exposure(pd.DataFrame([_child((2015, 3), (2019, 6))]), met, {"T2": (-6, -4)})
    assert out[["T2_spi3", "T2_tmax_anom", "T2_hw_days"]].isna().all(axis=None)


def test_window_past_interview_is_nan():
    # born Jan 2019, interviewed Jun 2020 (age 17 m): inf12_23 ends Dec 2020 > interview
    kids = pd.DataFrame([_child((2019, 1), (2020, 6))])
    out = window_exposure(kids, _metrics(), WINDOWS)
    assert out.loc[0, ["inf12_23_spi3", "inf12_23_tmax_anom", "inf12_23_hw_days"]].isna().all()
    assert np.isfinite(out.loc[0, "inf0_11_spi3"])
    assert np.isfinite(out.loc[0, "placebo_spi3"])


def test_window_ending_at_interview_month_is_kept():
    kids = pd.DataFrame([_child((2019, 1), (2019, 12))])        # inf0_11 = Jan..Dec 2019
    out = window_exposure(kids, _metrics(), {"inf0_11": (0, 11)})
    assert np.isfinite(out.loc[0, "inf0_11_spi3"])
    kids = pd.DataFrame([_child((2019, 1), (2019, 11))])
    assert np.isnan(window_exposure(kids, _metrics(), {"inf0_11": (0, 11)}).loc[0, "inf0_11_spi3"])


def test_window_outside_climate_series_is_nan_not_partial():
    met = _metrics(y0=2010, y1=2020)
    kids = pd.DataFrame([_child((2010, 6), (2015, 1)),    # T1 = Sep-Nov 2009: before the series
                         _child((2020, 6), (2021, 12))])  # inf0_11 runs into 2021: after the series
    out = window_exposure(kids, met, WINDOWS)
    assert np.isnan(out.loc[0, "T1_spi3"]) and np.isnan(out.loc[0, "T1_hw_days"])
    assert np.isfinite(out.loc[0, "T3_spi3"])
    assert np.isnan(out.loc[1, "inf0_11_spi3"]) and np.isnan(out.loc[1, "inf0_11_hw_days"])
    assert np.isfinite(out.loc[1, "T3_spi3"])


def test_recent_windows_are_interview_relative():
    kids = pd.DataFrame([_child((2016, 1), (2019, 6))])
    out = window_exposure(kids, _metrics(), WINDOWS)
    assert out.loc[0, "recent1_spi3"] == pytest.approx(cmc(2019, 5))
    assert out.loc[0, "recent3_spi3"] == pytest.approx(np.mean([cmc(2019, m) for m in (3, 4, 5)]))
    assert out.loc[0, "recent12_hw_days"] == pytest.approx(12.0)


def test_unknown_or_missing_region_is_nan_and_index_preserved():
    kids = pd.DataFrame([_child((2015, 3), (2019, 6)), _child((2015, 3), (2019, 6), region="ZZ"),
                         _child((2015, 3), (2019, 6), region=None)], index=[10, 20, 30])
    out = window_exposure(kids, _metrics(), {"T3": (-3, -1)})
    assert list(out.index) == [10, 20, 30]
    assert np.isfinite(out.loc[10, "T3_spi3"])
    assert out.loc[[20, 30], "T3_spi3"].isna().all()


def test_regions_are_not_mixed():
    met = _metrics(regions=("R1", "R2"))
    met.loc[met.region_id == "R2", "spi3"] += 1000
    kids = pd.DataFrame([_child((2015, 3), (2019, 6), "R1"), _child((2015, 3), (2019, 6), "R2")])
    out = window_exposure(kids, met, {"T3": (-3, -1)})
    assert out.loc[1, "T3_spi3"] - out.loc[0, "T3_spi3"] == pytest.approx(1000)


def test_non_default_cluster_key():
    """Task 4b reuses window_exposure with key ["round", "v001"] on cluster-cell metrics."""
    a = _metrics(regions=(1, 2), key="v001").assign(round="nfhs4")
    b = _metrics(regions=(1,), key="v001").assign(round="nfhs5")
    b["spi3"] += 500                                       # same v001, other round -> other values
    met = pd.concat([a, b], ignore_index=True)
    met.loc[(met.v001 == 2) & (met["round"] == "nfhs4"), "spi3"] -= 7
    kids = pd.DataFrame({"round": ["nfhs4", "nfhs4", "nfhs5", "nfhs5"], "v001": [1, 2, 1, 2],
                         "birth_cmc": cmc(2015, 3), "interview_cmc": cmc(2019, 6)})
    out = window_exposure(kids, met, {"T3": (-3, -1)}, key=["round", "v001"])
    base = np.mean([cmc(2014, 12), cmc(2015, 1), cmc(2015, 2)])
    assert out["T3_spi3"].iloc[0] == pytest.approx(base)
    assert out["T3_spi3"].iloc[1] == pytest.approx(base - 7)
    assert out["T3_spi3"].iloc[2] == pytest.approx(base + 500)
    assert np.isnan(out["T3_spi3"].iloc[3])                # (nfhs5, 2) not in metrics
    with pytest.raises(KeyError):
        window_exposure(kids, met, {"T3": (-3, -1)})       # default key region_id absent


def test_duplicate_metric_rows_raise():
    met = _metrics()
    with pytest.raises(ValueError, match="duplicate"):
        window_exposure(pd.DataFrame([_child((2015, 3), (2019, 6))]), pd.concat([met, met.head(1)]),
                        {"T3": (-3, -1)})


def test_bad_window_raises():
    with pytest.raises(ValueError):
        window_exposure(pd.DataFrame([_child((2015, 3), (2019, 6))]), _metrics(), {"X": (2, 1)})
    with pytest.raises(ValueError):
        window_exposure(pd.DataFrame([_child((2015, 3), (2019, 6))]), _metrics(), {"X": (0, 1, "conception")})


# ---------------------------------------------------------------- district_monthly_metrics

def _never(y, m):
    return False


def _always(y, m):
    return True


def _month_table(regions=("R1", "R2"), y0=1981, y1=2016, seed=0, states=None,
                 fallback_rain=None, fallback_temp=None):
    """states: {region: state_name} (default "Kerala"); fallback_*: {region: predicate(year, month)}."""
    rng = np.random.default_rng(seed)
    states, fallback_rain, fallback_temp = states or {}, fallback_rain or {}, fallback_temp or {}
    rows = []
    for r in regions:
        for y in range(y0, y1 + 1):
            for m in range(1, 13):
                rows.append({"region_id": r, "year": y, "month": m,
                             "rain_mm": float(rng.gamma(2.0, 40.0)), "rainy_days": 5.0,
                             "tmax_mean": 30.0 + m + (y - 1981) * 0.01, "tmin_mean": 20.0,
                             "tmax_ge40_days": float(m == 5), "heatwave_days": float(m == 5 and y % 2 == 0),
                             "fallback_rain": bool(fallback_rain.get(r, _never)(y, m)),
                             "fallback_temp": bool(fallback_temp.get(r, _never)(y, m)),
                             "state_name": states.get(r, "Kerala")})
    return pd.DataFrame(rows)


def _qc_file(tmp_path, rows):
    q = pd.DataFrame(rows, columns=["level", "region_id", "region_name", "year", "variable", "ratio",
                                    "disposition", "reason", "evidence"])
    p = tmp_path / "qc_reviewed.csv"
    q.to_csv(p, index=False)
    return p


def test_district_monthly_metrics_spi_tmax_and_mask(tmp_path):
    month = _month_table()
    qc = _qc_file(tmp_path, [["district", "R1", "x", 2014, "rain_jjas_mm", 9.9, "mask", "test", "test"]])
    met = district_monthly_metrics(month, qc_path=qc)
    assert set(["region_id", "year", "month", "spi3", "tmax_anom", "heatwave_days", "tmax_ge40_days",
                "rain_mm", "qc_masked"]) <= set(met.columns)
    assert not met.duplicated(["region_id", "year", "month"]).any()
    r1 = met[met.region_id == "R1"].set_index(["year", "month"])
    # masked Jun-Sep 2014: SPI-3 NaN from Jun to Nov (3-month accumulations touching a masked month)
    for mo in range(6, 12):
        assert np.isnan(r1.loc[(2014, mo), "spi3"]), mo
    assert np.isfinite(r1.loc[(2014, 5), "spi3"]) and np.isfinite(r1.loc[(2014, 12), "spi3"])
    assert r1.loc[(2014, 7), "qc_masked"] and np.isnan(r1.loc[(2014, 7), "rain_mm"])
    # R2 is untouched by R1's mask
    assert met[(met.region_id == "R2") & (met.year == 2014)].spi3.notna().all()
    # SPI equals climate_annual.spi on the region's own series
    from pipeline.climate_annual import spi
    g = month[month.region_id == "R2"]
    ref = spi(pd.Series(g.rain_mm.to_numpy(), index=pd.PeriodIndex.from_fields(
        year=g.year, month=g.month, freq="M")), 3)
    r2 = met[met.region_id == "R2"].sort_values(["year", "month"])
    np.testing.assert_allclose(r2.spi3.to_numpy(), ref.to_numpy(), equal_nan=True)
    # tmax anomaly vs the 1981-2010 calendar-month mean: tmax = 30 + m + 0.01 (y-1981)
    assert r1.loc[(2015, 5), "tmax_anom"] == pytest.approx(0.01 * (2015 - 1981) - 0.01 * 14.5)
    assert r1.loc[(2014, 7), "tmax_anom"] == pytest.approx(0.01 * (2014 - 1981) - 0.01 * 14.5)  # temp unmasked


def test_metrics_feed_window_exposure_masked_t2(tmp_path):
    month = _month_table(regions=("R1",))
    qc = _qc_file(tmp_path, [["district", "R1", "x", 2014, "rain_jjas_mm", 9.9, "mask", "test", "test"]])
    met = district_monthly_metrics(month, qc_path=qc)
    kids = pd.DataFrame([_child((2015, 3), (2016, 6)),    # T2 = Sep-Nov 2014: masked
                         _child((2014, 3), (2016, 6))])   # T2 = Sep-Nov 2013: clean
    out = window_exposure(kids, met, WINDOWS)
    assert np.isnan(out.loc[0, "T2_spi3"]) and np.isfinite(out.loc[1, "T2_spi3"])
    assert np.isfinite(out.loc[0, "T2_tmax_anom"])


def test_metrics_keep_gaps_out_of_rolling_windows():
    month = _month_table(regions=("R1",))
    month = month[~((month.year == 2012) & (month.month == 4))]   # missing row
    met = district_monthly_metrics(month, qc_path=None)
    r1 = met.set_index(["year", "month"])
    assert np.isnan(r1.loc[(2012, 5), "spi3"]) and np.isnan(r1.loc[(2012, 6), "spi3"])
    assert np.isfinite(r1.loc[(2012, 7), "spi3"])


def test_island_proxy_climate_is_nan():
    """Lakshadweep rain (Kerala-coast proxy) and A&N rain / fallback temperature are not exposures
    (ruling P7); mainland fallbacks (here a Mizoram district like Saiha) stay, flagged."""
    month = _month_table(
        regions=("R1", "LK", "AN", "SA"),
        states={"LK": "Lakshadweep", "AN": "Andaman and Nicobar", "SA": "Mizoram"},
        fallback_rain={"LK": _always, "AN": _always, "SA": lambda y, m: y <= 2012},
        fallback_temp={"AN": lambda y, m: y >= 1984 and y != 2006, "SA": _always})
    met = district_monthly_metrics(month, qc_path=None)
    by = {r: g.set_index(["year", "month"]) for r, g in met.groupby("region_id")}
    for r in ("LK", "AN"):
        assert by[r][["rain_mm", "spi3"]].isna().all(axis=None), r
        assert by[r].island_proxy.all()
    # A&N temperature NaN exactly where fallback_temp; real-cell months kept
    an = by["AN"]
    proxy_t = an.fallback_temp.astype(bool)
    assert an.loc[proxy_t, ["tmax_mean", "tmax_anom", "heatwave_days", "tmax_ge40_days"]].isna().all(axis=None)
    assert an.loc[~proxy_t, ["tmax_mean", "heatwave_days", "tmax_ge40_days"]].notna().all(axis=None)
    assert np.isfinite(an.loc[(2006, 5), "tmax_anom"])
    # Lakshadweep has real-cell temperature (no fallback_temp): kept
    assert by["LK"][["tmax_anom", "heatwave_days"]].notna().all(axis=None)
    # mainland: fallback rain and temperature kept, flagged
    sa = by["SA"]
    assert sa.spi3.notna().sum() > 0.9 * len(sa) and sa.tmax_anom.notna().all()
    assert sa.loc[(2010, 7), "fallback_rain"] and not sa.loc[(2014, 7), "fallback_rain"]
    assert not sa.island_proxy.any() and not by["R1"].island_proxy.any()
    assert by["R1"].spi3.notna().sum() > 0.9 * len(by["R1"])
    # windows and sensitivity follow
    kids = pd.DataFrame([_child((2015, 3), (2016, 6), r) for r in ("R1", "LK", "AN", "SA")])
    out = window_exposure(kids, met, WINDOWS)
    assert np.isfinite(out.loc[0, "T2_spi3"]) and np.isfinite(out.loc[3, "T2_spi3"])
    assert np.isnan(out.loc[1, "T2_spi3"]) and np.isnan(out.loc[2, "T2_spi3"])
    assert np.isfinite(out.loc[1, "T2_tmax_anom"])
    assert np.isnan(out.loc[2, "T2_tmax_anom"]) and np.isnan(out.loc[2, "T2_hw_days"])
    ann = pd.DataFrame({"year": range(1981, 2017), "oni_jjas": np.linspace(-1, 1, 36)})
    sens = district_enso_sensitivity(met, ann).set_index("region_id")
    assert sens.loc[["LK", "AN"], "sens_beta"].isna().all()
    assert sens.loc[["R1", "SA"], "sens_beta"].notna().all()


def test_fallback_rain_any_window_flag():
    met = _metrics()
    met["fallback_rain"] = (met.year == 2014) & (met.month == 10)      # inside T2 of a Mar 2015 birth
    kids = pd.DataFrame([_child((2015, 3), (2019, 6)), _child((2017, 3), (2019, 6)),
                         _child((2015, 3), (2019, 6), region="ZZ")])
    out = window_exposure(kids, met, WINDOWS)
    assert out.loc[0, "T2_fallback_rain_any"] and not out.loc[0, "T1_fallback_rain_any"]
    assert not out.loc[0, "T3_fallback_rain_any"] and not out.loc[1, "T2_fallback_rain_any"]
    assert np.isfinite(out.loc[0, "T2_spi3"])                         # flag only, exposure kept
    assert out["T2_fallback_rain_any"].dtype == bool
    assert not out.loc[2, "T2_fallback_rain_any"]                     # unlinked child: no fallback month
    assert not any(c.endswith("_fallback_rain_any") for c in window_exposure(kids, _metrics(), WINDOWS))


# ---------------------------------------------------------------- ENSO

def _enso(years=range(2005, 2022), dev=(2009, 2015)):
    mon = []
    for y in years:
        for m in range(1, 13):
            mon.append({"year": y, "month": m, "oni_anom": (y - 2000) + m / 100,
                        "roni": -((y - 2000) + m / 100), "nino34_anom": 2 * ((y - 2000) + m / 100)})
    ann = pd.DataFrame({"year": list(years), "oni_jjas": 0.0,
                        "el_nino_developing": [y in dev for y in years]})
    return pd.DataFrame(mon), ann


def test_preg_elnino_dev():
    mon, ann = _enso()
    kids = pd.DataFrame([
        {"birth_cmc": cmc(2015, 11), "interview_cmc": cmc(2019, 1)},   # pregnancy Feb-Oct 2015: Jul 2015 in
        {"birth_cmc": cmc(2016, 3), "interview_cmc": cmc(2019, 1)},    # Jun 2015-Feb 2016: overlaps Jul 2015
        {"birth_cmc": cmc(2017, 11), "interview_cmc": cmc(2019, 1)},   # Feb-Oct 2017, not developing
        {"birth_cmc": cmc(2015, 7), "interview_cmc": cmc(2019, 1)},    # Oct 2014-Jun 2015: Jun 2015 in JJAS
        {"birth_cmc": cmc(2015, 6), "interview_cmc": cmc(2019, 1)},    # Sep 2014-May 2015: no JJAS of 2015
    ])
    out = enso_exposure(kids, mon, ann)
    assert out["preg_elnino_dev"].tolist() == [True, True, False, True, False]


def test_preg_oni_means_and_parts():
    mon, ann = _enso()
    kids = pd.DataFrame([{"birth_cmc": cmc(2015, 10), "interview_cmc": cmc(2019, 1)},
                         {"birth_cmc": cmc(2004, 6), "interview_cmc": cmc(2008, 1)}])   # before ENSO table
    out = enso_exposure(kids, mon, ann)
    preg = [(2015 - 2000) + m / 100 for m in range(1, 10)]          # Jan-Sep 2015
    assert out.loc[0, "preg_oni"] == pytest.approx(np.mean(preg))
    assert out.loc[0, "T1_oni"] == pytest.approx(np.mean(preg[:3]))
    assert out.loc[0, "T2_oni"] == pytest.approx(np.mean(preg[3:6]))
    assert out.loc[0, "T3_oni"] == pytest.approx(np.mean(preg[6:]))
    assert out.loc[0, "preg_oni_pos"] == pytest.approx(np.mean(preg)) and out.loc[0, "preg_oni_neg"] == 0
    assert out.loc[0, "preg_roni"] == pytest.approx(-np.mean(preg))
    assert out.loc[0, "preg_nino34"] == pytest.approx(2 * np.mean(preg))
    assert out.loc[1, ["preg_oni", "preg_oni_pos", "preg_oni_neg", "T1_oni", "preg_roni"]].isna().all()
    assert pd.isna(out.loc[1, "preg_elnino_dev"])


def test_preg_oni_neg_part():
    mon, ann = _enso()
    mon["oni_anom"] = -0.5
    out = enso_exposure(pd.DataFrame([{"birth_cmc": cmc(2015, 10), "interview_cmc": cmc(2019, 1)}]), mon, ann)
    assert out.loc[0, "preg_oni_neg"] == pytest.approx(-0.5) and out.loc[0, "preg_oni_pos"] == 0


# ---------------------------------------------------------------- sensitivity

def _sens_inputs(post_slope=None):
    """District rain whose JJAS % departure from its 1981-2009 mean is exactly -10 x oni_jjas."""
    years = list(range(1981, 2021))
    rng = np.random.default_rng(3)
    oni = rng.normal(0, 1, len(years))
    pre = np.array([y <= 2009 for y in years])
    oni[pre] -= oni[pre].mean()                                   # zero mean over 1981-2009
    rows = []
    for y, o in zip(years, oni):
        dep = -10 * o if (y <= 2009 or post_slope is None) else post_slope * o
        jjas = 1000.0 * (1 + dep / 100)
        for m in range(1, 13):
            rows.append({"region_id": "D1", "year": y, "month": m,
                         "rain_mm": jjas / 4 if m in (6, 7, 8, 9) else 7.0})
    ann = pd.DataFrame({"year": years, "oni_jjas": oni})
    return pd.DataFrame(rows), ann


def test_sensitivity_recovers_slope():
    met, ann = _sens_inputs()
    s = district_enso_sensitivity(met, ann)
    assert list(s.columns[:3]) == ["region_id", "sens_beta", "sens_se"]
    assert s.loc[0, "sens_beta"] == pytest.approx(-10, abs=1e-9)
    assert s.loc[0, "sens_se"] == pytest.approx(0, abs=1e-9)
    assert s.loc[0, "sens_n"] == 29


def test_sensitivity_uses_no_year_from_2010():
    met, ann = _sens_inputs()
    met2, ann2 = _sens_inputs(post_slope=+50)
    met2.loc[met2.year >= 2010, "rain_mm"] *= 3                    # also shifts any normal that used them
    a = district_enso_sensitivity(met, ann)
    b = district_enso_sensitivity(met2, ann2)
    assert b.loc[0, "sens_beta"] == pytest.approx(a.loc[0, "sens_beta"], abs=1e-12)
    assert b.loc[0, "sens_beta"] == pytest.approx(-10, abs=1e-9)


def test_sensitivity_skips_masked_years_and_short_series():
    met, ann = _sens_inputs()
    met.loc[(met.year == 1990) & (met.month == 7), "rain_mm"] = np.nan
    s = district_enso_sensitivity(met, ann)
    assert s.loc[0, "sens_n"] == 28
    short = met[met.year >= 2000]
    s2 = district_enso_sensitivity(short, ann)
    assert np.isnan(s2.loc[0, "sens_beta"]) and np.isnan(s2.loc[0, "sens_se"])


# ---------------------------------------------------------------- real data (controlled; aggregates only)

@pytest.mark.nfhs
def test_real_t2_t3_coverage():
    from pipeline.nfhs.outcomes import SUBSET_PARQUET
    if not SUBSET_PARQUET.exists():
        pytest.skip("NFHS subset not built")
    cov = exposure.coverage_check()
    print("\n" + cov.to_string())
    main = cov[cov.group == "outside A&N/Lakshadweep"]
    assert len(main) == 2
    assert (main[["T2_spi3", "T3_spi3", "T2_tmax_anom", "T3_tmax_anom"]] >= 0.90).all(axis=None)
    # these would catch a mask failure: reviewed Tamil Nadu 2011-12 masks; island proxies (ruling P7)
    st = cov[cov.level == "state"].set_index(["round", "group"])
    spi_cols = [f"{w}_spi3" for w in exposure.WINDOWS]
    tmax_cols = [f"{w}_tmax_anom" for w in exposure.WINDOWS]
    assert st.loc[("nfhs4", "Tamil Nadu"), "T2_spi3"] < 0.75
    for rnd in ("nfhs4", "nfhs5"):
        assert (st.loc[(rnd, "Andaman and Nicobar"), spi_cols] == 0).all()
        assert (st.loc[(rnd, "Lakshadweep"), spi_cols] == 0).all()
        assert (st.loc[(rnd, "Andaman and Nicobar"), tmax_cols] == 0).all()
