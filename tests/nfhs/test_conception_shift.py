"""Tests for models/B/py/conception_shift.py (Paper B task 7, ruling P13 Q5).

Synthetic cluster-month metrics and children only (no NFHS records)."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pipeline.nfhs import gps

ROOT = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("conception_shift", ROOT / "models" / "B" / "py" / "conception_shift.py")
cs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(cs)


def _cmc(y, m):
    return (y - 1900) * 12 + m


def _metrics():
    months = pd.period_range("2013-01", "2016-12", freq="M")
    rows = []
    for v001 in (1, 2):
        for k, p in enumerate(months):
            rows.append({"round": "nfhs4", "v001": v001, "year": p.year, "month": p.month,
                         "spi3": 0.1 * k + v001, "spi3_buf10": 0.0, "tmax_anom": -0.05 * k,
                         "heatwave_days": 0.0, "tmax_ge40_days": 0.0,
                         "fallback_rain": False, "fallback_temp": False, "mis_district": False,
                         "qc_masked": False, "qc_cell_masked": False})
    return pd.DataFrame(rows)


def _kids():
    return pd.DataFrame({"round": ["nfhs4", "nfhs4", "nfhs4"], "v001": [1, 2, 1],
                         "birth_cmc": [_cmc(2015, 6), _cmc(2014, 12), _cmc(2015, 1)],
                         "interview_cmc": [_cmc(2016, 3), _cmc(2016, 3), _cmc(2015, 1)]})


def test_shifted_windows_move_every_trimester():
    assert cs.shifted_windows(0) == {"T1": (-9, -7), "T2": (-6, -4), "T3": (-3, -1)}
    assert cs.shifted_windows(-1) == {"T1": (-10, -8), "T2": (-7, -5), "T3": (-4, -2)}
    assert cs.shifted_windows(+1) == {"T1": (-8, -6), "T2": (-5, -3), "T3": (-2, 0)}


def test_shift_zero_reproduces_primary_windows():
    m, k = _metrics(), _kids()
    ref = gps.cluster_window_exposure(k, m, cs.shifted_windows(0))
    out = cs.shift_exposures(k, m, shifts=(0,))
    for t in ("T1", "T2", "T3"):
        np.testing.assert_allclose(out[f"{t}_spi3_s0"], ref[f"{t}_spi3"])
        np.testing.assert_allclose(out[f"{t}_tmax_anom_s0"], ref[f"{t}_tmax_anom"])


def test_shift_minus_one_is_one_month_earlier():
    m, k = _metrics(), _kids()
    out = cs.shift_exposures(k, m, shifts=(-1, 1))
    # child 0, v001 1: spi3 at month index k is 0.1k + 1; birth Jun 2015 = index 29
    b = 29
    np.testing.assert_allclose(out.loc[0, "T2_spi3_sm1"], np.mean([0.1 * (b + j) + 1 for j in (-7, -6, -5)]))
    np.testing.assert_allclose(out.loc[0, "T3_spi3_sp1"], np.mean([0.1 * (b + j) + 1 for j in (-2, -1, 0)]))
    # child 2 interviewed in its birth month: T3 shift +1 includes the birth month (<= interview) -> finite
    assert np.isfinite(out.loc[2, "T3_spi3_sp1"])
    assert list(out.columns) == [f"{t}_{v}_s{s}" for s in ("m1", "p1") for t in ("T1", "T2", "T3")
                                 for v in ("spi3", "tmax_anom")]


def test_post_birth_oni_is_mean_of_months_plus1_to_plus12():
    months = pd.period_range("2014-01", "2017-12", freq="M")
    em = pd.DataFrame({"year": months.year, "month": months.month,
                       "oni_anom": np.arange(len(months), dtype=float),
                       "roni": 0.0, "nino34_anom": 0.0})
    k = _kids()
    out = cs.post_birth_oni(k, em)
    b = 17  # Jun 2015 is index 17 from Jan 2014
    assert out.iloc[0] == pytest.approx(np.mean([b + j for j in range(1, 13)]))


# ---------------------------------------------------------------- fingerprint (ruling P14, M5-M7)

def _write_source(path, values):
    pd.DataFrame({"round": ["nfhs4"] * len(values), "T2_spi3": values}).to_parquet(path, index=False)


def test_source_fingerprint_hashes_bytes_and_counts_rows(tmp_path):
    import hashlib
    src = tmp_path / "analysis_child.parquet"
    _write_source(src, [0.1, 0.2, 0.3])
    fp = cs.source_fingerprint(src)
    raw = src.read_bytes()
    assert fp == {"sha256": hashlib.sha256(raw).hexdigest(), "md5": hashlib.md5(raw).hexdigest(),
                  "size": len(raw), "rows": 3}


def test_write_output_stores_fingerprint_and_is_current_detects_a_changed_source(tmp_path):
    src, out = tmp_path / "analysis_child.parquet", tmp_path / "conception_shift.parquet"
    _write_source(src, [0.1, 0.2])
    assert not cs.is_current(out, src)                       # missing output
    cs.write_output(pd.DataFrame({"row": [0, 1], "x": [1.0, 2.0]}), out, src)
    meta = cs.read_fingerprint(out)
    assert meta == {k: str(v) for k, v in cs.source_fingerprint(src).items()}
    assert cs.is_current(out, src)
    _write_source(src, [0.1, 0.25])                          # same rows, changed exposure
    assert not cs.is_current(out, src)


def test_ensure_builds_only_when_missing_or_stale(tmp_path, monkeypatch):
    src, out = tmp_path / "analysis_child.parquet", tmp_path / "conception_shift.parquet"
    _write_source(src, [0.1])
    calls = []
    monkeypatch.setattr(cs, "main", lambda analysis=None, out=None: calls.append(1) or
                        cs.write_output(pd.DataFrame({"row": [0]}), out, analysis))
    assert cs.ensure(out, src) is True and calls == [1]
    assert cs.ensure(out, src) is False and calls == [1]
    _write_source(src, [0.2])
    assert cs.ensure(out, src) is True and calls == [1, 1]
