"""Tests for NFHS KR extraction. Synthetic data only; built inside the test."""
import zipfile

import numpy as np
import pandas as pd
import pytest

from pipeline.nfhs import extract
from pipeline.nfhs.extract import HARMONISED, district_labels, extract_all, read_kr

LABELS = {1: "Nicobars", 2: "North & Middle Andaman"}
STATE_LABELS = {1: "andaman & nicobar islands", 2: "andhra pradesh"}


def _synthetic_dta(path, rnd, labelled=True):
    """Write a small DTA with every source column of the round plus extras.

    Each column gets distinct values (offset 100 * column index) so a swapped
    or mis-renamed column is detectable.
    """
    n = 6
    data = {}
    for i, (out, (s4, s5)) in enumerate(HARMONISED.items()):
        src = s4 if rnd == "nfhs4" else s5
        if src is None:
            continue
        data[src] = np.arange(n, dtype="int32") + 1 + 100 * i
    data["h11"] = np.array([2, 0, 8, 2, 0, 8], dtype="int32")
    data["extra_col"] = np.arange(n, dtype="int32")
    df = pd.DataFrame(data)
    df["caseid"] = ["  1  1  1"] * n
    dist = extract.DISTRICT_VAR[rnd]
    df[dist] = np.array([1, 2, 1, 2, 1, 2], dtype="int32")
    df["v024"] = np.array([1, 2, 1, 2, 1, 2], dtype="int32")
    labels = {dist: LABELS, "v024": STATE_LABELS} if labelled else None
    df.to_stata(path, write_index=False, value_labels=labels, version=118)
    return df


@pytest.fixture
def zips(tmp_path):
    out = {}
    for rnd, (_, member) in extract.KR_SOURCES.items():
        dta = tmp_path / member
        _synthetic_dta(dta, rnd)
        z = tmp_path / f"{rnd}.zip"
        with zipfile.ZipFile(z, "w") as zf:
            zf.write(dta, member)
        out[rnd] = z
    return out


def test_read_kr_returns_numeric_codes(zips):
    df = read_kr("nfhs5", ["h11"], zips["nfhs5"])
    assert pd.api.types.is_numeric_dtype(df["h11"])
    assert sorted(df["h11"].unique()) == [0, 2, 8]


def test_read_kr_selects_columns_only(zips):
    df = read_kr("nfhs4", ["h11", "v005"], zips["nfhs4"])
    assert set(df.columns) == {"h11", "v005"}
    assert len(df) == 6


def test_district_labels_roundtrip(zips):
    assert district_labels("nfhs4", zips["nfhs4"]) == LABELS
    assert district_labels("nfhs5", zips["nfhs5"]) == LABELS


def test_district_labels_missing_label_set_raises(tmp_path):
    member = extract.KR_SOURCES["nfhs5"][1]
    dta = tmp_path / member
    _synthetic_dta(dta, "nfhs5", labelled=False)
    z = tmp_path / "nolabels.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.write(dta, member)
    with pytest.raises(KeyError):
        district_labels("nfhs5", z)


def test_state_labels_roundtrip(zips):
    assert extract.state_labels("nfhs4", zips["nfhs4"]) == STATE_LABELS


def test_extract_all_harmonises_district_name_and_round(zips, tmp_path, monkeypatch):
    monkeypatch.setattr(
        extract, "KR_SOURCES",
        {r: (zips[r], m) for r, (_, m) in extract.KR_SOURCES.items()},
    )
    path = extract_all(tmp_path / "out")
    df = pd.read_parquet(path)
    assert path.name == "kr_harmonised.parquet"
    assert set(df["round"]) == {"nfhs4", "nfhs5"}
    assert len(df) == 12
    assert "district" in df.columns
    assert "sdistri" not in df.columns and "sdist" not in df.columns
    assert df.loc[df["round"] == "nfhs4", "sphase"].isna().all()
    assert df.loc[df["round"] == "nfhs5", "sphase"].notna().all()
    assert set(df["district"]) == {1, 2}
    for rnd in ("nfhs4", "nfhs5"):
        src = pd.read_stata(tmp_path / extract.KR_SOURCES[rnd][1],
                            convert_categoricals=False)
        got = df.loc[df["round"] == rnd]
        assert got["v005"].tolist() == src["v005"].tolist()
        assert got["h11"].tolist() == src["h11"].tolist()
        assert got["v012"].tolist() == src["v012"].tolist()      # mother's age (ruling P11)


@pytest.mark.nfhs
def test_real_district_label_counts():
    assert len(district_labels("nfhs4")) == 640
    assert len(district_labels("nfhs5")) == 707


@pytest.mark.nfhs
def test_real_row_counts_and_size():
    path = extract_all()
    df = pd.read_parquet(path, columns=["round"])
    counts = df["round"].value_counts()
    assert counts["nfhs4"] == 259_627
    assert counts["nfhs5"] == 232_920
    assert path.stat().st_size < 150 * 1024 * 1024
