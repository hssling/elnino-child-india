"""Tests for the NFHS district -> GADM district crosswalk.

Synthetic names/codes only (no NFHS records); the real-data checks are marked
`nfhs`.
"""
import pandas as pd
import pytest

from pipeline.nfhs import districts
from pipeline.nfhs.districts import (
    CROSSWALK_CSV,
    REVIEW_CSV,
    boundary_flags,
    build_crosswalk,
    crosswalk,
    review_rows,
)

R = "nfhs5"


def _gadm(rows):
    return pd.DataFrame(rows, columns=["region_id", "region_name", "state_id", "state_name"])


GADM = _gadm([
    ("IND.16.1_1", "Bengaluru Urban", "IND.16_1", "Karnataka"),
    ("IND.16.2_1", "Mysore", "IND.16_1", "Karnataka"),
    ("IND.16.3_1", "Chamrajnagar", "IND.16_1", "Karnataka"),
    ("IND.5.1_1", "Aurangabad", "IND.5_1", "Bihar"),
    ("IND.5.2_1", "Patna", "IND.5_1", "Bihar"),
    ("IND.20.1_1", "Aurangabad", "IND.20_1", "Maharashtra"),
    ("IND.20.2_1", "Pune", "IND.20_1", "Maharashtra"),
    ("IND.32.1_1", "Khammam", "IND.32_1", "Telangana"),
    ("IND.32.2_1", "Karimnagar", "IND.32_1", "Telangana"),
    ("IND.2.1_1", "Krishna", "IND.2_1", "Andhra Pradesh"),
    ("IND.25.1_1", "West", "IND.25_1", "NCT of Delhi"),
])
STATES = {R: {29: "karnataka", 10: "bihar", 27: "maharashtra", 36: "telangana",
              28: "andhra pradesh", 7: "nct of delhi"}}
LABELS = {R: {
    1: "mysore",                 # exact
    2: "bangalore",              # alias -> Bengaluru Urban
    3: "aurangabad",             # Bihar code: must not go to Maharashtra
    4: "bhadradri kothagudem",   # split (2016) -> Khammam
    5: "atlantis",               # unmatched -> review
    6: "pune",                   # Bihar code, name only exists in Maharashtra
    7: "chamarajanagar",         # fuzzy (ratio >= 0.9)
    8: "shahdara",               # single-polygon state -> manual
    9: "khammam",                # Telangana code: must not match AP
    10: "west",                  # Delhi: GADM's one polygon is named "West"
}}
DSTATE = {R: {1: 29, 2: 29, 3: 10, 4: 36, 5: 27, 6: 10, 7: 29, 8: 7, 9: 36, 10: 7}}


@pytest.fixture
def xw():
    return build_crosswalk(labels=LABELS, state_labels=STATES,
                           district_state=DSTATE, gadm=GADM)


def _row(xw, code):
    return xw.set_index(["round", "nfhs_code"]).loc[(R, code)]


def test_columns(xw):
    for c in ["round", "nfhs_code", "nfhs_name", "state_code", "region_id", "match_type"]:
        assert c in xw.columns
    assert len(xw) == len(LABELS[R])


def test_exact_match(xw):
    r = _row(xw, 1)
    assert r.region_id == "IND.16.2_1" and r.match_type == "exact"


def test_alias_bangalore_to_bengaluru_urban(xw):
    r = _row(xw, 2)
    assert r.region_id == "IND.16.1_1" and r.match_type == "normalised"


def test_cross_state_collision_not_matched(xw):
    r = _row(xw, 3)
    assert r.region_id == "IND.5.1_1"           # Bihar's Aurangabad, not Maharashtra's
    r = _row(xw, 6)
    assert r.match_type == "unmatched" and pd.isna(r.region_id)


def test_split_district_maps_to_parent(xw):
    r = _row(xw, 4)
    assert r.region_id == "IND.32.1_1" and r.match_type == "split"
    assert "Khammam" in r.note


def test_fuzzy_accepted_and_logged(xw):
    r = _row(xw, 7)
    assert r.region_id == "IND.16.3_1" and r.match_type == "fuzzy"
    rev = review_rows(xw)
    assert 7 in set(rev.loc[rev["round"] == R, "nfhs_code"])


def test_single_polygon_state_is_manual_with_reason(xw):
    for code in (8, 10):
        r = _row(xw, code)
        assert r.region_id == "IND.25.1_1" and r.match_type == "manual"
        assert "one polygon" in r.note


def test_telangana_does_not_match_andhra(xw):
    assert _row(xw, 9).region_id == "IND.32.1_1"


def test_unmatched_goes_to_review(xw):
    r = _row(xw, 5)
    assert r.match_type == "unmatched" and pd.isna(r.region_id)
    rev = review_rows(xw)
    hit = rev[(rev["round"] == R) & (rev["nfhs_code"] == 5)]
    assert len(hit) == 1 and hit["reason"].iloc[0]
    # exact/normalised rows are not in review; every review row has a reason
    assert 1 not in set(rev["nfhs_code"]) and rev["reason"].str.len().gt(0).all()


def test_crosswalk_dict_excludes_unmatched(xw, tmp_path):
    p = tmp_path / "xw.csv"
    xw.to_csv(p, index=False)
    m = crosswalk(R, path=p)
    assert m[1] == "IND.16.2_1" and m[4] == "IND.32.1_1"
    assert 5 not in m and all(isinstance(k, int) for k in m)


def test_committed_csvs_hold_no_records():
    for path in (CROSSWALK_CSV, REVIEW_CSV):
        df = pd.read_csv(path)
        assert not {"v001", "caseid", "v005", "bidx"} & set(df.columns)
    xw = pd.read_csv(CROSSWALK_CSV)
    assert not xw.duplicated(["round", "nfhs_code"]).any()
    assert (xw["match_type"] != "unmatched").all()
    assert "boundary_flag" in xw.columns
    assert set(xw["boundary_flag"].fillna("")) == {"", "remnant", "split", "single_polygon"}


@pytest.mark.nfhs
def test_real_crosswalk_covers_every_kr_district():
    xw = build_crosswalk()
    assert (xw["round"] == "nfhs4").sum() == 640
    assert (xw["round"] == "nfhs5").sum() == 707
    kr = pd.read_parquet(districts.KR_PARQUET, columns=["round", "district"])
    used = kr.drop_duplicates()
    mapped = xw[xw["match_type"] != "unmatched"]
    keys = set(zip(mapped["round"], mapped["nfhs_code"]))
    missing = [k for k in zip(used["round"], used["district"].astype(int)) if k not in keys]
    assert missing == []
    flags = xw.groupby(["round", "boundary_flag"]).size().to_dict()
    assert flags[("nfhs4", "remnant")] == 30
    assert flags[("nfhs5", "remnant")] == 1
    assert flags[("nfhs5", "split")] == 35
    assert flags[("nfhs4", "single_polygon")] == 9
    assert flags[("nfhs5", "single_polygon")] == 11
    assert ("nfhs4", "split") not in flags
    committed = pd.read_csv(CROSSWALK_CSV)
    pd.testing.assert_frame_equal(
        committed[["round", "nfhs_code", "region_id", "match_type"]].reset_index(drop=True),
        xw[["round", "nfhs_code", "region_id", "match_type"]].reset_index(drop=True),
        check_dtype=False,
    )
    assert boundary_flags("nfhs4") == dict(zip(
        xw.loc[xw["round"] == "nfhs4", "nfhs_code"], xw.loc[xw["round"] == "nfhs4", "boundary_flag"]))
    assert len(  # GADM_CARVED size as documented
        [k for t in districts.GADM_CARVED.values() for kids in t.values() for k in kids]) == 34
    assert sum(len(t) for t in districts.GADM_CARVED.values()) == 30
    pd.testing.assert_frame_equal(
        committed[["round", "nfhs_code"]].reset_index(drop=True),
        xw[["round", "nfhs_code"]].reset_index(drop=True),
        check_dtype=False,
    )


def test_remnant_note_when_gadm_child_not_reached():
    gadm = _gadm([
        ("IND.20.1_1", "Thane", "IND.20_1", "Maharashtra"),
        ("IND.20.2_1", "Palghar", "IND.20_1", "Maharashtra"),
    ])
    xw = build_crosswalk(labels={R: {1: "thane"}}, state_labels={R: {27: "maharashtra"}},
                         district_state={R: {1: 27}}, gadm=gadm)
    r = xw.iloc[0]
    assert r.match_type == "exact" and "Palghar" in r.note
    assert r.boundary_flag == "remnant"


def test_remnant_note_names_secondary_parent():
    gadm = _gadm([
        ("IND.34.1_1", "Sultanpur", "IND.34_1", "Uttar Pradesh"),
        ("IND.34.2_1", "Amethi", "IND.34_1", "Uttar Pradesh"),
    ])
    xw = build_crosswalk(labels={R: {1: "sultanpur"}}, state_labels={R: {9: "uttar pradesh"}},
                         district_state={R: {1: 9}}, gadm=gadm)
    r = xw.iloc[0]
    assert r.region_id == "IND.34.1_1" and r.boundary_flag == "remnant"
    assert "Amethi (+part of Rae Bareli)" in r.note


def test_boundary_flag_values(xw):
    assert _row(xw, 1).boundary_flag == ""
    assert _row(xw, 4).boundary_flag == "split"
    assert _row(xw, 8).boundary_flag == "single_polygon"
    assert _row(xw, 7).boundary_flag == ""      # fuzzy: spelling only
    assert set(xw.boundary_flag) <= {"", "remnant", "split", "single_polygon"}


def test_boundary_flags_helper(xw, tmp_path):
    p = tmp_path / "xw.csv"
    xw.to_csv(p, index=False)
    f = boundary_flags(R, path=p)
    assert f[4] == "split" and f[8] == "single_polygon" and f[1] == ""
    assert set(f) == set(LABELS[R])
    with pytest.raises(KeyError):
        boundary_flags("nfhs9", path=p)


def test_single_polygon_rule_only_for_allow_listed_units():
    gadm = _gadm([("IND.6.1_1", "Chandigarh", "IND.6_1", "Chandigarh")])
    xw = build_crosswalk(labels={R: {1: "chandigarh", 2: "atlantis"}},
                         state_labels={R: {4: "chandigarh"}},
                         district_state={R: {1: 4, 2: 4}}, gadm=gadm)
    assert "manual" not in set(xw.match_type)
    assert xw.set_index("nfhs_code").loc[1, "match_type"] == "exact"
    assert xw.set_index("nfhs_code").loc[2, "match_type"] == "unmatched"


def test_curated_tables_name_real_gadm_districts():
    gadm = pd.read_csv(districts.GADM_CSV,
                       usecols=["region_name", "state_id"]).drop_duplicates()
    gadm["unit"] = gadm["state_id"].map(districts.states.unit_of_gid1)
    have = set(zip(gadm["unit"], gadm["region_name"]))
    for unit, table in districts.SPLITS.items():
        for _, (parent, reason) in table.items():
            assert (unit, parent) in have and reason
    for unit, table in districts.GADM_CARVED.items():
        for parent, kids in table.items():
            assert (unit, parent) in have
            assert all((unit, k) in have for k in kids)
