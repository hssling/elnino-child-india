"""Tests for NFHS KR outcome cleaning. Synthetic rows only; built inside the test."""
import numpy as np
import pandas as pd
import pytest

from pipeline.nfhs.outcomes import (
    REPORT_CSV,
    check_against_report,
    clean_outcomes,
    national_estimates,
)


def _kr(n, **cols):
    """n synthetic KR rows with neutral defaults; keyword args override columns."""
    base = {
        "v005": np.full(n, 1_000_000),
        "v008": np.full(n, 1400),
        "b0": np.zeros(n, dtype=int),
        "b3": np.full(n, 1380),
        "b5": np.ones(n, dtype=int),
        "hw1": np.full(n, 20.0),
        "hw56": np.full(n, 120.0),
        "hw70": np.full(n, 0.0),
        "hw71": np.full(n, 0.0),
        "hw72": np.full(n, 0.0),
        "h11": np.zeros(n),
        "m19": np.full(n, 3000.0),
        "m19a": np.ones(n, dtype=int),
        "round": ["nfhs4"] * n,
    }
    base.update(cols)
    return pd.DataFrame(base)


def test_birth_weight_grams_and_missing_codes():
    out = clean_outcomes(_kr(5, m19=[3000.0, 9996.0, 9998.0, 7000.0, 2400.0]))
    assert out["bw_g"].iloc[0] == 3000.0
    assert out["bw_g"].iloc[1:4].isna().all()
    assert out["bw_g"].iloc[4] == 2400.0
    assert out["lbw"].iloc[0] == 0 and out["lbw"].iloc[4] == 1
    assert out["lbw"].iloc[1:4].isna().all()  # missing is NaN, not 0


def test_birth_weight_card_flag():
    out = clean_outcomes(_kr(3, m19a=[1, 2, 0]))
    assert out["bw_card"].dtype == bool
    assert out["bw_card"].tolist() == [True, False, False]


def test_haemoglobin_and_anaemia():
    out = clean_outcomes(_kr(3, hw56=[105.0, 95.0, 997.0]))
    assert out["hb_adj"].iloc[0] == pytest.approx(10.5)
    assert out["anaemia_any"].iloc[0] == 1 and out["anaemia_modsev"].iloc[0] == 0
    assert out["anaemia_modsev"].iloc[1] == 1 and out["anaemia_any"].iloc[1] == 1
    assert np.isnan(out["hb_adj"].iloc[2])
    assert np.isnan(out["anaemia_any"].iloc[2]) and np.isnan(out["anaemia_modsev"].iloc[2])


def test_anthropometry_z_scores_and_flags():
    out = clean_outcomes(_kr(4, hw70=[-250.0, 9998.0, -650.0, 100.0],
                             hw71=[-250.0, 9996.0, 550.0, -100.0],
                             hw72=[-150.0, 9998.0, -550.0, 600.0]))
    assert out["haz"].iloc[0] == pytest.approx(-2.5) and out["stunted"].iloc[0] == 1
    assert out["waz"].iloc[0] == pytest.approx(-2.5) and out["underweight"].iloc[0] == 1
    assert out["whz"].iloc[0] == pytest.approx(-1.5) and out["wasted"].iloc[0] == 0
    assert out.loc[1, ["haz", "waz", "whz", "stunted", "underweight", "wasted"]].isna().all()
    # WHO flags: HAZ outside [-6, 6], WAZ outside [-6, 5], WHZ outside [-5, 5]
    assert np.isnan(out["haz"].iloc[2]) and np.isnan(out["waz"].iloc[2]) and np.isnan(out["whz"].iloc[2])
    assert out["haz"].iloc[3] == pytest.approx(1.0) and out["stunted"].iloc[3] == 0
    assert np.isnan(out["whz"].iloc[3]) and np.isnan(out["wasted"].iloc[3])


def test_waz_upper_flag_is_five():
    out = clean_outcomes(_kr(2, hw71=[500.0, 501.0]))
    assert out["waz"].iloc[0] == pytest.approx(5.0)
    assert np.isnan(out["waz"].iloc[1])


def test_diarrhoea_codes():
    out = clean_outcomes(_kr(4, h11=[2.0, 0.0, 8.0, np.nan]))
    assert out["diarrhoea"].iloc[0] == 1 and out["diarrhoea"].iloc[1] == 0
    assert np.isnan(out["diarrhoea"].iloc[2]) and np.isnan(out["diarrhoea"].iloc[3])


def test_twins_alive_age_and_cmc():
    out = clean_outcomes(_kr(3, b0=[0, 1, 2], b5=[1, 0, 1], b3=[1380, 1390, 1399], v008=[1400, 1400, 1400]))
    assert out["twin"].tolist() == [False, True, True]
    assert out["alive"].tolist() == [True, False, True]
    assert out["age_m"].tolist() == [20, 10, 1]
    assert out["birth_cmc"].tolist() == [1380, 1390, 1399]
    assert out["interview_cmc"].tolist() == [1400, 1400, 1400]


def test_weights_rescale_to_each_round_row_count():
    df = _kr(5, v005=[1_000_000, 3_000_000, 2_000_000, 500_000, 1_500_000],
             round=["nfhs4", "nfhs4", "nfhs4", "nfhs5", "nfhs5"])
    out = clean_outcomes(df)
    sums = out.groupby("round")["weight"].sum()
    assert sums["nfhs4"] == pytest.approx(3.0)
    assert sums["nfhs5"] == pytest.approx(2.0)
    assert out["weight"].iloc[1] / out["weight"].iloc[0] == pytest.approx(3.0)


def test_clean_outcomes_does_not_mutate_input():
    df = _kr(2, m19=[9996.0, 3000.0])
    before = df.copy()
    clean_outcomes(df)
    pd.testing.assert_frame_equal(df, before)


def test_birth_weight_boundaries():
    out = clean_outcomes(_kr(3, m19=[6000.0, 6001.0, 2500.0]))
    assert out["bw_g"].iloc[0] == 6000.0
    assert np.isnan(out["bw_g"].iloc[1]) and np.isnan(out["lbw"].iloc[1])
    assert out["lbw"].iloc[2] == 0  # 2500 g is not low birth weight


def test_haemoglobin_boundaries():
    out = clean_outcomes(_kr(5, hw56=[110.0, 100.0, 993.0, 994.0, 999.0]))
    assert out["anaemia_any"].iloc[0] == 0
    assert out["anaemia_modsev"].iloc[1] == 0 and out["anaemia_any"].iloc[1] == 1
    assert out["hb_adj"].iloc[2] == pytest.approx(99.3)  # 993 kept by the primary rule
    assert out[["hb_adj", "anaemia_any", "anaemia_modsev"]].iloc[3:].isna().all().all()


def test_haz_boundaries():
    out = clean_outcomes(_kr(5, hw70=[-200.0, -600.0, 600.0, 9996.0, 9997.0]))
    assert out["stunted"].iloc[0] == 0  # exactly -2 SD is not stunted
    assert out["haz"].iloc[1] == pytest.approx(-6.0) and out["haz"].iloc[2] == pytest.approx(6.0)
    assert out[["haz", "stunted"]].iloc[3:].isna().all().all()


def test_diarrhoea_last_24_hours_code_is_yes():
    out = clean_outcomes(_kr(1, h11=[1.0]))
    assert out["diarrhoea"].iloc[0] == 1


def test_hb_plausible_flag():
    out = clean_outcomes(_kr(6, hw56=[39.0, 40.0, 180.0, 181.0, 993.0, 997.0]))
    assert out["hb_plausible"].iloc[:5].tolist() == [0, 1, 1, 0, 0]
    assert np.isnan(out["hb_plausible"].iloc[5])
    # primary rule unchanged: implausible values still classified
    assert out["anaemia_any"].iloc[0] == 1


def test_check_against_report_exactly_one_point_passes(tmp_path):
    report = tmp_path / "report.csv"
    pd.DataFrame({"round": ["nfhs4"], "indicator": ["lbw"], "denominator_note": ["x"],
                  "value_pct": [18.0], "source": ["s"], "page": [1]}).to_csv(report, index=False)
    est = pd.DataFrame({"round": ["nfhs4"], "indicator": ["lbw"], "value_pct": [19.0]})
    assert check_against_report(est, report) == []


def test_national_estimates_weighted_with_report_denominators():
    df = _kr(
        4,
        v005=[1_000_000, 3_000_000, 1_000_000, 1_000_000],
        hw70=[-300.0, 0.0, 9998.0, 0.0],
        hw1=[10.0, 20.0, 30.0, 3.0],
        hw56=[100.0, 120.0, 120.0, 100.0],  # last child is <6 months: excluded from anaemia
        h11=[2.0, 0.0, 0.0, 8.0],
        b5=[1, 1, 0, 1],                    # dead child excluded from diarrhoea
        m19=[2000.0, 3000.0, 2000.0, 9996.0],
    )
    est = national_estimates(clean_outcomes(df)).set_index(["round", "indicator"])["value_pct"]
    # stunting: rows 0,1,3 valid; weights 1,3,1 -> 1/5
    assert est[("nfhs4", "stunting")] == pytest.approx(20.0)
    # anaemia 6-59 months: rows 0,1,2; anaemic row 0 -> 1/5
    assert est[("nfhs4", "anaemia_any")] == pytest.approx(20.0)
    # diarrhoea among living children with a valid answer: rows 0,1 -> 1/4
    assert est[("nfhs4", "diarrhoea")] == pytest.approx(25.0)
    # LBW among all births with a reported weight: rows 0,1,2 -> 2/5
    assert est[("nfhs4", "lbw")] == pytest.approx(40.0)


def test_check_against_report_tolerance(tmp_path):
    report = tmp_path / "report.csv"
    pd.DataFrame({
        "round": ["nfhs4", "nfhs4"],
        "indicator": ["stunting", "wasting"],
        "denominator_note": ["x", "y"],
        "value_pct": [38.0, 21.0],
        "source": ["s", "s"],
        "page": [1, 2],
    }).to_csv(report, index=False)
    est = pd.DataFrame({"round": ["nfhs4", "nfhs4"], "indicator": ["stunting", "wasting"],
                        "value_pct": [39.5, 21.5]})
    failures = check_against_report(est, report)
    assert len(failures) == 1 and "stunting" in failures[0]


def test_check_against_report_missing_estimate_is_a_failure(tmp_path):
    report = tmp_path / "report.csv"
    pd.DataFrame({"round": ["nfhs5"], "indicator": ["lbw"], "denominator_note": ["x"],
                  "value_pct": [18.2], "source": ["s"], "page": [1]}).to_csv(report, index=False)
    est = pd.DataFrame({"round": ["nfhs4"], "indicator": ["lbw"], "value_pct": [18.2]})
    assert len(check_against_report(est, report)) == 1


def test_report_csv_covers_every_indicator_and_round():
    ref = pd.read_csv(REPORT_CSV)
    assert list(ref.columns) == ["round", "indicator", "denominator_note", "value_pct", "source", "page"]
    want = {"stunting", "underweight", "wasting", "anaemia_any", "diarrhoea", "lbw"}
    for rnd in ("nfhs4", "nfhs5"):
        assert set(ref.loc[ref["round"] == rnd, "indicator"]) == want
    assert ref["value_pct"].between(0, 100).all()


@pytest.mark.nfhs
def test_real_kr_matches_report_within_one_point():
    from pipeline.nfhs.outcomes import SUBSET_PARQUET

    kr = pd.read_parquet(SUBSET_PARQUET)
    est = national_estimates(clean_outcomes(kr))
    failures = check_against_report(est, REPORT_CSV)
    assert failures == [], failures
