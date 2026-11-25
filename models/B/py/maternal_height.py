"""Paper B additional analysis A2 (papers/B/envres/analysis_plan_addendum.md): mother's height.

Reads caseid, bidx and v438 (respondent's height, cm x 10) from the raw NFHS-4 and NFHS-5
children's recode files, and writes them in the row order of kr_harmonised.parquet, which is
row-aligned with analysis_child.parquet (checked again on the R side, as for hw57).

Output: data/external/nfhs_subset/maternal_height.parquet. This is controlled data, git-ignored
and never committed. Columns: round, caseid, v001, bidx, b3, mheight_cm. mheight_cm is v438 / 10
for 130.0-200.0 cm; other values (DHS special codes, implausible heights) are missing.

Run from <repo>:  python -I models/B/py/maternal_height.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from pipeline.nfhs import extract as ex  # noqa: E402
from pipeline.paths import EXTERNAL  # noqa: E402

SUBSET = EXTERNAL / "nfhs_subset"
KR = SUBSET / "kr_harmonised.parquet"
OUT = SUBSET / "maternal_height.parquet"


def main() -> None:
    keys = pd.read_parquet(KR, columns=["round", "caseid", "v001", "bidx", "b3"])
    if keys.duplicated(["round", "caseid", "bidx"]).any():
        raise SystemExit("kr_harmonised: (round, caseid, bidx) is not unique")
    parts = []
    for rnd in sorted(keys["round"].unique()):
        raw = ex.read_kr(rnd, ["caseid", "bidx", "v438"], ex.KR_SOURCES[rnd][0])
        raw["round"] = rnd
        parts.append(raw)
    h = pd.concat(parts, ignore_index=True)
    h["caseid"] = h["caseid"].astype(str)
    keys_m = keys.assign(caseid=keys["caseid"].astype(str))
    if h.duplicated(["round", "caseid", "bidx"]).any():
        raise SystemExit("raw KR: (round, caseid, bidx) is not unique")
    m = keys_m.merge(h, on=["round", "caseid", "bidx"], how="left", validate="one_to_one", indicator=True)
    if len(m) != len(keys) or (m["_merge"] != "both").any():
        raise SystemExit(f"{(m['_merge'] != 'both').sum()} kr_harmonised rows not found in the raw KR")
    v = pd.to_numeric(m["v438"], errors="coerce")
    m["mheight_cm"] = np.where((v >= 1300) & (v <= 2000), v / 10.0, np.nan)
    out = m[["round", "caseid", "v001", "bidx", "b3", "mheight_cm"]]
    out.to_parquet(OUT, index=False)
    ok = out["mheight_cm"].notna()
    print(f"wrote {OUT.name}: {len(out):,} rows; valid height {ok.sum():,} ({100 * ok.mean():.1f}%); "
          f"mean {out.loc[ok, 'mheight_cm'].mean():.1f} cm")


if __name__ == "__main__":
    main()
