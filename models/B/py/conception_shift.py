"""Paper B task 7 extra exposures (ruling P13 Q5): conception-shift windows and post-birth ONI.

Run from <repo>:  python models/B/py/conception_shift.py

* Conception shift +/-1 month (spec §6, results_keys.md `shiftm1` / `shiftp1`): the three
  trimester windows move together. Shift -1 = conception at birth - 10: T1 (-10,-8), T2 (-7,-5),
  T3 (-4,-2). Shift +1 = conception at birth - 8: T1 (-8,-6), T2 (-5,-3), T3 (-2,0). SPI-3 and
  Tmax-anomaly window means on the cluster-GPS cell series (gps.cluster_window_exposure, same NaN
  rules as the primary windows).
* Post-birth ONI (`ensoplacebo`): mean oni_anom over birth +1..+12 months by central month
  (exposure.window_exposure on the ENSO series, as exposure.enso_exposure does for pregnancy).

Pipeline modules are imported read-only. The output, aligned to analysis_child.parquet's row
order (columns `row` and `birth_cmc` for the alignment check), goes to the git-ignored
data/external/nfhs_subset/conception_shift.parquet. Only aggregate coverage is printed.

Fingerprint (ruling P14, M5-M7): the output's parquet metadata records the source
analysis_child.parquet's sha256, md5, size in bytes and row count (keys b_source_*).
`is_current()` compares them with the source on disk; `ensure()` rebuilds when the output is
missing or stale (called by the `nfhs` step of run_all.py). models/B/R/io.R::load_b checks the
same fingerprint (md5 + size) before using the shifted windows.
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pipeline.nfhs import gps  # noqa: E402
from pipeline.nfhs.exposure import PREG_WINDOWS, cmc, window_exposure  # noqa: E402
from pipeline.paths import DERIVED, EXTERNAL  # noqa: E402

SUBSET = EXTERNAL / "nfhs_subset"
ANALYSIS = SUBSET / "analysis_child.parquet"
OUT = SUBSET / "conception_shift.parquet"
SHIFTS = (-1, 1)
SHIFT_VARIABLES = {"spi3": ("spi3", "mean"), "tmax_anom": ("tmax_anom", "mean")}
POST_WINDOW = (1, 12)
FP_PREFIX = "b_source_"
FP_KEYS = ("sha256", "md5", "size", "rows")


def source_fingerprint(path=ANALYSIS) -> dict:
    """sha256 / md5 of the file bytes, its size and its parquet row count."""
    path = Path(path)
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            sha.update(chunk)
            md5.update(chunk)
    return {"sha256": sha.hexdigest(), "md5": md5.hexdigest(), "size": path.stat().st_size,
            "rows": pq.ParquetFile(path).metadata.num_rows}


def write_output(df: pd.DataFrame, out=OUT, analysis=ANALYSIS) -> Path:
    """Write `df` with the source fingerprint in the parquet metadata."""
    fp = source_fingerprint(analysis)
    table = pa.Table.from_pandas(df, preserve_index=False)
    meta = dict(table.schema.metadata or {})
    meta.update({(FP_PREFIX + k).encode(): str(fp[k]).encode() for k in FP_KEYS})
    pq.write_table(table.replace_schema_metadata(meta), out)
    return Path(out)


def read_fingerprint(out=OUT) -> dict:
    meta = pq.read_schema(out).metadata or {}
    return {k: meta[(FP_PREFIX + k).encode()].decode() for k in FP_KEYS if (FP_PREFIX + k).encode() in meta}


def is_current(out=OUT, analysis=ANALYSIS) -> bool:
    """True if `out` exists and was built from the current bytes of `analysis`."""
    if not Path(out).exists():
        return False
    fp = {k: str(v) for k, v in source_fingerprint(analysis).items()}
    return read_fingerprint(out) == fp


def ensure(out=OUT, analysis=ANALYSIS) -> bool:
    """Build the output if missing or stale; True if it was (re)built."""
    if is_current(out, analysis):
        return False
    main(analysis=analysis, out=out)
    return True


def _tag(shift: int) -> str:
    return "s0" if shift == 0 else ("sm" if shift < 0 else "sp") + str(abs(shift))


def shifted_windows(shift: int) -> dict:
    """Trimester windows with conception moved by `shift` months (-1: earlier, +1: later)."""
    return {t: (lo + shift, hi + shift) for t, (lo, hi) in PREG_WINDOWS.items()}


def shift_exposures(children: pd.DataFrame, cluster_monthly: pd.DataFrame, shifts=SHIFTS) -> pd.DataFrame:
    """<T>_<spi3|tmax_anom>_<s0|sm1|sp1> per child (indexed like `children`)."""
    parts = []
    for s in shifts:
        w = gps.cluster_window_exposure(children, cluster_monthly, shifted_windows(s),
                                        variables=SHIFT_VARIABLES, flags={})
        tag = _tag(s)
        cols = [f"{t}_{v}" for t in PREG_WINDOWS for v in SHIFT_VARIABLES]
        parts.append(w[cols].add_suffix(f"_{tag}"))
    return pd.concat(parts, axis=1)


def post_birth_oni(children: pd.DataFrame, enso_monthly: pd.DataFrame) -> pd.Series:
    """Mean oni_anom over birth +1..+12 (central months); NaN outside the ENSO series."""
    em = enso_monthly[["year", "month", "oni_anom"]].copy()
    em["series"] = "enso"
    kids = pd.DataFrame({"series": "enso", "birth_cmc": children["birth_cmc"].to_numpy(),
                         "_far": cmc(9999, 12)}, index=children.index)
    w = window_exposure(kids, em, {"post": POST_WINDOW}, key="series",
                        variables={"oni_anom": ("oni", "mean")}, flags={},
                        birth_col="birth_cmc", interview_col="_far")
    return w["post_oni"]


def main(analysis=ANALYSIS, out=OUT) -> Path:
    kids = pd.read_parquet(analysis, columns=["round", "v001", "birth_cmc", "interview_cmc"])
    kids["v001"] = kids["v001"].astype(np.int64)
    cm = pd.read_parquet(gps.CLUSTER_MONTHLY_PARQUET)
    cm["v001"] = cm["v001"].astype(np.int64)
    target = out
    out = shift_exposures(kids, cm)
    # sanity: shift 0 reproduces the primary windows stored in the analysis parquet
    ref = pd.read_parquet(analysis, columns=["T2_spi3", "T2_tmax_anom"])
    zero = shift_exposures(kids, cm, shifts=(0,))
    for c in ("T2_spi3", "T2_tmax_anom"):
        a, b = zero[f"{c}_s0"].to_numpy(float), ref[c].to_numpy(float)
        same = (np.isnan(a) & np.isnan(b)) | (np.abs(a - b) < 1e-4)
        if not same.all():
            raise ValueError(f"shift-0 {c} does not reproduce the analysis parquet ({(~same).sum()} rows)")
    em = pd.read_csv(DERIVED / "enso_monthly.csv")
    out["post_oni"] = post_birth_oni(kids, em).astype(float)
    out.insert(0, "birth_cmc", kids["birth_cmc"].to_numpy())
    out.insert(0, "row", np.arange(len(kids), dtype=np.int64))
    for c in out.columns:
        if out[c].dtype == float:
            out[c] = out[c].astype(np.float32)
    write_output(out, target, analysis)
    print(f"rows {len(out):,}; non-missing share:")
    for c in [c for c in out.columns if c not in ("row", "birth_cmc")]:
        print(f"  {c}: {out[c].notna().mean():.4f}")
    return Path(target)


if __name__ == "__main__":
    main()
