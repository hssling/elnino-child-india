"""Stream and harmonise NFHS-4 / NFHS-5 KR (children under 5) columns.

The zipped Stata files are controlled data (DHS data-use agreement): they are
streamed from the zip member and never extracted to disk. Only the listed
columns are kept, with categorical codes left numeric
(convert_categoricals=False) so the numeric missing-code rules hold.
"""
from __future__ import annotations

import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline.paths import EXTERNAL

_NFHS4_ZIP = Path(EXTERNAL / "nfhs" / r"IAKR74DT.zip")
_NFHS5_ZIP = Path(EXTERNAL / "nfhs" / r"IAKR7EDT.zip")

KR_SOURCES = {
    "nfhs4": (_NFHS4_ZIP, "IAKR74FL.DTA"),
    "nfhs5": (_NFHS5_ZIP, "IAKR7EFL.DTA"),
}

DISTRICT_VAR = {"nfhs4": "sdistri", "nfhs5": "sdist"}

# output name -> (NFHS-4 source variable, NFHS-5 source variable); None = absent
# in that round (output is NaN). sphase: the NFHS-4 file has a variable of the
# same name but it is "questionnaire version", not the NFHS-5 survey phase, so
# it is deliberately not mapped for NFHS-4.
HARMONISED = {
    "caseid": ("caseid", "caseid"),
    "v001": ("v001", "v001"),
    "v005": ("v005", "v005"),
    "v008": ("v008", "v008"),
    "v012": ("v012", "v012"),        # mother's age at interview (spec §12 covariate; ruling P11)
    "v016": ("v016", "v016"),
    "v024": ("v024", "v024"),
    "district": ("sdistri", "sdist"),
    "v025": ("v025", "v025"),
    "v104": ("v104", "v104"),
    "v106": ("v106", "v106"),
    "v130": ("v130", "v130"),
    "v190": ("v190", "v190"),
    "s116": ("s116", "s116"),
    "b0": ("b0", "b0"),
    "b1": ("b1", "b1"),
    "b2": ("b2", "b2"),
    "b3": ("b3", "b3"),
    "b4": ("b4", "b4"),
    "b5": ("b5", "b5"),
    "bord": ("bord", "bord"),
    "bidx": ("bidx", "bidx"),
    "hw1": ("hw1", "hw1"),
    "hw53": ("hw53", "hw53"),
    "hw56": ("hw56", "hw56"),
    "hw57": ("hw57", "hw57"),
    "hw70": ("hw70", "hw70"),
    "hw71": ("hw71", "hw71"),
    "hw72": ("hw72", "hw72"),
    "h11": ("h11", "h11"),
    "m19": ("m19", "m19"),
    "m19a": ("m19a", "m19a"),
    "sphase": (None, "sphase"),
}

_CHUNK = 50_000


def _source(round_id: str, zip_path=None):
    if round_id not in KR_SOURCES:
        raise ValueError(f"unknown round {round_id!r}")
    default_zip, member = KR_SOURCES[round_id]
    return Path(zip_path or default_zip), member


def read_kr(round_id: str, columns: list[str], zip_path=None) -> pd.DataFrame:
    """Stream a KR zip member keeping only `columns` (numeric codes, no labels)."""
    path, member = _source(round_id, zip_path)
    columns = list(columns)
    with zipfile.ZipFile(path) as zf, zf.open(member) as fh:
        if fh.seekable():
            return pd.read_stata(fh, columns=columns, convert_categoricals=False)
        parts = []
        with pd.read_stata(
            fh, columns=columns, convert_categoricals=False,
            iterator=True, chunksize=_CHUNK,
        ) as reader:
            for chunk in reader:
                parts.append(chunk[columns])
        return pd.concat(parts, ignore_index=True)


def _value_labels(round_id: str, var: str, zip_path=None) -> dict[int, str]:
    """Value labels attached to `var` ({code: label}); KeyError if it has none."""
    path, member = _source(round_id, zip_path)
    with zipfile.ZipFile(path) as zf, zf.open(member) as fh:
        with pd.read_stata(fh, iterator=True, convert_categoricals=False) as reader:
            labels = reader.value_labels()
            # variable -> label-set name is held in the reader's label list
            names = list(reader._varlist)  # noqa: SLF001 (pandas has no public map)
            lblset = reader._lbllist[names.index(var)]  # noqa: SLF001
    out = {int(k): str(v) for k, v in labels[lblset].items()}
    if not out:
        raise ValueError(f"empty value-label set {lblset!r} for {round_id} {var}")
    return out


def district_labels(round_id: str, zip_path=None) -> dict[int, str]:
    """Value labels of the district variable ({code: name})."""
    return _value_labels(round_id, DISTRICT_VAR[round_id], zip_path)


def state_labels(round_id: str, zip_path=None) -> dict[int, str]:
    """Value labels of the state variable v024 ({code: name})."""
    return _value_labels(round_id, "v024", zip_path)


def extract_all(out_dir=EXTERNAL / "nfhs_subset") -> Path:
    """Write kr_harmonised.parquet (both rounds, harmonised names) and return it."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    frames = []
    for rnd in KR_SOURCES:
        idx = 0 if rnd == "nfhs4" else 1
        src = {out: pair[idx] for out, pair in HARMONISED.items() if pair[idx]}
        raw = read_kr(rnd, list(src.values()), KR_SOURCES[rnd][0])
        df = raw.rename(columns={v: k for k, v in src.items()})
        for out, pair in HARMONISED.items():
            if pair[idx] is None:
                df[out] = np.nan
        df["round"] = rnd
        frames.append(df[[*HARMONISED, "round"]])
    full = pd.concat(frames, ignore_index=True)
    path = out_dir / "kr_harmonised.parquet"
    full.to_parquet(path, index=False)
    return path
