"""Crosswalk NFHS-4/5 district codes to GADM districts (climate `region_id`).

Matching is done only within the same state (NFHS `v024` -> GADM state via the
harmonised state frame) in this order:

1. exact       normalised names equal (case, punctuation, "&"/"and")
2. normalised  compact keys overlap after dropping spaces, a bracketed part, a
               "district"/"the" affix, or through a curated alias set (keys are
               built for both the NFHS and the GADM name; unique hit required)
3. manual      the GADM layer has a single polygon for a state that NFHS splits
               into several districts; only for units in SINGLE_POLYGON_UNITS
               (NCT of Delhi); reason recorded
4. split       district carved after 2011 that GADM 4.1 does not have; mapped to
               its parent GADM district (curated table with year + reason)
5. fuzzy       best difflib ratio >= 0.9 on compact keys, unique best, target not
               already taken by an exact/normalised match; logged for review
6. unmatched   everything else; listed in the review file with a reason

GADM is version 4.1 (`gadm41_IND.gpkg`; district names via
climate_district_year.csv). `boundary_flag` marks rows whose climate unit does
not coincide with the NFHS district: remnant | split | single_polygon | "".

The CSV outputs hold district codes and names only (no NFHS records).
"""
from __future__ import annotations

import re
from difflib import SequenceMatcher
from pathlib import Path

import pandas as pd

from pipeline.outcomes import states
from pipeline.paths import DERIVED, EXTERNAL

REFERENCE = Path(__file__).resolve().parents[1] / "reference"
CROSSWALK_CSV = REFERENCE / "nfhs_district_crosswalk.csv"
REVIEW_CSV = REFERENCE / "nfhs_district_review.csv"
KR_PARQUET = EXTERNAL / "nfhs_subset" / "kr_harmonised.parquet"
GADM_CSV = DERIVED / "climate_district_year.csv"

FUZZY_MIN = 0.9
COLUMNS = ["round", "nfhs_code", "nfhs_name", "state_code", "state_name",
           "region_id", "gadm_name", "match_type", "boundary_flag", "note"]

# Units where GADM 4.1 has one polygon for a state NFHS splits into districts.
SINGLE_POLYGON_UNITS = {"DL"}

# Spelling / renaming equivalence sets, per state-frame unit (names as written;
# they are normalised before use). Symmetric: any member matches any other.
ALIASES: dict[str, list[set[str]]] = {
    "AN": [{"nicobars", "nicobar islands", "nicobar"}],
    "AP_TG": [{"sri potti sriramulu nellore", "nellore"},
              {"y.s.r.", "ysr kadapa", "kadapa", "cuddapah"}],
    "CG": [{"dakshin bastar dantewada", "dantewada"},
           {"kabirdham", "kabeerdham", "kawardha"},
           {"kodagaon", "kondagaon"}],
    "GJ": [{"dohad", "dahod"}],
    "JK": [{"punch", "poonch"}],
    "KA": [{"bangalore", "bengaluru urban", "bengaluru", "bangalore urban"},
           {"bangalore rural", "bengaluru rural"}],
    "MH": [{"gadchiroli", "garhchiroli"}, {"mumbai", "mumbai city"}],
    "OD": [{"baudh", "bauda", "boudh"}],
    "SK": [{"north district", "north sikkim"}, {"south district", "south sikkim"},
           {"east district", "east sikkim"}, {"west district", "west sikkim"}],
    "TR": [{"unakoti", "unokoti"}],
    "UP": [{"jyotiba phule nagar", "amroha"}, {"mahamaya nagar", "hathras"},
           {"kanshiram nagar", "kasganj"}, {"kheri", "lakhimpur kheri"},
           {"shrawasti", "shravasti"}],
    "WB": [{"north twenty four parganas", "north 24 parganas"},
           {"south twenty four parganas", "south 24 parganas"}],
}

# NFHS districts created after 2011 that GADM 4.1 does not separate:
# unit -> {nfhs name: (GADM parent name, reason)}. Parent = district the new
# one was (mostly) carved from.
SPLITS: dict[str, dict[str, tuple[str, str]]] = {
    "AP_TG": {
        "bhadradri kothagudem": ("Khammam", "2016 Telangana reorganisation; carved from Khammam"),
        "jagitial": ("Karimnagar", "2016; carved from Karimnagar"),
        "jangoan": ("Warangal", "2016; mainly from Warangal (small parts Nalgonda, Karimnagar)"),
        "jayashankar bhupalapally": ("Warangal", "2016; mainly from Warangal (part Karimnagar)"),
        "jogulamba gadwal": ("Mahbubnagar", "2016; carved from Mahbubnagar"),
        "kamareddy": ("Nizamabad", "2016; carved from Nizamabad"),
        "komaram bheem asifabad": ("Adilabad", "2016; carved from Adilabad"),
        "mahabubabad": ("Warangal", "2016; mainly from Warangal (part Khammam)"),
        "mancherial": ("Adilabad", "2016; carved from Adilabad"),
        "medchal-malkajgiri": ("Ranga Reddy", "2016; carved from Ranga Reddy"),
        "nagarkurnool": ("Mahbubnagar", "2016; carved from Mahbubnagar"),
        "nirmal": ("Adilabad", "2016; carved from Adilabad"),
        "peddapalli": ("Karimnagar", "2016; carved from Karimnagar"),
        "rajanna sircilla": ("Karimnagar", "2016; carved from Karimnagar"),
        "sangareddy": ("Medak", "2016; carved from Medak (Sangareddy was the old Medak district HQ)"),
        "siddipet": ("Medak", "2016; mainly from Medak (parts Karimnagar, Warangal)"),
        "suryapet": ("Nalgonda", "2016; carved from Nalgonda"),
        "vikarabad": ("Ranga Reddy", "2016; carved from Ranga Reddy (part Mahbubnagar)"),
        "wanaparthy": ("Mahbubnagar", "2016; carved from Mahbubnagar"),
        "warangal rural": ("Warangal", "2016; Warangal split into Rural and Urban"),
        "warangal urban": ("Warangal", "2016; Warangal split into Rural and Urban (part Karimnagar)"),
        "yadadri bhuvanagiri": ("Nalgonda", "2016; carved from Nalgonda"),
    },
    "AR": {
        "kra daadi": ("Kurung Kumey", "2015; carved from Kurung Kumey"),
        "siang": ("West Siang", "2015; carved from East and West Siang; parent chosen = West Siang (ambiguous)"),
    },
    "AS": {
        "biswanath": ("Sonitpur", "2015; carved from Sonitpur"),
        "charaideo": ("Sivasagar", "2015; carved from Sivasagar"),
        "hojai": ("Nagaon", "2015; carved from Nagaon"),
        "majuli": ("Jorhat", "2016; carved from Jorhat"),
        "south salmara mancachar": ("Dhubri", "2015; carved from Dhubri"),
        "west karbi anglong": ("Karbi Anglong", "2016; carved from Karbi Anglong"),
    },
    "HR": {
        "charkhi dadri": ("Bhiwani", "2016; carved from Bhiwani"),
    },
    "ML": {
        "east jantia hills": ("Jaintia Hills", "2012; Jaintia Hills split into East and West"),
        "west jaintia hills": ("Jaintia Hills", "2012; Jaintia Hills split into East and West"),
    },
    "WB": {
        "paschim barddhaman": ("Barddhaman", "2017; Barddhaman split into Paschim and Purba"),
        "purba barddhaman": ("Barddhaman", "2017; Barddhaman split into Paschim and Purba"),
    },
}

# GADM 4.1 districts carved after the 2011 Census frame (34 children under 30
# parents): unit -> {GADM parent: [GADM children]}. Only used to annotate (and
# flag `remnant`) an NFHS district that maps to the parent while the children
# are not reached by any code in that round (the NFHS district then also covers
# them; its climate is the parent remnant's). region_id is never changed.
GADM_CARVED: dict[str, dict[str, list[str]]] = {
    "AR": {"Lohit": ["Namsai"], "Tirap": ["Longding"]},
    "CG": {"Durg": ["Balod", "Bemetara"], "Raipur": ["Baloda Bazar", "Gariaband"],
           "Bilaspur": ["Mungeli"], "Surguja": ["Balrampur", "Surajpur"],
           "Bastar": ["Kondagaon"], "Dantewada": ["Sukma"]},
    "GJ": {"Bhavnagar": ["Botad"], "Junagadh": ["Gir Somnath"],
           "Jamnagar": ["Devbhumi Dwarka"], "Rajkot": ["Morbi"],
           "Sabar Kantha": ["Aravalli"], "Panch Mahals": ["Mahisagar"],
           "Vadodara": ["Chhota Udaipur"]},
    "MH": {"Thane": ["Palghar"]},
    "ML": {"East Garo Hills": ["North Garo Hills"],
           "West Garo Hills": ["South West Garo Hills"],
           "West Khasi Hills": ["South West Khasi Hills"]},
    "MP": {"Shajapur": ["Agar Malwa"]},
    "PB": {"Firozpur": ["Fazilka"], "Gurdaspur": ["Pathankot"]},
    "TR": {"West Tripura": ["Khowai", "Sipahijala"], "South Tripura": ["Gomati"],
           "North Tripura": ["Unokoti"]},
    "UP": {"Moradabad": ["Sambhal"], "Ghaziabad": ["Hapur"],
           "Muzaffarnagar": ["Shamli"], "Sultanpur": ["Amethi"]},
    "WB": {"Jalpaiguri": ["Alipurduar"]},
}

# Secondary parents of carved GADM children (annotation only): (unit, child) ->
# other 2011 districts that also contributed territory.
GADM_SECONDARY: dict[tuple[str, str], list[str]] = {
    ("UP", "Amethi"): ["Rae Bareli"],
    ("GJ", "Botad"): ["Ahmadabad"],
    ("GJ", "Morbi"): ["Surendranagar", "Jamnagar"],
    ("GJ", "Mahisagar"): ["Kheda"],
    ("UP", "Sambhal"): ["Budaun"],
}


def norm(name: str) -> str:
    """Lower-case, '&' -> 'and', punctuation -> space, collapse whitespace."""
    s = str(name).lower().replace("&", " and ")
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def _compact(name: str) -> str:
    return norm(name).replace(" ", "")


def _variants(name: str, unit: str) -> set[str]:
    """Compact keys for a name: whole, bracket removed, bracket content, affixes, aliases."""
    raw = str(name).lower()
    forms = {raw}
    inner = re.findall(r"\(([^)]*)\)", raw)
    forms.add(re.sub(r"\([^)]*\)", " ", raw))
    forms.update(inner)
    for f in list(forms):
        n = norm(f)
        forms.add(re.sub(r"^the ", "", n))
        forms.add(re.sub(r" district$", "", n))
    keys = {_compact(f) for f in forms if _compact(f)}
    for group in ALIASES.get(unit, []):
        gkeys = {_compact(g) for g in group}
        if keys & gkeys:
            keys |= gkeys
    return keys


def _state_targets(state_label: str, gadm: pd.DataFrame) -> tuple[str | None, list[str]]:
    """(state-frame unit, GADM state_ids) for an NFHS v024 label."""
    unit = states.resolve(state_label)
    if unit is None or unit == "TOTAL":
        return None, []
    ids = [s for s in gadm["state_id"].unique() if states.unit_of_gid1(s) == unit]
    named = [s for s in ids
             if states.normalise(gadm.loc[gadm.state_id == s, "state_name"].iloc[0])
             == states.normalise(state_label)]
    return unit, (named or ids)


def _load_defaults():
    from pipeline.nfhs import extract

    rounds = list(extract.KR_SOURCES)
    labels = {r: extract.district_labels(r) for r in rounds}
    state_labels = {r: extract.state_labels(r) for r in rounds}
    kr = pd.read_parquet(KR_PARQUET, columns=["round", "v024", "district"]).drop_duplicates()
    n_states = kr.groupby(["round", "district"])["v024"].nunique()
    if not (n_states == 1).all():
        raise ValueError("a district code occurs in more than one state")
    district_state = {r: {int(d): int(s) for d, s in zip(g["district"], g["v024"])}
                      for r, g in kr.groupby("round")}
    gadm = pd.read_csv(GADM_CSV, usecols=["region_id", "region_name", "state_id", "state_name"])
    return labels, state_labels, district_state, gadm.drop_duplicates().reset_index(drop=True)


def _match_round(rnd, labels, state_labels, district_state, gadm) -> pd.DataFrame:
    rows = []
    for code, name in sorted(labels.items()):
        sc = district_state.get(code)
        row = {"round": rnd, "nfhs_code": int(code), "nfhs_name": name,
               "state_code": sc, "state_name": state_labels.get(sc),
               "region_id": None, "gadm_name": None, "match_type": "unmatched",
               "boundary_flag": "", "note": ""}
        if sc is None:
            row["note"] = "code absent from KR: state unknown"
            rows.append(row)
            continue
        unit, sids = _state_targets(row["state_name"], gadm)
        row["_unit"], row["_sids"] = unit, sids
        if not sids:
            row["note"] = f"state label {row['state_name']!r} not resolved to a GADM state"
            rows.append(row)
            continue
        cand = gadm[gadm.state_id.isin(sids)]
        n_nfhs_state = sum(1 for c in labels if district_state.get(c) == sc)

        def hit(gname, mtype, note=""):
            g = cand[cand.region_name == gname]
            if len(g) != 1:
                raise ValueError(f"{rnd} {name}: GADM {gname!r} not unique in {sids}")
            row.update(region_id=g.region_id.iloc[0], gadm_name=gname,
                       match_type=mtype, note=note)

        exact = cand[cand.region_name.map(norm) == norm(name)]
        keys = _variants(name, unit)
        loose = cand[cand.region_name.map(lambda g: bool(_variants(g, unit) & keys))]
        split = {norm(k): v for k, v in SPLITS.get(unit, {}).items()}.get(norm(name))
        if unit in SINGLE_POLYGON_UNITS and len(cand) == 1 and n_nfhs_state > 1:
            hit(cand.region_name.iloc[0], "manual",
                f"GADM has one polygon ({cand.region_name.iloc[0]!r}) for the whole state; "
                f"NFHS has {n_nfhs_state} districts")
            row["boundary_flag"] = "single_polygon"
        elif len(exact) == 1:
            hit(exact.region_name.iloc[0], "exact")
        elif len(loose) == 1:
            hit(loose.region_name.iloc[0], "normalised")
        elif split is not None:
            hit(split[0], "split", f"parent {split[0]}: {split[1]}")
            row["boundary_flag"] = "split"
        rows.append(row)

    out = pd.DataFrame(rows)
    # fuzzy pass: only for still-unmatched rows with a resolved state
    taken = set(out.loc[out.match_type.isin(["exact", "normalised"]), "region_id"])
    for i, row in out[out.match_type == "unmatched"].iterrows():
        sids = row.get("_sids")
        if not isinstance(sids, list) or not sids:
            continue
        cand = gadm[gadm.state_id.isin(sids) & ~gadm.region_id.isin(taken)]
        keys = _variants(row.nfhs_name, row._unit)
        scored = sorted(
            ((max(SequenceMatcher(None, k, _compact(g)).ratio() for k in keys), g, rid)
             for g, rid in zip(cand.region_name, cand.region_id)),
            reverse=True,
        )
        if not scored:
            out.at[i, "note"] = "no GADM district left in state"
            continue
        best = scored[0]
        second = scored[1][0] if len(scored) > 1 else 0.0
        if best[0] >= FUZZY_MIN and best[0] > second:
            out.loc[i, ["region_id", "gadm_name", "match_type", "note"]] = [
                best[2], best[1], "fuzzy", f"fuzzy ratio {best[0]:.3f} to {best[1]}"]
        else:
            out.at[i, "note"] = (f"no exact/alias/split match; best fuzzy {best[1]} "
                                 f"ratio {best[0]:.3f} (< {FUZZY_MIN} or tied)")
    return out


def _remnant_notes(out: pd.DataFrame) -> pd.DataFrame:
    """Flag NFHS districts mapped to a GADM parent whose carved children no code reaches."""
    reached = set(zip(out["_unit"], out["gadm_name"]))
    for i, row in out.iterrows():
        kids = GADM_CARVED.get(row.get("_unit"), {}).get(row.gadm_name)
        if not kids:
            continue
        missing = [k for k in kids if (row._unit, k) not in reached]
        if missing:
            kids_txt = ", ".join(
                k + (f" (+part of {', '.join(GADM_SECONDARY[(row._unit, k)])})"
                     if (row._unit, k) in GADM_SECONDARY else "")
                for k in missing)
            add = (f"GADM {row.gadm_name} is the post-2011 remnant; this NFHS district "
                   f"also covers GADM {kids_txt}")
            out.at[i, "note"] = f"{row.note}; {add}" if row.note else add
            out.at[i, "boundary_flag"] = "remnant"
    return out


def build_crosswalk(labels=None, state_labels=None, district_state=None, gadm=None) -> pd.DataFrame:
    """Crosswalk for every round: one row per NFHS district code.

    All inputs are dicts keyed by round (labels {code: name}, state_labels
    {v024: name}, district_state {code: v024}) plus the GADM district table;
    when omitted they are loaded from the NFHS zips, the KR subset and
    climate_district_year.csv.
    """
    if labels is None:
        labels, state_labels, district_state, gadm = _load_defaults()
    frames = []
    for rnd in labels:
        out = _match_round(rnd, labels[rnd], state_labels[rnd], district_state.get(rnd, {}), gadm)
        frames.append(_remnant_notes(out))
    xw = pd.concat(frames, ignore_index=True)
    xw["state_code"] = xw["state_code"].astype("Int64")
    return xw[COLUMNS]


def review_rows(xw: pd.DataFrame) -> pd.DataFrame:
    """Rows a reviewer must check: fuzzy, manual, split (accepted), unmatched, and
    any other row carrying a note (e.g. a GADM post-2011 remnant)."""
    flagged = xw.match_type.isin(["fuzzy", "manual", "split", "unmatched"])
    rev = xw[flagged | xw.note.fillna("").ne("")].copy()
    rev["status"] = rev.match_type.map(lambda m: "unmatched" if m == "unmatched" else "accepted")
    rev = rev.rename(columns={"note": "reason"})
    rev["reason"] = rev["reason"].fillna("").replace("", "no match")
    return rev[["round", "nfhs_code", "nfhs_name", "state_code", "state_name",
                "match_type", "boundary_flag", "status", "region_id", "gadm_name", "reason"]]


def write_crosswalk(xw: pd.DataFrame | None = None) -> pd.DataFrame:
    xw = build_crosswalk() if xw is None else xw
    xw.to_csv(CROSSWALK_CSV, index=False)
    review_rows(xw).to_csv(REVIEW_CSV, index=False)
    return xw


def crosswalk(round_id: str, path=CROSSWALK_CSV) -> dict[int, str]:
    """{NFHS district code: GADM region_id} for one round (unmatched codes omitted)."""
    xw = pd.read_csv(path)
    xw = xw[(xw["round"] == round_id) & xw["region_id"].notna()]
    if xw.empty:
        raise KeyError(f"no crosswalk rows for round {round_id!r}")
    return {int(c): str(r) for c, r in zip(xw["nfhs_code"], xw["region_id"])}


def boundary_flags(round_id: str, path=CROSSWALK_CSV) -> dict[int, str]:
    """{NFHS district code: boundary_flag} for one round ("" = boundaries coincide)."""
    xw = pd.read_csv(path, keep_default_na=False, dtype={"boundary_flag": str})
    xw = xw[xw["round"] == round_id]
    if xw.empty:
        raise KeyError(f"no crosswalk rows for round {round_id!r}")
    return {int(c): str(f) for c, f in zip(xw["nfhs_code"], xw["boundary_flag"])}


if __name__ == "__main__":
    x = write_crosswalk()
    print(x.groupby(["round", "match_type"]).size().to_string())
