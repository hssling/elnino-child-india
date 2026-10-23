"""Cluster-GPS exposure linkage for NFHS children (Paper B task 4b; primary linkage, ruling P4).

Each NFHS cluster (DHS GE file, displaced GPS point) is linked to the IMD grid:

* rain: the 0.25 deg cell that contains the point (Phase 1 R19-masked monthly grids,
  data/derived/grids/rain_monthly_<year>.npz). SPI-3 is fitted per cell with
  pipeline.climate_annual.spi on 1981-2010. A point whose own cell has no IMD data (sea /
  outside the land mask) takes the nearest valid cell within RAIN_FALLBACK_KM, flagged
  `fallback_rain`; beyond that its rain is NaN.
* temperature: the 1 deg cell that contains the point (= the nearest 1 deg cell centre).
  Monthly tmax_mean, heatwave_days and tmax_ge40_days use the Phase 1 code
  (pipeline.climate_daily: tmax_climatology on 1981-2010, heatwave_mask, monthly_counts),
  computed from the raw daily .grd files. tmax_anom = tmax_mean minus the cell's 1981-2010
  calendar-month mean. A cell with no data falls back to the nearest valid 1 deg cell,
  flagged `fallback_temp`.

Cells: a point belongs to the cell whose half-open interval [edge, edge + step) contains it,
so a point exactly on a cell edge goes to the cell north / east of the edge.

Regions: `region_id` is the GADM 4.1 district containing the point (point in polygon; a
point on a shared border takes the lowest region_id). A point outside every polygon takes
the nearest polygon within NEAREST_MAX_KM (`nearest`), else the Task 3 crosswalk district of
its NFHS district code (`district_fallback`). MIS clusters (no coordinates; (0,0) in GE)
take the crosswalk district of the GE district code (DHSREGCO = the KR district code) and
no point (`district`).

Masking (no partial values, missing is NaN):
* QC (ruling R22: artefacts never enter the gamma fits): before the SPI fit a cell's rain is
  NaN in the months masked (pipeline/reference/qc_reviewed.csv) for ANY of: the cell's
  dominant district (largest overlap area) and the own (point-in-polygon) districts of every
  cluster using the cell as its point cell or as a buffer cell. A cluster's SPI-3 is then NaN
  in every month whose 3-month accumulation touches a masked month. Flags: `qc_masked` = the
  cluster's own district is masked that month (raw month); `qc_cell_masked` = the point cell's
  SPI-3 accumulation (t-2..t) touches a pre-fit-masked cell month (any of the reasons above).
  Temperature is never QC-masked.
* islands (ruling P7): clusters in Andaman & Nicobar and Lakshadweep have NaN rain/SPI;
  their temperature is NaN where the 1 deg cell is a fallback.
* MIS clusters (and any (0,0) point, whatever its SOURCE) take the district-month values of
  exposure.district_monthly_metrics (which applies the same QC and island rules), flagged
  `mis_district` (window output: <window>_mis_district_any).

Sensitivity: `spi3_buf10` is the area-weighted mean of cell SPI-3 over the 0.25 deg cells
intersecting a 10 km geodesic buffer around rural points (weights renormalised over cells
with IMD data; NaN if any such cell is NaN that month). It averages standardised values (a
mean of cell SPIs, variance slightly below 1), not the SPI of buffer-mean rain; and its
uniform-area weighting is a neighbourhood sensitivity, not an expectation over the DHS
displacement kernel (whose density concentrates near the reported point). Urban clusters
(displacement <= 2 km) use the point value; so do points whose buffer has no IMD cell, and
MIS clusters carry their district value.

Distances: the nearest-polygon tolerance (NEAREST_MAX_KM) is measured in the Albers
equal-area projection of grid_weights.EQUAL_AREA (not equidistant: a few % scale error near
India's extremes, immaterial at 10 km); the nearest-cell fallback uses an equirectangular km
distance; buffers are geodesic (pyproj Geod, WGS84).

CONTROLLED DATA: GE files are opened read-only straight from the zip (/vsizip/); cluster
coordinates and cluster-level exposures are written only to the git-ignored
data/external/nfhs_subset/ and never to committed files or reports.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd

from pipeline import climate_daily as cd
from pipeline import imd_grid
from pipeline.climate_annual import BASE, QC_MONTHS, QC_REVIEWED, qc_mask_for_level, spi
from pipeline.grid_weights import EQUAL_AREA, compute_weights
from pipeline.nfhs.exposure import EXCLUDED_STATES, WINDOW_VARIABLES, WINDOWS, window_exposure
from pipeline.paths import DERIVED, EXTERNAL, RAW

GE_SOURCES = {
    "nfhs4": (Path(EXTERNAL / "nfhs" / r"IAGE71FL.zip"), "IAGE71FL.shp"),
    "nfhs5": (Path(EXTERNAL / "nfhs" / r"IAGE7AFL.zip"), "IAGE7AFL.shp"),
}
GE_COLUMNS = ["DHSCLUST", "LATNUM", "LONGNUM", "URBAN_RURA", "SOURCE", "DHSREGCO"]
SUBSET_DIR = EXTERNAL / "nfhs_subset"
CLUSTERS_PARQUET = SUBSET_DIR / "clusters.parquet"
CLUSTER_MONTHLY_PARQUET = SUBSET_DIR / "cluster_cell_monthly.parquet"
GRIDS_DIR = DERIVED / "grids"
WEIGHTS_R025_DISTRICT = DERIVED / "_parts" / "weights_r025_district.csv"
GADM_DISTRICTS_CSV = DERIVED / "climate_district_year.csv"
KR_PARQUET = SUBSET_DIR / "kr_harmonised.parquet"

YEARS = (1981, 2021)          # normals (1981-2010) + every NFHS-4/5 birth window
NEAREST_MAX_KM = 10.0         # point outside every GADM polygon: nearest polygon within this
RAIN_FALLBACK_KM = 30.0       # own 0.25 deg cell without data: nearest valid cell within this
BUFFER_KM = 10.0
BUFFER_SEGMENTS = 64
EARTH_KM_PER_DEG = 111.195

CLUSTER_COLUMNS = ["round", "v001", "lat", "lon", "urban", "source", "region_id", "district",
                   "region_method", "region_xw", "boundary_flag", "state_name"]
VALUE_COLUMNS = ["spi3", "spi3_buf10", "tmax_anom", "heatwave_days", "tmax_ge40_days"]
FLAG_COLUMNS = ["qc_masked", "qc_cell_masked", "fallback_rain", "fallback_temp", "island_proxy", "mis_district"]
CLUSTER_KEY = ["round", "v001"]
CLUSTER_WINDOW_VARIABLES = {**WINDOW_VARIABLES, "spi3_buf10": ("spi3_buf10", "mean")}
CLUSTER_WINDOW_FLAGS = {"fallback_rain": "fallback_rain_any", "fallback_temp": "fallback_temp_any",
                        "mis_district": "mis_district_any"}
OUTPUT_PAD_MONTHS = 12        # cluster-month output starts this many months before the earliest birth


# ---------------------------------------------------------------- grid cells

def cell_index(lon, lat, grid) -> tuple[np.ndarray, np.ndarray]:
    """(ilat, ilon) of the cell containing each point; -1 if off the grid or NaN.

    Cells are half-open [centre - step/2, centre + step/2): a point on an edge belongs to the
    cell north / east of it."""
    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)

    def idx(x, x0, n):
        q = (x - (x0 - grid.step / 2)) / grid.step
        q = np.round(q, 9)                         # exact edges stay exact despite float error
        with np.errstate(invalid="ignore"):
            k = np.floor(q)
        ok = np.isfinite(k) & (k >= 0) & (k < n)
        return np.where(ok, k, -1).astype(int)

    ilat, ilon = idx(lat, grid.lat0, grid.nlat), idx(lon, grid.lon0, grid.nlon)
    off = (ilat < 0) | (ilon < 0)
    return np.where(off, -1, ilat), np.where(off, -1, ilon)


def _km(lon1, lat1, lon2, lat2):
    """Equirectangular distance in km (accurate to < 0.1% at these scales)."""
    dlon = (np.asarray(lon2) - lon1) * np.cos(np.radians((np.asarray(lat2) + lat1) / 2))
    return EARTH_KM_PER_DEG * np.hypot(dlon, np.asarray(lat2) - lat1)


def nearest_valid_cell(lon, lat, ilat, ilon, valid, grid, max_km=None):
    """(ilat, ilon, fallback): the point's own cell if `valid` there, else the nearest valid cell
    centre (within `max_km` if given; -1 beyond it). fallback is True where a neighbour is used."""
    lon, lat = np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)
    ilat, ilon = np.asarray(ilat).copy(), np.asarray(ilon).copy()
    own = (ilat >= 0) & (ilon >= 0)
    own[own] = valid[ilat[own], ilon[own]]
    fb = np.zeros(len(lon), dtype=bool)
    vi, vj = np.nonzero(valid)
    clat, clon = grid.lats[vi], grid.lons[vj]
    for k in np.flatnonzero(~own):
        ilat[k] = ilon[k] = -1
        if not np.isfinite(lon[k]) or not np.isfinite(lat[k]) or not len(vi):
            continue
        d = _km(lon[k], lat[k], clon, clat)
        b = int(np.argmin(d))
        if max_km is None or d[b] <= max_km:
            ilat[k], ilon[k], fb[k] = vi[b], vj[b], True
    return ilat, ilon, fb


# ---------------------------------------------------------------- GE clusters and regions

def read_ge(round_id: str, zip_path=None, member=None) -> pd.DataFrame:
    """The GE attribute table (GE_COLUMNS), read-only straight from the zip; never extracted."""
    import pyogrio
    if round_id not in GE_SOURCES:
        raise ValueError(f"unknown round {round_id!r}")
    z, m = GE_SOURCES[round_id]
    path = "/vsizip/" + str(Path(zip_path or z)).replace("\\", "/") + "/" + (member or m)
    return pyogrio.read_dataframe(path, columns=GE_COLUMNS, read_geometry=False)


def normalise_ge(ge: pd.DataFrame, round_id: str) -> pd.DataFrame:
    """round, v001, lat, lon, urban, source, district; MIS / (0,0) points get NaN coordinates."""
    src = ge["SOURCE"].astype(str).str.strip().str.upper()
    lat = pd.to_numeric(ge["LATNUM"], errors="coerce").astype(float)
    lon = pd.to_numeric(ge["LONGNUM"], errors="coerce").astype(float)
    nopoint = (src == "MIS") | ((lat == 0) & (lon == 0))
    out = pd.DataFrame({
        "round": round_id,
        "v001": pd.to_numeric(ge["DHSCLUST"]).astype(int).to_numpy(),
        "lat": lat.where(~nopoint).to_numpy(),
        "lon": lon.where(~nopoint).to_numpy(),
        "urban": (ge["URBAN_RURA"].astype(str).str.strip().str.upper() == "U").to_numpy(),
        "source": src.to_numpy(),
        "district": pd.to_numeric(ge["DHSREGCO"], errors="coerce").astype("Int64").to_numpy(),
    })
    if out.v001.duplicated().any():
        raise ValueError(f"{round_id}: duplicate DHSCLUST in GE file")
    return out


def assign_regions(clusters: pd.DataFrame, districts: gpd.GeoDataFrame, xw: dict,
                   max_km=NEAREST_MAX_KM) -> pd.DataFrame:
    """Copy of `clusters` with region_id, region_method (pip | nearest | district |
    district_fallback) and region_xw (Task 3 crosswalk district of the cluster's district code).
    `xw` is {round: {NFHS district code: region_id}}."""
    out = clusters.reset_index(drop=True).copy()
    out["region_xw"] = [xw.get(r, {}).get(int(d)) if pd.notna(d) else None
                        for r, d in zip(out["round"], out["district"])]
    out["region_id"] = None
    out["region_method"] = None
    has_pt = out["lat"].notna() & out["lon"].notna()
    reg = districts[["region_id", "geometry"]].to_crs("EPSG:4326").reset_index(drop=True)
    if has_pt.any():
        pts = gpd.GeoDataFrame({"_row": np.flatnonzero(has_pt)},
                               geometry=gpd.points_from_xy(out.loc[has_pt, "lon"], out.loc[has_pt, "lat"]),
                               crs="EPSG:4326")
        j = gpd.sjoin(pts, reg, how="inner", predicate="intersects")
        pip = j.groupby("_row")["region_id"].min()          # shared border: lowest region_id
        out.loc[pip.index, "region_id"] = pip.to_numpy()
        out.loc[pip.index, "region_method"] = "pip"
        rest = pts[~pts["_row"].isin(pip.index)]
        if len(rest):
            nn = gpd.sjoin_nearest(rest.to_crs(EQUAL_AREA), reg.to_crs(EQUAL_AREA), how="inner",
                                   distance_col="_dist")
            nn = nn.sort_values(["_row", "_dist", "region_id"]).drop_duplicates("_row")
            nn = nn[nn["_dist"] <= max_km * 1000]
            out.loc[nn["_row"].to_numpy(), "region_id"] = nn["region_id"].to_numpy()
            out.loc[nn["_row"].to_numpy(), "region_method"] = "nearest"
    left = out["region_id"].isna()
    out.loc[left & ~has_pt, "region_method"] = "district"
    out.loc[left & has_pt, "region_method"] = "district_fallback"
    out.loc[left, "region_id"] = out.loc[left, "region_xw"]
    return out


def _gadm_table() -> pd.DataFrame:
    return pd.read_csv(GADM_DISTRICTS_CSV, usecols=["region_id", "region_name", "state_id", "state_name"]
                       ).drop_duplicates("region_id").reset_index(drop=True)


def cluster_table(rounds=("nfhs4", "nfhs5"), districts: gpd.GeoDataFrame | None = None,
                  save=True, path=CLUSTERS_PARQUET) -> pd.DataFrame:
    """One row per GE cluster: CLUSTER_COLUMNS (coordinates included: controlled data, saved only
    to the git-ignored `path`)."""
    from pipeline import boundaries
    from pipeline.nfhs import districts as dmod
    if districts is None:
        districts = boundaries.load("district")
    frames = []
    for r in rounds:
        ge = normalise_ge(read_ge(r), r)
        out = assign_regions(ge, districts, {r: dmod.crosswalk(r)})
        flags = dmod.boundary_flags(r)
        out["boundary_flag"] = [flags.get(int(d), "") if pd.notna(d) else "" for d in out["district"]]
        frames.append(out)
    cl = pd.concat(frames, ignore_index=True)
    gadm = _gadm_table()
    cl["state_name"] = cl["region_id"].map(dict(zip(gadm["region_id"], gadm["state_name"])))
    cl = cl[CLUSTER_COLUMNS]
    if save:
        path.parent.mkdir(parents=True, exist_ok=True)
        cl.to_parquet(path, index=False)
    return cl


def load_or_build_clusters(path=CLUSTERS_PARQUET) -> pd.DataFrame:
    return pd.read_parquet(path) if Path(path).exists() else cluster_table(path=path)


# ---------------------------------------------------------------- buffers

def buffer_polygons(lon, lat, km=BUFFER_KM, n=BUFFER_SEGMENTS):
    """Geodesic circles (WGS84) of radius `km` as shapely polygons in lon/lat."""
    import shapely
    from pyproj import Geod
    lon, lat = np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)
    az = np.linspace(0.0, 360.0, n, endpoint=False)
    lo, la, _ = Geod(ellps="WGS84").fwd(np.repeat(lon, n), np.repeat(lat, n), np.tile(az, len(lon)),
                                        np.full(len(lon) * n, km * 1000.0))
    ring = np.stack([lo, la], axis=1).reshape(len(lon), n, 2)
    ring = np.concatenate([ring, ring[:, :1]], axis=1)
    return shapely.polygons(ring)


def buffer_weights(points: pd.DataFrame, grid, km=BUFFER_KM) -> pd.DataFrame:
    """Area weights of the `grid` cells intersecting a `km` buffer around each point.

    Returns row (position in `points`), ilat, ilon, area_km2, weight (sums to 1 per row).
    Uses grid_weights.compute_weights (intersection in lon/lat, area in equal-area projection)."""
    ok = np.flatnonzero(points["lat"].notna().to_numpy() & points["lon"].notna().to_numpy())
    if not len(ok):
        return pd.DataFrame(columns=["row", "ilat", "ilon", "area_km2", "weight"])
    polys = buffer_polygons(points["lon"].to_numpy()[ok], points["lat"].to_numpy()[ok], km)
    g = gpd.GeoDataFrame({"region_id": ok}, geometry=polys, crs="EPSG:4326")
    w = compute_weights(grid, g).rename(columns={"region_id": "row"})
    w["row"] = w["row"].astype(int)
    return w.sort_values(["row", "ilat", "ilon"]).reset_index(drop=True)


def buffer_mean(values: np.ndarray, weights: pd.DataFrame, n_points: int, valid: np.ndarray) -> np.ndarray:
    """(T, n_points) area-weighted means of values (T, nlat, nlon) over each point's buffer cells.

    Cells where `valid` is False (no IMD data at all) are dropped and the weights renormalised; if
    any remaining cell is NaN at t the mean is NaN (no partial means); no valid cell -> NaN."""
    from scipy import sparse
    t, nlat, nlon = values.shape
    w = weights[valid[weights["ilat"].to_numpy(), weights["ilon"].to_numpy()]]
    tot = w.groupby("row")["weight"].transform("sum")
    wn = (w["weight"] / tot).to_numpy()
    cols = w["ilat"].to_numpy() * nlon + w["ilon"].to_numpy()
    W = sparse.csr_matrix((wn, (w["row"].to_numpy(), cols)), shape=(n_points, nlat * nlon))
    B = sparse.csr_matrix((np.ones(len(w)), (w["row"].to_numpy(), cols)), shape=(n_points, nlat * nlon))
    v = values.reshape(t, -1)
    nan = np.isnan(v)
    num = (W @ np.where(nan, 0.0, v).T).T
    bad = (B @ nan.T.astype(float)).T > 0
    has = np.zeros(n_points, dtype=bool)
    has[np.unique(w["row"].to_numpy())] = True
    return np.where(bad | ~has[None, :], np.nan, num)


# ---------------------------------------------------------------- climate inputs

def load_rain(years=YEARS, grids_dir=GRIDS_DIR) -> tuple[np.ndarray, pd.PeriodIndex]:
    """(T, nlat, nlon) monthly 0.25 deg rain (R19-masked) and its month index."""
    parts = []
    for y in range(years[0], years[1] + 1):
        p = Path(grids_dir) / f"rain_monthly_{y}.npz"
        if not p.exists():
            raise FileNotFoundError(p)
        with np.load(p) as z:
            parts.append(z["rain"].astype(np.float32))
    return np.concatenate(parts), pd.period_range(f"{years[0]}-01", f"{years[1]}-12", freq="M")


def temperature_monthly(years=YEARS, raw_dir=RAW / "imd", base=BASE) -> dict[str, np.ndarray]:
    """1 deg monthly tmax_mean, heatwave_days, tmax_ge40_days (T, nlat, nlon), Phase 1 logic
    (build_climate.build_temperature): counts NaN where tmax_mean is NaN."""
    def read(y):
        p = Path(raw_dir) / f"tmax_{y}.grd"
        if not p.exists():
            raise FileNotFoundError(p)
        return imd_grid.read("tmax", p, y)

    clim = cd.tmax_climatology({y: read(y) for y in range(base[0], base[1] + 1)})
    out = {"tmax_mean": [], "heatwave_days": [], "tmax_ge40_days": []}
    for y in range(years[0], years[1] + 1):
        tmax = read(y)
        mean = cd.daily_to_monthly(tmax, y, "mean")
        with np.errstate(invalid="ignore"):
            ge40 = cd.monthly_counts(tmax >= 40, y).astype(float)
        hw = cd.monthly_counts(cd.heatwave_mask(tmax, y, clim), y).astype(float)
        bad = np.isnan(mean)
        ge40[bad] = np.nan
        hw[bad] = np.nan
        out["tmax_mean"].append(mean)
        out["heatwave_days"].append(hw)
        out["tmax_ge40_days"].append(ge40)
    return {k: np.concatenate(v).astype(np.float32) for k, v in out.items()}


def qc_months(path=QC_REVIEWED) -> pd.DataFrame:
    """Reviewed district masks expanded to (region_id, year, month)."""
    q = qc_mask_for_level("district", path)
    rows = [(str(m.region_id), int(m.year), mo) for m in q.itertuples() for mo in QC_MONTHS[m.variable]]
    return pd.DataFrame(rows, columns=["region_id", "year", "month"]).drop_duplicates().reset_index(drop=True)


def dominant_region(grid=imd_grid.RAIN, weights_path=WEIGHTS_R025_DISTRICT) -> np.ndarray:
    """(nlat, nlon) object array: district with the largest overlap area of each cell (None if none)."""
    if Path(weights_path).exists():
        w = pd.read_csv(weights_path, dtype={"region_id": str})
    else:
        from pipeline import boundaries
        w = compute_weights(grid, boundaries.load("district"))
    w = w.sort_values(["ilat", "ilon", "area_km2", "region_id"], ascending=[True, True, False, True])
    w = w.drop_duplicates(["ilat", "ilon"])
    out = np.full((grid.nlat, grid.nlon), None, dtype=object)
    out[w.ilat.to_numpy(), w.ilon.to_numpy()] = w.region_id.to_numpy()
    return out


def output_years(kr_path=KR_PARQUET) -> dict[str, tuple[int, int]]:
    """Per round: (year of earliest birth - OUTPUT_PAD_MONTHS, year of latest interview).

    The pad covers every birth-relative window down to -OUTPUT_PAD_MONTHS (T1 = -9, or -10 with
    the conception shift); cluster_window_exposure refuses windows reaching further back."""
    kr = pd.read_parquet(kr_path, columns=["round", "b3", "v008"])
    out = {}
    for r, g in kr.groupby("round"):
        lo, hi = int(g["b3"].min()) - OUTPUT_PAD_MONTHS, int(g["v008"].max())
        out[r] = ((lo - 1) // 12 + 1900, (hi - 1) // 12 + 1900)
    return out


# ---------------------------------------------------------------- cluster-month series

def _region_mask_lookup(qc: pd.DataFrame, months: pd.PeriodIndex) -> dict[str, np.ndarray]:
    """{region_id: (T,) bool month mask} for the regions in `qc`."""
    out: dict[str, np.ndarray] = {}
    if qc is None or not len(qc):
        return out
    t = pd.PeriodIndex.from_fields(year=qc["year"].astype(int), month=qc["month"].astype(int), freq="M")
    pos = months.get_indexer(t)
    for r, p in zip(qc["region_id"].astype(str), pos):
        if p >= 0:
            out.setdefault(r, np.zeros(len(months), dtype=bool))[p] = True
    return out


def _region_mask_matrix(regions, qc: pd.DataFrame, months: pd.PeriodIndex) -> np.ndarray:
    """(n, T) bool: month masked for each entry of `regions`."""
    regions = pd.Series(regions, dtype=object).astype("string").fillna("")
    out = np.zeros((len(regions), len(months)), dtype=bool)
    if qc is None or not len(qc):
        return out
    t = pd.PeriodIndex.from_fields(year=qc["year"].astype(int), month=qc["month"].astype(int), freq="M")
    pos = months.get_indexer(t)
    q = pd.DataFrame({"region_id": qc["region_id"].astype(str).to_numpy(), "pos": pos})
    q = q[q.pos >= 0]
    by = {r: g.pos.to_numpy() for r, g in q.groupby("region_id")}
    for r, rows in pd.Series(np.arange(len(regions))).groupby(regions.to_numpy()):
        if r in by:
            out[np.ix_(rows.to_numpy(), by[r])] = True
    return out


def _tail3(mask: np.ndarray) -> np.ndarray:
    """True at t if mask at t, t-1 or t-2 (an SPI-3 accumulation touching a masked month)."""
    out = mask.copy()
    out[..., 1:] |= mask[..., :-1]
    out[..., 2:] |= mask[..., :-2]
    return out


def build_cluster_monthly(clusters: pd.DataFrame, rain: np.ndarray, temp: dict, months: pd.PeriodIndex,
                          district_metrics: pd.DataFrame, qc_months: pd.DataFrame, cell_region: np.ndarray,
                          island_regions: set, years_by_round: dict, buffer_w: pd.DataFrame | None = None,
                          rain_grid=imd_grid.RAIN, temp_grid=imd_grid.TEMP, base=BASE,
                          rain_fallback_km=RAIN_FALLBACK_KM) -> pd.DataFrame:
    """Cluster-month exposures keyed round, v001, year, month (see module doc).

    `rain` (T, nlat, nlon) and `temp` {tmax_mean, heatwave_days, tmax_ge40_days} (T, nlat, nlon) are
    on `months`; `cell_region` gives each rain cell's dominant district (QC masks before the SPI fit);
    `buffer_w` (row = position in `clusters`, ilat, ilon, weight) gives the rural 10 km buffers;
    `years_by_round` limits the output months per round."""
    cl = clusters.reset_index(drop=True)
    n, T = len(cl), len(months)
    is_mis = (cl["source"].astype(str) == "MIS").to_numpy() | cl["lat"].isna().to_numpy() | cl["lon"].isna().to_numpy()
    lon = np.where(is_mis, np.nan, cl["lon"].to_numpy(dtype=float))
    lat = np.where(is_mis, np.nan, cl["lat"].to_numpy(dtype=float))
    region = cl["region_id"].astype(object).to_numpy()
    island = np.array([r in island_regions for r in region], dtype=bool)

    # rain cells and SPI-3 per cell on QC-masked rain
    valid_r = np.isfinite(rain).any(axis=0)
    ri, rj = cell_index(lon, lat, rain_grid)
    ri, rj, fb_r = nearest_valid_cell(lon, lat, ri, rj, valid_r, rain_grid, max_km=rain_fallback_km)
    need = {(i, j) for i, j in zip(ri, rj) if i >= 0}
    if buffer_w is not None and len(buffer_w):
        bw = buffer_w[valid_r[buffer_w.ilat.to_numpy(), buffer_w.ilon.to_numpy()]]
        need |= set(zip(bw.ilat.to_numpy(), bw.ilon.to_numpy()))
    need = sorted(need)
    # pre-fit masks (R22): union of the cell's dominant district and the own districts of every
    # cluster using the cell (point cell or buffer cell)
    users: dict = {c: {cell_region[c]} for c in need}
    for k in np.flatnonzero((ri >= 0) & ~is_mis):
        users[(ri[k], rj[k])].add(region[k])
    if buffer_w is not None and len(buffer_w):
        for r, i, j in zip(bw.row.to_numpy(), bw.ilat.to_numpy(), bw.ilon.to_numpy()):
            if not is_mis[r]:
                users[(i, j)].add(region[r])
    rmask = _region_mask_lookup(qc_months, months)
    none = np.zeros(T, dtype=bool)
    cell_mask = np.zeros((len(need), T), dtype=bool)
    for k, c in enumerate(need):
        for r in users[c]:
            if r is not None:
                cell_mask[k] |= rmask.get(str(r), none)
    cell_pos = {c: k for k, c in enumerate(need)}
    spi_grid = np.full(rain.shape, np.nan, dtype=np.float32)
    for k, (i, j) in enumerate(need):
        s = rain[:, i, j].astype(float)
        s[cell_mask[k]] = np.nan
        spi_grid[:, i, j] = spi(pd.Series(s, index=months), 3, base).to_numpy()
    qc_cell = np.zeros((n, T), dtype=bool)
    for k in np.flatnonzero(ri >= 0):
        qc_cell[k] = cell_mask[cell_pos[(ri[k], rj[k])]]
    qc_cell = _tail3(qc_cell) & ~is_mis[:, None]

    spi3 = np.full((n, T), np.nan)
    ok = ri >= 0
    spi3[ok] = spi_grid[:, ri[ok], rj[ok]].T
    buf = spi3.copy()
    if buffer_w is not None and len(buffer_w):
        bm = buffer_mean(spi_grid, buffer_w, n, valid_r).T            # (n, T)
        has_buf = np.zeros(n, dtype=bool)
        has_buf[np.unique(buffer_w.row[valid_r[buffer_w.ilat.to_numpy(), buffer_w.ilon.to_numpy()]])] = True
        use = has_buf & ~cl["urban"].astype(bool).to_numpy() & ~is_mis
        buf[use] = bm[use]

    qc_cl = _region_mask_matrix(region, qc_months, months)
    blank = _tail3(qc_cl) | island[:, None]
    spi3[blank] = np.nan
    buf[blank] = np.nan

    # temperature
    tm = temp["tmax_mean"]
    valid_t = np.isfinite(tm).any(axis=0)
    ti, tj = cell_index(lon, lat, temp_grid)
    ti, tj, fb_t = nearest_valid_cell(lon, lat, ti, tj, valid_t, temp_grid, max_km=None)
    in_base = np.asarray((months.year >= base[0]) & (months.year <= base[1]))
    normal = np.full((12,) + tm.shape[1:], np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for m in range(12):
            normal[m] = np.nanmean(tm[in_base & (np.asarray(months.month) == m + 1)], axis=0)
    anom = tm - normal[np.asarray(months.month) - 1]
    tvals = {}
    okt = ti >= 0
    for name, arr in (("tmax_anom", anom), ("heatwave_days", temp["heatwave_days"]),
                      ("tmax_ge40_days", temp["tmax_ge40_days"])):
        v = np.full((n, T), np.nan)
        v[okt] = arr[:, ti[okt], tj[okt]].T
        v[island & fb_t] = np.nan
        tvals[name] = v

    # MIS: district-month values
    dm = district_metrics.copy()
    dm["region_id"] = dm["region_id"].astype(str)
    dm_idx = pd.MultiIndex.from_arrays([dm.region_id, pd.PeriodIndex.from_fields(
        year=dm.year.astype(int), month=dm.month.astype(int), freq="M")])
    dmi = pd.Series(np.arange(len(dm)), index=dm_idx)

    frames = []
    for rnd, rows in pd.Series(np.arange(n)).groupby(cl["round"].to_numpy()):
        rows = rows.to_numpy()
        y0, y1 = years_by_round[rnd]
        tsel = np.flatnonzero((months.year >= y0) & (months.year <= y1))
        nt = len(tsel)
        rep = np.repeat(rows, nt)
        tt = np.tile(tsel, len(rows))
        df = pd.DataFrame({"round": rnd, "v001": cl["v001"].to_numpy()[rep].astype(np.int32),
                           "year": np.asarray(months.year)[tt].astype(np.int16),
                           "month": np.asarray(months.month)[tt].astype(np.int8)})
        df["spi3"] = spi3[rep, tt]
        df["spi3_buf10"] = buf[rep, tt]
        for k, v in tvals.items():
            df[k] = v[rep, tt]
        df["qc_masked"] = qc_cl[rep, tt]
        df["qc_cell_masked"] = qc_cell[rep, tt]
        df["fallback_rain"] = fb_r[rep]
        df["fallback_temp"] = fb_t[rep]
        df["island_proxy"] = island[rep]
        df["mis_district"] = is_mis[rep]
        mis = df["mis_district"].to_numpy()
        if mis.any():
            key = pd.MultiIndex.from_arrays([pd.Series(region[rep[mis]]).astype(str).to_numpy(),
                                             months[tt[mis]]])
            pos = dmi.reindex(key).to_numpy()
            hit = ~np.isnan(pos)
            src = dm.iloc[pos[hit].astype(int)]
            idx = np.flatnonzero(mis)
            for c in ("spi3", "tmax_anom", "heatwave_days", "tmax_ge40_days"):
                vals = np.full(mis.sum(), np.nan)
                vals[hit] = src[c].to_numpy(dtype=float)
                df.loc[idx, c] = vals
            df.loc[idx, "spi3_buf10"] = df.loc[idx, "spi3"].to_numpy()
            for c in ("qc_masked", "fallback_rain", "fallback_temp", "island_proxy"):
                vals = np.zeros(mis.sum(), dtype=bool)
                if c in src:
                    vals[hit] = src[c].astype("boolean").fillna(False).to_numpy(dtype=bool)
                df.loc[idx, c] = vals
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    for c in VALUE_COLUMNS:
        out[c] = out[c].astype(np.float32)
    for c in FLAG_COLUMNS:
        out[c] = out[c].astype(bool)
    return out[CLUSTER_KEY + ["year", "month"] + VALUE_COLUMNS + FLAG_COLUMNS].sort_values(
        CLUSTER_KEY + ["year", "month"]).reset_index(drop=True)


def cluster_cell_monthly(clusters: pd.DataFrame | None = None, years=YEARS, save=True,
                         path=CLUSTER_MONTHLY_PARQUET, buffer=True) -> pd.DataFrame:
    """Real-data cluster-month exposures (controlled; saved only to the git-ignored `path`)."""
    from pipeline.nfhs.exposure import district_monthly_metrics
    cl = load_or_build_clusters() if clusters is None else clusters
    cl = cl.reset_index(drop=True)
    rain, months = load_rain(years)
    temp = temperature_monthly(years)
    if len(temp["tmax_mean"]) != len(months):
        raise ValueError("rain and temperature month axes differ")
    gadm = _gadm_table()
    islands = set(gadm.loc[gadm.state_name.isin(EXCLUDED_STATES), "region_id"])
    bw = None
    if buffer:
        rural = (~cl["urban"].astype(bool)) & (cl["source"] != "MIS") & cl["lat"].notna()
        pos = np.flatnonzero(rural.to_numpy())
        bw = buffer_weights(cl.iloc[pos].reset_index(drop=True), imd_grid.RAIN)
        bw["row"] = pos[bw["row"].to_numpy()]
    out = build_cluster_monthly(cl, rain, temp, months, district_monthly_metrics(), qc_months(),
                                dominant_region(), islands, output_years(), buffer_w=bw)
    if save:
        path.parent.mkdir(parents=True, exist_ok=True)
        out.to_parquet(path, index=False, compression="zstd")
    return out


def load_or_build_monthly(clusters=None, path=CLUSTER_MONTHLY_PARQUET) -> pd.DataFrame:
    return pd.read_parquet(path) if Path(path).exists() else cluster_cell_monthly(clusters, path=path)


# ---------------------------------------------------------------- windows, agreement, coverage

def cluster_window_exposure(children: pd.DataFrame, metrics: pd.DataFrame, windows: dict = WINDOWS,
                            variables: dict = CLUSTER_WINDOW_VARIABLES, flags: dict = CLUSTER_WINDOW_FLAGS
                            ) -> pd.DataFrame:
    """Task 4 window_exposure on cluster-month metrics, key ["round", "v001"].

    Raises ValueError for a birth-relative window starting before -OUTPUT_PAD_MONTHS (the saved
    cluster-month table does not reach that far back, so it would silently be NaN)."""
    for name, w in windows.items():
        anchor = w[2] if len(w) == 3 else "birth"
        if anchor == "birth" and int(w[0]) < -OUTPUT_PAD_MONTHS:
            raise ValueError(f"window {name!r} starts {w[0]} months from birth; cluster-month output "
                             f"is padded only {OUTPUT_PAD_MONTHS} months (raise OUTPUT_PAD_MONTHS and rebuild)")
    return window_exposure(children, metrics, windows, key=CLUSTER_KEY, variables=variables, flags=flags)


def _carved_regions() -> dict[str, set[str]]:
    """{GADM parent region_id: region_ids of its post-2011 carved children} (districts.GADM_CARVED)."""
    from pipeline.nfhs.districts import GADM_CARVED
    from pipeline.outcomes import states
    g = _gadm_table()
    g["unit"] = g["state_id"].map(states.unit_of_gid1)
    look = {(u, n): r for u, n, r in zip(g.unit, g.region_name, g.region_id)}
    out = {}
    for unit, parents in GADM_CARVED.items():
        for parent, kids in parents.items():
            p = look.get((unit, parent))
            if p is not None:
                out[p] = {look[(unit, k)] for k in kids if (unit, k) in look}
    return out


def agreement(clusters: pd.DataFrame, accepted: dict | None = None) -> pd.DataFrame:
    """Share of clusters with a point whose GADM district (point in polygon / nearest) agrees with
    the Task 3 crosswalk district of their NFHS district code.

    Agreement: region_id == region_xw, or, for a crosswalk row flagged `remnant` (2011-boundary
    NFHS district larger than the GADM 4.1 polygon), region_id in the parent's carved children.
    `accepted` maps (round, region_xw) -> extra region_ids accepted for remnant rows (default:
    districts.GADM_CARVED). Groups: all, non_remnant, remnant. Aggregates only."""
    c = clusters[clusters["source"] != "MIS"]
    if "region_method" in c:
        c = c[c["region_method"].isin(["pip", "nearest"])]
    if accepted is None:
        carved = _carved_regions()
        accepted = {(r, p): kids for r in c["round"].unique() for p, kids in carved.items()}
    rem = c["boundary_flag"].fillna("").eq("remnant").to_numpy()
    exact = (c["region_id"].astype(str) == c["region_xw"].astype(str)).to_numpy()
    extra = np.array([rid in accepted.get((r, x), set())
                      for r, x, rid in zip(c["round"], c["region_xw"], c["region_id"])], dtype=bool)
    agree = exact | (rem & extra)
    d = pd.DataFrame({"round": c["round"].to_numpy(), "rem": rem, "agree": agree, "exact": exact})
    rows = []
    for rnd, g in d.groupby("round"):
        for name, sel in (("all", np.ones(len(g), bool)), ("non_remnant", ~g.rem.to_numpy()),
                          ("remnant", g.rem.to_numpy())):
            s = g[sel]
            rows.append({"round": rnd, "group": name, "n": len(s),
                         "agree": s.agree.mean() if len(s) else np.nan,
                         "agree_exact": s.exact.mean() if len(s) else np.nan})
    return pd.DataFrame(rows, columns=["round", "group", "n", "agree", "agree_exact"])


def coverage_check(clusters: pd.DataFrame | None = None, metrics: pd.DataFrame | None = None,
                   kids: pd.DataFrame | None = None) -> pd.DataFrame:
    """Share of children with non-missing cluster-linked window exposures, by round and
    A&N/Lakshadweep group (cluster's GADM state). Aggregates only."""
    from pipeline.nfhs.exposure import load_children
    cl = load_or_build_clusters() if clusters is None else clusters
    metrics = load_or_build_monthly(cl) if metrics is None else metrics
    kids = load_children() if kids is None else kids
    w = cluster_window_exposure(kids, metrics, {k: WINDOWS[k] for k in ("T1", "T2", "T3", "recent3")})
    w = w[[c for c in w if not c.endswith("_any")]]
    state = pd.MultiIndex.from_frame(kids[CLUSTER_KEY]).map(
        pd.Series(cl["state_name"].to_numpy(), index=pd.MultiIndex.from_frame(cl[CLUSTER_KEY])).to_dict())
    grp = np.where(pd.Series(state).isin(EXCLUDED_STATES), "A&N/Lakshadweep", "outside A&N/Lakshadweep")
    df = w.notna().assign(round=kids["round"].to_numpy(), group=grp)
    agg = df.groupby(["round", "group"]).mean()
    agg.insert(0, "n_children", df.groupby(["round", "group"]).size())
    return agg.reset_index()
