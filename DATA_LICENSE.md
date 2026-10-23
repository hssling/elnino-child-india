# Data licence

The aggregate data and figure data in this repository, and any figures drawn from them, are licensed under the Creative Commons Attribution 4.0 International licence (CC BY 4.0), https://creativecommons.org/licenses/by/4.0/.

Copyright (c) 2026 Siddalingaiah H S.

This covers:

- `papers/B/results.json`, `papers/B/output/` and `papers/B/figures/data/` (results, aggregate tables and figure data);
- `data/derived/` (NFHS aggregates in `data/derived/nfhs/`, ENSO and climate tables, QC-mask tables);
- `pipeline/reference/` (reviewed QC masks, state frame, NFHS district crosswalk and review, NFHS report estimates).

Attribution: cite the paper (see README.md) or this repository (CITATION.cff), and the original sources below.

**GADM exception.** GADM-derived identifiers, names and boundaries (the GADM district and state codes such as `region_id` and `district_gid`, the district and state names taken from GADM, and the district and national outlines drawn in Figure 2, including the disputed-area units from the GADM Pakistan and China files) remain subject to the GADM licence (free for academic and other non-commercial use; redistribution or commercial use needs GADM's permission, https://gadm.org/license.html), not to CC BY 4.0.

## Third-party sources

CC BY 4.0 applies to the aggregation and processing in this repository. The original data remain subject to their producers' terms, and users must also attribute them.

| Source | Producer | Terms |
|---|---|---|
| National Family Health Survey, NFHS-4 (2015–16) and NFHS-5 (2019–21): children's recode (KR) and cluster GPS (GE) datasets | International Institute for Population Sciences (IIPS), Mumbai, and ICF; distributed by the DHS Program | Available from the DHS Program on registration; **not redistributed here**. Only aggregates are published (cells with fewer than 30 children suppressed); no record-level or cluster-level data and no coordinates. Cite IIPS and ICF (2017, 2021; see README.md). |
| Oceanic Niño Index, Relative ONI, Niño-3.4 | NOAA Climate Prediction Center | US Government public data; attribute NOAA CPC |
| Dipole Mode Index | NOAA Physical Sciences Laboratory | US Government public data; attribute NOAA PSL |
| Gridded daily rainfall (0.25°) and temperature (1°) | India Meteorological Department (IMD), Pune | IMD gridded-data terms; raw grids are not redistributed, only area-weighted aggregates. Cite Pai DS et al. 2014, MAUSAM 65(1):1–18 (rainfall) and Srivastava AK et al. 2009, Atmos Sci Lett 10:249–254 (temperature) |
| Administrative boundaries | GADM 4.1 (https://gadm.org) | Used for area weights and district assignment only; not redistributed (GADM licence) |
