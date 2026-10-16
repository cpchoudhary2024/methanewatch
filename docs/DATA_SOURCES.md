# Data sources

| Dataset | Use | Access |
|---|---|---|
| Sentinel-5P TROPOMI L2 `CH4` (OFFL, processors 02.06 in 2024 and 02.09 in 2026), bias-corrected XCH4, `qa_value >= 0.5` | Detection, mass balance | [Copernicus Data Space](https://dataspace.copernicus.eu/). Catalogue search is open; download needs a free account. |
| ERA5 hourly 10 m and 100 m wind | Point-source rates (10 m), basin transport wind (100 m), detection floor | [Open-Meteo archive API](https://open-meteo.com/en/docs/historical-weather-api), no key. Responses cached in `data/cache/wind/`. |
| EPA GHGRP Subpart W, CH4, TX and NM, 2011 to 2023 | Reported emissions | [Envirofacts](https://data.epa.gov/efservice), table `V_GHG_EMITTER_SUBPART`. |
| Carbon Mapper annotated CH4 plumes, query through 2026-09-01 | Source locations and rates for comparison | [Carbon Mapper API](https://api.carbonmapper.org/api/v1/catalog/plumes/annotated), no key. |
| Flaring Monitor detailed observations, 2019 to 2021 | Persistent flares | [flaringmonitor/viirs-flare-data](https://github.com/flaringmonitor/viirs-flare-data), built on VIIRS Nightfire (Elvidge et al., 2016). |
| OpenStreetMap nodes tagged `man_made=petroleum_well` or `oil_well` | Infrastructure-conditioned null and QA-pixel weights | Overpass API. Retrieved 2026-10-14, 52,546 nodes in the study box. |

## Basin attribution of GHGRP reports

Onshore production and gathering and boosting facilities report per AAPG basin, and the
basin code usually appears in the facility name. A report is counted as Permian if its name
contains "Permian" or basin code 430, regardless of its nominal coordinate. Reports whose
name carries another basin code are excluded even when their coordinate falls inside the
box. Other reports (processing plants, compressor stations) are counted if their coordinate
is inside the box.

## Optional local inputs

A Texas RRC or New Mexico OCD well export saved as `data/wellheads/wellheads_permian.csv`
(columns `lon`/`longitude` and `lat`/`latitude`) replaces the OpenStreetMap layer.
