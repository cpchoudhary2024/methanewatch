# Design notes

What each part of the pipeline does, why it is built that way, and which earlier mistakes it
guards against. Read this before changing the methods.

## Layout

| Path | Contents |
|---|---|
| `src/methane/read.py` | Reads a TROPOMI L2 CH4 granule: bias-corrected XCH4, QA >= 0.5, surface pressure, overpass time |
| `src/methane/detect.py` | Single-scene detector (0.05 degree grid, 0.6 degree median background, 3 sigma MAD, >= 3 cells) |
| `src/methane/flux.py` | Integrated mass enhancement rates and the column conversion |
| `src/methane/wind.py` | ERA5 hourly wind from Open-Meteo, cached per point and day |
| `src/methane/completeness.py` | Detection floor implied by each scene's threshold |
| `src/methane/spatial_stats.py` | Uniform and well-density nulls, day bootstrap, Moran's I, spatial-error model |
| `src/methane/opportunity.py` | Sensitivity analysis conditioning both nulls on QA-passing pixels in candidate scenes |
| `src/methane/validate.py` | Carbon Mapper download, proximity tests, rate comparison |
| `src/methane/flaring.py` | Flare proximity tests and the grid-scale correlation |
| `src/methane/gap.py` | Box mass balance and GHGRP Subpart W totals |
| `src/methane/wells.py` | Well locations (local export or OpenStreetMap) |
| `src/methane/report.py` | Builds `outputs/permian/index.html` from every `outputs/<basin>/` folder |
| `src/methane/cli.py` | `methane` command: catalog, run, validate, gap, flaring, report, provenance |
| `outputs/<basin>[-year]/` | JSON results of each run; `outputs/permian/` is the 2026 main run and the website root |
| `paper/manuscript.md` | Manuscript; `make figures` and `make paper` rebuild its figures and HTML |

Run order for one basin and season: `run`, `validate`, `flaring`, then `gap` (Permian only),
`report`, `provenance`. `make all` does this for the default window.

## Decisions that the results depend on

1. **Two nulls, always side by side.** Every proximity statistic is reported against points spread
   uniformly in the box and against points drawn from a kernel density of mapped wells. The uniform
   null alone is misleading in a production basin, and that contrast is the paper's main result. Do
   not report a proximity statistic against one null only.
2. **Sign-inverted control.** `cli._process_granule` runs the detector on `-XCH4` for
   every scene and stores the clusters in `summary.json` (`overpasses[].inverted`). Keep this when
   changing the detector. The count diagnoses negative excursions; it is not a calibrated
   false-positive rate.
3. **Days, not candidates, are the bootstrap unit.** Candidates from one overpass are not
   independent (`spatial_stats.enrichment_ci`).
4. **The headline gap uses the published inversion.** The ratio of 7.6 is Zhang et al. (2020)
   against GHGRP 2023. The box mass balance in `gap.py` is a consistency check and must not replace
   it: its value moves from 0.33 to 4.8 Tg/yr with the background choice.
5. **GHGRP attribution is by basin code, not coordinate.** Production and gathering report per
   AAPG basin, so names containing 430 or "Permian" are kept wherever their coordinate lies.
6. **Fixed study windows.** Results are for June to August of a named year. There is no rolling
   or daily update, and none of the outputs describe a "live" state.
7. **Observation opportunity is a paired sensitivity check.** `opportunity.py` retains the same
   positive candidates and flare definition, evaluates both nulls on QA-passing pixels from each
   candidate's own granule, and resamples candidate days with their control fractions. It does not
   overwrite the box-wide results.

## Mistakes already fixed, which must not return

| Mistake | Effect | Guard |
|---|---|---|
| Mass balance background from each scene's 10th percentile, enhancements clipped at zero | Measured pixel spread, gave 4.2 Tg/yr | `gap.scene_balance` uses a terrain-corrected ring |
| Catalogue query capped at 100 results | Silently dropped the last two weeks of August | `search_granules(limit=1000)` |
| Wind taken at the granule start time | Up to an hour off the Permian overpass | `read._box_time` uses scanlines over the box |
| First 5,000 OpenStreetMap wells by node ID | Spatially biased null | `wells.fetch_wellheads_osm` takes all nodes |
| Carbon Mapper plumes summed without merging repeat visits | One source counted many times | `validate.carbon_mapper_sources` |
| Toroidal null computed on a different cell mask from the observed rho | Invalid permutation p | Shared `observed` mask in `flaring.grid_correlation` |
| New CDSE token per granule | Account locked by the identity server | `download.get_access_token` caches the token |
| Detection floor from a hand-tuned albedo and sun-angle formula | Claimed 2 t/h leaks were visible half the time | `completeness.scene_floor_t_hr` uses the detector's own threshold |

## Tests

`pytest` runs without network access in about two seconds. `tests/test_core.py` covers the mass
balance geometry, terrain correction, window filtering and both nulls on synthetic data. Add a test
with any change to a statistic that appears in the paper.

## Data that changes over time

The Carbon Mapper catalogue, OpenStreetMap wells and the GHGRP database are live services. Runs
record the Carbon Mapper end date (`validation.json`), the well retrieval date (cache files in
`data/wellheads/`) and the reporting year (`gap.json`). Re-running later can shift numbers
slightly. The granule SHA-256 hashes in `provenance.json` fix the satellite input exactly.

## Next validation steps

1. Replace the OpenStreetMap wells with dated Texas RRC and New Mexico OCD records, and match flare
   activity to the satellite windows. Report how far both controls move.
2. Find same-day high-resolution plume observations or documented facility events for candidate
   scenes. The sign-inverted count shows ambiguity but cannot say which candidates are real.
3. Weight the QA-matched references by each scene's detection floor rather than treating every
   passing pixel as equally observable.
4. Test the fixed method on further seasons and on basins with at least ten candidates.
5. Compare reported and observed totals for matching years once a newer Subpart W year is published.
