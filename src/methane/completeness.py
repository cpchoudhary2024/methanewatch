# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Where and how often the detector could have seen a source of a given size.

For each overpass the smallest detectable source follows from the detector's own
threshold T (ppb) and minimum cluster size: a cluster of `min_pixels` cells of
area A_cell with mean enhancement T gives, through Q = U * IME / L,

    Q_min = U * c(p_s) * T * sqrt(min_pixels * A_cell)

with U the 10 m ERA5 wind at the overpass. A grid cell observed on an overpass
counts as covered for a reference source Q_ref when Q_min <= Q_ref. Completeness
is the fraction of the window's overpasses on which a cell was both observed and
covered. Cells never observed are null.
"""

from __future__ import annotations

import math

import numpy as np

from . import wind
from .detect import cell_area_km2, detect
from .flux import column_excess_kg_per_m2
from .read import read_scene
from .regions import Region


def scene_floor_t_hr(
    threshold_ppb: float, wind_ms: float, p_surf: float, cell_km2: float, min_pixels: int
) -> float:
    area_m2 = min_pixels * cell_km2 * 1e6
    q_kg_s = wind_ms * column_excess_kg_per_m2(threshold_ppb, p_surf) * math.sqrt(area_m2)
    return q_kg_s * 3.6


def compute_completeness(
    paths: list[str],
    region: Region,
    *,
    reference_t_hr: float = 25.0,
    grid_res_deg: float = 0.05,
    min_pixels: int = 3,
) -> dict:
    x0, y0, x1, y1 = region.bbox
    lon_e = np.arange(x0, x1 + grid_res_deg / 2, grid_res_deg)
    lat_e = np.arange(y0, y1 + grid_res_deg / 2, grid_res_deg)
    shape = (len(lat_e) - 1, len(lon_e) - 1)
    n_total = 0
    n_seen = np.zeros(shape)
    n_covered = np.zeros(shape)
    floors = []

    cell_km2 = cell_area_km2(grid_res_deg, (y0 + y1) / 2)
    for path in sorted(paths):
        try:
            s = read_scene(path, region.bbox)
        except (OSError, KeyError):
            continue
        n_total += 1
        res = detect(s, grid_res_deg=grid_res_deg, min_pixels=min_pixels)
        if not math.isfinite(res.threshold_ppb) or s.time_utc == "unknown":
            continue
        u10, _ = wind.area_mean(region.bbox, s.time_utc, level=10)
        floor = scene_floor_t_hr(res.threshold_ppb, u10, s.p_surf_median, cell_km2, min_pixels)
        floors.append(floor)
        seen = np.histogram2d(s.lon, s.lat, [lon_e, lat_e])[0].T > 0
        n_seen += seen
        if floor <= reference_t_hr:
            n_covered += seen

    with np.errstate(invalid="ignore", divide="ignore"):
        pct = np.where(n_seen > 0, 100.0 * n_covered / max(n_total, 1), np.nan)

    return {
        "grid_res_deg": grid_res_deg,
        "lon_edges": [round(float(v), 3) for v in lon_e],
        "lat_edges": [round(float(v), 3) for v in lat_e],
        "reference_t_hr": reference_t_hr,
        "n_overpasses": n_total,
        "n_scenes_with_floor": len(floors),
        "floor_t_hr_median": round(float(np.median(floors)), 1) if floors else None,
        "floor_t_hr_p10_p90": (
            [round(float(v), 1) for v in np.percentile(floors, [10, 90])] if floors else None
        ),
        "completeness_pct": [
            [None if not np.isfinite(v) else round(float(v), 1) for v in row] for row in pct
        ],
        "n_seen": n_seen.astype(int).tolist(),
    }
