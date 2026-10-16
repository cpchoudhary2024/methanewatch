# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Reader for Sentinel-5P TROPOMI L2 CH4 granules.

Uses the bias-corrected XCH4 field and the qa_value >= 0.5 filter recommended in
the product user manual. Pixels failing QA are dropped, never filled.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

QA_THRESHOLD = 0.5


@dataclass
class Scene:
    lon: np.ndarray
    lat: np.ndarray
    xch4: np.ndarray  # ppb
    precision: np.ndarray  # ppb, 1-sigma retrieval precision
    qa: np.ndarray
    p_surf: np.ndarray  # Pa
    source: str
    n_total: int  # pixels in the box before QA filtering
    time_utc: str  # median time of the scanlines that cross the box

    @property
    def n(self) -> int:
        return int(self.xch4.size)

    @property
    def qa_pass_fraction(self) -> float:
        return self.n / self.n_total if self.n_total else 0.0

    @property
    def p_surf_median(self) -> float:
        p = self.p_surf[np.isfinite(self.p_surf)]
        return float(np.median(p)) if p.size else float("nan")


def _var(group, path: str):
    *groups, name = path.split("/")
    for g in groups:
        if g not in group.groups:
            return None
        group = group.groups[g]
    return group.variables.get(name)


def _filled(var, shape) -> np.ndarray:
    if var is None:
        return np.full(shape, np.nan)
    return np.ma.filled(var[:].astype(float), np.nan).squeeze()


def _box_time(prod, in_box: np.ndarray) -> str:
    """Median timestamp of the scanlines that have at least one pixel in the box."""
    try:
        times = np.asarray(prod.variables["time_utc"][:]).squeeze()
    except KeyError:
        return "unknown"
    rows = np.where(in_box.any(axis=1))[0]
    if rows.size == 0 or times.ndim != 1:
        return str(times.flat[0]) if times.size else "unknown"
    return str(times[rows[rows.size // 2]])


def read_scene(
    path: str | Path,
    bbox: tuple[float, float, float, float],
    *,
    qa_threshold: float = QA_THRESHOLD,
) -> Scene:
    lon_min, lat_min, lon_max, lat_max = bbox
    path = Path(path)
    with Dataset(path, "r") as ds:
        prod = ds.groups["PRODUCT"]
        lon = _filled(prod.variables["longitude"], None)
        lat = _filled(prod.variables["latitude"], None)
        xch4 = _filled(prod.variables["methane_mixing_ratio_bias_corrected"], None)
        qa = _filled(prod.variables["qa_value"], None)
        precision = _filled(_var(prod, "methane_mixing_ratio_precision"), xch4.shape)
        p_surf = _filled(_var(prod, "SUPPORT_DATA/INPUT_DATA/surface_pressure"), xch4.shape)

        in_box = (
            (lon >= lon_min)
            & (lon <= lon_max)
            & (lat >= lat_min)
            & (lat <= lat_max)
            & np.isfinite(xch4)
        )
        time_utc = _box_time(prod, in_box) if in_box.ndim == 2 else "unknown"

    keep = in_box & (qa >= qa_threshold)
    return Scene(
        lon=lon[keep],
        lat=lat[keep],
        xch4=xch4[keep],
        precision=precision[keep],
        qa=qa[keep],
        p_surf=p_surf[keep],
        source=path.name,
        n_total=int(in_box.sum()),
        time_utc=time_utc,
    )
