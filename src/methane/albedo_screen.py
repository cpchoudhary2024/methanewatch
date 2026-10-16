# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Screen candidates for surface-albedo artifacts.

Sharp SWIR albedo boundaries can produce spurious XCH4 enhancements. Around each
candidate the enhancement field is correlated with the albedo-gradient magnitude
on the detection grid, and the plume's principal axis is compared with the wind
axis. A candidate is flagged only if both hold: r above R_THRESHOLD and the axis
more than MISALIGN_DEG off the wind. The p-value on r ignores spatial
autocorrelation between cells and is optimistic.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np

R_THRESHOLD = 0.65
MISALIGN_DEG = 60.0
MIN_CELLS = 30
QA_THRESHOLD = 0.5


def albedo_gradient(albedo: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Gradient magnitude and bearing of steepest increase (degrees from north).

    Missing cells are median-filled before differencing and masked again after.
    """
    a = np.asarray(albedo, dtype=float)
    invalid = ~np.isfinite(a)
    if invalid.all():
        return np.full_like(a, np.nan), np.full_like(a, np.nan)
    filled = np.where(invalid, np.nanmedian(a), a)
    gy, gx = np.gradient(filled)
    magnitude = np.hypot(gx, gy)
    direction = (np.degrees(np.arctan2(gx, gy))) % 360.0
    magnitude[invalid] = np.nan
    direction[invalid] = np.nan
    return magnitude, direction


def pearson_with_p(a: np.ndarray, b: np.ndarray) -> tuple[float | None, float | None, int]:
    """Pearson r, its naive two-sided p-value and the number of cells used."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    n = int(ok.sum())
    if n < 4:
        return None, None, n
    x, y = a[ok], b[ok]
    if np.std(x) == 0 or np.std(y) == 0:
        return None, None, n
    from scipy import stats

    r, p = stats.pearsonr(x, y)
    return float(r), float(p), n


def principal_axis_deg(field: np.ndarray, lat_deg: float) -> float | None:
    """Bearing in [0, 180) of the major axis of the positive part of `field`.

    Uses intensity-weighted second moments with x scaled by cos(lat).
    """
    f = np.asarray(field, dtype=float)
    w = np.where(np.isfinite(f), np.clip(f, 0, None), 0.0)
    total = w.sum()
    if total <= 0 or np.count_nonzero(w) < 3:
        return None
    ny, nx = w.shape
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    xx = xx * math.cos(math.radians(lat_deg))
    cx = float((w * xx).sum() / total)
    cy = float((w * yy).sum() / total)
    dx, dy = xx - cx, yy - cy
    sxx = float((w * dx * dx).sum() / total)
    syy = float((w * dy * dy).sum() / total)
    sxy = float((w * dx * dy).sum() / total)
    cov = np.array([[sxx, sxy], [sxy, syy]])
    if not np.all(np.isfinite(cov)):
        return None
    vals, vecs = np.linalg.eigh(cov)
    vx, vy = vecs[:, int(np.argmax(vals))]
    return float(math.degrees(math.atan2(vx, vy)) % 180.0)


def axis_misalignment_deg(axis_deg: float | None, wind_from_deg: float | None) -> float | None:
    """Angle in [0, 90] between a plume axis and the wind axis (direction ignored)."""
    if axis_deg is None or wind_from_deg is None or not math.isfinite(wind_from_deg):
        return None
    diff = abs((axis_deg % 180.0) - (wind_from_deg % 180.0)) % 180.0
    return float(min(diff, 180.0 - diff))


def _candidate_patch(
    granule_path: str | Path,
    lon_c: float,
    lat_c: float,
    *,
    half_width_deg: float = 0.5,
    res_deg: float = 0.05,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Enhancement and albedo around a candidate on a common grid, or None if too sparse."""
    from netCDF4 import Dataset

    try:
        with Dataset(granule_path, "r") as ds:
            prod = ds.groups["PRODUCT"]
            lon = np.ma.filled(prod.variables["longitude"][:], np.nan).squeeze()
            lat = np.ma.filled(prod.variables["latitude"][:], np.nan).squeeze()
            xch4 = np.ma.filled(
                prod.variables["methane_mixing_ratio_bias_corrected"][:], np.nan
            ).squeeze()
            qa = np.ma.filled(prod.variables["qa_value"][:], np.nan).squeeze()
            albedo = np.ma.filled(
                prod.groups["SUPPORT_DATA"]
                .groups["DETAILED_RESULTS"]
                .variables["surface_albedo_SWIR"][:],
                np.nan,
            ).squeeze()
    except Exception:
        return None

    lon_e = np.arange(lon_c - half_width_deg, lon_c + half_width_deg + res_deg, res_deg)
    lat_e = np.arange(lat_c - half_width_deg, lat_c + half_width_deg + res_deg, res_deg)
    keep = (
        np.isfinite(qa)
        & (qa >= QA_THRESHOLD)
        & np.isfinite(xch4)
        & np.isfinite(albedo)
        & (albedo > 0.01)
        & (albedo < 1.5)
        & (lon >= lon_e[0])
        & (lon <= lon_e[-1])
        & (lat >= lat_e[0])
        & (lat <= lat_e[-1])
    )
    if keep.sum() < 10:
        return None

    def _grid(vals):
        s, _, _ = np.histogram2d(lon[keep], lat[keep], [lon_e, lat_e], weights=vals[keep])
        c, _, _ = np.histogram2d(lon[keep], lat[keep], [lon_e, lat_e])
        with np.errstate(invalid="ignore"):
            return np.where(c > 0, s / c, np.nan).T

    ch4_grid = _grid(xch4)
    albedo_grid = _grid(albedo)
    enhancement = ch4_grid - np.nanmedian(ch4_grid)
    return enhancement, albedo_grid


def analyze_candidate(
    granule_path: str | Path,
    lon_c: float,
    lat_c: float,
    wind_from_deg: float | None,
    *,
    half_width_deg: float = 0.5,
    res_deg: float = 0.05,
) -> dict[str, Any]:
    patch = _candidate_patch(
        granule_path, lon_c, lat_c, half_width_deg=half_width_deg, res_deg=res_deg
    )
    if patch is None:
        return {"available": False, "note": "no QA-passing pixels around the candidate"}

    enhancement, albedo_grid = patch
    grad_mag, _ = albedo_gradient(albedo_grid)
    r, p, n_cells = pearson_with_p(enhancement, grad_mag)
    if r is None or n_cells < MIN_CELLS:
        return {
            "available": False,
            "n_cells": n_cells,
            "note": f"{n_cells} valid cells, fewer than the {MIN_CELLS} required",
        }

    axis = principal_axis_deg(enhancement, lat_c)
    misalign = axis_misalignment_deg(axis, wind_from_deg)
    albedo_flag = r > R_THRESHOLD
    wind_flag = misalign is not None and misalign > MISALIGN_DEG
    return {
        "available": True,
        "n_cells": n_cells,
        "albedo_gradient_r": round(r, 3),
        "albedo_gradient_p_naive": None if p is None else round(p, 4),
        "plume_axis_deg": None if axis is None else round(axis, 1),
        "wind_from_deg": wind_from_deg,
        "axis_misalignment_deg": None if misalign is None else round(misalign, 1),
        "albedo_flag": bool(albedo_flag),
        "wind_flag": bool(wind_flag),
        "flagged": bool(albedo_flag and wind_flag),
    }
