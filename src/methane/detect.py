# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Single-overpass enhancement detection.

QA-passing pixels are mean-binned to a 0.05 degree grid (about the native pixel
size), a 0.6 degree median filter supplies the regional background, and cells
whose residual exceeds max(k * 1.4826 * MAD, median retrieval precision) are
flagged. Connected clusters of at least `min_pixels` cells are candidates.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy import ndimage

from .read import Scene


@dataclass
class Candidate:
    lon: float
    lat: float
    delta_xch4_mean: float  # ppb
    delta_xch4_peak: float  # ppb
    n_pixels: int
    area_km2: float
    confidence_sigma: float  # cluster-mean enhancement / robust residual sigma

    def as_dict(self) -> dict[str, Any]:
        return {
            "lon": round(self.lon, 4),
            "lat": round(self.lat, 4),
            "delta_xch4_mean_ppb": round(self.delta_xch4_mean, 1),
            "delta_xch4_peak_ppb": round(self.delta_xch4_peak, 1),
            "n_pixels": self.n_pixels,
            "area_km2": round(self.area_km2, 1),
            "confidence_sigma": round(self.confidence_sigma, 1),
        }


@dataclass
class DetectionResult:
    candidates: list[Candidate]
    background_ppb: float
    residual_sigma_ppb: float
    noise_floor_ppb: float
    threshold_ppb: float
    grid_res_deg: float
    scene_source: str
    diagnostics: dict[str, Any] = field(default_factory=dict)

    def summary(self) -> str:
        return (
            f"{len(self.candidates)} candidate(s) in {self.scene_source}\n"
            f"  background={self.background_ppb:.0f} ppb, "
            f"robust sigma={self.residual_sigma_ppb:.1f} ppb, "
            f"retrieval noise floor={self.noise_floor_ppb:.1f} ppb\n"
            f"  detection threshold={self.threshold_ppb:.1f} ppb enhancement"
        )


def _bin_to_grid(scene: Scene, res_deg: float):
    """Mean-bin pixels onto a regular grid; returns (grid[lat, lon], lon centres, lat centres)."""
    lon_edges = np.arange(scene.lon.min(), scene.lon.max() + res_deg, res_deg)
    lat_edges = np.arange(scene.lat.min(), scene.lat.max() + res_deg, res_deg)
    if lon_edges.size < 2 or lat_edges.size < 2:
        raise ValueError("Region too small for the chosen grid resolution.")

    sum_, _, _ = np.histogram2d(
        scene.lon, scene.lat, bins=[lon_edges, lat_edges], weights=scene.xch4
    )
    cnt, _, _ = np.histogram2d(scene.lon, scene.lat, bins=[lon_edges, lat_edges])
    with np.errstate(invalid="ignore", divide="ignore"):
        grid = np.where(cnt > 0, sum_ / cnt, np.nan)
    lon_c = 0.5 * (lon_edges[:-1] + lon_edges[1:])
    lat_c = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    return grid.T, lon_c, lat_c


def cell_area_km2(res_deg: float, lat: float) -> float:
    km_lat = 110.574 * res_deg
    km_lon = 111.320 * np.cos(np.radians(lat)) * res_deg
    return km_lat * km_lon


def detect(
    scene: Scene,
    *,
    grid_res_deg: float = 0.05,
    background_window_deg: float = 0.6,
    k_sigma: float = 3.0,
    min_pixels: int = 3,
) -> DetectionResult:
    if scene.n < 50:
        return DetectionResult(
            candidates=[],
            background_ppb=float("nan"),
            residual_sigma_ppb=float("nan"),
            noise_floor_ppb=float("nan"),
            threshold_ppb=float("nan"),
            grid_res_deg=grid_res_deg,
            scene_source=scene.source,
            diagnostics={"reason": "insufficient QA-passing pixels", "n": scene.n},
        )

    grid, lon_c, lat_c = _bin_to_grid(scene, grid_res_deg)
    valid = np.isfinite(grid)

    fill = np.where(valid, grid, np.nanmedian(grid))
    win = max(3, int(round(background_window_deg / grid_res_deg)))
    background_field = ndimage.median_filter(fill, size=win, mode="nearest")

    residual = np.where(valid, grid - background_field, np.nan)
    r = residual[valid]
    med = float(np.nanmedian(r))
    mad = float(np.nanmedian(np.abs(r - med)))
    sigma = 1.4826 * mad if mad > 0 else float(np.nanstd(r))

    noise_floor = (
        float(np.nanmedian(scene.precision)) if np.isfinite(scene.precision).any() else 0.0
    )

    threshold = max(k_sigma * sigma, noise_floor)

    flagged = valid & (residual > threshold)
    labels, n = ndimage.label(flagged)

    candidates: list[Candidate] = []
    for lab in range(1, n + 1):
        mask = labels == lab
        npix = int(mask.sum())
        if npix < min_pixels:
            continue
        dvals = residual[mask]
        iy, ix = np.where(mask)
        weights = np.clip(dvals - threshold, 1e-6, None)
        clon = float(np.average(lon_c[ix], weights=weights))
        clat = float(np.average(lat_c[iy], weights=weights))
        conf = float(dvals.mean() / sigma) if sigma > 0 else float("nan")
        area = npix * cell_area_km2(grid_res_deg, clat)
        candidates.append(
            Candidate(
                lon=clon,
                lat=clat,
                delta_xch4_mean=float(dvals.mean()),
                delta_xch4_peak=float(dvals.max()),
                n_pixels=npix,
                area_km2=area,
                confidence_sigma=conf,
            )
        )

    candidates.sort(key=lambda c: c.confidence_sigma, reverse=True)
    return DetectionResult(
        candidates=candidates,
        background_ppb=float(np.nanmedian(grid)),
        residual_sigma_ppb=sigma,
        noise_floor_ppb=noise_floor,
        threshold_ppb=threshold,
        grid_res_deg=grid_res_deg,
        scene_source=scene.source,
        diagnostics={
            "n_pixels_qa": scene.n,
            "qa_pass_fraction": round(scene.qa_pass_fraction, 3),
            "grid_cells_valid": int(valid.sum()),
            "k_sigma": k_sigma,
            "min_pixels": min_pixels,
        },
    )
