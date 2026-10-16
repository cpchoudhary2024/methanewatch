# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Sensitivity analysis conditioning spatial controls on observed pixels."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

from .flaring import load_flares
from .read import read_scene
from .regions import REGIONS

EARTH_RADIUS_KM = 6371.0


def _xyz(lon, lat):
    lon = np.radians(np.asarray(lon, dtype=float))
    lat = np.radians(np.asarray(lat, dtype=float))
    return np.column_stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def within_km(tree: cKDTree, lon, lat, distance_km: float) -> np.ndarray:
    """Exact great-circle threshold via monotone 3-D chord distance."""
    chord = 2 * np.sin(distance_km / (2 * EARTH_RADIUS_KM))
    distance, _ = tree.query(_xyz(lon, lat), workers=-1)
    return distance <= chord


def well_density_grid(wells: list[dict], bbox, *, bandwidth_km=8.0, res=0.02):
    """Reproduce the original well-density surface for an isolated mask check."""
    x0, y0, x1, y1 = bbox
    lon_edges = np.arange(x0, x1 + res, res)
    lat_edges = np.arange(y0, y1 + res, res)
    hist, _, _ = np.histogram2d(
        [w["lon"] for w in wells], [w["lat"] for w in wells], [lon_edges, lat_edges]
    )
    km_per_deg = 111.0 * np.cos(np.radians((y0 + y1) / 2))
    sigma_cells = max(bandwidth_km / (km_per_deg * res), 0.5)
    density = np.clip(gaussian_filter(hist, sigma=sigma_cells, mode="constant"), 0, None)
    return lon_edges, lat_edges, density


def pixel_weights(lon, lat, surface) -> np.ndarray:
    lon_edges, lat_edges, density = surface
    ix = np.clip(np.searchsorted(lon_edges, lon, side="right") - 1, 0, density.shape[0] - 1)
    iy = np.clip(np.searchsorted(lat_edges, lat, side="right") - 1, 0, density.shape[1] - 1)
    weights = density[ix, iy]
    if weights.sum() <= 0:
        raise ValueError("No mapped well density in observed pixels")
    return weights / weights.sum()


def paired_enrichment_ci(rows: list[dict], key: str, *, n_boot=2000, seed=0):
    """Resample candidate days, carrying each day's observed and null fractions together."""
    days = sorted({row["day"] for row in rows})
    groups = [[row for row in rows if row["day"] == day] for day in days]
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(n_boot):
        selected = [row for i in rng.integers(0, len(groups), len(groups)) for row in groups[i]]
        numerator = np.mean([row["near_flare"] for row in selected])
        denominator = np.mean([row[key] for row in selected])
        if denominator:
            values.append(numerator / denominator)
    return [round(float(v), 2) for v in np.percentile(values, [2.5, 97.5])]


def analyze(out: Path, *, cache: Path = Path("data/cache"), tol_km: float = 5.0) -> dict:
    region = REGIONS[out.name.split("-")[0]]
    detections = json.loads((out / "detections.json").read_text())
    wells = json.loads(Path(f"data/wellheads/osm_wells_{region.key}.json").read_text())["points"]
    sites = load_flares("data/flaring/flares_permian.csv", region)
    sites = sites[sites.rh > sites.rh.quantile(0.5)]
    flare_tree = cKDTree(_xyz(sites.klon, sites.klat))
    surface = well_density_grid(wells, region.bbox)
    scenes = {}
    rows = []
    for detection in detections:
        name = detection["granule"]
        if name not in scenes:
            scene = read_scene(cache / name, region.bbox)
            near = within_km(flare_tree, scene.lon, scene.lat, tol_km).astype(float)
            weights = pixel_weights(scene.lon, scene.lat, surface)
            scenes[name] = {
                "n_pixels": scene.n,
                "uniform_fraction": float(np.mean(near)),
                "wells_fraction": float(np.dot(near, weights)),
            }
        scene_result = scenes[name]
        candidate = detection["detection"]
        near_candidate = bool(
            within_km(flare_tree, [candidate["lon"]], [candidate["lat"]], tol_km)[0]
        )
        rows.append(
            {
                "day": detection["overpass_utc"][:10],
                "granule": name,
                "near_flare": near_candidate,
                **scene_result,
            }
        )
    observed = float(np.mean([row["near_flare"] for row in rows]))
    result = {
        "region": region.key,
        "year": json.loads((out / "summary.json").read_text())["window"][0][:4],
        "n_candidates": len(rows),
        "n_candidate_days": len({row["day"] for row in rows}),
        "n_candidate_scenes": len(scenes),
        "n_wells": len(wells),
        "n_flares": len(sites),
        "tolerance_km": tol_km,
        "candidate_fraction": round(observed, 4),
        "method": (
            "QA-passing pixel centres in each candidate scene; well-density weighting "
            "uses the original 8 km kernel surface; day-paired bootstrap"
        ),
        "rows": rows,
    }
    for key, label in [("uniform_fraction", "uniform"), ("wells_fraction", "wells")]:
        null = float(np.mean([row[key] for row in rows]))
        result[label] = {
            "null_fraction": round(null, 4),
            "enrichment": round(observed / null, 3),
            "ci_days": paired_enrichment_ci(rows, key),
        }
    return result
