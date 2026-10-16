# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Spatial comparison of TROPOMI candidates against Carbon Mapper CH4 plumes.

Carbon Mapper resolves ~30 m point sources and rarely observes on the same day as
a TROPOMI overpass, so candidates are compared with the full-mission set of
plume locations rather than matched event by event. Each statistic is computed
for two control populations: points drawn uniformly in the box, and points drawn
from the density of mapped wellheads. Only the second controls for the fact that
candidates and sources both sit where the industry is.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .spatial_stats import clustered_null_points, enrichment_ci, uniform_null_points

EARTH_RADIUS_KM = 6371.0
CM_API = "https://api.carbonmapper.org/api/v1/catalog/plumes/annotated"


@dataclass
class ReferenceDetection:
    lon: float
    lat: float
    source_id: str = ""
    emission_kg_hr: float | None = None
    date: str = ""


def nearest_km(qlon, qlat, ref_lon, ref_lat) -> np.ndarray:
    """Great-circle distance from each query point to its nearest reference point."""
    qlon = np.radians(np.atleast_1d(qlon))[:, None]
    qlat = np.radians(np.atleast_1d(qlat))[:, None]
    rlon = np.radians(np.asarray(ref_lon))[None, :]
    rlat = np.radians(np.asarray(ref_lat))[None, :]
    a = (
        np.sin((rlat - qlat) / 2) ** 2
        + np.cos(qlat) * np.cos(rlat) * np.sin((rlon - qlon) / 2) ** 2
    )
    return (2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))).min(axis=1)


def fetch_carbon_mapper(
    bbox: tuple[float, float, float, float],
    start: str | None = None,
    end: str | None = None,
    *,
    page: int = 500,
    session=None,
) -> list[ReferenceDetection]:
    import requests

    sess = session or requests.Session()
    x0, y0, x1, y1 = bbox
    out: list[ReferenceDetection] = []
    offset = 0
    while True:
        params = [("bbox", v) for v in (x0, y0, x1, y1)]
        params += [("plume_gas", "CH4"), ("limit", page), ("offset", offset)]
        if start and end:
            params.append(("datetime", f"{start}T00:00:00Z/{end}T00:00:00Z"))
        r = sess.get(CM_API, params=params, timeout=90)
        r.raise_for_status()
        j = r.json()
        items = j.get("items", [])
        for it in items:
            g = it.get("geometry_json") or {}
            if g.get("type") != "Point":
                continue
            lon, lat = g["coordinates"][:2]
            out.append(
                ReferenceDetection(
                    lon=float(lon),
                    lat=float(lat),
                    source_id=str(it.get("plume_id", "")),
                    emission_kg_hr=it.get("emission_auto"),
                    date=str(it.get("scene_timestamp", ""))[:10],
                )
            )
        offset += page
        if not items or offset >= j.get("bbox_count", len(items)):
            break
    return out


def _null_block(dist_det, days, null_lon, null_lat, method, ref_lon, ref_lat, tol_km):
    dist_null = nearest_km(null_lon, null_lat, ref_lon, ref_lat)
    frac_det = float(np.mean(dist_det <= tol_km))
    frac_null = float(np.mean(dist_null <= tol_km))
    return {
        "method": method,
        "median_km_null": round(float(np.median(dist_null)), 2),
        "frac_within_tol_null": round(frac_null, 3),
        "enrichment": round(frac_det / frac_null, 2) if frac_null else None,
        "enrichment_ci_days": enrichment_ci(dist_det <= tol_km, days, frac_null),
    }


def build_validation(
    detections: list[dict],
    references: list[ReferenceDetection],
    *,
    bbox: tuple[float, float, float, float],
    infra_lon,
    infra_lat,
    tol_km: float = 10.0,
    n_null: int = 2000,
    seed: int = 0,
    references_in_window: int | None = None,
    reference_end: str | None = None,
) -> dict:
    ref_lon = np.array([r.lon for r in references])
    ref_lat = np.array([r.lat for r in references])
    det_lon = np.array([d["detection"]["lon"] for d in detections])
    det_lat = np.array([d["detection"]["lat"] for d in detections])
    dist = nearest_km(det_lon, det_lat, ref_lon, ref_lat) if len(det_lon) else np.array([])

    rows = [
        {
            "lon": d["detection"]["lon"],
            "lat": d["detection"]["lat"],
            "day": d["overpass_utc"][:10],
            "sigma": d["detection"].get("confidence_sigma"),
            "q_kg_hr": round(d["emission"]["q_kg_hr"]),
            "nearest_cm_km": round(float(km), 1),
            "within_tol": bool(km <= tol_km),
        }
        for d, km in zip(detections, dist)
    ]

    days = [d["overpass_utc"][:10] for d in detections]
    ulon, ulat, umethod = uniform_null_points(n_null, bbox, seed=seed)
    clon, clat, cmethod = clustered_null_points(n_null, bbox, infra_lon, infra_lat, seed=seed)
    by_bandwidth = {}
    for h in (4.0, 16.0):
        blon, blat, _ = clustered_null_points(
            n_null, bbox, infra_lon, infra_lat, bandwidth_km=h, seed=seed
        )
        frac = float(np.mean(nearest_km(blon, blat, ref_lon, ref_lat) <= tol_km))
        by_bandwidth[f"{h:.0f}"] = round(float(np.mean(dist <= tol_km)) / frac, 2) if frac else None
    em = [r.emission_kg_hr for r in references if r.emission_kg_hr]
    return {
        "reference": "Carbon Mapper annotated CH4 plumes, full mission",
        "reference_end": reference_end,
        "n_detections": len(detections),
        "n_cm_plumes": len(references),
        "n_cm_in_window": references_in_window,
        "tolerance_km": tol_km,
        "median_km_detections": round(float(np.median(dist)), 2) if dist.size else None,
        "frac_within_tol_detections": (
            round(float(np.mean(dist <= tol_km)), 3) if dist.size else None
        ),
        "n_within_tol": int(np.sum(dist <= tol_km)),
        "uniform_null": _null_block(dist, days, ulon, ulat, umethod, ref_lon, ref_lat, tol_km),
        "clustered_null": _null_block(dist, days, clon, clat, cmethod, ref_lon, ref_lat, tol_km),
        "clustered_enrichment_by_bandwidth_km": by_bandwidth,
        "detections": rows,
        "farthest": max(rows, key=lambda r: r["nearest_cm_km"]) if rows else None,
        "cm_emission_kg_hr": {
            "min": round(min(em)) if em else None,
            "median": round(float(np.median(em))) if em else None,
            "max": round(max(em)) if em else None,
        },
    }


def carbon_mapper_sources(cm_points: list[dict], cell_deg: float = 0.002) -> list[dict]:
    """Collapse repeat plume detections into sources on a ~200 m grid.

    Each source carries the mean of its plume rates and the number of plumes, so
    a site imaged many times is counted once.
    """
    sources: dict[tuple[int, int], dict] = {}
    for p in cm_points:
        em = p.get("emission_kg_hr")
        if p.get("lon") is None or p.get("lat") is None:
            continue
        key = (round(p["lon"] / cell_deg), round(p["lat"] / cell_deg))
        s = sources.setdefault(key, {"lon": 0.0, "lat": 0.0, "n": 0, "rates": []})
        s["lon"] += p["lon"]
        s["lat"] += p["lat"]
        s["n"] += 1
        if em is not None and math.isfinite(em) and em > 0:
            s["rates"].append(em)
    return [
        {
            "lon": s["lon"] / s["n"],
            "lat": s["lat"] / s["n"],
            "n_plumes": s["n"],
            "mean_kg_hr": sum(s["rates"]) / len(s["rates"]) if s["rates"] else None,
        }
        for s in sources.values()
    ]


def carbon_mapper_cross_check(
    detections: list[dict], cm_points: list[dict], *, radius_km: float = 15.0
) -> list[dict]:
    """Ratio of each candidate's rate to the summed mean rate of Carbon Mapper sources nearby.

    A source's rate is the mean over the plumes in which it was detected. Sources
    are intermittent and go undetected when quiet, so the sum is an upper estimate
    of their time-averaged contribution.
    """
    out = [dict(d, carbon_mapper=None) for d in detections]
    sources = [s for s in carbon_mapper_sources(cm_points) if s["mean_kg_hr"]]
    if not sources or not out:
        return out
    slon = np.array([s["lon"] for s in sources])
    slat = np.array([s["lat"] for s in sources])
    rate = np.array([s["mean_kg_hr"] for s in sources])
    for d in out:
        lon, lat = np.radians(d["detection"]["lon"]), np.radians(d["detection"]["lat"])
        a = (
            np.sin((np.radians(slat) - lat) / 2) ** 2
            + np.cos(lat) * np.cos(np.radians(slat)) * np.sin((np.radians(slon) - lon) / 2) ** 2
        )
        dist = 2 * EARTH_RADIUS_KM * np.arcsin(np.sqrt(np.clip(a, 0, 1)))
        near = dist <= radius_km
        q = d.get("emission", {}).get("q_kg_hr")
        info = {"n_sources": int(near.sum()), "nearest_km": round(float(dist.min()), 2)}
        if near.any():
            summed = float(rate[near].sum())
            info.update(summed_kg_hr=round(summed, 1), ratio=round(q / summed, 2) if q else None)
        d["carbon_mapper"] = info
    return out
