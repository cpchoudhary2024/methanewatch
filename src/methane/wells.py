# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Oil and gas well locations: a local state export when present, otherwise OpenStreetMap."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

from .regions import Region

OVERPASS_MIRRORS = [
    "https://overpass-api.de/api/interpreter",
    "https://overpass.kumi.systems/api/interpreter",
]


def load_wellheads(csv_path: str, region: Region) -> list[dict] | None:
    """Wells from a local CSV export (columns lon/longitude/x and lat/latitude/y)."""
    p = Path(csv_path)
    if not p.exists():
        return None
    import pandas as pd

    df = pd.read_csv(p, low_memory=False)
    cols = {c.lower(): c for c in df.columns}
    lon_c = next((cols[k] for k in ("lon", "longitude", "x") if k in cols), None)
    lat_c = next((cols[k] for k in ("lat", "latitude", "y") if k in cols), None)
    if lon_c is None or lat_c is None:
        raise ValueError(f"{csv_path}: no longitude/latitude columns in {list(df.columns)}")
    name_c = next(
        (cols[k] for k in ("operator", "lease_name", "well_name", "name") if k in cols), None
    )

    lon = pd.to_numeric(df[lon_c], errors="coerce")
    lat = pd.to_numeric(df[lat_c], errors="coerce")
    x0, y0, x1, y1 = region.bbox
    keep = (lon.between(x0, x1) & lat.between(y0, y1)).fillna(False)
    names = df.loc[keep, name_c].astype(str) if name_c else [""] * int(keep.sum())
    return [
        {"lon": float(a), "lat": float(b), "name": n}
        for a, b, n in zip(lon[keep], lat[keep], names)
    ]


def fetch_wellheads_osm(
    region: Region, *, timeout: int = 180, mirrors: list[str] | None = None
) -> list[dict] | None:
    """All OSM nodes tagged man_made=petroleum_well or oil_well in the region."""
    import requests

    x0, y0, x1, y1 = region.bbox
    query = (
        f"[out:json][timeout:{timeout}];"
        f'node["man_made"~"^(petroleum_well|oil_well)$"]({y0},{x0},{y1},{x1});'
        "out skel qt;"
    )
    headers = {"User-Agent": "methanewatch (research; github.com/cpchoudhary2024/methanewatch)"}
    for url in mirrors or OVERPASS_MIRRORS:
        try:
            r = requests.post(url, data={"data": query}, timeout=timeout + 30, headers=headers)
            r.raise_for_status()
            elements = r.json().get("elements", [])
        except (requests.RequestException, ValueError) as e:
            print(f"[wells] {url}: {e}", file=sys.stderr)
            continue
        pts = [{"lon": el["lon"], "lat": el["lat"], "name": ""} for el in elements if "lon" in el]
        if pts:
            return pts
    return None


def get_wellheads(
    region: Region,
    csv_path: str = "data/wellheads/wellheads_permian.csv",
    cache_path: str | None = None,
    max_age_days: float = 30.0,
) -> tuple[list[dict] | None, str]:
    """Local export if present, else a cached or live OSM query."""
    pts = load_wellheads(csv_path, region)
    if pts:
        return pts, f"local export ({csv_path})"

    cache = Path(cache_path or f"data/wellheads/osm_wells_{region.key}.json")
    if cache.exists() and (time.time() - cache.stat().st_mtime) < max_age_days * 86400:
        cached = json.loads(cache.read_text())
        if cached.get("region") == region.key and cached.get("points"):
            return cached["points"], f"OpenStreetMap (cached {cached.get('retrieved', '')})"

    pts = fetch_wellheads_osm(region)
    if pts:
        cache.parent.mkdir(parents=True, exist_ok=True)
        stamp = time.strftime("%Y-%m-%d")
        cache.write_text(json.dumps({"region": region.key, "retrieved": stamp, "points": pts}))
        return pts, f"OpenStreetMap (retrieved {stamp})"
    return None, "unavailable"
