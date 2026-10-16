# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""ERA5 hourly wind from the Open-Meteo archive API, cached on disk.

Responses are cached per (lat, lon, day) so a re-run reproduces the same winds
without network access.
"""

from __future__ import annotations

import json
import math
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import requests

ENDPOINT = "https://archive-api.open-meteo.com/v1/archive"
SOURCE = "ERA5 via Open-Meteo"


def _cache_path(lat: float, lon: float, day: str) -> Path:
    root = Path(os.environ.get("METHANE_CACHE_DIR", "data/cache")) / "wind"
    root.mkdir(parents=True, exist_ok=True)
    return root / f"{day}_{lat:.3f}_{lon:.3f}.json"


def hourly(lat: float, lon: float, day: str) -> dict:
    path = _cache_path(lat, lon, day)
    if path.exists():
        return json.loads(path.read_text())
    params = {
        "latitude": round(lat, 3),
        "longitude": round(lon, 3),
        "start_date": day,
        "end_date": day,
        "hourly": "wind_speed_10m,wind_direction_10m,wind_speed_100m,wind_direction_100m",
        "wind_speed_unit": "ms",
        "timezone": "UTC",
    }
    r = requests.get(ENDPOINT, params=params, timeout=30)
    r.raise_for_status()
    data = r.json()["hourly"]
    path.write_text(json.dumps(data))
    return data


def _parse(when: str | datetime) -> datetime:
    if isinstance(when, datetime):
        return when
    return datetime.fromisoformat(when.replace("Z", "+00:00"))


def at(lat: float, lon: float, when: str | datetime, level: int = 10) -> tuple[float, float]:
    """(speed m/s, meteorological direction deg) linearly interpolated in time."""
    t = _parse(when)
    data = hourly(lat, lon, t.strftime("%Y-%m-%d"))
    h = t.hour + t.minute / 60 + t.second / 3600
    i0 = min(int(h), 22)
    f = h - i0
    spd = data[f"wind_speed_{level}m"]
    drn = data[f"wind_direction_{level}m"]
    u = [-spd[i] * math.sin(math.radians(drn[i])) for i in (i0, i0 + 1)]
    v = [-spd[i] * math.cos(math.radians(drn[i])) for i in (i0, i0 + 1)]
    ui = (1 - f) * u[0] + f * u[1]
    vi = (1 - f) * v[0] + f * v[1]
    return _to_speed_dir(ui, vi)


def _to_speed_dir(u: float, v: float) -> tuple[float, float]:
    return math.hypot(u, v), math.degrees(math.atan2(-u, -v)) % 360.0


def area_mean(
    bbox: tuple[float, float, float, float], when: str | datetime, level: int = 10, n: int = 3
) -> tuple[float, float]:
    """Vector-mean wind over an n x n grid of points inside bbox."""
    x0, y0, x1, y1 = bbox
    us, vs = [], []
    for lat in np.linspace(y0, y1, n + 2)[1:-1]:
        for lon in np.linspace(x0, x1, n + 2)[1:-1]:
            spd, drn = at(lat, lon, when, level)
            us.append(-spd * math.sin(math.radians(drn)))
            vs.append(-spd * math.cos(math.radians(drn)))
    # scalar mean speed: averaging vectors would cancel opposing flow and bias U low
    speed = float(np.mean(np.hypot(us, vs)))
    _, direction = _to_speed_dir(float(np.mean(us)), float(np.mean(vs)))
    return speed, direction
