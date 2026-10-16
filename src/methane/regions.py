# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Study regions as (lon_min, lat_min, lon_max, lat_max) in WGS84."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Region:
    key: str
    name: str
    bbox: tuple[float, float, float, float]

    @property
    def lon_min(self) -> float:
        return self.bbox[0]

    @property
    def lat_min(self) -> float:
        return self.bbox[1]

    @property
    def lon_max(self) -> float:
        return self.bbox[2]

    @property
    def lat_max(self) -> float:
        return self.bbox[3]

    def wkt_polygon(self) -> str:
        x0, y0, x1, y1 = self.bbox
        return f"POLYGON(({x0} {y0},{x1} {y0},{x1} {y1},{x0} {y1},{x0} {y0}))"


REGIONS: dict[str, Region] = {
    "permian": Region(
        key="permian",
        name="Permian Basin (Delaware + Midland)",
        bbox=(-105.0, 30.8, -101.0, 33.5),
    ),
    "bakken": Region(
        key="bakken",
        name="Bakken (Williston Basin, North Dakota)",
        bbox=(-104.5, 47.0, -102.0, 48.8),
    ),
    "eagleford": Region(
        key="eagleford",
        name="Eagle Ford (South Texas)",
        bbox=(-100.0, 27.8, -96.5, 30.0),
    ),
    "anadarko": Region(
        key="anadarko",
        name="Anadarko Basin (Oklahoma and Texas Panhandle)",
        bbox=(-100.5, 34.5, -97.0, 36.8),
    ),
    "dj": Region(
        key="dj",
        name="Denver-Julesburg Basin (Colorado)",
        bbox=(-105.2, 39.8, -103.8, 41.0),
    ),
    "delaware": Region(
        key="delaware",
        name="Delaware sub-basin (W. Texas / SE New Mexico)",
        bbox=(-104.6, 31.3, -103.2, 32.5),
    ),
    "midland": Region(
        key="midland",
        name="Midland sub-basin (W. Texas)",
        bbox=(-102.6, 31.5, -101.4, 32.7),
    ),
}


def get_region(key: str) -> Region:
    try:
        return REGIONS[key.lower()]
    except KeyError:
        raise SystemExit(f"Unknown region '{key}'. Known regions: {', '.join(REGIONS)}")
