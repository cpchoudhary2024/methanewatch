# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""TROPOMI L2 CH4 granule search on the Copernicus Data Space OData catalogue.

The catalogue is open; only downloads need an account.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Iterable

import requests

from .regions import Region

CATALOGUE_ODATA = "https://catalogue.dataspace.copernicus.eu/odata/v1/Products"
_TIMEOUT = 60


@dataclass(frozen=True)
class Granule:
    product_id: str
    name: str
    start: datetime
    size_bytes: int
    online: bool

    @property
    def _tokens(self) -> list[str]:
        # S5P_OFFL_L2__CH4____<start>_<stop>_<orbit>_<coll>_<proc>_<prod>.nc
        return [t for t in self.name.split("_") if t]

    @property
    def mode(self) -> str:
        t = self._tokens
        return t[1] if len(t) > 1 else "?"

    @property
    def orbit(self) -> str:
        t = self._tokens
        return t[6] if len(t) > 6 else "?"

    @property
    def size_mb(self) -> float:
        return round(self.size_bytes / 1e6, 1)


def _parse_dt(s: str) -> datetime:
    return datetime.fromisoformat(s.replace("Z", "+00:00")).astimezone(timezone.utc)


def _iso(d: datetime | str) -> str:
    if isinstance(d, str):
        d = datetime.fromisoformat(d).replace(tzinfo=timezone.utc)
    return d.astimezone(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S.000Z")


def search_granules(
    region: Region,
    start: datetime | str,
    end: datetime | str,
    *,
    prefer_offl: bool = True,
    limit: int = 100,
    session: requests.Session | None = None,
) -> list[Granule]:
    """Granules intersecting `region` with sensing start in [start, end), one per orbit."""
    sess = session or requests.Session()
    poly = region.wkt_polygon()
    filt = (
        "Collection/Name eq 'SENTINEL-5P' "
        "and contains(Name,'L2__CH4') "
        f"and OData.CSC.Intersects(area=geography'SRID=4326;{poly}') "
        f"and ContentDate/Start gt {_iso(start)} "
        f"and ContentDate/Start lt {_iso(end)}"
    )
    params = {
        "$filter": filt,
        "$orderby": "ContentDate/Start",
        "$top": str(limit),
    }
    resp = sess.get(CATALOGUE_ODATA, params=params, timeout=_TIMEOUT)
    resp.raise_for_status()
    values = resp.json().get("value", [])

    granules: list[Granule] = []
    for v in values:
        granules.append(
            Granule(
                product_id=v["Id"],
                name=v["Name"],
                start=_parse_dt(v["ContentDate"]["Start"]),
                size_bytes=int(v.get("ContentLength", 0)),
                online=bool(v.get("Online", True)),
            )
        )

    if prefer_offl:
        granules = _dedupe_prefer_offl(granules)
    return granules


def _dedupe_prefer_offl(granules: Iterable[Granule]) -> list[Granule]:
    rank = {"OFFL": 0, "RPRO": 1, "NRTI": 2}
    best: dict[str, Granule] = {}
    for g in granules:
        cur = best.get(g.orbit)
        if cur is None or rank.get(g.mode, 9) < rank.get(cur.mode, 9):
            best[g.orbit] = g
    return sorted(best.values(), key=lambda g: g.start)


def summarize(granules: list[Granule]) -> str:
    if not granules:
        return "No granules found for that region/time window."
    lines = [f"{len(granules)} granule(s):"]
    total = 0.0
    for g in granules:
        total += g.size_mb
        flag = "" if g.online else "  [offline, needs LTA retrieval]"
        lines.append(
            f"  {g.start:%Y-%m-%d %H:%M} UTC  {g.mode}  orbit {g.orbit}  "
            f"{g.size_mb:6.1f} MB  {g.name}{flag}"
        )
    lines.append(f"  total download: {total:.0f} MB")
    return "\n".join(lines)
