# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Basin top-down estimate (TROPOMI mass balance) against GHGRP Subpart W self-reports.

Top-down, per overpass:

    Q = U * W * c(p_s) * (X_basin - X_background)

X_basin is the median XCH4 of QA-passing pixels inside the basin box and
X_background the median over the pixels in a 1 degree ring around it, after both
are corrected for the dependence of XCH4 on surface pressure and latitude fitted
on the ring.
W is the box width perpendicular to the wind, U the area-mean ERA5 wind at the
overpass hour and c converts ppb to kg m^-2 using the scene surface pressure.
This is the well-mixed-box form; if sources are spread uniformly along the wind
the outflow enhancement is twice the box mean and Q doubles, so that factor is
carried as a systematic term. The basin value is the median over overpasses.

The result depends strongly on how the background is defined (see `sensitivity`).
In summer the flow is mostly southerly and the upwind strip lies over the Eagle
Ford region, so it is not clean inflow. The estimate is reported for that reason
as a consistency check on the published inversion, not as an independent one.
"""

from __future__ import annotations

import math
import re
from dataclasses import asdict, dataclass

import numpy as np

from . import wind
from .flux import column_excess_kg_per_m2
from .read import read_scene
from .regions import Region

SECONDS_PER_YEAR = 365.25 * 86400
KM_PER_DEG = 111.32

ZHANG2020 = {
    "tg_yr": 2.7,
    "tg_yr_sigma": 0.5,
    "period": "2018-05 to 2019-03",
    "cite": "Zhang et al. 2020, Sci. Adv. 6, eaaz5120, doi:10.1126/sciadv.aaz5120",
}

# AAPG basin codes that appear in GHGRP facility names; 430 is the Permian.
_OTHER_BASINS = {"395", "360", "345", "355", "375", "160", "400", "420", "585", "515"}


@dataclass
class Overpass:
    granule: str
    time_utc: str
    n_basin: int
    n_upwind: int
    delta_ppb: float
    p_surf_pa: float
    wind_dir_deg: float
    u10_ms: float
    u100_ms: float
    width_km: float
    q_tg_yr: float


@dataclass
class TopDown:
    tg_yr_central: float
    ci_low: float
    ci_high: float
    delta_ppb_median: float
    u100_ms_median: float
    n_overpasses: int
    method: str = (
        "box mass balance against a terrain-corrected background ring, ERA5 100 m wind "
        "at overpass time; "
        "median over overpasses, 95% CI by bootstrap over overpasses"
    )


def _local_km(lon, lat, lon0, lat0):
    x = (np.asarray(lon) - lon0) * KM_PER_DEG * math.cos(math.radians(lat0))
    y = (np.asarray(lat) - lat0) * KM_PER_DEG
    return x, y


def _wind_frame(x, y, wind_from_deg):
    """Along-wind (positive downwind) and cross-wind coordinates."""
    t = math.radians(wind_from_deg)
    dx, dy = -math.sin(t), -math.cos(t)
    return x * dx + y * dy, -x * dy + y * dx


def _detrend(s, bg: np.ndarray) -> np.ndarray:
    """Remove the XCH4 dependence on surface pressure and latitude fitted on background pixels.

    Over higher terrain a larger share of the column is CH4-poor stratospheric air,
    so XCH4 falls with surface elevation independently of any surface source.
    """
    A = np.c_[np.ones(bg.sum()), s.p_surf[bg] / 100.0, s.lat[bg]]
    coef, *_ = np.linalg.lstsq(A, s.xch4[bg], rcond=None)
    return s.xch4 - coef[1] * (s.p_surf / 100.0) - coef[2] * s.lat


def scene_balance(
    path: str,
    region: Region,
    *,
    buffer_deg: float = 1.0,
    upwind_km: float = 100.0,
    background: str = "ring",
    terrain_correction: bool = True,
    min_basin: int = 300,
    min_background: int = 100,
) -> Overpass | None:
    x0, y0, x1, y1 = region.bbox
    big = (x0 - buffer_deg, y0 - buffer_deg, x1 + buffer_deg, y1 + buffer_deg)
    s = read_scene(path, big)
    inside = (s.lon >= x0) & (s.lon <= x1) & (s.lat >= y0) & (s.lat <= y1)
    if inside.sum() < min_basin or s.time_utc == "unknown":
        return None

    u100, wdir = wind.area_mean(region.bbox, s.time_utc, level=100)
    u10, _ = wind.area_mean(region.bbox, s.time_utc, level=10)

    lon0, lat0 = (x0 + x1) / 2, (y0 + y1) / 2
    px, py = _local_km(s.lon, s.lat, lon0, lat0)
    cx, cy = _local_km([x0, x1, x1, x0], [y0, y0, y1, y1], lon0, lat0)
    pa, pc = _wind_frame(px, py, wdir)
    ca, cc = _wind_frame(cx, cy, wdir)
    width_km = float(cc.max() - cc.min())

    if background == "upwind":
        bg = (~inside) & (pa < ca.min()) & (pa > ca.min() - upwind_km)
        bg &= (pc > cc.min()) & (pc < cc.max())
    else:
        bg = ~inside
    ok = np.isfinite(s.p_surf)
    bg &= ok
    if bg.sum() < min_background:
        return None

    xch4 = s.xch4
    if terrain_correction:
        xch4 = _detrend(s, bg)
    delta = float(np.median(xch4[inside & ok]) - np.median(xch4[bg]))
    p_s = s.p_surf_median
    q_kg_s = u100 * width_km * 1e3 * column_excess_kg_per_m2(delta, p_s)
    return Overpass(
        granule=s.source,
        time_utc=s.time_utc,
        n_basin=int(inside.sum()),
        n_upwind=int(bg.sum()),
        delta_ppb=round(delta, 2),
        p_surf_pa=round(p_s),
        wind_dir_deg=round(wdir, 1),
        u10_ms=round(u10, 2),
        u100_ms=round(u100, 2),
        width_km=round(width_km, 1),
        q_tg_yr=round(q_kg_s * SECONDS_PER_YEAR / 1e9, 3),
    )


def overpass_balances(paths: list[str], region: Region, **kw) -> list[Overpass]:
    out = []
    for f in sorted(paths):
        try:
            ob = scene_balance(f, region, **kw)
        except (OSError, KeyError, ValueError):
            continue
        if ob is not None:
            out.append(ob)
    return out


def _median_ci(values: np.ndarray, n_boot: int = 2000, seed: int = 0) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    boots = [np.median(rng.choice(values, values.size)) for _ in range(n_boot)]
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(lo), float(hi)


def top_down(obs: list[Overpass]) -> TopDown:
    q = np.array([o.q_tg_yr for o in obs])
    lo, hi = _median_ci(q)
    return TopDown(
        tg_yr_central=round(float(np.median(q)), 2),
        ci_low=round(lo, 2),
        ci_high=round(hi, 2),
        delta_ppb_median=round(float(np.median([o.delta_ppb for o in obs])), 2),
        u100_ms_median=round(float(np.median([o.u100_ms for o in obs])), 2),
        n_overpasses=len(obs),
    )


def sensitivity(paths: list[str], region: Region, obs: list[Overpass], bottom_up_tg: float):
    def row(case, values):
        tg = round(float(np.median(values)), 2)
        return {"case": case, "top_down_tg": tg, "ratio": round(tg / bottom_up_tg, 1)}

    q = np.array([o.q_tg_yr for o in obs])
    rows = [
        row("central: background ring, terrain-corrected, 100 m wind", q),
        row("10 m wind", [o.q_tg_yr * o.u10_ms / o.u100_ms for o in obs]),
        row("linear accumulation along wind (x2)", 2 * q),
    ]
    for label, kw in (
        ("background ring, no terrain correction", {"terrain_correction": False}),
        ("upwind strip background, terrain-corrected", {"background": "upwind"}),
        (
            "upwind strip background, no terrain correction",
            {"background": "upwind", "terrain_correction": False},
        ),
    ):
        alt = overpass_balances(paths, region, **kw)
        if alt:
            rows.append(row(label, [o.q_tg_yr for o in alt]))
    return rows


@dataclass
class BottomUp:
    tg_yr: float
    tonnes_yr: float
    year: int
    n_reports: int
    top: list
    source: str = "EPA GHGRP Subpart W, attributed by AAPG basin code 430"


def _basin_code(name: str) -> str | None:
    m = re.search(r"\((\d{3})\)|\b(\d{3})\b", name or "")
    return (m.group(1) or m.group(2)) if m else None


def _ghgrp_permian_year(region: Region, year: int, states=("TX", "NM")) -> list[dict]:
    """Subpart W CH4 reports attributable to the Permian for one reporting year.

    Onshore production and gathering are reported per basin, so reports named
    for basin 430 are kept wherever their nominal coordinate lies, and reports
    tagged to another basin are dropped even when their coordinate is inside
    the box.
    """
    import requests

    base = "https://data.epa.gov/efservice"
    x0, y0, x1, y1 = region.bbox
    rows: list[dict] = []
    for st in states:
        off = 0
        while True:
            url = (
                f"{base}/V_GHG_EMITTER_SUBPART/SUBPART_NAME/W/GAS_CODE/CH4/"
                f"STATE/{st}/YEAR/{year}/ROWS/{off}:{off + 499}/JSON"
            )
            page = requests.get(url, timeout=60).json()
            rows += page
            if len(page) < 500:
                break
            off += 500

    kept = []
    for x in rows:
        name = x.get("facility_name") or ""
        bc = _basin_code(name)
        lon, lat = x.get("longitude"), x.get("latitude")
        permian = "permian" in name.lower() or bc == "430"
        other = bc in _OTHER_BASINS and not permian
        in_box = lon is not None and lat is not None and x0 <= lon <= x1 and y0 <= lat <= y1
        if permian or (in_box and not other):
            kept.append(
                {
                    "name": name,
                    "ch4_t": float(x.get("ghg_quantity") or 0.0),
                    "lon": lon,
                    "lat": lat,
                    "county": x.get("county"),
                }
            )
    return kept


def bottom_up_ghgrp(region: Region, year: int = 2023) -> tuple[BottomUp, list[dict]]:
    kept = _ghgrp_permian_year(region, year)
    total = sum(k["ch4_t"] for k in kept)
    top = sorted(kept, key=lambda k: -k["ch4_t"])[:8]
    bu = BottomUp(
        tg_yr=round(total / 1e6, 3),
        tonnes_yr=round(total),
        year=year,
        n_reports=len(kept),
        top=[{"name": t["name"], "ch4_t": round(t["ch4_t"])} for t in top],
    )
    return bu, kept


def ghgrp_trend(region: Region, years=range(2011, 2024)) -> list[dict]:
    out = []
    for y in years:
        kept = _ghgrp_permian_year(region, y)
        out.append(
            {
                "year": y,
                "tg": round(sum(k["ch4_t"] for k in kept) / 1e6, 3),
                "n_reports": len(kept),
            }
        )
    return out


def compute_gap(
    region: Region,
    paths: list[str],
    *,
    with_trend: bool = True,
) -> dict:
    obs = overpass_balances(paths, region)
    td = top_down(obs)
    bu, points = bottom_up_ghgrp(region)
    sens = sensitivity(paths, region, obs, bu.tg_yr)
    trend = ghgrp_trend(region) if with_trend else []

    anchor = ZHANG2020["tg_yr"]
    return {
        "region": region.name,
        "top_down": asdict(td),
        "overpasses": [asdict(o) for o in obs],
        "anchor": ZHANG2020,
        "bottom_up": asdict(bu),
        "ratio_vs_anchor": round(anchor / bu.tg_yr, 1),
        "ratio_vs_this_study": round(td.tg_yr_central / bu.tg_yr, 1),
        "sensitivity": sens,
        "range_over_cases_tg_yr": [
            min(r["top_down_tg"] for r in sens),
            max(r["top_down_tg"] for r in sens),
        ],
        "ghgrp_trend": trend,
        "ghgrp_points": points,
    }
