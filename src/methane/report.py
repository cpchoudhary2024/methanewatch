# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Static results page built from the JSON outputs of every analysed basin."""

from __future__ import annotations

import json
from html import escape
from pathlib import Path

from .regions import REGIONS
from .spatial_stats import clustered_null_points, uniform_null_points
from .wells import get_wellheads

REPO = "https://github.com/cpchoudhary2024/methanewatch"
N_CONTROLS = 300
MIN_TESTABLE = 10


def _load(out: Path, name: str, default):
    p = out / name
    return json.loads(p.read_text()) if p.exists() else default


def _r(v, nd=2):
    return None if v is None or v != v else round(v, nd)


def _basin(out: Path) -> dict | None:
    summary = _load(out, "summary.json", None)
    if not summary:
        return None
    region = REGIONS[out.name.split("-")[0]]
    det = _load(out, "detections.json", [])
    val = _load(out, "validation.json", {})
    prox = _load(out, "flaring.json", {}).get("proximity", {})

    def block(src, key):
        b = src.get(key, {})
        return {
            "e": _r(b.get("enrichment")),
            "ci": b.get("enrichment_ci_days"),
            "p": _r(b.get("mannwhitney_p"), 3),
            "frac": _r(b.get("frac_within_tol_null"), 3),
        }

    basin = {
        "id": out.name,
        "name": region.name.split(" (")[0],
        "year": summary["window"][0][:4],
        "bbox": list(region.bbox),
        "n_overpasses": summary["n_overpasses"],
        "n_candidates": summary["n_candidates"],
        "n_noise": summary.get("n_inverted_candidates"),
        "frac_det": _r(prox.get("frac_within_tol_detections"), 3),
        "flare": {"uniform": block(prox, "uniform_null"), "wells": block(prox, "clustered_null")},
        "testable": summary["n_candidates"] >= MIN_TESTABLE,
    }
    if not basin["testable"]:
        return basin

    opportunity = _load(out, "opportunity.json", {})
    basin["opportunity"] = {
        "uniform": opportunity.get("uniform", {}),
        "wells": opportunity.get("wells", {}),
    }

    wells, _ = get_wellheads(region)
    ilon = [w["lon"] for w in wells or []]
    ilat = [w["lat"] for w in wells or []]
    ulon, ulat, _ = uniform_null_points(N_CONTROLS, region.bbox, seed=1)
    clon, clat, _ = clustered_null_points(N_CONTROLS, region.bbox, ilon, ilat, seed=1)
    x0, y0, x1, y1 = region.bbox
    basin["candidates"] = [
        [
            d["detection"]["lat"],
            d["detection"]["lon"],
            round(d["emission"]["q_kg_hr"] / 1000),
            d["overpass_utc"][:10],
        ]
        for d in det
    ]
    basin["controls"] = {
        "uniform": [[round(a, 3), round(b, 3)] for a, b in zip(ulat, ulon)],
        "wells": [[round(a, 3), round(b, 3)] for a, b in zip(clat, clon)],
    }
    flare_points = [
        [p[1], p[0]]
        for p in _load(out, "flare_points.json", [])
        if x0 <= p[0] <= x1 and y0 <= p[1] <= y1
    ]
    stride = max(1, len(flare_points) // 800)
    basin["flares"] = flare_points[::stride][:800]
    basin["cm_wells"] = _r(val.get("clustered_null", {}).get("enrichment"))
    return basin


def build(region, out: Path, paths: list[str]) -> Path:
    roots = sorted(p for p in out.parent.iterdir() if (p / "summary.json").exists())
    basins = [b for b in (_basin(p) for p in roots) if b]
    basins.sort(key=lambda b: (not b["testable"], b["id"] != out.name, -int(b["year"])))
    (out / "index.html").write_text(_page(basins, _load(out, "gap.json", {})))
    return out / "index.html"


def _summary_rows(basins) -> str:
    rows = []
    for b in basins:
        fu, fw = b["flare"]["uniform"], b["flare"]["wells"]
        if b["testable"]:
            cells = (
                f"<td>{fu['e']:.2f}&times; <span class=mut>p = {fu['p']}</span></td>"
                f"<td>{fw['e']:.2f}&times; <span class=mut>p = {fw['p']}</span></td>"
            )
        else:
            cells = "<td colspan=2 class=mut>too few plumes to test</td>"
        rows.append(
            f"<tr><td>{escape(b['name'])}</td><td>{b['year']}</td><td>{b['n_overpasses']}</td>"
            f"<td>{b['n_candidates']}</td><td>{b['n_noise']}</td>{cells}</tr>"
        )
    return "".join(rows)


def _page(basins: list[dict], gap: dict) -> str:
    """Render a fixed-window research page from the checked JSON outputs."""
    base = Path(__file__).parent
    bu = gap.get("bottom_up", {})
    anchor = gap.get("anchor", {})
    reported = bu.get("tg_yr", 0.0)
    observed = anchor.get("tg_yr", 0.0)
    sigma = anchor.get("tg_yr_sigma", 0.0)
    explorer = [b for b in basins if b["testable"]]
    template = (base / "site_template.html").read_text()
    replacements = {
        "@@CSS@@": (base / "site.css").read_text(),
        "@@JS@@": (base / "site.js").read_text(),
        "@@DATA@@": json.dumps(explorer, separators=(",", ":")),
        "@@ROWS@@": _summary_rows(basins),
        "@@REPO@@": REPO,
        "@@REPORTED@@": f"{reported:.3f}",
        "@@OBSERVED@@": f"{observed:.1f}",
        "@@SIGMA@@": f"{sigma:.1f}",
        "@@RATIO@@": f"{observed / reported:.1f}" if reported else "n/a",
    }
    for token, value in replacements.items():
        template = template.replace(token, value)
    return template
