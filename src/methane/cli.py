# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Command-line entry point.

    methane catalog   --start 2026-06-01 --end 2026-09-01
    methane run       --start 2026-06-01 --end 2026-09-01     # download, detect, quantify
    methane validate | gap | flaring | report | provenance

Every subcommand reads and writes outputs/<region>/ by default.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import date, datetime
from pathlib import Path

from . import __version__
from .regions import get_region

STUDY_START, STUDY_END = "2026-06-01", "2026-09-01"


def granules_in_window(cache: str | Path, start: str, end: str) -> list[str]:
    """Cached OFFL granules whose sensing start date lies in [start, end)."""
    d0, d1 = date.fromisoformat(start), date.fromisoformat(end)
    out = []
    for p in sorted(Path(cache).glob("S5P_OFFL_L2__CH4____*.nc")):
        day = datetime.strptime(p.name[20:28], "%Y%m%d").date()
        if d0 <= day < d1:
            out.append(str(p))
    return out


def _cmd_catalog(args) -> int:
    from .catalog import search_granules, summarize

    region = get_region(args.region)
    print(summarize(search_granules(region, args.start, args.end)))
    return 0


def _process_granule(path: Path, region, granule_name: str, orbit: str):
    import dataclasses
    import math

    from . import wind
    from .albedo_screen import analyze_candidate
    from .detect import detect
    from .flux import P_SURF_DEFAULT, estimate_emission
    from .read import read_scene

    scene = read_scene(path, region.bbox)
    result = detect(scene)
    inverted = detect(dataclasses.replace(scene, xch4=-scene.xch4))
    overpass = {
        "day": scene.time_utc[:10],
        "time_utc": scene.time_utc,
        "orbit": orbit,
        "n_qa": scene.n,
        "qa_pass_fraction": round(scene.qa_pass_fraction, 3),
        "threshold_ppb": (
            round(result.threshold_ppb, 1) if math.isfinite(result.threshold_ppb) else None
        ),
        "n_candidates": len(result.candidates),
        "n_inverted": len(inverted.candidates),
        "inverted": [{"lon": round(c.lon, 4), "lat": round(c.lat, 4)} for c in inverted.candidates],
    }
    p_surf = scene.p_surf_median if math.isfinite(scene.p_surf_median) else P_SURF_DEFAULT

    records = []
    for cand in result.candidates:
        spd, wdir = wind.at(cand.lat, cand.lon, scene.time_utc, level=10)
        emis = estimate_emission(
            cand, spd, wind_source=f"{wind.SOURCE} 10 m", wind_direction_deg=wdir, p_surf=p_surf
        )
        records.append(
            {
                "granule": granule_name,
                "overpass_utc": scene.time_utc,
                "detection": cand.as_dict(),
                "emission": emis.as_dict(),
                "albedo_screen": analyze_candidate(path, cand.lon, cand.lat, wdir),
            }
        )
    print(f"  {granule_name[20:35]}  {scene.n:5d} px  {len(records)} candidate(s)")
    return overpass, records


def _download_with_retry(granule, attempts: int = 3):
    import requests

    from .download import download_granule, get_access_token

    for i in range(attempts):
        try:
            return download_granule(granule, get_access_token())
        except requests.RequestException as e:
            if i == attempts - 1:
                raise
            print(f"  retrying {granule.name[20:35]}: {e}", file=sys.stderr)


class _LocalGranule:
    def __init__(self, path: Path):
        self.path = path
        self.name = path.name
        self.orbit = path.name.split("_")[-4]


def _cmd_run(args) -> int:
    from .catalog import search_granules
    from .download import CredentialsError, cache_dir

    region = get_region(args.region)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    if args.cached:
        granules = [
            _LocalGranule(Path(p)) for p in granules_in_window(args.cache, args.start, args.end)
        ]
    else:
        granules = search_granules(region, args.start, args.end, limit=1000)
        granules = [g for g in granules if g.mode == "OFFL"]
    print(f"{len(granules)} OFFL granules for {region.name}, {args.start} to {args.end}")

    overpasses, records = [], []
    for g in granules:
        path = g.path if isinstance(g, _LocalGranule) else cache_dir() / g.name
        if not path.exists():
            try:
                path = _download_with_retry(g)
            except CredentialsError as e:
                print(e, file=sys.stderr)
                return 2
        ov, recs = _process_granule(path, region, g.name, g.orbit)
        overpasses.append(ov)
        records += recs

    (out / "detections.json").write_text(json.dumps(records, indent=2))
    (out / "summary.json").write_text(
        json.dumps(
            {
                "region": region.name,
                "bbox": list(region.bbox),
                "window": [args.start, args.end],
                "n_overpasses": len(overpasses),
                "n_usable": sum(1 for o in overpasses if o["threshold_ppb"] is not None),
                "n_candidates": len(records),
                "n_inverted_candidates": sum(o["n_inverted"] for o in overpasses),
                "overpasses": overpasses,
            },
            indent=2,
        )
    )
    print(f"{len(records)} candidates written to {out / 'detections.json'}")
    return 0


def _infra_points(region) -> tuple[list[float], list[float], str]:
    from .wells import get_wellheads

    pts, source = get_wellheads(region)
    if not pts:
        return [], [], source
    return [p["lon"] for p in pts], [p["lat"] for p in pts], source


def _noise_clusters(out: Path) -> list[dict]:
    """Clusters found in the sign-inverted field, shaped like detection records."""
    summary = json.loads((out / "summary.json").read_text())
    return [
        {"detection": c, "overpass_utc": o["time_utc"], "emission": {"q_kg_hr": 0.0}}
        for o in summary["overpasses"]
        for c in o.get("inverted", [])
    ]


def _cmd_validate(args) -> int:
    from .validate import build_validation, carbon_mapper_cross_check, fetch_carbon_mapper

    region = get_region(args.region)
    out = Path(args.out)
    records = json.loads((out / "detections.json").read_text())

    refs = fetch_carbon_mapper(region.bbox, "2016-01-01", args.reference_end)
    in_window = fetch_carbon_mapper(region.bbox, args.start, args.end)
    print(f"{len(refs)} Carbon Mapper plumes to {args.reference_end}, {len(in_window)} in window")
    (out / "carbon_mapper_plumes.json").write_text(
        json.dumps(
            [
                {
                    "lon": round(r.lon, 5),
                    "lat": round(r.lat, 5),
                    "emission_kg_hr": r.emission_kg_hr,
                    "date": r.date,
                }
                for r in refs
            ]
        )
    )

    ilon, ilat, source = _infra_points(region)
    print(f"clustered null from {len(ilon)} wells ({source})")
    report = build_validation(
        records,
        refs,
        bbox=region.bbox,
        infra_lon=ilon,
        infra_lat=ilat,
        references_in_window=len(in_window),
        reference_end=args.reference_end,
    )
    noise = _noise_clusters(out)
    if noise:
        report["noise_clusters"] = build_validation(
            noise, refs, bbox=region.bbox, infra_lon=ilon, infra_lat=ilat
        )
        report["noise_clusters"].pop("detections")
    (out / "validation.json").write_text(json.dumps(report, indent=2))
    cm = [{"lon": r.lon, "lat": r.lat, "emission_kg_hr": r.emission_kg_hr} for r in refs]
    records = carbon_mapper_cross_check(records, cm)
    (out / "detections.json").write_text(json.dumps(records, indent=2))
    for key in ("uniform_null", "clustered_null"):
        print(key, report[key])
    return 0


def _cmd_gap(args) -> int:
    from .gap import compute_gap

    region = get_region(args.region)
    out = Path(args.out)
    result = compute_gap(region, granules_in_window(args.cache, args.start, args.end))
    points = result.pop("ghgrp_points")
    (out / "gap.json").write_text(json.dumps(result, indent=2))
    (out / "ghgrp_points.json").write_text(
        json.dumps(
            [
                {"lon": p["lon"], "lat": p["lat"], "ch4_t": round(p["ch4_t"]), "name": p["name"]}
                for p in points
                if p.get("lon") and p.get("lat")
            ]
        )
    )
    td = result["top_down"]
    print(
        f"top-down {td['tg_yr_central']} Tg/yr (95% CI {td['ci_low']}-{td['ci_high']}, "
        f"{td['n_overpasses']} overpasses); GHGRP {result['bottom_up']['tg_yr']} Tg/yr"
    )
    return 0


def _cmd_flaring(args) -> int:
    from .flaring import detection_proximity, grid_correlation, load_flares, top_flare_points

    region = get_region(args.region)
    out = Path(args.out)
    det = json.loads((out / "detections.json").read_text())
    noise = _noise_clusters(out)
    paths = granules_in_window(args.cache, args.start, args.end)
    sites = load_flares(args.flares, region)
    ilon, ilat, _ = _infra_points(region)
    result = {
        "flare_data": "Flaring Monitor (Colorado School of Mines), VIIRS Nightfire 2019-2021",
        "n_flare_sites": int(len(sites)),
        "grid": grid_correlation(args.flares, paths, region),
        "proximity": detection_proximity(det, sites, region, infra_lon=ilon, infra_lat=ilat),
        "proximity_noise_clusters": detection_proximity(
            noise, sites, region, infra_lon=ilon, infra_lat=ilat
        ),
    }
    (out / "flaring.json").write_text(json.dumps(result, indent=2))
    (out / "flare_points.json").write_text(json.dumps(top_flare_points(sites)))
    g, p = result["grid"], result["proximity"]
    print(f"grid rho {g['rho']} (permutation p {g['p_permutation']}, {g['n_cells']} cells)")
    for key in ("uniform_null", "clustered_null"):
        print(key, p[key])
    return 0


def _cmd_report(args) -> int:
    from .report import build

    region = get_region(args.region)
    paths = granules_in_window(args.cache, args.start, args.end)
    print(f"wrote {build(region, Path(args.out), paths)}")
    return 0


def _cmd_provenance(args) -> int:
    from .provenance import build_provenance

    rec = build_provenance(args.out, granules_in_window(args.cache, args.start, args.end))
    print(f"hashed {len(rec['granules'])} granules and {len(rec['outputs'])} outputs")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="methane", description=__doc__.split("\n")[0])
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    commands = {
        "catalog": _cmd_catalog,
        "run": _cmd_run,
        "validate": _cmd_validate,
        "gap": _cmd_gap,
        "flaring": _cmd_flaring,
        "report": _cmd_report,
        "provenance": _cmd_provenance,
    }
    for name, fn in commands.items():
        sp = sub.add_parser(name)
        sp.add_argument("--region", default="permian")
        sp.add_argument("--start", default=STUDY_START)
        sp.add_argument("--end", default=STUDY_END)
        sp.add_argument("--out", default="outputs/permian")
        sp.add_argument("--cache", default="data/cache")
        if name == "flaring":
            sp.add_argument("--flares", default="data/flaring/flares_permian.csv")
        if name == "run":
            sp.add_argument(
                "--cached", action="store_true", help="use cached granules, no catalogue query"
            )
        if name == "validate":
            sp.add_argument("--reference-end", default=STUDY_END)
        sp.set_defaults(func=fn)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
