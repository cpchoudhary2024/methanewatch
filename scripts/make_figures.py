"""Manuscript figures from the run outputs in outputs/permian, written to paper/figures."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from methane.wells import get_wellheads
from methane.flaring import load_flares
from methane.regions import get_region
from methane.spatial_stats import clustered_null_points, uniform_null_points
from methane.validate import nearest_km

mpl.rcParams.update(
    {
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 300,
    }
)

OUT = Path("outputs/permian")
FIG = Path("paper/figures")
REGION = get_region("permian")


def _load(name):
    return json.loads((OUT / name).read_text())


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.png", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def fig_region(det, cm):
    x0, y0, x1, y1 = REGION.bbox
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    hb = ax.hexbin(
        [c["lon"] for c in cm],
        [c["lat"] for c in cm],
        gridsize=45,
        extent=(x0, x1, y0, y1),
        cmap="Greys",
        mincnt=1,
        linewidths=0,
    )
    fig.colorbar(hb, ax=ax, fraction=0.045, pad=0.02, label="Carbon Mapper plumes per cell")
    ax.plot([-103.06, -103.06, x0], [y1, 32.0, 32.0], color="#334155", lw=0.8, ls="--")
    for key in ("delaware", "midland"):
        bx0, by0, bx1, by1 = get_region(key).bbox
        ax.add_patch(
            plt.Rectangle((bx0, by0), bx1 - bx0, by1 - by0, fill=False, ec="#1d4ed8", lw=1)
        )
        ax.text(bx0 + 0.05, by1 + 0.05, key.capitalize(), color="#1d4ed8", fontsize=8)
    near = np.array([(d.get("carbon_mapper") or {}).get("nearest_km", 99) <= 10 for d in det])
    lon = np.array([d["detection"]["lon"] for d in det])
    lat = np.array([d["detection"]["lat"] for d in det])
    size = np.array([12 + d["emission"]["q_kg_hr"] / 1500 for d in det])
    ax.scatter(
        lon[near],
        lat[near],
        s=size[near],
        c="#1d4ed8",
        ec="w",
        lw=0.5,
        label="within 10 km of a CM plume",
    )
    ax.scatter(
        lon[~near],
        lat[~near],
        s=size[~near],
        c="#ea580c",
        ec="w",
        lw=0.5,
        label="farther than 10 km",
    )
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    ax.legend(loc="lower left", fontsize=7)
    save(fig, "fig1_region")


def fig_gap(gap):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw={"wspace": 0.45})
    td, anchor = gap["top_down"], gap["anchor"]
    rows = [
        ("GHGRP 2023", gap["bottom_up"]["tg_yr"], None),
        ("Zhang et al. 2020", anchor["tg_yr"], (anchor["tg_yr_sigma"],) * 2),
        (
            "Box mass balance",
            td["tg_yr_central"],
            (td["tg_yr_central"] - td["ci_low"], td["ci_high"] - td["tg_yr_central"]),
        ),
    ]
    y = np.arange(len(rows))[::-1]
    for yi, (label, v, err) in zip(y, rows):
        a.barh(yi, v, color="#94a3b8" if err is None else "#1d4ed8", height=0.6)
        if err:
            a.errorbar(v, yi, xerr=[[err[0]], [err[1]]], color="k", capsize=3, lw=1)
    lo, hi = gap["range_over_cases_tg_yr"]
    a.plot(
        [lo, hi],
        [y[-1] - 0.42] * 2,
        color="#ea580c",
        lw=2,
        label="range over background and wind cases",
    )
    a.set_yticks(y)
    a.set_yticklabels([r[0] for r in rows])
    a.axvline(0, color="k", lw=0.6)
    a.set_xlabel("Tg CH$_4$ yr$^{-1}$")
    a.legend(fontsize=7, loc="lower right")
    a.set_title("(a)", loc="left")

    trend = [(t["year"], t["tg"]) for t in gap["ghgrp_trend"]]
    b.plot(*zip(*trend), "-o", color="#1d4ed8", ms=3)
    b.set_ylim(0, None)
    b.set_xlabel("Reporting year")
    b.set_ylabel("Subpart W CH$_4$ (Tg)")
    b.set_title("(b)", loc="left")
    save(fig, "fig2_gap")


def fig_nulls(det, flaring):
    sites = load_flares("data/flaring/flares_permian.csv", REGION)
    sites = sites[sites.rh > sites.rh.quantile(0.5)]
    flon, flat = sites.klon.values, sites.klat.values
    wells, _ = get_wellheads(REGION)
    dist = {
        "candidates": nearest_km(
            [d["detection"]["lon"] for d in det], [d["detection"]["lat"] for d in det], flon, flat
        ),
        "uniform null": nearest_km(*uniform_null_points(3000, REGION.bbox)[:2], flon, flat),
        "well-density null": nearest_km(
            *clustered_null_points(
                3000, REGION.bbox, [w["lon"] for w in wells], [w["lat"] for w in wells]
            )[:2],
            flon,
            flat,
        ),
    }
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw={"wspace": 0.35})
    style = {
        "candidates": ("#ea580c", "-"),
        "uniform null": ("#94a3b8", "-"),
        "well-density null": ("#1d4ed8", "--"),
    }
    for name, d in dist.items():
        c, ls = style[name]
        a.step(np.sort(d), np.arange(1, d.size + 1) / d.size, color=c, ls=ls, lw=1.5, label=name)
    a.axvline(5, color="k", ls=":", lw=0.8)
    a.set_xlim(0, 30)
    a.set_xlabel("Distance to nearest flare (km)")
    a.set_ylabel("Cumulative fraction")
    a.legend(fontsize=7, loc="lower right")
    a.set_title("(a)", loc="left")

    sar = flaring["grid"].get("sar_check", {})
    vals = [sar.get("ols_beta", np.nan), sar.get("sar_beta", np.nan)]
    b.bar(["OLS", "Spatial error"], vals, color=["#94a3b8", "#1d4ed8"], width=0.55)
    for i, v in enumerate(vals):
        b.text(i, v, f"{v:.3f}", ha="center", va="bottom", fontsize=8)
    b.set_ylabel(r"$\beta$, enhancement on log radiant heat")
    b.set_title("(b)", loc="left")
    save(fig, "fig3_flaring")


if __name__ == "__main__":
    det = _load("detections.json")
    fig_region(det, _load("carbon_mapper_plumes.json"))
    fig_gap(_load("gap.json"))
    fig_nulls(det, _load("flaring.json"))
