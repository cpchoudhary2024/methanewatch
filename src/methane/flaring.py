# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Flaring against methane, at two scales.

Grid: Spearman correlation between log flaring radiant heat (persistent flares,
Flaring Monitor 2019-2021) and mean TROPOMI enhancement on a 0.05 degree grid.
Significance comes from a toroidal-shift permutation of the flaring field, which
keeps its spatial structure and breaks its alignment with methane. Moran's I and
a spatial-error regression are reported alongside.

Points: distance from each candidate to the nearest high-radiance flare, against
a uniform null and an infrastructure-conditioned null (one-sided Mann-Whitney).
"""

from __future__ import annotations


import numpy as np

from .read import read_scene
from .regions import Region
from .spatial_stats import clustered_null_points, enrichment_ci, uniform_null_points
from .validate import nearest_km


def _read_flares(csv_path: str, region: Region):
    import pandas as pd

    x0, y0, x1, y1 = region.bbox
    df = pd.read_csv(csv_path, low_memory=False)
    keep = df.longitude.between(x0, x1) & df.latitude.between(y0, y1)
    return df[keep].copy()


def load_flares(csv_path: str, region: Region):
    """Flare sites (rounded to 0.001 degree) with summed radiant heat."""
    df = _read_flares(csv_path, region)
    df["klon"] = df.longitude.round(3)
    df["klat"] = df.latitude.round(3)
    return df.groupby(["klon", "klat"]).agg(rh=("sum_rh", "sum")).reset_index()


def top_flare_points(sites, n: int = 1500) -> list[list[float]]:
    s = sites.sort_values("rh", ascending=False).head(n)
    return [[round(float(r.klon), 4), round(float(r.klat), 4)] for r in s.itertuples()]


def methane_enhancement_grid(paths: list[str], region: Region, lon_e, lat_e):
    """Mean per-cell enhancement over the scene median, and per-cell counts."""
    ny, nx = len(lat_e) - 1, len(lon_e) - 1
    msum = np.zeros((ny, nx))
    mcnt = np.zeros((ny, nx))
    for path in sorted(paths):
        try:
            s = read_scene(path, region.bbox)
        except (OSError, KeyError):
            continue
        if s.n < 300:
            continue
        enh = s.xch4 - np.median(s.xch4)
        sh, _, _ = np.histogram2d(s.lon, s.lat, [lon_e, lat_e], weights=enh)
        ch, _, _ = np.histogram2d(s.lon, s.lat, [lon_e, lat_e])
        msum += sh.T
        mcnt += ch.T
    with np.errstate(invalid="ignore"):
        return np.where(mcnt > 0, msum / mcnt, np.nan), mcnt


def grid_correlation(
    csv_path: str,
    paths: list[str],
    region: Region,
    *,
    res: float = 0.05,
    n_perm: int = 500,
    seed: int = 0,
    min_obs: int = 10,
    min_months: int = 6,
) -> dict:
    from scipy import stats

    from .spatial_stats import moran_geary, sar_error_check

    x0, y0, x1, y1 = region.bbox
    lon_e = np.arange(x0, x1 + res, res)
    lat_e = np.arange(y0, y1 + res, res)

    df = _read_flares(csv_path, region)
    key = [df.longitude.round(2), df.latitude.round(2)]
    df = df[df.groupby(key)["month"].transform("nunique") >= min_months]
    fh, _, _ = np.histogram2d(
        df.longitude, df.latitude, [lon_e, lat_e], weights=df.sum_rh.fillna(0)
    )
    flare_log = np.log1p(fh.T)
    meth, mcnt = methane_enhancement_grid(paths, region, lon_e, lat_e)

    observed = (mcnt >= min_obs) & np.isfinite(meth)
    mask = observed & (flare_log > 0)
    rho, p_naive = stats.spearmanr(flare_log[mask], meth[mask])

    rng = np.random.default_rng(seed)
    ny, nx = flare_log.shape
    null = []
    for _ in range(n_perm):
        shifted = np.roll(flare_log, (int(rng.integers(1, ny)), int(rng.integers(1, nx))), (0, 1))
        m = observed & (shifted > 0)
        if m.sum() > 30:
            null.append(stats.spearmanr(shifted[m], meth[m])[0])
    null = np.array(null)
    p_perm = (np.sum(null >= rho) + 1) / (null.size + 1)

    fx, my = flare_log[mask], meth[mask]
    idx = np.linspace(0, fx.size - 1, min(400, fx.size)).astype(int)
    result = {
        "rho": round(float(rho), 3),
        "p_naive": float(p_naive),
        "p_permutation": round(float(p_perm), 3),
        "n_cells": int(mask.sum()),
        "null_rho_sd": round(float(null.std()), 3),
        "scatter": [[round(float(fx[i]), 3), round(float(my[i]), 1)] for i in idx],
    }

    mg = moran_geary(meth, mask)
    result["moran_geary"] = {k: v for k, v in mg.items() if k not in ("coords", "vals")}
    if mg.get("significant_clustering"):
        result["sar_check"] = sar_error_check(mg["vals"], flare_log[mask], mg["coords"])
    return result


def _proximity_block(dd, days, null_lon, null_lat, method, flon, flat, tol_km):
    from scipy.stats import mannwhitneyu

    rd = nearest_km(null_lon, null_lat, flon, flat)
    frac_det = float(np.mean(dd <= tol_km))
    frac_null = float(np.mean(rd <= tol_km))
    return {
        "method": method,
        "median_km_null": round(float(np.median(rd)), 2),
        "frac_within_tol_null": round(frac_null, 3),
        "enrichment": round(frac_det / frac_null, 2) if frac_null else None,
        "enrichment_ci_days": enrichment_ci(dd <= tol_km, days, frac_null),
        "mannwhitney_p": round(float(mannwhitneyu(dd, rd, alternative="less").pvalue), 4),
    }


def detection_proximity(
    detections: list[dict],
    sites,
    region: Region,
    *,
    infra_lon,
    infra_lat,
    top_quantile: float = 0.5,
    tol_km: float = 5.0,
    n_null: int = 3000,
    seed: int = 0,
) -> dict:
    sites = sites[sites.rh > sites.rh.quantile(top_quantile)]
    flon, flat = sites.klon.values, sites.klat.values
    dlon = np.array([d["detection"]["lon"] for d in detections])
    dlat = np.array([d["detection"]["lat"] for d in detections])
    dd = nearest_km(dlon, dlat, flon, flat)

    days = [d.get("overpass_utc", "")[:10] for d in detections]
    ulon, ulat, umethod = uniform_null_points(n_null, region.bbox, seed=seed)
    clon, clat, cmethod = clustered_null_points(
        n_null, region.bbox, infra_lon, infra_lat, seed=seed
    )
    by_bandwidth = {}
    for h in (4.0, 16.0):
        blon, blat, _ = clustered_null_points(
            n_null, region.bbox, infra_lon, infra_lat, bandwidth_km=h, seed=seed
        )
        frac = float(np.mean(nearest_km(blon, blat, flon, flat) <= tol_km))
        by_bandwidth[f"{h:.0f}"] = round(float(np.mean(dd <= tol_km)) / frac, 2) if frac else None
    return {
        "n_flare_sites": int(len(sites)),
        "n_detections": int(dd.size),
        "tolerance_km": tol_km,
        "median_km_detections": round(float(np.median(dd)), 2),
        "frac_within_tol_detections": round(float(np.mean(dd <= tol_km)), 3),
        "uniform_null": _proximity_block(dd, days, ulon, ulat, umethod, flon, flat, tol_km),
        "clustered_null": _proximity_block(dd, days, clon, clat, cmethod, flon, flat, tol_km),
        "clustered_enrichment_by_bandwidth_km": by_bandwidth,
    }
