# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Null distributions and spatial-autocorrelation diagnostics.

Detections, flares, wells and Carbon Mapper sources all concentrate in the same
productive sub-basins. A control point drawn uniformly in the bounding box mostly
lands where none of them can be, so "closer than random" is true for anything
industrial. `clustered_null_points` draws controls from a Gaussian kernel density
of mapped infrastructure instead, which asks whether a detection is closer to a
feature than a typical industrial location is.
"""

from __future__ import annotations

import numpy as np


def uniform_null_points(n: int, bbox, *, seed: int = 0):
    x0, y0, x1, y1 = bbox
    rng = np.random.default_rng(seed)
    return rng.uniform(x0, x1, n), rng.uniform(y0, y1, n), "uniform in bounding box"


def clustered_null_points(
    n: int,
    bbox: tuple[float, float, float, float],
    infra_lon,
    infra_lat,
    *,
    bandwidth_km: float = 8.0,
    grid_res_deg: float = 0.02,
    seed: int = 0,
):
    """Draw `n` points from a kernel density surface over the infrastructure points.

    Returns (lon, lat, method). Falls back to a uniform draw, and says so in
    `method`, when fewer than 5 infrastructure points are supplied.
    """
    from scipy.ndimage import gaussian_filter

    x0, y0, x1, y1 = bbox
    rng = np.random.default_rng(seed)
    infra_lon = np.asarray(infra_lon, dtype=float)
    infra_lat = np.asarray(infra_lat, dtype=float)
    finite = np.isfinite(infra_lon) & np.isfinite(infra_lat)
    infra_lon, infra_lat = infra_lon[finite], infra_lat[finite]
    if infra_lon.size < 5:
        lon, lat, _ = uniform_null_points(n, bbox, seed=seed)
        return lon, lat, "uniform in bounding box (fewer than 5 infrastructure points)"

    lon_e = np.arange(x0, x1 + grid_res_deg, grid_res_deg)
    lat_e = np.arange(y0, y1 + grid_res_deg, grid_res_deg)
    hist, _, _ = np.histogram2d(infra_lon, infra_lat, bins=[lon_e, lat_e])

    km_per_deg = 111.0 * np.cos(np.radians((y0 + y1) / 2))
    sigma_cells = max(bandwidth_km / (km_per_deg * grid_res_deg), 0.5)
    density = np.clip(gaussian_filter(hist, sigma=sigma_cells, mode="constant"), 0, None)

    probs = (density / density.sum()).ravel()
    ix, iy = np.unravel_index(rng.choice(density.size, size=n, p=probs), density.shape)
    lon = lon_e[ix] + rng.uniform(0, grid_res_deg, n)
    lat = lat_e[iy] + rng.uniform(0, grid_res_deg, n)
    return lon, lat, f"kernel density of {infra_lon.size} wells, h = {bandwidth_km:.0f} km"


def enrichment_ci(
    within: np.ndarray, days: list[str], null_frac: float, *, n_boot: int = 2000, seed: int = 0
) -> list[float] | None:
    """95% interval for enrichment from a bootstrap that resamples whole days.

    Candidates from one overpass share wind, retrieval conditions and often a plume,
    so days rather than candidates are treated as the independent units.
    """
    if not null_frac or len(within) == 0:
        return None
    within = np.asarray(within, dtype=float)
    groups = {}
    for w, d in zip(within, days):
        groups.setdefault(d, []).append(w)
    blocks = list(groups.values())
    rng = np.random.default_rng(seed)
    stats = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(blocks), len(blocks))
        sample = np.concatenate([blocks[i] for i in pick])
        stats.append(sample.mean() / null_frac)
    return [round(float(v), 2) for v in np.percentile(stats, [2.5, 97.5])]


def moran_geary(
    grid: np.ndarray,
    valid_mask: np.ndarray,
    *,
    k: int = 8,
    permutations: int = 499,
    seed: int = 0,
    alpha: float = 0.05,
) -> dict:
    """Global Moran's I, Geary's C and local Moran hot/cold spots on a gridded field.

    Weights are row-standardised k-nearest neighbours in grid-cell units.
    """
    from esda.geary import Geary
    from esda.moran import Moran, Moran_Local
    from libpysal.weights import KNN

    iy, ix = np.where(valid_mask)
    vals = np.asarray(grid)[iy, ix].astype(float)
    if vals.size < k + 2:
        return {"n_cells": int(vals.size), "sufficient_data": False}

    coords = np.column_stack([ix, iy]).astype(float)
    w = KNN.from_array(coords, k=k)
    w.transform = "r"
    mi = Moran(vals, w, permutations=permutations)
    gc = Geary(vals, w, permutations=permutations)
    ml = Moran_Local(vals, w, permutations=permutations, seed=seed, alternative="two-sided")
    sig = ml.p_sim < alpha
    return {
        "n_cells": int(vals.size),
        "sufficient_data": True,
        "k_neighbors": k,
        "morans_i": round(float(mi.I), 4),
        "morans_i_p_sim": round(float(mi.p_sim), 4),
        "gearys_c": round(float(gc.C), 4),
        "gearys_c_p_sim": round(float(gc.p_sim), 4),
        "significant_clustering": bool(mi.p_sim < alpha and mi.I > 0),
        "n_hotspots_local": int(np.sum(sig & (ml.q == 1))),
        "n_coldspots_local": int(np.sum(sig & (ml.q == 3))),
        "alpha": alpha,
        "coords": coords,
        "vals": vals,
    }


def sar_error_check(y, x, coords, *, k: int = 8) -> dict:
    """OLS against maximum-likelihood spatial-error estimates of y ~ x.

    The weights encode contiguity, not a physical decay, so lambda is read
    qualitatively; the point is how far beta moves once the error term is spatial.
    """
    import warnings

    from libpysal.weights import KNN
    from spreg import OLS, ML_Error

    y = np.asarray(y, dtype=float).reshape(-1, 1)
    x = np.asarray(x, dtype=float).reshape(-1, 1)
    w = KNN.from_array(np.asarray(coords, dtype=float), k=k)
    w.transform = "r"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ols = OLS(y, x, w=w, spat_diag=True, moran=True)
        sar = ML_Error(y, x, w=w)

    moran_i, _z, moran_p = ols.moran_res
    b_ols, b_sar = float(ols.betas[1, 0]), float(sar.betas[1, 0])
    return {
        "n": int(y.size),
        "ols_beta": round(b_ols, 4),
        "ols_p": round(float(ols.t_stat[1][1]), 4),
        "ols_residual_morans_i": round(float(moran_i), 4),
        "ols_residual_morans_i_p": round(float(moran_p), 4),
        "sar_beta": round(b_sar, 4),
        "sar_beta_p": round(float(sar.z_stat[1][1]), 4),
        "sar_lambda": round(float(sar.lam), 4),
        "sar_lambda_p": round(float(sar.z_stat[-1][1]), 4),
        "beta_shrinkage_pct": round(100.0 * (1 - abs(b_sar) / abs(b_ols)), 1) if b_ols else None,
    }
