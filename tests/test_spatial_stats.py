import numpy as np

from methane.spatial_stats import clustered_null_points, moran_geary, sar_error_check

PERMIAN_BBOX = (-105.0, 30.8, -101.0, 33.5)


def test_clustered_null_falls_back_to_uniform_with_no_infra_points():
    lon, lat, method = clustered_null_points(500, PERMIAN_BBOX, [], [], seed=0)
    assert len(lon) == len(lat) == 500
    assert method.startswith("uniform")
    assert (lon >= PERMIAN_BBOX[0]).all() and (lon <= PERMIAN_BBOX[2]).all()
    assert (lat >= PERMIAN_BBOX[1]).all() and (lat <= PERMIAN_BBOX[3]).all()


def test_clustered_null_falls_back_with_fewer_than_five_points():
    lon, lat, method = clustered_null_points(
        100, PERMIAN_BBOX, [-103.5, -103.4], [32.0, 32.1], seed=0
    )
    assert method.startswith("uniform")
    assert len(lon) == 100


def test_clustered_null_concentrates_near_dense_infrastructure():
    rng = np.random.default_rng(1)
    infra_lon = rng.normal(-103.5, 0.02, 200)
    infra_lat = rng.normal(32.0, 0.02, 200)

    lon, lat, method = clustered_null_points(
        1000, PERMIAN_BBOX, infra_lon, infra_lat, bandwidth_km=5.0, seed=0
    )
    assert method.startswith("kernel density")
    dist = np.sqrt((lon - (-103.5)) ** 2 + (lat - 32.0) ** 2)
    assert np.mean(dist < 0.1) > 0.5


def test_clustered_null_drops_nan_infra_points():
    infra_lon = [np.nan, -103.5, -103.4, -103.6, -103.3, -103.55]
    infra_lat = [32.0, np.nan, 32.1, 31.9, 32.05, 32.02]
    lon, lat, method = clustered_null_points(50, PERMIAN_BBOX, infra_lon, infra_lat, seed=0)
    assert len(lon) == 50
    assert np.isfinite(lon).all() and np.isfinite(lat).all()


def test_clustered_null_reproducible_with_same_seed():
    rng = np.random.default_rng(2)
    infra_lon = rng.uniform(-104, -103, 50)
    infra_lat = rng.uniform(31.5, 32.5, 50)
    lon1, lat1, _ = clustered_null_points(20, PERMIAN_BBOX, infra_lon, infra_lat, seed=7)
    lon2, lat2, _ = clustered_null_points(20, PERMIAN_BBOX, infra_lon, infra_lat, seed=7)
    assert np.array_equal(lon1, lon2) and np.array_equal(lat1, lat2)


def _clustered_grid(n=20, seed=0):
    rng = np.random.default_rng(seed)
    x = np.linspace(0, 1, n)
    y = np.linspace(0, 1, n)
    X, Y = np.meshgrid(x, y)
    grid = X + rng.normal(0, 0.01, X.shape)
    mask = np.ones_like(grid, dtype=bool)
    return grid, mask


def _random_grid(n=20, seed=0):
    rng = np.random.default_rng(seed)
    grid = rng.normal(0, 1, (n, n))
    mask = np.ones_like(grid, dtype=bool)
    return grid, mask


def test_moran_geary_detects_strong_clustering():
    grid, mask = _clustered_grid()
    result = moran_geary(grid, mask, k=8, permutations=199, seed=0)
    assert result["sufficient_data"] is True
    assert result["morans_i"] > 0.5
    assert result["significant_clustering"] is True
    assert result["n_hotspots_local"] + result["n_coldspots_local"] > 0


def test_moran_geary_insufficient_data_reports_gracefully():
    grid = np.array([[1.0, 2.0], [3.0, 4.0]])
    mask = np.array([[True, False], [False, False]])
    result = moran_geary(grid, mask, k=8)
    assert result["sufficient_data"] is False
    assert result["n_cells"] == 1


def test_moran_geary_all_invalid_mask_is_handled():
    grid = np.zeros((5, 5))
    mask = np.zeros((5, 5), dtype=bool)
    result = moran_geary(grid, mask, k=4)
    assert result["sufficient_data"] is False
    assert result["n_cells"] == 0


def test_sar_error_check_runs_on_synthetic_correlated_field():
    grid, mask = _clustered_grid(seed=3)
    mg = moran_geary(grid, mask, k=8, permutations=199, seed=0)
    assert mg["significant_clustering"] is True

    rng = np.random.default_rng(0)
    y = mg["vals"] + rng.normal(0, 0.05, mg["vals"].shape)
    result = sar_error_check(y, mg["vals"], mg["coords"], k=8)

    assert "ols_beta" in result and "sar_beta" in result
    assert "sar_lambda" in result
    assert isinstance(result["sar_beta_p"], float)


def test_sar_error_check_with_no_relationship_gives_small_beta():
    grid, mask = _random_grid(seed=5)
    mg = moran_geary(grid, mask, k=8, permutations=199, seed=0)
    rng = np.random.default_rng(9)
    unrelated_x = rng.normal(0, 1, mg["vals"].shape)
    result = sar_error_check(mg["vals"], unrelated_x, mg["coords"], k=8)
    assert abs(result["ols_beta"]) < 1.0


def test_enrichment_ci_resamples_days_not_candidates():
    from methane.spatial_stats import enrichment_ci

    assert enrichment_ci([True] * 6, ["a"] * 3 + ["b"] * 3, 0.5) == [2.0, 2.0]
    lo, hi = enrichment_ci([True, True, True, False, False, False], ["a"] * 3 + ["b"] * 3, 0.5)
    assert lo == 0.0 and hi == 2.0
    assert enrichment_ci([], [], 0.5) is None
