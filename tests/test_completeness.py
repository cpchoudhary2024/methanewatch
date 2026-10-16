import math

import numpy as np
import pytest

from methane import completeness
from methane.detect import Candidate
from methane.flux import estimate_emission
from methane.read import Scene
from methane.regions import get_region


def test_floor_equals_ime_rate_of_a_minimum_cluster_at_threshold():
    cell_km2, min_pixels, threshold, wind, p = 26.0, 3, 30.0, 4.0, 91000.0
    floor = completeness.scene_floor_t_hr(threshold, wind, p, cell_km2, min_pixels)
    cand = Candidate(0, 0, threshold, threshold, min_pixels, min_pixels * cell_km2, 3.0)
    q = estimate_emission(cand, wind, p_surf=p).q_kg_hr / 1000
    assert floor == pytest.approx(q)


def test_floor_scales_linearly_with_wind_and_threshold():
    base = completeness.scene_floor_t_hr(30.0, 4.0, 91000.0, 26.0, 3)
    assert completeness.scene_floor_t_hr(60.0, 4.0, 91000.0, 26.0, 3) == pytest.approx(2 * base)
    assert completeness.scene_floor_t_hr(30.0, 8.0, 91000.0, 26.0, 3) == pytest.approx(2 * base)


def test_completeness_counts_only_covered_overpasses(mocker):
    region = get_region("delaware")
    x0, y0, x1, y1 = region.bbox
    lon, lat = np.meshgrid(np.linspace(x0, x1, 60), np.linspace(y0, y1, 60))
    rng = np.random.default_rng(0)
    scene = Scene(
        lon=lon.ravel(),
        lat=lat.ravel(),
        xch4=1900 + rng.normal(0, 5, lon.size),
        precision=np.full(lon.size, 8.0),
        qa=np.full(lon.size, 0.9),
        p_surf=np.full(lon.size, 91000.0),
        source="s.nc",
        n_total=lon.size,
        time_utc="2024-06-05T19:30:00Z",
    )
    mocker.patch.object(completeness, "read_scene", return_value=scene)
    mocker.patch.object(completeness.wind, "area_mean", side_effect=[(1.0, 0.0), (50.0, 0.0)])

    res = completeness.compute_completeness(["calm.nc", "windy.nc"], region, reference_t_hr=25.0)
    pct = np.array(res["completeness_pct"], dtype=float)
    assert res["n_overpasses"] == 2
    assert np.nanmax(pct) == pytest.approx(50.0)
    assert math.isfinite(res["floor_t_hr_median"])
