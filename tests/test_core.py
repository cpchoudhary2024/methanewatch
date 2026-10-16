import math

import numpy as np
import pandas as pd
import pytest

from methane import gap
from methane.cli import _cmd_run, granules_in_window
from methane.detect import Candidate, detect
from methane.flaring import detection_proximity
from methane.flux import column_excess_kg_per_m2, estimate_emission
from methane.read import Scene
from methane.regions import get_region
from methane.validate import ReferenceDetection, build_validation, nearest_km


def _scene(lon, lat, xch4, p_surf=91000.0, time_utc="2024-06-05T19:30:00Z"):
    lon, lat, xch4 = (np.asarray(v, dtype=float).ravel() for v in (lon, lat, xch4))
    return Scene(
        lon=lon,
        lat=lat,
        xch4=xch4,
        precision=np.full(lon.size, 8.0),
        qa=np.full(lon.size, 0.8),
        p_surf=np.full(lon.size, p_surf),
        source="synthetic.nc",
        n_total=lon.size,
        time_utc=time_utc,
    )


def _plume_scene(n_side=120, plume_ppb=80.0):
    lon, lat = np.meshgrid(np.linspace(-104, -103, n_side), np.linspace(31.5, 32.5, n_side))
    rng = np.random.default_rng(0)
    field = 1900.0 + rng.normal(0, 3.0, lon.shape)
    field += plume_ppb * np.exp(-(((lon + 103.5) / 0.06) ** 2 + ((lat - 32.0) / 0.06) ** 2))
    return _scene(lon, lat, field)


def test_detect_finds_planted_plume():
    res = detect(_plume_scene(), k_sigma=3.0, min_pixels=3)
    top = res.candidates[0]
    assert abs(top.lon + 103.5) < 0.1 and abs(top.lat - 32.0) < 0.1
    assert top.confidence_sigma > 3.0


def test_detect_empty_when_no_plume():
    assert detect(_plume_scene(plume_ppb=0.0), k_sigma=4.0, min_pixels=4).candidates == []


def test_column_conversion_is_linear_in_enhancement_and_pressure():
    a = column_excess_kg_per_m2(10.0, 90000.0)
    assert column_excess_kg_per_m2(20.0, 90000.0) == pytest.approx(2 * a)
    assert column_excess_kg_per_m2(10.0, 45000.0) == pytest.approx(a / 2)


def test_emission_scales_with_wind():
    c = Candidate(-103.5, 32.0, 30.0, 70.0, 6, 150.0, 5.0)
    slow, fast = estimate_emission(c, 2.0), estimate_emission(c, 4.0)
    assert fast.q_kg_hr == pytest.approx(2 * slow.q_kg_hr)
    assert slow.q_kg_hr_lo < slow.q_kg_hr < slow.q_kg_hr_hi


def test_basin_code_parsing():
    assert gap._basin_code("ConocoPhillips Permian (430)") == "430"
    assert gap._basin_code("Endeavor Energy Resources LP 430 Permian Basin") == "430"
    assert gap._basin_code("ConocoPhillips' Williston (395)") == "395"
    assert gap._basin_code("Old Ocean Gas Plant") is None


@pytest.mark.parametrize("wind_from", [0.0, 90.0, 225.0])
def test_wind_frame_puts_downwind_positive(wind_from):
    t = math.radians(wind_from)
    downwind = (-math.sin(t) * 10, -math.cos(t) * 10)
    along, cross = gap._wind_frame(np.array([downwind[0]]), np.array([downwind[1]]), wind_from)
    assert along[0] == pytest.approx(10.0)
    assert cross[0] == pytest.approx(0.0, abs=1e-9)


def test_detrend_removes_a_pressure_gradient():
    rng = np.random.default_rng(1)
    p = rng.uniform(85000, 98000, 2000)
    lat = rng.uniform(30, 34, 2000)
    xch4 = 1880 + 0.05 * (p / 100) + 2.0 * lat
    s = _scene(np.zeros(2000), lat, xch4, p_surf=0.0)
    s.p_surf = p
    out = gap._detrend(s, np.ones(2000, dtype=bool))
    assert np.std(out) < 1e-6


def test_scene_balance_recovers_a_uniform_enhancement(mocker):
    region = get_region("permian")
    lon, lat = np.meshgrid(np.arange(-106, -100, 0.05), np.arange(29.8, 34.5, 0.05))
    inside = (lon >= -105) & (lon <= -101) & (lat >= 30.8) & (lat <= 33.5)
    xch4 = np.where(inside, 1910.0, 1900.0)
    mocker.patch.object(gap, "read_scene", return_value=_scene(lon, lat, xch4))
    mocker.patch.object(gap.wind, "area_mean", side_effect=[(5.0, 180.0), (3.0, 180.0)])

    ob = gap.scene_balance("x.nc", region)
    assert ob.delta_ppb == pytest.approx(10.0)
    expected = 5.0 * ob.width_km * 1e3 * column_excess_kg_per_m2(10.0, 91000.0)
    assert ob.q_tg_yr == pytest.approx(expected * gap.SECONDS_PER_YEAR / 1e9, rel=1e-3)
    # southerly wind: the cross-wind width is the east-west extent of the box
    assert ob.width_km == pytest.approx(4 * 111.32 * math.cos(math.radians(32.15)), rel=1e-3)


def test_granules_in_window_filters_by_date_and_product(tmp_path):
    names = [
        "S5P_OFFL_L2__CH4____20240531T190000_20240531T200000_34368_03_020600_x.nc",
        "S5P_OFFL_L2__CH4____20240601T190000_20240601T200000_34382_03_020600_x.nc",
        "S5P_OFFL_L2__CH4____20240831T190000_20240831T200000_35700_03_020600_x.nc",
        "S5P_OFFL_L2__CH4____20240901T190000_20240901T200000_35714_03_020600_x.nc",
        "S5P_NRTI_L2__CH4____20240615T190000_20240615T200000_34580_03_020600_x.nc",
    ]
    for n in names:
        (tmp_path / n).touch()
    got = [p.split("/")[-1] for p in granules_in_window(tmp_path, "2024-06-01", "2024-09-01")]
    assert got == names[1:3]


def test_cached_run_reads_requested_cache(tmp_path, mocker):
    cache = tmp_path / "granules"
    cache.mkdir()
    name = "S5P_OFFL_L2__CH4____20240601T190000_20240601T200000_34382_03_020600_x.nc"
    path = cache / name
    path.touch()
    mocker.patch("methane.download.cache_dir", return_value=tmp_path / "wrong-cache")
    process = mocker.patch(
        "methane.cli._process_granule",
        return_value=(
            {"threshold_ppb": 20.0, "n_inverted": 0, "n_candidates": 0},
            [],
        ),
    )
    args = type(
        "Args",
        (),
        {
            "region": "permian",
            "out": str(tmp_path / "results"),
            "cache": str(cache),
            "start": "2024-06-01",
            "end": "2024-09-01",
            "cached": True,
        },
    )()
    assert _cmd_run(args) == 0
    assert process.call_args.args[0] == path


def test_nearest_km_matches_a_known_distance():
    d = nearest_km([-103.0], [32.0], np.array([-103.0]), np.array([33.0]))
    assert d[0] == pytest.approx(111.2, abs=0.2)


def _det(lon, lat):
    return {
        "detection": {"lon": lon, "lat": lat, "confidence_sigma": 4.0},
        "emission": {"q_kg_hr": 20000.0},
        "overpass_utc": "2024-06-05T19:30:00Z",
    }


def test_clustered_null_removes_proximity_that_uniform_null_rewards():
    region = get_region("permian")
    rng = np.random.default_rng(2)
    wells_lon = -103.5 + rng.normal(0, 0.1, 3000)
    wells_lat = 32.0 + rng.normal(0, 0.1, 3000)
    refs = [ReferenceDetection(lon, lat) for lon, lat in zip(wells_lon[:400], wells_lat[:400])]
    dets = [_det(lon, lat) for lon, lat in zip(wells_lon[500:520], wells_lat[500:520])]

    v = build_validation(dets, refs, bbox=region.bbox, infra_lon=wells_lon, infra_lat=wells_lat)
    assert v["uniform_null"]["enrichment"] > 3
    assert v["clustered_null"]["enrichment"] == pytest.approx(1.0, abs=0.15)


def test_flare_proximity_reports_both_nulls():
    region = get_region("permian")
    sites = pd.DataFrame(
        {
            "klon": [-103.5, -102.0, -103.51, -101.99],
            "klat": [32.0, 31.6, 32.01, 31.61],
            "rh": [10, 8, 9, 7],
        }
    )
    dets = [_det(-103.5, 32.0), _det(-102.0, 31.6)]
    res = detection_proximity(dets, sites, region, infra_lon=[], infra_lat=[], top_quantile=0.0)
    assert res["uniform_null"]["mannwhitney_p"] < 0.05
    assert res["clustered_null"]["method"].startswith("uniform")


def test_inverted_field_has_no_planted_plume():
    import dataclasses

    scene = _plume_scene()
    inverted = detect(dataclasses.replace(scene, xch4=-scene.xch4), k_sigma=3.0, min_pixels=3)
    assert all(abs(c.lon + 103.5) > 0.1 or abs(c.lat - 32.0) > 0.1 for c in inverted.candidates)
