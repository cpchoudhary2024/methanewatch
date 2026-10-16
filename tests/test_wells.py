import json
from unittest.mock import MagicMock

import pytest
import requests

from methane import validate
from methane import wells as attr
from methane.regions import get_region

REGION = get_region("delaware")


def test_load_wellheads_missing_file_returns_none(tmp_path):
    assert attr.load_wellheads(str(tmp_path / "nope.csv"), REGION) is None


def test_load_wellheads_filters_to_bbox_and_skips_bad_rows(tmp_path):
    csv = tmp_path / "wells.csv"
    csv.write_text(
        "longitude,latitude,operator\n-104.0,32.0,Inside\n-50.0,10.0,Outside\nx,32.1,Bad\n"
    )
    pts = attr.load_wellheads(str(csv), REGION)
    assert [p["name"] for p in pts] == ["Inside"]


def test_load_wellheads_rejects_missing_columns(tmp_path):
    csv = tmp_path / "bad.csv"
    csv.write_text("foo,bar\n1,2\n")
    with pytest.raises(ValueError, match="longitude/latitude"):
        attr.load_wellheads(str(csv), REGION)


def _overpass(elements):
    resp = MagicMock()
    resp.raise_for_status.return_value = None
    resp.json.return_value = {"elements": elements}
    return resp


def test_fetch_wellheads_osm_parses_nodes(mocker):
    els = [{"lat": 32.0, "lon": -104.0}, {"lat": 32.1, "lon": -104.1}, {"id": 3}]
    post = mocker.patch("requests.post", return_value=_overpass(els))
    pts = attr.fetch_wellheads_osm(REGION, mirrors=["https://a"])
    assert len(pts) == 2
    assert "User-Agent" in post.call_args.kwargs["headers"]


def test_fetch_wellheads_osm_tries_next_mirror(mocker):
    mocker.patch(
        "requests.post",
        side_effect=[requests.ConnectionError("down"), _overpass([{"lat": 32.0, "lon": -104.0}])],
    )
    assert len(attr.fetch_wellheads_osm(REGION, mirrors=["https://a", "https://b"])) == 1


def test_fetch_wellheads_osm_returns_none_when_all_fail(mocker):
    mocker.patch("requests.post", side_effect=requests.ConnectionError("down"))
    assert attr.fetch_wellheads_osm(REGION, mirrors=["https://a"]) is None


def test_get_wellheads_prefers_local_csv(tmp_path, mocker):
    csv = tmp_path / "w.csv"
    csv.write_text("lon,lat\n-104.0,32.0\n")
    fetch = mocker.patch.object(attr, "fetch_wellheads_osm")
    pts, source = attr.get_wellheads(REGION, str(csv), str(tmp_path / "c.json"))
    assert len(pts) == 1 and source.startswith("local")
    fetch.assert_not_called()


def test_get_wellheads_uses_cache_then_writes_fresh_fetch(tmp_path, mocker):
    cache = tmp_path / "c.json"
    cache.write_text(
        json.dumps({"region": "delaware", "retrieved": "x", "points": [{"lon": 1, "lat": 2}]})
    )
    fetch = mocker.patch.object(attr, "fetch_wellheads_osm", return_value=[{"lon": 3, "lat": 4}])
    pts, _ = attr.get_wellheads(REGION, str(tmp_path / "none.csv"), str(cache))
    assert pts == [{"lon": 1, "lat": 2}]
    fetch.assert_not_called()

    cache.unlink()
    pts, _ = attr.get_wellheads(REGION, str(tmp_path / "none.csv"), str(cache))
    assert pts == [{"lon": 3, "lat": 4}] and cache.exists()


def _det(lon=-104.0, lat=32.0, q=10000.0):
    return {"detection": {"lon": lon, "lat": lat}, "emission": {"q_kg_hr": q}}


def test_carbon_mapper_sources_collapse_repeat_plumes():
    pts = [
        {"lon": -104.0, "lat": 32.0, "emission_kg_hr": 100.0},
        {"lon": -104.0001, "lat": 32.0001, "emission_kg_hr": 300.0},
        {"lon": -103.0, "lat": 32.0, "emission_kg_hr": float("nan")},
    ]
    src = sorted(validate.carbon_mapper_sources(pts), key=lambda s: s["lon"])
    assert len(src) == 2
    assert src[0]["n_plumes"] == 2 and src[0]["mean_kg_hr"] == pytest.approx(200.0)
    assert src[1]["mean_kg_hr"] is None


def test_carbon_mapper_cross_check_sums_sources_in_radius():
    cm = [
        {"lon": -104.0, "lat": 32.01, "emission_kg_hr": 1000.0},
        {"lon": -104.0, "lat": 32.01, "emission_kg_hr": 3000.0},
        {"lon": -104.0, "lat": 32.05, "emission_kg_hr": 2000.0},
        {"lon": -102.0, "lat": 32.0, "emission_kg_hr": 9999.0},
    ]
    res = validate.carbon_mapper_cross_check([_det(q=8000.0)], cm)[0]["carbon_mapper"]
    assert res["n_sources"] == 2
    assert res["summed_kg_hr"] == pytest.approx(4000.0)
    assert res["ratio"] == pytest.approx(2.0)


def test_carbon_mapper_cross_check_without_nearby_sources():
    cm = [{"lon": -102.0, "lat": 32.0, "emission_kg_hr": 500.0}]
    res = validate.carbon_mapper_cross_check([_det()], cm)[0]["carbon_mapper"]
    assert res["n_sources"] == 0 and res["nearest_km"] > 15
