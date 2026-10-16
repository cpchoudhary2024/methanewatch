from datetime import datetime, timezone
from unittest.mock import MagicMock

import pytest
import requests

from methane.catalog import Granule, search_granules, summarize, _dedupe_prefer_offl
from methane.regions import get_region


def _fake_product(name, pid="00000000-0000-0000-0000-000000000000", size=123456789, online=True):
    return {
        "Id": pid,
        "Name": name,
        "ContentDate": {"Start": "2024-06-01T18:57:46.000Z"},
        "ContentLength": size,
        "Online": online,
    }


def test_granule_name_parsing_extracts_mode_and_orbit():
    g = Granule(
        product_id="x",
        name="S5P_OFFL_L2__CH4____20240601T185746_20240601T203916_34382_03_020600_20240608T044218.nc",
        start=datetime(2024, 6, 1, tzinfo=timezone.utc),
        size_bytes=10_800_000,
        online=True,
    )
    assert g.mode == "OFFL"
    assert g.orbit == "34382"
    assert g.size_mb == pytest.approx(10.8, abs=0.05)


def test_granule_name_parsing_handles_short_malformed_name():
    # A name with fewer underscore-tokens than expected must not raise,
    # mode/orbit degrade to "?" rather than an IndexError.
    g = Granule(
        product_id="x",
        name="onlyonetoken",
        start=datetime.now(timezone.utc),
        size_bytes=0,
        online=True,
    )
    assert g.mode == "?"
    assert g.orbit == "?"


def test_search_granules_parses_and_dedupes(mocker):
    region = get_region("permian")
    payload = {
        "value": [
            _fake_product(
                "S5P_NRTI_L2__CH4____20240601T185746_20240601T203916_34382_03_020600_20240608T044218.nc",
                pid="nrti-1",
            ),
            _fake_product(
                "S5P_OFFL_L2__CH4____20240601T185746_20240601T203916_34382_03_020600_20240608T044218.nc",
                pid="offl-1",
            ),
        ]
    }
    mock_resp = MagicMock()
    mock_resp.json.return_value = payload
    mock_resp.raise_for_status.return_value = None

    mock_session = MagicMock(spec=requests.Session)
    mock_session.get.return_value = mock_resp

    granules = search_granules(region, "2024-06-01", "2024-06-02", session=mock_session)

    # same orbit, OFFL + NRTI both returned by the (fake) API -> OFFL wins
    assert len(granules) == 1
    assert granules[0].mode == "OFFL"
    assert granules[0].product_id == "offl-1"
    mock_session.get.assert_called_once()
    called_kwargs = mock_session.get.call_args.kwargs
    assert "params" in called_kwargs
    assert called_kwargs["timeout"] == 60


def test_search_granules_empty_response_returns_empty_list(mocker):
    """Empty dataset response: CDSE returns zero granules for the window,
    a normal outcome (e.g. no overpass), not an error."""
    region = get_region("permian")
    mock_resp = MagicMock()
    mock_resp.json.return_value = {"value": []}
    mock_resp.raise_for_status.return_value = None
    mock_session = MagicMock(spec=requests.Session)
    mock_session.get.return_value = mock_resp

    granules = search_granules(region, "2024-06-01", "2024-06-02", session=mock_session)
    assert granules == []


def test_search_granules_propagates_http_error(mocker):
    """A live CDSE outage / 5xx must surface, not be silently swallowed."""
    region = get_region("permian")
    mock_resp = MagicMock()
    mock_resp.raise_for_status.side_effect = requests.exceptions.HTTPError("503 Server Error")
    mock_session = MagicMock(spec=requests.Session)
    mock_session.get.return_value = mock_resp

    with pytest.raises(requests.exceptions.HTTPError):
        search_granules(region, "2024-06-01", "2024-06-02", session=mock_session)


def test_search_granules_sudden_network_dropout(mocker):
    """Simulates a mid-request connection drop (not an HTTP error status,
    the TCP connection itself fails)."""
    region = get_region("permian")
    mock_session = MagicMock(spec=requests.Session)
    mock_session.get.side_effect = requests.exceptions.ConnectionError("Connection reset by peer")

    with pytest.raises(requests.exceptions.ConnectionError):
        search_granules(region, "2024-06-01", "2024-06-02", session=mock_session)


def test_dedupe_prefers_offl_over_rpro_over_nrti():
    mk = lambda name, orbit_suffix: Granule(  # noqa: E731
        product_id=name,
        name=f"S5P_{name}_L2__CH4____20240601T000000_20240601T000000_{orbit_suffix}_03_020600_20240608T000000.nc",
        start=datetime(2024, 6, 1, tzinfo=timezone.utc),
        size_bytes=1,
        online=True,
    )
    granules = [mk("NRTI", "11111"), mk("RPRO", "11111"), mk("OFFL", "11111")]
    result = _dedupe_prefer_offl(granules)
    assert len(result) == 1
    assert result[0].mode == "OFFL"


def test_summarize_empty_dataset():
    assert summarize([]) == "No granules found for that region/time window."


def test_summarize_flags_offline_granule():
    g = Granule(
        product_id="x",
        name="S5P_OFFL_L2__CH4____20240601T185746_20240601T203916_34382_03_020600_20240608T044218.nc",
        start=datetime(2024, 6, 1, tzinfo=timezone.utc),
        size_bytes=10_800_000,
        online=False,
    )
    out = summarize([g])
    assert "offline" in out
    assert "1 granule(s)" in out
