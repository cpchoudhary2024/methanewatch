from unittest.mock import MagicMock

import numpy as np

from methane.read import QA_THRESHOLD, read_scene


class _FakeVar:
    """Mimics a netCDF4 Variable: subscriptable, returns a numpy array."""

    def __init__(self, data):
        self._data = np.asarray(data)

    def __getitem__(self, item):
        return self._data


class _FakeGroup:
    def __init__(self, variables=None, groups=None):
        self.variables = variables or {}
        self.groups = groups or {}


def _fake_dataset(n=6, with_precision=True):
    lon = np.array([-104.0, -104.1, -104.2, -103.9, -50.0, -104.05][:n])
    lat = np.array([32.0, 32.1, 32.2, 31.9, 10.0, 32.05][:n])
    xch4 = np.array([1900.0, 1950.0, np.nan, 1920.0, 1900.0, 1980.0][:n])
    qa = np.array([0.9, 0.3, 0.9, 0.6, 0.9, 0.55][:n])  # index 1 fails QA (0.3<0.5)

    prod_vars = {
        "longitude": _FakeVar(lon),
        "latitude": _FakeVar(lat),
        "methane_mixing_ratio_bias_corrected": _FakeVar(xch4),
        "qa_value": _FakeVar(qa),
        "time_utc": _FakeVar(np.array(["2024-06-01T19:19:20Z"])),
    }
    detailed = _FakeGroup(
        variables=(
            {"methane_mixing_ratio_precision": _FakeVar(np.full(n, 5.0))} if with_precision else {}
        )
    )
    prod = _FakeGroup(
        variables=prod_vars,
        groups={"SUPPORT_DATA": _FakeGroup(groups={"DETAILED_RESULTS": detailed})},
    )
    ds = MagicMock()
    ds.groups = {"PRODUCT": prod}
    ds.__enter__.return_value = ds
    ds.__exit__.return_value = False
    return ds


def test_read_scene_applies_qa_filter_and_bbox(mocker, tmp_path):
    fake_ds = _fake_dataset()
    mocker.patch("methane.read.Dataset", return_value=fake_ds)

    p = tmp_path / "fake_granule.nc"
    p.write_bytes(b"not a real netcdf file")  # Dataset() is mocked, content unused

    bbox = (-105.0, 30.8, -101.0, 33.5)  # Permian
    scene = read_scene(p, bbox, qa_threshold=QA_THRESHOLD)

    # index 1 fails QA, index 2 has NaN xch4, index 4 is outside bbox
    # (lon -50, lat 10) -> only indices 0, 3, 5 should survive.
    assert scene.n == 3
    assert scene.n_total == 4  # in-box with finite xch4: 0,1,3,5 (2 is NaN xch4)
    assert set(np.round(scene.lon, 2)) == {-104.0, -103.9, -104.05}
    assert scene.source == "fake_granule.nc"


def test_read_scene_missing_precision_variable_fills_nan(mocker, tmp_path):
    fake_ds = _fake_dataset(with_precision=False)
    mocker.patch("methane.read.Dataset", return_value=fake_ds)
    p = tmp_path / "no_precision.nc"
    p.write_bytes(b"stub")

    scene = read_scene(p, (-105.0, 30.8, -101.0, 33.5))
    assert np.isnan(scene.precision).all()


def test_read_scene_empty_dataset_response():
    """All pixels fail QA / are outside the bbox -> an empty Scene, not a
    crash. This is the 'insufficient coverage' case detect.py already
    handles gracefully downstream."""
    import methane.read as read_mod

    n = 4
    lon = np.array([-50.0, -51.0, -52.0, -53.0])  # all outside Permian
    lat = np.array([10.0, 11.0, 12.0, 13.0])
    xch4 = np.full(n, 1900.0)
    qa = np.full(n, 0.9)
    prod = _FakeGroup(
        variables={
            "longitude": _FakeVar(lon),
            "latitude": _FakeVar(lat),
            "methane_mixing_ratio_bias_corrected": _FakeVar(xch4),
            "qa_value": _FakeVar(qa),
            "time_utc": _FakeVar(np.array(["2024-06-01T19:19:20Z"])),
        }
    )
    ds = MagicMock()
    ds.groups = {"PRODUCT": prod}
    ds.__enter__.return_value = ds
    ds.__exit__.return_value = False

    orig = read_mod.Dataset
    read_mod.Dataset = MagicMock(return_value=ds)
    try:
        scene = read_scene("unused_path.nc", (-105.0, 30.8, -101.0, 33.5))
    finally:
        read_mod.Dataset = orig

    assert scene.n == 0
    assert scene.n_total == 0
    assert scene.qa_pass_fraction == 0.0


def test_read_scene_without_time_variable_marks_time_unknown(mocker, tmp_path):
    prod = _FakeGroup(
        variables={
            "longitude": _FakeVar([-104.0, -104.1]),
            "latitude": _FakeVar([32.0, 32.1]),
            "methane_mixing_ratio_bias_corrected": _FakeVar([1900.0, 1910.0]),
            "qa_value": _FakeVar([0.9, 0.9]),
        }
    )
    ds = MagicMock()
    ds.groups = {"PRODUCT": prod}
    ds.__enter__.return_value = ds
    ds.__exit__.return_value = False
    mocker.patch("methane.read.Dataset", return_value=ds)

    scene = read_scene(tmp_path / "x.nc", (-105.0, 30.8, -101.0, 33.5))
    assert scene.time_utc == "unknown"
    assert scene.n == 2
