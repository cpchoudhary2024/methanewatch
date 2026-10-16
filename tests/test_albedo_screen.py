from unittest.mock import MagicMock

import numpy as np
import pytest

from methane.albedo_screen import (
    MIN_CELLS,
    MISALIGN_DEG,
    R_THRESHOLD,
    albedo_gradient,
    analyze_candidate,
    axis_misalignment_deg,
    pearson_with_p,
    principal_axis_deg,
)


def test_gradient_of_eastward_ramp_points_east():
    a = np.tile(np.linspace(0, 1, 10), (10, 1))
    mag, direction = albedo_gradient(a)
    assert np.nanmean(mag) > 0
    assert np.nanmedian(direction) == pytest.approx(90.0, abs=1.0)


def test_gradient_of_northward_ramp_points_north():
    a = np.tile(np.linspace(0, 1, 10).reshape(-1, 1), (1, 10))
    _mag, direction = albedo_gradient(a)
    assert np.nanmedian(direction) == pytest.approx(0.0, abs=1.0)


def test_gradient_of_flat_field_is_zero():
    mag, _ = albedo_gradient(np.full((8, 8), 0.2))
    assert np.allclose(mag, 0.0)


def test_gradient_masks_invalid_cells_but_does_not_propagate_nan():
    a = np.tile(np.linspace(0, 1, 10), (10, 1))
    a[5, 5] = np.nan
    mag, _ = albedo_gradient(a)
    assert np.isnan(mag[5, 5])
    assert np.isfinite(mag[5, 3]) and np.isfinite(mag[5, 7])


def test_gradient_all_nan_returns_all_nan():
    mag, direction = albedo_gradient(np.full((4, 4), np.nan))
    assert np.isnan(mag).all() and np.isnan(direction).all()


def test_east_west_blob_axis_is_90():
    f = np.zeros((11, 11))
    f[5, 2:9] = 1.0
    assert principal_axis_deg(f, 32.0) == pytest.approx(90.0, abs=1.0)


def test_north_south_blob_axis_is_zero():
    f = np.zeros((11, 11))
    f[2:9, 5] = 1.0
    axis = principal_axis_deg(f, 32.0)
    assert min(axis, 180.0 - axis) == pytest.approx(0.0, abs=1.0)


def test_principal_axis_none_when_no_positive_signal():
    assert principal_axis_deg(np.zeros((6, 6)), 32.0) is None
    assert principal_axis_deg(np.full((6, 6), -1.0), 32.0) is None


def test_misalignment_is_zero_when_aligned():
    assert axis_misalignment_deg(90.0, 90.0) == pytest.approx(0.0)


def test_misalignment_is_90_when_perpendicular():
    assert axis_misalignment_deg(0.0, 90.0) == pytest.approx(90.0)


def test_antiparallel_wind_is_the_same_axis_not_a_conflict():
    assert axis_misalignment_deg(90.0, 270.0) == pytest.approx(0.0)


def test_misalignment_none_without_wind():
    assert axis_misalignment_deg(90.0, None) is None
    assert axis_misalignment_deg(None, 90.0) is None


def test_pearson_detects_perfect_correlation():
    a = np.arange(50, dtype=float)
    r, p, n = pearson_with_p(a, 2 * a + 1)
    assert r == pytest.approx(1.0)
    assert n == 50


def test_pearson_none_on_constant_input():
    r, p, n = pearson_with_p(np.ones(30), np.arange(30, dtype=float))
    assert r is None and p is None


def test_pearson_ignores_nan_pairs():
    a = np.arange(20, dtype=float)
    b = a.copy()
    b[:5] = np.nan
    r, _p, n = pearson_with_p(a, b)
    assert n == 15 and r == pytest.approx(1.0)


class _FakeVar:
    def __init__(self, data):
        self._data = np.ma.masked_invalid(np.asarray(data, dtype=float))

    def __getitem__(self, item):
        return self._data


class _FakeGroup:
    def __init__(self, variables=None, groups=None):
        self.variables = variables or {}
        self.groups = groups or {}


def _fake_granule(clon=-104.0, clat=32.0, n=40, *, artifact=False, span=0.5):
    """Synthetic patch with a N-S albedo edge, albedo = 0.2 + 0.1 tanh(dlon / a).

    With artifact=True methane follows the gradient magnitude (sech^2) of that edge;
    otherwise it is a compact blob away from the edge.
    """
    a = 0.06
    lin = np.linspace(-span, span, n)
    dlon, dlat = np.meshgrid(lin, lin)
    lon = clon + dlon
    lat = clat + dlat

    albedo = 0.2 + 0.1 * np.tanh(dlon / a)

    if artifact:
        xch4 = 1850.0 + 100.0 / np.cosh(dlon / a) ** 2
    else:
        xch4 = 1850.0 + 60.0 * np.exp(-((dlon - 0.3) ** 2 + dlat**2) / 0.004)

    qa = np.full((n, n), 0.9)
    prod = _FakeGroup(
        variables={
            "longitude": _FakeVar(lon),
            "latitude": _FakeVar(lat),
            "methane_mixing_ratio_bias_corrected": _FakeVar(xch4),
            "qa_value": _FakeVar(qa),
        },
        groups={
            "SUPPORT_DATA": _FakeGroup(
                groups={
                    "DETAILED_RESULTS": _FakeGroup(
                        variables={"surface_albedo_SWIR": _FakeVar(albedo)}
                    )
                }
            )
        },
    )
    ds = MagicMock()
    ds.groups = {"PRODUCT": prod}
    ds.__enter__.return_value = ds
    ds.__exit__.return_value = False
    return ds


def test_real_plume_is_not_flagged(mocker):
    mocker.patch("netCDF4.Dataset", return_value=_fake_granule(artifact=False))
    d = analyze_candidate("fake.nc", -104.0, 32.0, 90.0)
    assert d["available"] is True
    assert d["albedo_gradient_r"] < R_THRESHOLD
    assert d["flagged"] is False


def test_albedo_tracking_and_cross_wind_is_flagged(mocker):
    mocker.patch("netCDF4.Dataset", return_value=_fake_granule(artifact=True))
    d = analyze_candidate("fake.nc", -104.0, 32.0, 90.0)
    assert d["available"] is True
    assert d["albedo_gradient_r"] > R_THRESHOLD
    assert d["axis_misalignment_deg"] > MISALIGN_DEG
    assert d["flagged"] is True


def test_albedo_tracking_but_wind_aligned_is_retained(mocker):
    mocker.patch("netCDF4.Dataset", return_value=_fake_granule(artifact=True))
    d = analyze_candidate("fake.nc", -104.0, 32.0, 0.0)
    assert d["albedo_gradient_r"] > R_THRESHOLD
    assert d["albedo_flag"] is True
    assert d["flagged"] is False


def test_unreadable_granule_is_graceful(mocker):
    mocker.patch("netCDF4.Dataset", side_effect=OSError("no such file"))
    d = analyze_candidate("missing.nc", -104.0, 32.0, 90.0)
    assert d["available"] is False and "note" in d


def test_sparse_patch_issues_no_verdict(mocker):
    mocker.patch("netCDF4.Dataset", return_value=_fake_granule(artifact=True, n=5, span=0.1))
    d = analyze_candidate("fake.nc", -104.0, 32.0, 90.0)
    assert d["available"] is False
    assert d["n_cells"] < MIN_CELLS
    assert "flagged" not in d


def test_missing_wind_direction_cannot_flag(mocker):
    mocker.patch("netCDF4.Dataset", return_value=_fake_granule(artifact=True))
    d = analyze_candidate("fake.nc", -104.0, 32.0, None)
    assert d["axis_misalignment_deg"] is None
    assert d["wind_flag"] is False
    assert d["flagged"] is False
