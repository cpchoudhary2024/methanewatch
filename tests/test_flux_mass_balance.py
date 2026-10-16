import math

import pytest

from methane.detect import Candidate
from methane.flux import G, M_AIR, M_CH4, UEFF_REL_UNC, column_excess_kg_per_m2, estimate_emission


def _candidate(delta_ppb=100.0, area_km2=25.0):
    return Candidate(-102.0, 31.9, delta_ppb, delta_ppb * 1.5, 10, area_km2, 5.0)


def test_dry_air_column_at_one_bar():
    # 1e5 / (9.80665 * 0.0289644)
    assert 1e5 / (G * M_AIR) == pytest.approx(3.5206e5, rel=1e-4)


def test_column_excess_hand_computed():
    # 100e-9 * 3.52064e5 mol m-2 * 0.0160425 kg mol-1
    assert column_excess_kg_per_m2(100.0, 1e5) == pytest.approx(5.6478e-4, rel=1e-4)
    assert column_excess_kg_per_m2(0.0, 1e5) == 0.0


def test_permian_surface_pressure_lowers_the_conversion():
    assert column_excess_kg_per_m2(100.0, 9.0e4) == pytest.approx(
        0.9 * column_excess_kg_per_m2(100.0, 1.0e5)
    )


def test_ime_hand_computed():
    # A = 2.5e7 m2, L = 5000 m, IME = 5.6478e-4 * 2.5e7 = 1.41194e4 kg,
    # Q = 5 * 1.41194e4 / 5000 * 3600 = 5.0830e4 kg/h
    est = estimate_emission(_candidate(100.0, 25.0), wind_speed_ms=5.0, p_surf=1e5)
    assert est.length_scale_m == pytest.approx(5000.0)
    assert est.ime_kg == pytest.approx(1.41194e4, rel=1e-4)
    assert est.q_kg_hr == pytest.approx(5.0830e4, rel=1e-4)


def test_rate_scales_as_sqrt_area_at_fixed_enhancement():
    small = estimate_emission(_candidate(area_km2=25.0), 5.0).q_kg_hr
    large = estimate_emission(_candidate(area_km2=100.0), 5.0).q_kg_hr
    assert large == pytest.approx(2.0 * small)


def test_uncertainty_is_quadrature_sum():
    est = estimate_emission(_candidate(), wind_speed_ms=5.0)
    assert est.rel_uncertainty == pytest.approx(math.sqrt(0.50**2 + 0.15**2 + 0.25**2))


def test_wind_uncertainty_is_floored_at_ueff_spread():
    tight = estimate_emission(_candidate(), 5.0, wind_speed_unc_ms=0.01).rel_uncertainty
    loose = estimate_emission(_candidate(), 5.0, wind_speed_unc_ms=5.0).rel_uncertainty
    default = estimate_emission(_candidate(), 5.0).rel_uncertainty
    assert UEFF_REL_UNC == 0.50
    assert tight == pytest.approx(default)
    assert loose > default


def test_lower_bound_is_not_negative():
    est = estimate_emission(_candidate(), wind_speed_ms=1.0, wind_speed_unc_ms=100.0)
    assert est.q_kg_hr_lo == 0.0
    assert est.q_kg_hr_lo <= est.q_kg_hr <= est.q_kg_hr_hi


def test_constants():
    assert (G, M_AIR, M_CH4) == (9.80665, 0.0289644, 0.0160425)
