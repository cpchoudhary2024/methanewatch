# Copyright (c) 2026 Chandra Prakash Choudhary. Licensed under the MIT License.
"""Point-source emission rate by integrated mass enhancement (Varon et al., 2018).

    Q = U_eff * IME / L,   IME = sum(dOmega * A),   L = sqrt(A)

U_eff is taken as the ERA5 10 m wind with a 50% relative uncertainty standing in
for the U_eff(U10) calibration spread. The column conversion uses the scene
surface pressure: N_air = p_s / (g * M_air).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .detect import Candidate

G = 9.80665  # m s-2
M_AIR = 0.0289644  # kg mol-1
M_CH4 = 0.0160425  # kg mol-1
P_SURF_DEFAULT = 101325.0  # Pa, used only when the product has no surface pressure

UEFF_REL_UNC = 0.50
RETRIEVAL_REL_UNC = 0.15
LENGTH_REL_UNC = 0.25
HR_PER_YR = 24 * 365.25


def column_excess_kg_per_m2(delta_xch4_ppb: float, p_surf: float = P_SURF_DEFAULT) -> float:
    n_air = p_surf / (G * M_AIR)
    return delta_xch4_ppb * 1e-9 * n_air * M_CH4


@dataclass
class EmissionEstimate:
    q_kg_hr: float
    q_kg_hr_lo: float
    q_kg_hr_hi: float
    ime_kg: float
    length_scale_m: float
    wind_speed_ms: float
    wind_direction_deg: float | None
    wind_source: str
    p_surf_pa: float
    rel_uncertainty: float

    def as_dict(self) -> dict:
        return {
            "q_kg_hr": round(self.q_kg_hr, 1),
            "q_kg_hr_range": [round(self.q_kg_hr_lo, 1), round(self.q_kg_hr_hi, 1)],
            "ime_kg": round(self.ime_kg, 1),
            "length_scale_m": round(self.length_scale_m),
            "wind_speed_ms": round(self.wind_speed_ms, 2),
            "wind_direction_deg": (
                None if self.wind_direction_deg is None else round(self.wind_direction_deg, 1)
            ),
            "wind_source": self.wind_source,
            "p_surf_pa": round(self.p_surf_pa),
            "rel_uncertainty_1sigma": round(self.rel_uncertainty, 2),
        }


def estimate_emission(
    cand: Candidate,
    wind_speed_ms: float,
    *,
    wind_speed_unc_ms: float | None = None,
    wind_source: str = "user-supplied",
    wind_direction_deg: float | None = None,
    p_surf: float = P_SURF_DEFAULT,
    retrieval_rel_unc: float = RETRIEVAL_REL_UNC,
) -> EmissionEstimate:
    area_m2 = cand.area_km2 * 1e6
    length_scale = float(np.sqrt(area_m2))
    ime_kg = column_excess_kg_per_m2(cand.delta_xch4_mean, p_surf) * area_m2
    q_kg_hr = wind_speed_ms * ime_kg / length_scale * 3600.0

    wind_rel = UEFF_REL_UNC
    if wind_speed_unc_ms is not None and wind_speed_ms > 0:
        wind_rel = max(wind_speed_unc_ms / wind_speed_ms, UEFF_REL_UNC)
    rel = float(np.sqrt(wind_rel**2 + retrieval_rel_unc**2 + LENGTH_REL_UNC**2))

    return EmissionEstimate(
        q_kg_hr=q_kg_hr,
        q_kg_hr_lo=q_kg_hr * max(0.0, 1 - rel),
        q_kg_hr_hi=q_kg_hr * (1 + rel),
        ime_kg=ime_kg,
        length_scale_m=length_scale,
        wind_speed_ms=wind_speed_ms,
        wind_direction_deg=wind_direction_deg,
        wind_source=wind_source,
        p_surf_pa=p_surf,
        rel_uncertainty=rel,
    )
