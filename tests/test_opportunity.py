import numpy as np
from scipy.spatial import cKDTree

from methane.opportunity import _xyz, paired_enrichment_ci, pixel_weights, within_km


def test_pixel_opportunity_weights_change_reference_fraction():
    surface = (
        np.array([0.0, 1.0, 2.0]),
        np.array([0.0, 1.0]),
        np.array([[1.0], [3.0]]),
    )
    weights = pixel_weights(np.array([0.5, 1.5]), np.array([0.5, 0.5]), surface)
    assert np.allclose(weights, [0.25, 0.75])
    assert np.dot([1.0, 0.0], weights) == 0.25


def test_great_circle_threshold_and_day_paired_bootstrap():
    tree = cKDTree(_xyz([0], [0]))
    assert within_km(tree, [0, 1], [0, 0], 5).tolist() == [True, False]
    rows = [
        {"day": "2024-06-01", "near_flare": True, "uniform_fraction": 0.5},
        {"day": "2024-06-01", "near_flare": False, "uniform_fraction": 0.5},
        {"day": "2024-06-02", "near_flare": True, "uniform_fraction": 0.5},
    ]
    low, high = paired_enrichment_ci(rows, "uniform_fraction", n_boot=200, seed=1)
    assert low <= 1.33 <= high
