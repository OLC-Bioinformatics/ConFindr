#!/usr/bin/env python3

"""
Tests for determine_cutoff tightening behavior in confindr_src.methods
"""

# Local imports
from confindr_src.methods import determine_cutoff


def test_tightening_respects_max_expected_positions():
    """
    Test that determine_cutoff tightens the cutoff to respect
    max_expected_positions when given very high-depth data.
    """
    # High-depth simulated qualities
    qualities = [30] * 50000

    # Short reference length to make per-pos depth large
    k, expected_positions, per_site_percent = determine_cutoff(
        qualities=qualities,
        reference_sequence='A' * 1000,
        base_cutoff=0,
        error_cutoff=1.0,
        max_expected_positions=0.001
    )

    # Expect the returned expected_positions to be <= threshold
    assert expected_positions <= 0.001
    assert isinstance(k, int)
    assert isinstance(per_site_percent, float)
