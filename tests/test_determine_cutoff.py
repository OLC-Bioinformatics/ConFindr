#!/usr/bin/env python3

"""
Tests for determine_cutoff in confindr_src.methods
"""

# Local imports
from confindr_src.methods import determine_cutoff


def test_determine_cutoff_empty_qualities_returns_one():
    """
    Test that determine_cutoff returns a cutoff of 1 and zero error
    when given an empty qualities list.
    """
    k, exp, err = determine_cutoff(
        qualities=[],
        reference_sequence='A' * 100,
        base_cutoff=0,
        error_cutoff=1.0
    )
    assert k == 1
    assert exp == 0.0
    assert err == 0.0


def test_determine_cutoff_dynamic_with_qualities_computes_reasonable_cutoff():
    """
    Test that determine_cutoff with a moderate qualities list computes a
    reasonable cutoff and returns expected positions and error percentage.
    """
    # Use a moderate number of high-quality bases; since we use a per-
    # position approximation, the computed cutoff should be 1 when average
    # per-position depth is small.
    qualities = [30] * 50  # Q30 -> p ~ 1e-3
    k, exp, err = determine_cutoff(
        qualities=qualities,
        reference_sequence='A' * 500,
        base_cutoff=0,
        error_cutoff=1.0
    )
    assert isinstance(k, int)
    assert k >= 1

    # Expected positions should be small due to low average depth
    assert isinstance(exp, float)
    assert exp == 0.5

    # With low average per-position depth, the dynamic cutoff should be 1
    assert k == 1
    assert isinstance(err, float)
    assert err >= 0.0


def test_determine_cutoff_high_depth_returns_larger_cutoff():
    """
    Test that determine_cutoff with very high-depth qualities computes a
    larger cutoff and returns expected positions and error percentage.
    """
    # Simulate high average per-position depth (e.g. 50x across 1000bp)
    qualities = [30] * 50000
    k, exp, err = determine_cutoff(
        qualities=qualities,
        reference_sequence='A' * 1000,
        base_cutoff=0,
        error_cutoff=1.0
    )
    assert isinstance(k, int) and k >= 1
    assert k > 1
    assert isinstance(exp, float) and exp >= 0.0
    assert isinstance(err, float) and err >= 0.0
