#!/usr/bin/env python3

"""
Tests for Benjamini-Hochberg and Fisher's method implementations in
confindr_src.methods
"""

# Third-party imports
from scipy.stats import combine_pvalues

# Local imports
from confindr_src.methods import (
    benjamini_hochberg,
    combine_pvalues_fisher
)


def test_bh_returns_monotonic_adjusted_qvalues():
    """
    Test that benjamini_hochberg returns q-values that are monotonic and
    between 0 and 1.
    """
    # Test input p-values
    pvals = [0.02, 0.5, 0.001, 0.04, 0.03]

    # Get q-values
    qvals = benjamini_hochberg(pvals=pvals)
    assert len(qvals) == len(pvals)

    # q-values ordered by original p should be non-decreasing
    order = sorted(range(len(pvals)), key=lambda i: pvals[i])

    # Check monotonicity
    q_in_order = [qvals[i] for i in order]
    assert q_in_order == sorted(q_in_order)

    # All q-values should be between 0 and 1
    assert all(0.0 <= q <= 1.0 for q in qvals)


def test_fisher_combination_matches_scipy():
    """
    Test that combine_pvalues_fisher matches SciPy's implementation.
    """
    # Test input p-values
    pvals = [0.01, 0.02, 0.05]
    combined = combine_pvalues_fisher(pvals=pvals)

    # Get expected value from SciPy
    expected = combine_pvalues(pvals, method='fisher')[1]

    assert abs(combined - expected) < 1e-9
