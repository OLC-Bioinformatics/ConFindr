#!/usr/bin/env python3
"""
Tests for Mann-Whitney U test with ties in confindr_src.methods
"""

# Third-party imports
from scipy.stats import mannwhitneyu

# Local imports
from confindr_src.methods import mann_whitney_u_p


def test_mann_whitney_with_ties_matches_scipy():
    """
    Test that mann_whitney_u_p handles ties correctly and matches SciPy's
    implementation.
    """
    # Sample data with ties
    x = [1, 1, 1, 2, 2]
    y = [1, 2, 3, 3, 3]

    # Get p-value from our implementation
    p_ours = mann_whitney_u_p(x=x, y=y)

    # Compare with SciPy's implementation
    p_scipy = mannwhitneyu(x, y, alternative='two-sided').pvalue

    # Allow a small tolerance for differences in approximations
    assert abs(p_ours - p_scipy) < 1e-6
