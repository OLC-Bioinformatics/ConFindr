#!/usr/bin/env python3

"""
Tests for Poisson Binomial distribution implementations in confindr_src.methods
"""

# Standard imports
import math

# Third-party imports
import numpy as np
import pytest

# Local imports
from confindr_src.methods import (
    _poisson_binomial_pmf_fft,
    poisson_binomial_tail
)


def test_fft_pmf_matches_convolution():
    """
    Test that the PMF computed via FFT matches that computed via brute-force
    convolution
    """
    # Example probability list
    p_list = [0.1, 0.2, 0.3]

    # PMF via FFT
    pmf_fft = _poisson_binomial_pmf_fft(p_list=p_list)

    # Brute-force convolution
    pmf_conv = np.array([1.0])

    # Convolve for each probability
    for p in p_list:
        pmf_conv = np.convolve(pmf_conv, [1.0 - p, p])

    assert pmf_fft.shape == pmf_conv.shape
    assert np.allclose(pmf_fft, pmf_conv, atol=1e-12)


def test_poisson_binomial_tail_exact_vs_pmf():
    """
    Test that the tail probability computed via PMF matches that computed
    directly
    """
    # Example probability list
    p_list = [0.1, 0.2, 0.25, 0.3]

    # Tail probability for k=2
    k = 2

    # Compute tail via PMF
    pmf = _poisson_binomial_pmf_fft(p_list=p_list)

    # Compute tail directly
    expected = pmf[k:].sum()

    # Compute tail via function
    actual = poisson_binomial_tail(k=k, p_list=p_list)

    assert math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-12)


def test_poisson_binomial_tail_fallback_to_normal():
    """
    Test that the tail probability falls back to normal approximation correctly
    """
    # Example probability list that may cause FFT issues
    p_list = [0.01] * 2000  # large-ish to make normal approx sensible
    k = 5

    # Compute tail via function
    tail = poisson_binomial_tail(k=k, p_list=p_list)

    # Manual normal approximation with continuity correction
    mu = sum(p_list)
    var = sum(p * (1 - p) for p in p_list)
    sigma = math.sqrt(var)

    # Survival function with continuity correction (k - 0.5)
    if sigma == 0:
        expected = 0.0 if k > mu else 1.0
    else:
        z = (k - 0.5 - mu) / (sigma * math.sqrt(2.0))
        expected = 0.5 * math.erfc(z)

    assert pytest.approx(expected, rel=1e-6) == tail
