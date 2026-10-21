"""Exact doubling does not force the normalized slope to converge."""
from __future__ import annotations

from minmax_quadratic import doubling_only_oscillator


def test_exact_doubling_keeps_two_slopes():
    amplitude = 0.05
    samples = (3.0, 10.0, 100.0, 2.0**20 + 3.0)
    for point in samples:
        doubled = doubling_only_oscillator(2.0 * point, amplitude)
        assert abs(doubled - 2.0 * doubling_only_oscillator(point, amplitude)) < 1e-9

    gaps = []
    for step in (8, 16, 32):
        low = doubling_only_oscillator(2.0 ** (step + 0.75), amplitude) / 2.0 ** (step + 0.75)
        high = doubling_only_oscillator(2.0 ** (step + 0.25), amplitude) / 2.0 ** (step + 0.25)
        gaps.append(high - low)
    assert gaps[0] > 2.0 * amplitude - 1e-9
    assert abs(gaps[0] - gaps[1]) < 1e-9
    assert abs(gaps[1] - gaps[2]) < 1e-9

    tripling_gaps = []
    for step in (8, 24):
        base = 2.0**step
        tripling_gaps.append(
            abs(
                doubling_only_oscillator(3.0 * base, amplitude) / (3.0 * base)
                - doubling_only_oscillator(base, amplitude) / base
            )
        )
    assert tripling_gaps[0] > 0.0
    assert abs(tripling_gaps[0] - tripling_gaps[1]) < 1e-9
