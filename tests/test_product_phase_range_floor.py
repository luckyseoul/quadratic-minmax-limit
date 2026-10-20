"""Product and phase-range payments do not lift a doubling above the Gaussian floor."""
from __future__ import annotations

from minmax_quadratic import (
    alpha,
    dual_gaussian_lower_bound,
    optimized_product_delta,
    phase_range_lower,
    product_forced_minimum,
    rademacher_absolute_moment,
    spherical_half_bound,
)


def _brute_absolute_moment(width: int) -> float:
    total = 0
    for mask in range(1 << width):
        signed = 0
        for bit in range(width):
            signed += 1 if mask & (1 << bit) else -1
        total += abs(signed)
    return total / float(1 << width)


def test_rademacher_moment_matches_the_cube_average():
    for width in range(1, 11):
        shipped = rademacher_absolute_moment(width)
        brute = _brute_absolute_moment(width)
        assert abs(shipped - brute) <= 1e-12 * max(1.0, brute)


def test_short_window_product_bound_stays_inside_the_sandwich():
    kept, nxt = 80, 88
    forced = product_forced_minimum(
        spherical_half_bound(kept),
        optimized_product_delta(nxt - kept),
    )
    assert alpha(dual_gaussian_lower_bound(nxt), nxt) < alpha(forced, nxt)
    assert alpha(forced, nxt) < alpha(spherical_half_bound(nxt), nxt)


def test_longer_windows_fall_through_the_gaussian_floor():
    kept = 80
    past = spherical_half_bound(kept)
    for nxt in (120, 160):
        deleted = nxt - kept
        product = product_forced_minimum(past, optimized_product_delta(deleted))
        phase = 0.5 * (past + phase_range_lower(deleted))
        floor = alpha(dual_gaussian_lower_bound(nxt), nxt)
        assert alpha(product, nxt) < floor
        assert alpha(phase, nxt) < floor
