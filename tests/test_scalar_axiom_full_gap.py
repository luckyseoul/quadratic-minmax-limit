"""The scalar axioms admit a gap equal to the whole sandwich."""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from scalar_axiom_full_gap import (
    axioms_hold,
    ends,
    first_n,
    frequency,
    full_gap_pair,
    high_slope,
)


def test_full_sandwich_gap_obeys_the_scalar_axioms():
    lam, lam_upper = full_gap_pair()
    assert lam_upper - lam == 0.5 - 1.0 / math.pi
    _a, b, margin = ends(lam, lam_upper)
    assert margin > 0.3
    assert high_slope(lam, lam_upper) < 1.0
    # The oscillation width is the whole proved sandwich.
    assert abs(2.0 * b - (lam_upper - lam)) < 1e-15
    n0 = first_n(lam, lam_upper)
    report = axioms_hold(lam, lam_upper, n0, n0 + 200_000)
    assert report["monotone"]
    assert report["step"]
    # One log-period is 2 pi.  The range is long enough to approach both ends.
    assert report["alpha_min"] < lam + 0.01
    assert report["alpha_max"] > lam_upper - 0.01


def test_every_inner_window_has_positive_margin():
    lam, lam_upper = full_gap_pair()
    for lower, upper in ((lam, lam), (lam, 0.4), (lam, lam_upper), (0.1, 0.5)):
        _a, b, margin = ends(lower, upper)
        assert margin >= 0.75 * lower - 1e-12
        assert b * frequency(lower, upper) <= 0.75 * lower + 1e-12
        assert first_n(lower, upper) < 10_000
