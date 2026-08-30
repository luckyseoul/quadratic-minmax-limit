"""Prop 15.250: R-path ⇔ E[s⁴]≤15n²−22n; odd moments; residual_i open."""
from __future__ import annotations

from fractions import Fraction

import pytest

from e1_gmin_m4_prop15170 import (
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15217 import m4_target_R
from e1_gmin_m4_prop15250 import (
    Es4_R_bound,
    Es4_from_m4_norm,
    hinge_status_250,
    residual_i_closed_via_250,
    theorem_Es4_expansion,
    theorem_R_path_iff_Es4,
    theorem_odd_moments_vanish,
)


def test_odd_moments_vanish():
    t = theorem_odd_moments_vanish()
    assert t["proved"] is True


def test_Es4_expansion_identity():
    t = theorem_Es4_expansion()
    assert t["proved"] is True
    for p in (5, 7, 11, 13, 17):
        n = n_of(p)
        # known part without m4
        known = (
            n
            + Fraction(4 * n * (n - 1), p * p)
            + 3 * n * (n - 1)
            + Fraction(6 * n * (n - 1) * (n - 2), p * p)
        )
        assert known == Fraction(9 * n * n - 10 * n)


def test_R_path_iff_Es4():
    t = theorem_R_path_iff_Es4()
    assert t["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        n = n_of(p)
        # at the R threshold, Es4 hits the bound
        es4 = Es4_from_m4_norm(p, m4_target_R(p))
        assert es4 == Es4_R_bound(p)
        assert Es4_R_bound(p) == Fraction(15 * n * n - 22 * n)
        # below threshold => below Es4 bound
        es4_lo = Es4_from_m4_norm(p, m4_target_R(p) - Fraction(1, 1000))
        assert es4_lo < Es4_R_bound(p)


def test_predicates_remain_false():
    assert residual_i_closed_via_250() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is True
    assert e1_closed_general() is False
    h = hinge_status_250()
    assert h["residual_i_closed"] is False
    assert h["Es4_bound_general"] is False
    assert h["R_path_iff_Es4_bound"] is True
    assert h["odd_moments_vanish"] is True


def test_census_p5_if_cache():
    from pathlib import Path

    import numpy as np

    cache = Path("/tmp/maxplus_p5.npy")
    if not cache.exists():
        pytest.skip("no Max+ p5 cache")
    Y = np.load(cache)
    n = Y.shape[1]
    M = Y.shape[0]
    # odd moment sample
    mu3 = float(np.mean(Y[:, 0] * Y[:, 1] * Y[:, 2]))
    assert abs(mu3) < 1e-12
    # E[s^4]
    total = 0.0
    for s in range(0, M, 50):
        G = Y[s : s + 50] @ Y.T
        total += float(np.sum(G.astype(np.float64) ** 4))
    Es4 = total / (M * M)
    bound = float(Es4_R_bound(5))
    assert Es4 <= bound + 1e-6
    # expansion vs m4
    from itertools import combinations

    m4n2 = 0.0
    for S in combinations(range(n), 4):
        m = float(np.mean(np.prod(Y[:, S], axis=1)))
        m4n2 += m * m
    Es4_pred = 9 * n * n - 10 * n + 24 * m4n2
    assert abs(Es4 - Es4_pred) < 1e-6
    assert m4n2 <= float(m4_target_R(5)) + 1e-6
