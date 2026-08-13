"""Tests for Prop 15.261 — Max± configuration multiplicities."""
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15261 import (
    hinge_status_261,
    mult_minus_formula,
    mult_plus_formula,
    residual_i_closed_via_261,
    theorem_mult_formula,
    theorem_mult_nonneg_bound,
    theorem_orientation_mult_equivalence,
)


def test_theorems_proved():
    assert theorem_mult_formula()["proved"] is True
    assert theorem_mult_nonneg_bound()["proved"] is True
    assert theorem_orientation_mult_equivalence()["proved"] is True
    assert theorem_mult_nonneg_bound()["bound_sufficient_for_residual_i"] is False
    assert theorem_orientation_mult_equivalence()["hypothesis_proved_general"] is False


def test_mult_formula_consistency():
    # sum over 16 alphas of mult formula = N when sum pr=0, sum qf pr=0
    N, p, m4 = 260, 5, 0.046153846153846156
    # use a simple A with known qf distribution not needed for sum
    total = 0.0
    from itertools import product
    import numpy as np

    A = np.array(
        [
            [0.0, 1, 1, 1],
            [1, 0, 1, -1],
            [1, 1, 0, 1],
            [1, -1, 1, 0],
        ]
    )
    for bits in product([-1.0, 1.0], repeat=4):
        alpha = np.array(bits)
        pr = int(np.prod(alpha))
        qf = int(round(float(alpha @ A @ alpha)))
        total += mult_plus_formula(N, p, qf, pr, m4)
    assert abs(total - N) < 1e-8


def test_bound_weaker_than_thr():
    for p in (5, 7, 11, 13):
        assert 1 - 4 / p > 1 / (2 * p)


def test_predicates_open():
    assert residual_i_closed_via_261() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_261()
    assert h["nu_zero_general"] is False
    assert h["residual_i_closed"] is False
