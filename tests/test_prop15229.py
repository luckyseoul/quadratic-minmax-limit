"""Tests for Prop 15.229 — size-3 vanishes; functional eq; residual open."""
from __future__ import annotations

import sys
from fractions import Fraction
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15229 import (  # noqa: E402
    R_bar_budget_for_mu_thr,
    hinge_status_229,
    inner_s_of_C6,
    kappa_of_C6,
    residual_i_closed_via_229,
    theorem_conditional_R_bound_closes,
    theorem_functional_equation,
    theorem_s_vanishes_on_kappa1,
    theorem_size3_collision_formula,
    theorem_size3_vanishes_kappa1,
)


def test_s_vanishes_proved():
    assert theorem_s_vanishes_on_kappa1()["proved"] is True


def test_size3_collision_proved():
    assert theorem_size3_collision_formula()["proved"] is True


def test_size3_vanishes_proved():
    assert theorem_size3_vanishes_kappa1()["proved"] is True


def test_functional_eq_proved():
    assert theorem_functional_equation()["proved"] is True


def test_conditional_implication_proved():
    t = theorem_conditional_R_bound_closes()
    assert t["proved_implication"] is True
    assert t["bound_proved_general"] is False


def test_s_zero_on_all_kappa1_labelings():
    for bits in product((-1, 1), repeat=6):
        k = kappa_of_C6(bits)
        s = inner_s_of_C6(bits)
        if abs(k) == 1:
            assert s == 0
        if abs(k) == 3:
            assert abs(s) == 12


def test_s_exhaustion_counts():
    counts = {1: 0, -1: 0, 3: 0, -3: 0}
    for bits in product((-1, 1), repeat=6):
        counts[kappa_of_C6(bits)] += 1
    assert counts[1] == 24
    assert counts[-1] == 24
    assert counts[3] == 8
    assert counts[-3] == 8


def test_R_bar_budget_positive():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 47):
        assert is_prime(p)
        bud = R_bar_budget_for_mu_thr(p)
        assert bud > 0
        # (4(p-2) + bud)/(p^4-1) == 1/(2p)
        lhs = (4 * (p - 2) + bud) / Fraction(p**4 - 1)
        assert lhs == Fraction(1, 2 * p)


def test_R_bar_budget_values():
    # p=5: (625-1)/10 - 12 = 62.4 - 12 = 50.4 = 252/5
    assert R_bar_budget_for_mu_thr(5) == Fraction(624, 10) - 12
    assert R_bar_budget_for_mu_thr(5) == Fraction(252, 5)


def test_residual_i_still_open():
    assert residual_i_closed_via_229() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True


def test_hinge_status():
    h = hinge_status_229()
    assert h["size3_vanishes"] is True
    assert h["functional_equation"] is True
    assert h["R_bound_general"] is False
    assert h["residual_i_closed_via_229"] is False
