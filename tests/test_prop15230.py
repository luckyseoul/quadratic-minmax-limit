"""Tests for Prop 15.230 — R_part_max ≤ budget; Cy≡δ unify; residual open."""
from __future__ import annotations

import sys
from fractions import Fraction
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
from e1_gmin_m4_prop15225 import room_delta  # noqa: E402
from e1_gmin_m4_prop15229 import R_bar_budget_for_mu_thr  # noqa: E402
from e1_gmin_m4_prop15230 import (  # noqa: E402
    R_part_coeffs,
    R_part_max,
    R_part_max_closed,
    budget_minus_R_part_max,
    hinge_status_230,
    residual_i_closed_via_230,
    theorem_R_part_coeffs,
    theorem_R_part_max_closed,
    theorem_R_part_max_le_budget,
    theorem_channel_unification,
)


def test_R_part_coeffs_proved():
    assert theorem_R_part_coeffs()["proved"] is True


def test_R_part_max_closed_proved():
    assert theorem_R_part_max_closed()["proved"] is True


def test_R_part_max_le_budget_proved():
    assert theorem_R_part_max_le_budget()["proved"] is True


def test_channel_unification_proved():
    assert theorem_channel_unification()["proved"] is True


def test_R_part_coeffs_signs():
    for p in (5, 7, 11, 13, 17, 19, 23, 31):
        ck, cph = R_part_coeffs(p)
        assert ck > 0
        assert cph < 0


def test_R_part_max_matches_closed():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 101):
        assert is_prime(p)
        assert R_part_max(p) == R_part_max_closed(p)


def test_R_part_max_le_budget_values():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 47):
        assert R_part_max_closed(p) <= R_bar_budget_for_mu_thr(p)
        slack = budget_minus_R_part_max(p)
        assert slack == R_bar_budget_for_mu_thr(p) - R_part_max_closed(p)
        assert slack > 0


def test_slack_equals_room_delta_times_p4m1():
    for p in (5, 7, 11, 13, 17, 19, 23):
        slack = budget_minus_R_part_max(p)
        assert slack / (p**4 - 1) == room_delta(p)


def test_R_part_max_p5():
    # 4116/125 = 32.928
    assert R_part_max_closed(5) == Fraction(4116, 125)
    assert R_bar_budget_for_mu_thr(5) == Fraction(252, 5)
    assert budget_minus_R_part_max(5) == Fraction(252, 5) - Fraction(4116, 125)


def test_cubic_positive():
    for p in range(5, 100):
        if is_prime(p):
            assert p**3 - 2 * p * p - 13 * p + 18 > 0


def test_residual_i_still_open():
    assert residual_i_closed_via_230() is False
    assert residual_i_dual_eq_empty_proved_general() is False
    assert type_I_k_3p_minus_2_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False


def test_hinge_status():
    h = hinge_status_230()
    assert h["R_part_max_le_budget"] is True
    assert h["channel_unification"] is True
    assert h["residual_i_closed_via_230"] is False
