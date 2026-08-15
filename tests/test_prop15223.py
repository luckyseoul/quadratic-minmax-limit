"""Tests for Prop 15.223 — Path-C ⇒ Wick residual-i; δ bound open."""
from __future__ import annotations

import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15117 import room_hyp_over_24  # noqa: E402
from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    is_prime,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)
from e1_gmin_m4_prop15198 import eta_sumsq_budget_for_wick_hi  # noqa: E402
from e1_gmin_m4_prop15215 import eta_part_norm_sq  # noqa: E402
from e1_gmin_m4_prop15217 import delta_room_for_R  # noqa: E402
from e1_gmin_m4_prop15223 import (  # noqa: E402
    gap_wick_closed,
    hinge_status_223,
    hyp_inside_wick_gap,
    residual_i_closed_via_223,
    theorem_CS_closes_if_k2_over_N_le_gap,
    theorem_delta_bound_open,
    theorem_gap_wick_closed_form,
    theorem_path_C_inside_wick_gap,
)


def test_gap_wick_proved():
    assert theorem_gap_wick_closed_form()["proved"] is True


def test_path_C_inside_proved():
    assert theorem_path_C_inside_wick_gap()["proved"] is True


def test_CS_conditional_proved():
    assert theorem_CS_closes_if_k2_over_N_le_gap()["proved"] is True


def test_delta_still_open():
    assert theorem_delta_bound_open()["proved"] is False


def test_gap_wick_matches_B_minus_eta():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 47):
        assert is_prime(p)
        assert gap_wick_closed(p) == (
            eta_sumsq_budget_for_wick_hi(p) - eta_part_norm_sq(p)
        )
        assert gap_wick_closed(p) > 0


def test_hyp_strictly_inside_gap():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 47, 101):
        assert is_prime(p)
        d = hyp_inside_wick_gap(p)
        assert d == gap_wick_closed(p) - room_hyp_over_24(p)
        assert d > 0
        # explicit formula
        assert d == Fraction(
            16 * p * p * (p * p - 9),
            3 * (p * p - 5) * (p * p + 1),
        )


def test_gap_wick_values():
    assert gap_wick_closed(5) == Fraction(416, 15)
    assert gap_wick_closed(7) == Fraction(2000, 33)
    assert room_hyp_over_24(5) == Fraction(1536, 65)
    assert room_hyp_over_24(5) < gap_wick_closed(5)


def test_certified_CS_beats_gaps_at_p5_p7():
    c = theorem_delta_bound_open()["certified_small_p"]
    assert c["p5"]["CS_le_gap"] is True
    assert c["p5"]["CS_le_R_room"] is True
    assert c["p7"]["CS_le_gap"] is True
    assert c["p7"]["CS_le_R_room"] is True
    # exact fractions
    assert Fraction(6912, 260) <= gap_wick_closed(5)
    assert Fraction(6912, 260) <= delta_room_for_R(5)
    assert Fraction(193705, 11452) <= gap_wick_closed(7)


def test_residual_i_still_open():
    assert residual_i_closed_via_223() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True


def test_hinge_status():
    h = hinge_status_223()
    assert h["gap_wick_closed"] is True
    assert h["path_C_inside_wick"] is True
    assert h["CS_conditional"] is True
    assert h["delta_bound"] is False
    assert h["residual_i_closed_via_223"] is False
