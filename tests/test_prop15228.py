"""Tests for Prop 15.228 — diamond; Tμ_part; |κ|=1 ν=δ₊−δ₋; residual open."""
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
from e1_gmin_m4_prop15217 import delta_room_for_R  # noqa: E402
from e1_gmin_m4_prop15228 import (  # noqa: E402
    T_mupart_norm_sq,
    T_mupart_star_coeff,
    T_mupart_star_coeff_via_ab,
    diamond_L2_delta_ub,
    hinge_status_228,
    residual_i_closed_via_228,
    star_norm_sq,
    theorem_Aut_SOS_diamond_hinge,
    theorem_L2_T_orthogonality,
    theorem_T_mupart_closed,
    theorem_boolean_diamond,
    theorem_coupling_masters,
    theorem_diamond_L2_delta_ub,
    theorem_kappa1_nu_delta,
)


def test_coupling_proved():
    assert theorem_coupling_masters()["proved"] is True


def test_diamond_proved():
    assert theorem_boolean_diamond()["proved"] is True


def test_T_mupart_proved():
    assert theorem_T_mupart_closed()["proved"] is True


def test_kappa1_proved():
    assert theorem_kappa1_nu_delta()["proved"] is True


def test_L2_orth_proved():
    assert theorem_L2_T_orthogonality()["proved"] is True


def test_diamond_L2_proved():
    assert theorem_diamond_L2_delta_ub()["proved"] is True


def test_Aut_SOS_structure_only():
    t = theorem_Aut_SOS_diamond_hinge()
    assert t["proved_structure"] is True
    assert t["bound_proved_general"] is False


def test_T_mupart_coeff_matches_ab():
    for p in (5, 7, 11, 13, 17, 19, 23, 31, 47):
        assert is_prime(p)
        assert T_mupart_star_coeff(p) == T_mupart_star_coeff_via_ab(p)
        assert T_mupart_star_coeff(p) == Fraction(-8, p * p - 5)


def test_T_mupart_coeff_values():
    assert T_mupart_star_coeff(5) == Fraction(-2, 5)  # -8/20
    assert T_mupart_star_coeff(7) == Fraction(-2, 11)  # -8/44


def test_star_norm_positive():
    for p in (5, 7, 11, 13):
        assert star_norm_sq(p) > 0
        assert T_mupart_norm_sq(p) > 0


def test_diamond_L2_ub_positive_but_weak():
    for p in (5, 7, 11, 13, 17):
        ub = diamond_L2_delta_ub(p)
        assert ub > 0
        # far above R-room — not a hinge
        assert ub > delta_room_for_R(p)


def test_diamond_L2_p5_value():
    # structural smoke: ub is a concrete Fraction
    ub = diamond_L2_delta_ub(5)
    assert isinstance(ub, Fraction)
    assert ub > Fraction(36)  # R-room at p=5 is 182/5=36.4


def test_residual_i_still_open():
    assert residual_i_closed_via_228() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False


def test_hinge_status():
    h = hinge_status_228()
    assert h["T_mupart_closed"] is True
    assert h["boolean_diamond"] is True
    assert h["kappa1_nu_delta"] is True
    assert h["Aut_SOS_bound"] is False
    assert h["residual_i_closed_via_228"] is False
