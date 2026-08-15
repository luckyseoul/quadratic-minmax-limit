"""Tests for Prop 15.218 — Max+⊥Max−; Gsum=2m₁μ₄; free-e dual structure."""
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
from e1_gmin_m4_prop15218 import (  # noqa: E402
    hinge_status_218,
    residual_i_closed_via_218,
    sum_ne_row_ones_dual,
    theorem_Gsum_eq_2_m1_mu4,
    theorem_Max_plus_perp_Max_minus,
    theorem_Tr_Gsum_sq,
    theorem_row_ones_dual_beats_need_if_Gmin,
)


def test_Max_plus_perp_proved():
    t = theorem_Max_plus_perp_Max_minus()
    assert t["proved"] is True


def test_Gsum_mu4_identity_proved():
    t = theorem_Gsum_eq_2_m1_mu4()
    assert t["proved"] is True
    assert "1/(2p)" in t["consequence"]


def test_Tr_Gsum_sq_proved():
    t = theorem_Tr_Gsum_sq()
    assert t["proved"] is True


def test_row_ones_dual_beats_need():
    t = theorem_row_ones_dual_beats_need_if_Gmin()
    assert t["proved"] is True
    # spot-check a=p/(2p+1), cost < need at sample primes
    for p in (5, 7, 11, 13, 17, 19, 23, 47, 101):
        assert is_prime(p)
        G_min = Fraction(-1, p)
        sne = sum_ne_row_ones_dual(p, G_min)
        al = Fraction(6, p * (p * p + 1))
        assert al * sne < 2 - al


def test_sum_ne_formula_p5():
    # a = 5/11, n=26
    # sum_ne = (26*25-2)/2 + (5/11)*26*(2-26)
    #        = 324 + (5/11)*26*(-24) = 324 - 3120/11 = (3564-3120)/11 = 444/11
    sne = sum_ne_row_ones_dual(5, Fraction(-1, 5))
    assert sne == Fraction(444, 11)


def test_residual_i_218_path_open_live_dual_eq_closed():
    assert residual_i_closed_via_218() is False
    assert residual_i_dual_eq_empty_proved_general() is True
    assert type_I_k_3p_minus_2_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True


def test_hinge_status():
    h = hinge_status_218()
    assert h["Max_plus_perp_Max_minus"] is True
    assert h["Gsum_eq_2_m1_mu4"] is True
    assert h["Tr_Gsum_sq"] is True
    assert h["row_ones_dual_if_Gmin"] is True
    assert h["residual_i_closed_via_218"] is False
    assert "OPEN" in h["open_hinge"] or "1/(2p)" in h["open_hinge"]
