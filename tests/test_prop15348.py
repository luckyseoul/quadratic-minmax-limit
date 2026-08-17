"""Tests for Prop 15.348 — n_4 not Jacobi-linear; rationals miss."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15346 import live_cross
from e1_gmin_m4_prop15348 import main, prove_A, prove_B, prove_C, prove_open


def test_linear_norm_features_miss():
    A = prove_A()
    assert A["proved"] is True
    assert A["err5"] > 0.5
    assert A["err7"] > 8.0
    assert A["n3_exact_subsets"] == 0


def test_weil_linear_misses_and_p_plus_W4_false():
    B = prove_B()
    assert B["proved"] is True
    assert B["err"] > 50
    # n_4 = p+W_4 is not an identity (some types may accidentally match)
    assert any(s["n4"] != s["p_plus_W4"] for s in B["samples"])


def test_no_deg31_and_neighbour_miss():
    C = prove_C()
    assert C["proved"] is True
    assert C["n_hits"] == 0
    assert Fraction(5 * 5 - 1, 2) != Fraction(735, 26)


def test_cross_still_those_rationals():
    assert live_cross(5)["E_n2n2"] == Fraction(735, 26)
    assert live_cross(7)["E_n2n2"] == Fraction(101585, 818)


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["n4_not_linear_in_norm_features"] is True
    assert out["proved"]["n4_not_linear_in_Weil"] is True
    assert out["proved"]["no_deg31_formula"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
