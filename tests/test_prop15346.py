"""Tests for Prop 15.346 — two-line moment; LP miss; p=5 2+2."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15344 import HPLUS
from e1_gmin_m4_prop15346 import (
    bivar_lp,
    live_cross,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    two_two_types,
)


def test_cross_moment_exact():
    A = prove_A()
    assert A["proved"] is True
    c5 = live_cross(5)
    c7 = live_cross(7)
    assert c5["E_n2n2"] == Fraction(735, 26)
    assert c7["E_n2n2"] == Fraction(101585, 818)
    assert c5["E_n2n2"] != Fraction(36, 1)
    assert c5["E_n2n2"].denominator == HPLUS[5] // 5
    assert c5["E_n2n2"].denominator != HPLUS[5] // 10


def test_bivar_lp_misses_and_contains_live():
    B = prove_B()
    assert B["proved"] is True
    rec = bivar_lp(5)
    live = float(Fraction(735, 26))
    assert rec["min"] < live < rec["max"]
    assert rec["max"] > live + 1.0


def test_p5_22_constant_p7_splits():
    C = prove_C()
    assert C["proved"] is True
    r5 = two_two_types(5)
    r7 = two_two_types(7)
    assert r5["n_split"] == 0
    assert r5["n_const"] == r5["n_types"]
    assert r7["n_split"] > 0


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["cross_moment_exact"] is True
    assert out["proved"]["bivar_lp_misses"] is True
    assert out["proved"]["p5_22_type_constant"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
