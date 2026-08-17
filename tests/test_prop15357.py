"""Tests for Prop 15.357 — w=1 at p=3; k named; p=11 hole-types."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15330 import HPLUS
from e1_gmin_m4_prop15355 import Q_AP_named
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15357 import (
    D_live,
    k_named,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_p3_no_NL():
    A = prove_A()
    assert A["proved"] is True
    assert A["n_Hplus"] == n_1d(3) == 6
    assert A["w"] == 1


def test_k_and_D():
    B = prove_B()
    assert B["proved"] is True
    assert k_named(5) == 3
    assert k_named(7) == 10
    assert k_named(11) == 126
    assert D_live(5) == 13
    assert D_live(7) == 409
    assert D_live(7) != (49 + 1) // 2
    assert Fraction(k_named(5), D_live(5)) == Fraction(n_1d(5), HPLUS[5])


def test_p11_holes_not_AP():
    C = prove_C()
    assert C["proved"] is True
    assert C["rows"]["d1"]["Q"] == "544/45"
    assert C["rows"]["d2"]["Q"] == "632/45"
    assert Fraction(C["rows"]["d1"]["Q"]) != Q_AP_named(11)
    assert Fraction(C["rows"]["d2"]["Q"]) != Q_1d_pp_named(11)


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_NL_named_in_p"] is False
    assert D["w_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["w_eq_1_at_p3"] is True
    assert out["proved"]["k_integer"] is True
    assert out["proved"]["p11_hole_types"] is True
    assert out["proved"]["Q_NL_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
