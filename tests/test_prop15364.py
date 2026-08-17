"""Tests for Prop 15.364 — 15.290 type Q; S0 Q_N* is a mixture."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15326 import Qn_from_Qpp
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15356 import Q_1d_pp_named
from e1_gmin_m4_prop15364 import (
    Q_types_of_cache,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_p5_types_and_p7_split():
    A = prove_A()
    assert A["proved"] is True
    q5 = Q_types_of_cache(5)
    assert q5[str((1, 1, "sub"))] == Fraction(48, 13)
    assert q5[str((-1, -1, "N*"))] == Fraction(32, 13)
    assert Fraction(A["p7_m1m1_N1"]) == Fraction(1376, 409)
    assert Fraction(A["p7_m1m1_Ns"]) == Fraction(1496, 409)
    assert A["p7_m1m1_N1"] != A["p7_m1m1_Ns"]


def test_paley11_is_ensemble_not_1d():
    B = prove_B()
    assert B["proved"] is True
    assert Fraction(B["p5_11_sub"]) == LIVE_QPP[5]
    assert Fraction(B["p5_11_sub"]) != Q_1d_pp_named(5)
    assert Fraction(B["p7_11_sub"]) == LIVE_QPP[7]
    assert Q_1d_pp_named(5) == Fraction(16, 9)


def test_S0_QN_is_not_a_type_value():
    C = prove_C()
    assert C["proved"] is True
    assert Qn_from_Qpp(7, LIVE_QPP[7]) == Fraction(1440, 409)
    assert Qn_from_Qpp(7, LIVE_QPP[7]) != Fraction(1496, 409)
    assert Qn_from_Qpp(7, LIVE_QPP[7]) != Fraction(1376, 409)


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
    assert out["proved"]["type_constant_and_p7_split"] is True
    assert out["proved"]["paley11_equals_ensemble_Qpp"] is True
    assert out["proved"]["S0_QN_is_mixture"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
