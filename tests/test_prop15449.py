"""Tests for Prop 15.449 — L_NL on 15.448 slots; 19 not GJ."""
from __future__ import annotations

from fractions import Fraction
from math import gcd

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import named_gj_ints
from e1_gmin_m4_prop15430 import Lpar_R_named
from e1_gmin_m4_prop15449 import (
    bucket_slot_N,
    main,
    nl_L,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
)


def test_not_constant_on_nsmix():
    F, L, nnl = nl_L(5)
    _by, by_slot = bucket_slot_N(F, 5, L)
    assert nnl == 100
    assert by_slot["ns_mix"] == {12, 13}
    A = prove_A()
    assert A["proved"] is True
    assert "ns_mix" in A["by_p"]["5"]["split_slots"]


def test_constant_slots_and_LparR_fail():
    B = prove_B()
    assert B["proved"] is True
    assert B["live"]["5"]["Sub"] == ["24"]
    assert B["live"]["7"]["Sub"] == ["2113/19"]
    assert Lpar_R_named(7) == 133
    assert Lpar_R_named(7) != Fraction(2113, 19)


def test_19_essential_not_gj():
    C = prove_C()
    assert C["proved"] is True
    assert gcd(2113, 19) == 1
    assert 19 not in named_gj_ints(7)
    assert 2113 % 19 != 0


def test_ystar_not_nl_at_p7():
    D = prove_D()
    assert D["proved"] is True
    assert D["ystar5"] == 100
    assert D["L_Sub_ystar7"] == "1321/12"
    assert Fraction(D["L_Sub_ystar7"]) != Fraction(2113, 19)


def test_flags_untouched():
    E = prove_open()
    assert E["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.449"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["const_on_slot_N"] is True
    assert out["proved"]["den19_not_GJ"] is True
    assert out["proved"]["ystar_not_NL"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["algebra"]["D"]["L_Sub_ystar7"] == "1321/12"
