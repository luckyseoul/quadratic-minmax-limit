"""Tests for Prop 15.467 — coarse L++ / ns L / named atoms."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15298 import mu_minus, mu_plus
from e1_gmin_m4_prop15414 import Lpar_1d_named
from e1_gmin_m4_prop15467 import (
    LIVE_COARSE,
    coarse_L,
    main,
    named_atom,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
)


def test_coarse_values():
    A = prove_A()
    assert A["proved"] is True
    C5, C7 = coarse_L(5), coarse_L(7)
    assert C5["++"][1] == LIVE_COARSE[5]["++"][1] == Fraction(985, 39)
    assert C7["++"][1] == LIVE_COARSE[7]["++"][1] == Fraction(91581, 818)
    assert C5["++"][1] != mu_plus(5) ** 2
    assert C7["++"][1] != mu_plus(7) ** 2


def test_ns_is_mupm():
    B = prove_B()
    assert B["proved"] is True
    assert coarse_L(5)["ns"][1] == mu_plus(5) * mu_minus(5) == Fraction(25, 2)
    assert coarse_L(7)["ns"][1] == mu_plus(7) * mu_minus(7) == Fraction(147, 2)
    assert Fraction(180, 13) != Fraction(25, 2)


def test_named_atoms_miss():
    C = prove_C()
    assert C["proved"] is True
    assert named_atom("Lpar_1d", 5) == Lpar_1d_named(5) != Fraction(985, 39)
    assert named_atom("mup_sq", 5) != coarse_L(5)["++"][1]
    assert named_atom("Lpar_1d", 7) != coarse_L(7)["++"][1]


def test_D_mod_p_dies_at_p3():
    D = prove_D()
    assert D["proved"] is True
    assert D["D"][3] % 3 != 0
    assert D["D"][5] % 5 == 3
    assert D["D"][7] % 7 == 3


def test_flags_untouched():
    E = prove_open()
    assert E["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.467"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["coarse_L_values"] is True
    assert out["proved"]["coarse_Lns_is_mupm"] is True
    assert out["proved"]["named_atoms_miss_Lpp"] is True
    assert out["proved"]["D_mod_p_not_3_at_p3"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
