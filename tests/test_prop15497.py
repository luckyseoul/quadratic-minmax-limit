"""Tests for Prop 15.497 — J_{N*} named; pairing reduced to D."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15497 import (
    D_live,
    D_need_named,
    J_Nstar_named,
    J_Nstar_swapped,
    main,
    max_pair_named,
    pair_named,
    pair_signflip,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_J_Nstar_divisibility():
    A = prove_A()
    assert A["proved"] is True
    assert J_Nstar_named(5, 4) == -2
    assert J_Nstar_swapped(5, 4) == -4
    assert J_Nstar_named(5, 6) == -(5 - 1)
    assert J_Nstar_named(5, 2) == 2
    assert J_Nstar_named(7, 8) == -(7 - 1)
    assert J_Nstar_named(7, 6) == -(7 - 3)
    assert A["by_p"][5]["values"] == [-4, -2, 2]


def test_pairing_identity_p5():
    B = prove_B()
    assert B["proved"] is True
    assert B["max_named"] == "24/13"
    D = D_live(5)
    assert abs(float(pair_named(5, 2, D)) - 24 / 13) < 1e-9
    assert abs(float(pair_signflip(5, 2, D)) - float(pair_named(5, 2, D))) > 0.1


def test_D_need_not_D_name():
    C = prove_C()
    assert C["proved"] is True
    assert D_need_named(5) == Fraction(64, 5)
    assert D_live(5) == 13
    assert max_pair_named(5, 13) == Fraction(24, 13)
    assert max_pair_named(5, 13) <= 2
    assert max_pair_named(5, 3) > 2


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["phi_F_imported"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["phi_F_ge_6"] is False
