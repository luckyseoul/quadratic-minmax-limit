"""Tests for Prop 15.501 — Hoffman window meets the p≡1 pairing bound."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15326 import Q_lo_named
from e1_gmin_m4_prop15501 import (
    Qub_expanded,
    ineq_num,
    ineq_num_drop_n1d,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
)


def test_Qub_expand():
    A = prove_A()
    assert A["proved"] is True
    assert Qub_expanded(5) == Fraction(26, 7)
    assert Qub_expanded(5) == Qpp_floor_ub(5)
    assert Qub_expanded(5) != Q_lo_named(5)


def test_inequality_holds_drop_n1d_fails():
    B = prove_B()
    assert B["proved"] is True
    assert ineq_num(5) == 8
    assert ineq_num(5) > 0
    assert ineq_num_drop_n1d(5) < 0
    for p in (5, 13, 17, 29, 37):
        assert ineq_num(p) > 0


def test_window_meets_need_does_not_name_D():
    C = prove_C()
    assert C["proved"] is True
    assert C["by_p"][5]["c_min"] == C["by_p"][5]["c_eq"] == 1
    assert C["by_p"][13]["nwin"] > 1
    assert C["by_p"][13]["D_cmin"] != C["by_p"][13]["D_ceq"]


def test_flags_untouched():
    D = prove_open()
    assert D["proved"] is False
    assert D["phi_F_imported"] is False
    assert D["p1_floor_unconditional"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is True
    assert out["phi_F_ge_6"] is False
