"""Tests for Prop 15.532 — T_v free on H_+ iff χ(v)=−1."""
from __future__ import annotations

from math import comb

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15522 import D_named
from e1_gmin_m4_prop15532 import (
    Fix_named,
    m_of,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_fix_named_is_binom_on_squares_else_zero():
    assert Fix_named(5, 1) == comb(5, 3) == 10
    assert Fix_named(7, 1) == comb(7, 4) == 35
    assert Fix_named(5, -1) == 0
    assert Fix_named(7, -1) == 0
    assert n_1d(5) == m_of(5) * comb(5, 3)
    assert n_1d(7) == m_of(7) * comb(7, 4)


def test_live_fix_matches_and_T1_not_free():
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"]["5"]["sq"]["live_fix"] == 10
    assert A["by_p"]["5"]["ns"]["live_fix"] == 0
    assert A["by_p"]["7"]["sq"]["live_fix"] == 35
    assert A["by_p"]["7"]["ns"]["live_fix"] == 0


def test_n_orb_is_two_D_not_D():
    B = prove_B()
    assert B["proved"] is True
    for p in (5, 7):
        assert B["by_p"][str(p)]["n_orb"] == 2 * D_named(p)
        assert B["by_p"][str(p)]["n_orb"] != D_named(p)


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["D_named_as_orbits"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["phi_F_imported"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False
    out = main()
    assert out["A"]["proved"] is True
    assert out["B"]["proved"] is True
    assert out["C"]["proved"] is False
    assert out["phi_F_ge_6"] is False
    assert out["L_status"] == "OPEN"
