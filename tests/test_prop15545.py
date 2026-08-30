"""Tests for Prop 15.545 — NL mass 16pA − n_1d Q_1d^{++}."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15545 import (
    M_NL_drop_1d,
    M_NL_drop_k,
    M_NL_live,
    M_NL_named,
    main,
    one_d_mass,
    one_d_mass_wrong_sub,
    prove_A,
    prove_B,
    prove_open,
)


def test_one_d_mass_uses_Q1d_pp_not_sub():
    A = prove_A()
    assert A["proved"] is True
    assert one_d_mass(5) == Fraction(160, 3)
    assert one_d_mass_wrong_sub(5) == 160
    assert one_d_mass(5) != one_d_mass_wrong_sub(5)
    assert one_d_mass(7) != one_d_mass_wrong_sub(7)


def test_M_NL_hits_live_and_fails_drops():
    B = prove_B()
    assert B["proved"] is True
    assert M_NL_named(5) == M_NL_live(5) == Fraction(1280, 3)
    assert M_NL_named(7) == M_NL_live(7) == Fraction(61936, 3)
    assert M_NL_drop_1d(5) == 480
    assert M_NL_drop_1d(5) != M_NL_live(5)
    assert M_NL_drop_k(5) != M_NL_live(5)
    assert M_NL_drop_k(7) != M_NL_live(7)


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert C["NUM_SUM_16pA_general"] is False
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
