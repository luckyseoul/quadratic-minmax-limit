"""Tests for Prop 15.478 — Paley++ type-indicator Fourier."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15478 import (
    LIVE_LIN,
    a_E_sq,
    fold_pack,
    main,
    n_pp_drop_aE,
    n_pp_interpolant,
    n_pp_named,
    prove_A,
    prove_B,
    prove_C,
    prove_D,
    prove_open,
)


def test_n_pp_formula():
    assert n_pp_named(3) == 0
    assert n_pp_named(5) == 2
    assert n_pp_named(7) == 6
    assert n_pp_named(11) == 16
    assert n_pp_drop_aE(5) == Fraction(5, 2)
    assert n_pp_interpolant(11) == 20
    assert n_pp_interpolant(11) != n_pp_named(11)
    assert a_E_sq(7) == 0
    assert a_E_sq(5) == 4
    A = prove_A()
    assert A["proved"] is True
    assert A["rows"]["11"]["live"] == 16


def test_K_lin_gauss():
    B = prove_B()
    assert B["proved"] is True
    assert B["rows"]["5"]["err"] < 1e-10
    assert B["rows"]["5"]["err_drop3"] > 0.3


def test_lin_misses_live_Qpp():
    C = prove_C()
    assert C["proved"] is True
    assert LIVE_LIN[5] == Fraction(87, 13)
    assert LIVE_LIN[5] != LIVE_QPP[5]
    rec = fold_pack(5)
    assert abs(rec["Q_lin"] - 87 / 13) < 1e-9
    assert abs(rec["Q_full"] - float(LIVE_QPP[5])) < 1e-9


def test_quartic_and_twolevel_dead():
    D = prove_D()
    assert D["proved"] is True
    assert D["n_twolevel"]["7"] == 0
    assert D["proj_err"]["5"] > 2.0


def test_flags_untouched():
    E = prove_open()
    assert E["proved"] is False
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.478"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["n_pp_CM_named"] is True
    assert out["proved"]["lin_misses_Qpp"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
