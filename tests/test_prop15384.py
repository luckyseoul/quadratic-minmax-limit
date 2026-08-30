"""Tests for Prop 15.384 — Var(M) inversion; E[M²] catalog empty."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15384 import (
    EM2_named,
    EM_named,
    LIVE_QPP,
    Var_M_drop_5mp,
    Var_M_named,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_Var_M_identity():
    A = prove_A()
    assert A["proved"] is True
    assert Var_M_named(5, LIVE_QPP[5]) == Fraction(600, 13)
    assert Var_M_named(7, LIVE_QPP[7]) == Fraction(53361, 409)
    assert Var_M_drop_5mp(7, LIVE_QPP[7]) != Fraction(53361, 409)


def test_EM2_catalog_empty():
    B = prove_B()
    assert B["proved"] is True
    assert EM2_named(5, LIVE_QPP[5]) == Fraction(12300, 13)
    assert EM2_named(7, LIVE_QPP[7]) == Fraction(2939265, 409)
    assert EM_named(5) ** 2 != EM2_named(5, LIVE_QPP[5])
    assert B["n_1term_hits"] == 0


def test_phi_F_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["Var_M_identity"] is True
    assert out["proved"]["EM2_catalog_empty"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
