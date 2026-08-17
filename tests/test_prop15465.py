"""Tests for Prop 15.465 — Paley-N*/Qn named ratio; drop (p-5)."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15326 import Qn_from_Qpp
from e1_gmin_m4_prop15330 import LIVE_QPP
from e1_gmin_m4_prop15464 import live_QN_paley
from e1_gmin_m4_prop15465 import (
    J_chi_chi,
    J_chi_chi_named,
    R_drop_pm5,
    R_named,
    four_from_J,
    live_ratio,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_ratio_hits_live():
    A = prove_A()
    assert A["proved"] is True
    assert R_named(5) == 1
    assert R_named(7) == Fraction(187, 180)
    assert live_ratio(5) == live_QN_paley(5) / Qn_from_Qpp(5, LIVE_QPP[5])
    assert live_ratio(5) == 1
    assert live_ratio(7) == Fraction(187, 180)
    assert R_named(5) == live_ratio(5)
    assert R_named(7) == live_ratio(7)


def test_drop_pm5_misses_p5():
    assert R_drop_pm5(5) == Fraction(77, 72)
    assert R_drop_pm5(5) != 1
    assert R_drop_pm5(5) != R_named(5)
    # 409 cancelled: live ratio is 1496/1440
    assert live_QN_paley(7) / Qn_from_Qpp(7, LIVE_QPP[7]) == Fraction(1496, 1440)


def test_Jacobi_gives_four():
    B = prove_B()
    assert B["proved"] is True
    assert J_chi_chi(5) == J_chi_chi_named(5) == -1
    assert J_chi_chi(7) == J_chi_chi_named(7) == 1
    assert four_from_J(5) == 4
    assert four_from_J(7) == 4
    assert J_chi_chi(5) != 0


def test_flags_untouched():
    C = prove_open()
    assert C["proved"] is False
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["prop"] == "15.465"
    assert out["L_status"] == "OPEN"
    assert out["proved"]["ratio_named_live"] is True
    assert out["proved"]["four_from_Jacobi"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
