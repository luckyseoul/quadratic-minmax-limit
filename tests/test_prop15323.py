"""Tests for Prop 15.323 — AP-orbit W; Var_2orb sufficient criterion."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general
from e1_gmin_m4_prop15323 import (
    EW_named,
    Var_2orb,
    main,
    prove_AP_hist,
    prove_open,
    prove_sufficient,
)


def test_AP_W_not_pm2p_at_p13():
    A = prove_AP_hist()
    assert A["proved"] is True
    assert A["APpm2p_all_p1_false"] is True
    assert A["by_p"]["5"]["hist"] == {20: 15, 40: 15}
    assert A["by_p"]["7"]["hist"] == {56: 28, 98: 56}
    assert A["by_p"]["13"]["hist"] == {312: 91, 520: 182, 624: 182, 676: 91}
    ew13 = EW_named(13)
    assert not set(A["by_p"]["13"]["hist"]).issubset({ew13 - 26, ew13 + 26})


def test_Var2orb_sufficient_not_equal_at_p7():
    B = prove_sufficient()
    assert B["proved"] is True
    assert B["Var2orb_eq_live_p7_false"] is True
    assert Var_2orb(5) == Fraction(660, 13)
    assert Var_2orb(7) == Fraction(280, 3)
    assert Var_2orb(7) != Fraction(21504, 409)
    assert B["by_p"]["5"]["le"] is True
    assert B["by_p"]["13"]["le"] is True


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["A_square_named"] is False
    assert C["Var_2orb_proved"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["AP_orbit_hists"] is True
    assert out["proved"]["Var2orb_sufficient"] is True
    assert out["proved"]["A_square_named"] is False
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["proved"]["lambda_ge_6"] is False
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["L_status"] == "OPEN"
