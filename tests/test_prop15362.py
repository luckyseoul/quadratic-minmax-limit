"""Tests for Prop 15.362 — N is the universal r-net 0-1 count."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15361 import n_1d_plus_hoff
from e1_gmin_m4_prop15362 import (
    T_A008300,
    count_net,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_net_count_independent_of_slopes():
    A = prove_A()
    assert A["proved"] is True
    assert A["p3"] == 6
    assert A["p5_inf01"] == 130
    assert A["p5_024"] == 130
    assert A["p7_mixed"] == 5726
    assert A["p7_QR"] == 5726
    assert A["p7_mixed"] != n_1d_plus_hoff(7)


def test_two_dir_is_A008300_not_N_at_p5():
    B = prove_B()
    assert B["proved"] is True
    assert B["T31"] == 6
    assert B["T52"] == 2040
    assert T_A008300[(5, 2)] == 2040
    assert T_A008300[(5, 2)] != 130
    assert count_net(5, ["inf", 0]) == 2040


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["D_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is False
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["net_count_independent"] is True
    assert out["proved"]["two_dir_is_A008300"] is True
    assert out["proved"]["D_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
