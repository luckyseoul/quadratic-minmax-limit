"""Tests for Prop 15.345 — AP n_3; single-line moments miss the floor."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15321 import Qpp_floor_ub
from e1_gmin_m4_prop15345 import (
    is_ap3,
    main,
    prove_A,
    prove_B,
    prove_C,
    prove_open,
    qpp_from_EM2,
)


def test_ap_classifier():
    assert is_ap3(0, 1, 2, 7) is True
    assert is_ap3(0, 1, 3, 7) is False
    # F_5: every 3-set is an AP
    assert is_ap3(0, 1, 2, 5) is True
    assert is_ap3(0, 1, 3, 5) is True


def test_n3_splits_on_AP_at_p7():
    A = prove_A()
    assert A["proved"] is True
    assert A["by_p"]["5"]["ap_unique"] == [14]
    assert A["by_p"]["7"]["ap_unique"] == [691]
    assert A["by_p"]["7"]["nap_unique"] == [699]
    assert 691 != 699


def test_moment3_lp_misses_floor():
    B = prove_B()
    assert B["proved"] is True
    assert B["by_p"]["5"]["closes"] is False
    assert B["by_p"]["5"]["Qpp_max"] > float(Qpp_floor_ub(5))


def test_moment4_lp_still_misses():
    C = prove_C()
    assert C["proved"] is True
    assert C["by_p"]["5"]["closes"] is False
    assert C["by_p"]["5"]["Qpp_max"] > float(Qpp_floor_ub(5))
    # fail-when-wrong neighbour: 26/7 is the floor ub, not the LP max
    assert C["by_p"]["5"]["Qpp_max"] != float(Qpp_floor_ub(5))


def test_qpp_from_EM2_p5_endpoints():
    # M-only extremes from 15.307 types
    assert abs(qpp_from_EM2(5, 900.0) - 0.0) < 1e-9
    assert abs(qpp_from_EM2(5, 1100.0) - 16.0) < 1e-9


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert D["Q_tau_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["n3_constant_on_AP"] is True
    assert out["proved"]["moment3_lp_misses_floor"] is True
    assert out["proved"]["moment4_lp_misses_floor"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
