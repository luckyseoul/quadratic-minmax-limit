"""Tests for Prop 15.367 — second NUM_SUM; Gal-classes unnamed."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15292 import n_1d
from e1_gmin_m4_prop15355 import n_pp_named
from e1_gmin_m4_prop15366 import NUM_SUM, neigh_16_n1d
from e1_gmin_m4_prop15367 import (
    A_num,
    NUM_SUM_drop_k,
    NUM_SUM_named,
    k_named,
    main,
    prove_A,
    prove_B,
    prove_open,
)


def test_named_NUM_hits_both():
    A = prove_A()
    assert A["proved"] is True
    assert NUM_SUM_named(5) == int(NUM_SUM(5)) == 480
    assert NUM_SUM_named(7) == int(NUM_SUM(7)) == 21616
    assert A_num(5) == 2 * n_1d(5) - 4 * n_pp_named(5) - 3 * k_named(5) == 6
    assert A_num(7) == 193


def test_fail_when_wrong_drop_k_and_16n1d():
    assert NUM_SUM_drop_k(5) != 480
    assert NUM_SUM_drop_k(5) == 2880
    assert neigh_16_n1d(7) == 2240
    assert neigh_16_n1d(7) != int(NUM_SUM(7))
    assert NUM_SUM_named(7) != neigh_16_n1d(7)


def test_no_3key_classifier():
    B = prove_B()
    assert B["proved"] is True
    assert B["Wchi_same"] is True
    assert B["Wchi_classes"][0] != B["Wchi_classes"][1]
    assert B["fiber_same"] is True
    assert B["fiber_classes"][0] != B["fiber_classes"][1]


def test_floor_still_open():
    C = prove_open()
    assert C["proved"] is False
    assert C["Q_tau_named_in_p"] is False
    assert C["D_named_in_p"] is False
    assert phi_F_ge_6_proved_general() is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["NUM_SUM_named_16pA"] is True
    assert out["proved"]["no_3key_gal_classifier"] is True
    assert out["proved"]["Q_tau_named_in_p"] is False
    assert out["L_status"] == "OPEN"
    assert out["series"].startswith("15.x")
