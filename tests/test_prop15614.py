"""Tests for Prop 15.614 — W1 for p≡3; named vectors miss W2."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15614 import (
    theorem_A_lift,
    theorem_B_W1_p_eq_3,
    theorem_C_two_fiber,
    theorem_D_W2_miss,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_lift_eps_is_weight():
    A = theorem_A_lift()
    assert A["proved"] is True
    assert A["rows"]["3"]["eps"] == 1
    assert A["rows"]["5"]["eps"] == 0
    assert A["rows"]["7"]["eps"] == 1
    for rec in A["rows"].values():
        assert rec["ker_S_aff"] is True
        assert rec["Dv_plus_v_is_diff"] is True
        assert rec["eps"] == rec["wt_field"] == rec["pred"] == rec["v0"]


def test_W1_p_eq_3():
    B = theorem_B_W1_p_eq_3()
    assert B["proved"] is True
    assert B["W1_p_eq_3"] is True
    assert B["W1_all_odd_p"] is False


def test_two_fiber_pattern():
    C = theorem_C_two_fiber()
    assert C["proved"] is False
    assert C["two_fiber_in_W0"] is True
    assert C["rows"]["3"]["eps"] == 0
    assert C["rows"]["5"]["eps"] == 1
    assert C["rows"]["7"]["eps"] == 0
    for rec in C["rows"].values():
        assert rec["in_W0"] is True


def test_W2_miss_and_open():
    D = theorem_D_W2_miss()
    assert D["proved"] is False
    assert D["named_Dspan_misses_g"] is True
    assert D["walsh_general_p"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15614 import main

    out = main()
    assert out["proved"]["W1_p_eq_3"] is True
    assert out["proved"]["W1_all_odd_p"] is False
    assert out["proved"]["walsh_general_p"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
