"""Tests for Prop 15.615 — two-fiber W1-1 false; L2 not closed."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15615 import (
    theorem_A_two_fiber_kill,
    theorem_B_stay_exists_p_eq_1,
    theorem_C_W2_named_pool,
    theorem_D_leftover2,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_two_fiber_kill():
    A = theorem_A_two_fiber_kill()
    assert A["proved"] is True
    assert A["biconditional_false"] is True
    assert A["rows"]["5"]["eps"] == 1
    assert A["rows"]["17"]["eps"] == 0


def test_stay_exists_not_p_law():
    B = theorem_B_stay_exists_p_eq_1()
    assert B["proved"] is False
    assert B["W1_p_eq_1"] is False
    assert B["rows"]["5"]["n_stay_eps1"] > 0


def test_W2_named_pool_miss():
    C = theorem_C_W2_named_pool()
    assert C["proved"] is False
    assert C["W2_p_law"] is False
    assert C["p11_named_Dspan_misses_g"] is True


def test_leftover2_not_closed():
    D = theorem_D_leftover2()
    assert D["proved"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False
    assert D["p5_k20_leftover_splus_empty"] is True


def test_main_does_not_flip():
    from e1_gmin_m4_prop15615 import main

    out = main()
    assert out["proved"]["two_fiber_biconditional"] is False
    assert out["proved"]["W1_p_eq_1"] is False
    assert out["proved"]["walsh_general_p"] is False
    assert out["proved"]["residual_ii"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
