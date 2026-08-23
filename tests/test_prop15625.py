"""Tests for Prop 15.625 — W1 eighth-interval / (2/p)_4."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15625 import (
    theorem_A_eighth,
    theorem_B_residual,
    theorem_C_class_exhaustive,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_eighth_interval():
    A = theorem_A_eighth()
    assert A["proved"] is True
    assert A["W1_all_odd_p"] is False
    assert A["rows"]["17"]["eps"] == 1
    assert A["rows"]["241"]["eps"] == 1
    assert A["rows"]["409"]["eps"] == 1
    assert A["rows"]["601"]["eps"] == 0


def test_class_and_open():
    C = theorem_C_class_exhaustive()
    assert C["proved"] is True
    assert C["generation_gap"] is False
    B = theorem_B_residual()
    assert B["p601_eps"] == 0
    D = theorem_D_open()
    assert D["walsh_general_p"] is False
    from e1_gmin_m4_prop15625 import main

    out = main()
    assert out["proved"]["W1_quartic2_minus"] is True
    assert out["proved"]["W1_all_odd_p"] is False
    assert out["L_status"] == "OPEN"
