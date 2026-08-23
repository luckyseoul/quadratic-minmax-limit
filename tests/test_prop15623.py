"""Tests for Prop 15.623 — W1 p≡73 or 97 (mod 120)."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15623 import (
    theorem_A_W1_73_97,
    theorem_B_open_class,
    theorem_C_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_W1_73_97():
    A = theorem_A_W1_73_97()
    assert A["proved"] is True
    assert A["W1_all_odd_p"] is False
    assert A["rows"]["73"]["eps"] == 1
    assert A["rows"]["97"]["eps"] == 1
    assert A["rows"]["241"]["eps"] == 0


def test_open_and_main():
    B = theorem_B_open_class()
    assert B["proved"] is False
    assert B["d_minus_3_eps_p241"] == 0
    C = theorem_C_open()
    assert C["walsh_general_p"] is False
    from e1_gmin_m4_prop15623 import main

    out = main()
    assert out["proved"]["W1_p_eq_73_or_97_mod_120"] is True
    assert out["proved"]["W1_all_odd_p"] is False
    assert out["L_status"] == "OPEN"
