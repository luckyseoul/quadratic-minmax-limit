"""Tests for Prop 15.607 — W G_aff-irred all odd p; Walsh still open."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15607 import (
    theorem_A_Fp_mixes,
    theorem_B_irred,
    theorem_C_maxminus_span,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_Fp_preserves_WH_not_simple_at_p7():
    A = theorem_A_Fp_mixes()
    assert A["proved"] is True
    assert A["rows"]["3"]["WH_simple_Cp"] is True
    assert A["rows"]["5"]["WH_simple_Cp"] is True
    assert A["rows"]["7"]["WH_simple_Cp"] is False
    assert A["rows"]["7"]["p7_cubic_Fp_fills_WH"] is True
    for rec in A["rows"].values():
        assert rec["Fp_preserves_WH"] is True
        assert rec["Fp_mixes_sq_nsq"] is False


def test_W_irred_all_odd_p():
    B = theorem_B_irred()
    assert B["proved"] is True
    assert B["W_irreducible_all_odd_p"] is True
    assert B["H0_quotient_irreducible"] is True


def test_Maxminus_span_not_walsh():
    C = theorem_C_maxminus_span()
    assert C["proved"] is True
    assert C["dir_affine_span_Maxminus_is_H0"] is True
    assert C["walsh_general_p"] is False


def test_walsh_open():
    D = theorem_D_open()
    assert D["proved"] is False
    assert D["walsh_general_p"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15607 import main

    out = main()
    assert out["proved"]["W_irreducible_all_odd_p"] is True
    assert out["proved"]["dir_affine_span_Maxminus_is_H0"] is True
    assert out["proved"]["walsh_general_p"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
