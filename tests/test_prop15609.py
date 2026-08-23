"""Tests for Prop 15.609 — I(H0)=H0; Walsh still open."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15609 import (
    theorem_A_dual_equals_nsq,
    theorem_B_no_mixed_tangency,
    theorem_C_I_preserves_H0,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_H0prime_is_dual():
    A = theorem_A_dual_equals_nsq()
    assert A["proved"] is True
    for rec in A["rows"].values():
        assert rec["rank_S"] == rec["n_over_2"]
        assert rec["rank_Sprime"] == rec["n_over_2"]
        assert rec["S_perp_Sprime"] is True
        assert all(c == 2 for c in rec["sample_caps"])


def test_no_mixed_tangency():
    B = theorem_B_no_mixed_tangency()
    assert B["proved"] is True
    for rec in B["rows"].values():
        assert rec["infty_opposite_caps"] == [2]
        assert 1 not in rec["inverted_vs_nsq_caps"]
        assert set(rec["inverted_vs_nsq_caps"]) <= {0, 2}


def test_I_H0_is_p_law():
    C = theorem_C_I_preserves_H0()
    assert C["proved"] is True
    assert C["H0_invariance_p_law"] is True
    for rec in C["rows"].values():
        assert rec["H0_preserved"] is True
        assert rec["off0_in_rowspan_S"] is True
        assert rec["off0_in_rowspan_Sprime"] is False
        assert rec["off0_is_row"] is False
        assert rec["dim_H0"] == rec["n_over_2"]


def test_walsh_open():
    D = theorem_D_open()
    assert D["proved"] is False
    assert D["walsh_general_p"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15609 import main

    out = main()
    assert out["proved"]["H0prime_equals_H0_perp"] is True
    assert out["proved"]["I_preserves_H0"] is True
    assert out["proved"]["walsh_general_p"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
