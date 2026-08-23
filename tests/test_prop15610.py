"""Tests for Prop 15.610 — pair-stabilizer uniqueness dead; Walsh open."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15610 import (
    theorem_A_W0_perp,
    theorem_B_I_preserves_flag,
    theorem_C_uniqueness_dead,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_W0_is_extra_perp():
    A = theorem_A_W0_perp()
    assert A["proved"] is True
    for rec in A["rows"].values():
        assert rec["dot_equals_ev0"] is True
        assert rec["extra_at_0"] == 0
        assert rec["dim_W0"] == rec["N_minus_1"]


def test_I_preserves_ker2():
    B = theorem_B_I_preserves_flag()
    assert B["proved"] is True
    for rec in B["rows"].values():
        assert rec["I_fixes_extra"] is True
        assert rec["I_preserves_ker2"] is True
        assert rec["dim_ker1_W0"] == 1
        assert rec["dim_ker2_W0"] == 2


def test_uniqueness_dead():
    C = theorem_C_uniqueness_dead()
    assert C["proved"] is True
    assert C["W0_quotient_irreducible"] is False
    assert C["pair_stabilizer_uniqueness_dead"] is True
    for rec in C["rows"].values():
        assert rec["four_divides_N"] is True
        assert rec["ker2_proper"] is True


def test_walsh_open():
    D = theorem_D_open()
    assert D["proved"] is False
    assert D["walsh_general_p"] is False
    assert D["pair_stabilizer_uniqueness_dead"] is True
    assert D["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15610 import main

    out = main()
    assert out["proved"]["W0_is_extra_perp"] is True
    assert out["proved"]["pair_stabilizer_uniqueness_dead"] is True
    assert out["proved"]["W0_quotient_irreducible"] is False
    assert out["proved"]["walsh_general_p"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
