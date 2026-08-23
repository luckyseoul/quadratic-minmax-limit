"""Tests for Prop 15.602 — G_aff^□ on H0; Walsh still open."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15602 import (
    theorem_A_affine_permutes,
    theorem_B_unique_line,
    theorem_C_inversion,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_affine_permutes_rows():
    A = theorem_A_affine_permutes()
    assert A["proved"] is True
    for p in ("3", "5", "7", "11"):
        assert A["rows"][p]["square_dil_permutes"]
        assert A["rows"][p]["trans_permutes"]
        assert A["rows"][p]["frob_permutes"]
        assert A["rows"][p]["nonsquare_image_is_row"] is False


def test_unique_1dim_invariant():
    B = theorem_B_unique_line()
    assert B["proved"] is True
    for p in ("3", "5", "7", "11"):
        assert B["rows"][p]["one_in_H0"]
        assert B["rows"][p]["dim_translation_invariants"] == 1


def test_inversion_pencil_not_all_rows():
    C = theorem_C_inversion()
    assert C["proved"] is True
    assert C["H0_invariance_p_law"] is False
    for p in ("3", "5", "7", "11"):
        assert C["rows"][p]["pencil_permutes"]
        assert C["rows"][p]["nonpencil_image_is_row"] is False
        assert C["rows"][p]["H0_preserved"]
        assert C["rows"][p]["n_pencil"] == C["rows"][p]["half"]


def test_walsh_and_irred_open():
    D = theorem_D_open()
    assert D["proved"] is False
    assert D["walsh_general_p"] is False
    assert D["H0_quotient_irreducible"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15602 import main

    out = main()
    assert out["proved"]["affine_permutes_rows"] is True
    assert out["proved"]["unique_1dim_invariant"] is True
    assert out["proved"]["H0_quotient_irreducible"] is False
    assert out["proved"]["walsh_general_p"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
