"""Tests for Prop 15.601 — QR in rowspan(S) or S+ℓ; Walsh still open."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15601 import theorem_A_pencil, theorem_B_on_H0, theorem_C_open


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_pencil_congruence():
    A = theorem_A_pencil()
    assert A["proved"] is True
    assert A["rows"]["5"]["p_mod_4"] == 1 and A["rows"]["5"]["match"]
    assert A["rows"]["7"]["p_mod_4"] == 3 and A["rows"]["7"]["match"]
    assert A["rows"]["11"]["wt_w"] == 6
    assert A["rows"]["13"]["wt_w"] == 7


def test_QR_on_H0_not_independent():
    B = theorem_B_on_H0()
    assert B["proved"] is True
    # p=5 ≡1: QR+ell in S, QR extra; p=7 ≡3: QR in S
    assert B["rows"]["5"]["rQE"] == B["rows"]["5"]["rS"]
    assert B["rows"]["5"]["rQ"] == B["rows"]["5"]["rS"] + 1
    assert B["rows"]["7"]["rQ"] == B["rows"]["7"]["rS"]
    assert B["rows"]["7"]["rQE"] == B["rows"]["7"]["rS"] + 1
    for p in ("5", "7", "11", "13"):
        assert B["rows"][p]["ell_extra"] is True


def test_walsh_not_closed():
    C = theorem_C_open()
    assert C["proved"] is False
    assert C["walsh_general_p"] is False
    assert C["single_orbit_spans_W"] is False
    assert C["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip_406E():
    from e1_gmin_m4_prop15601 import main

    out = main()
    assert out["proved"]["pencil_identity"] is True
    assert out["proved"]["QR_on_H0"] is True
    assert out["proved"]["walsh_general_p"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
