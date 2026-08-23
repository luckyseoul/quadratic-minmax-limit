"""Tests for Prop 15.620 — s_N not W1 p-law; χ_p-pullback Φ3-dead."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15620 import (
    theorem_A_sN_not_plaw,
    theorem_B_qr_chi,
    theorem_C_chi_miss,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_sN_killed():
    A = theorem_A_sN_not_plaw()
    assert A["proved"] is True
    assert A["sN_eps_p_law"] is False
    assert A["rows"]["5"]["phi"] == 1
    assert A["rows"]["29"]["phi"] == 0


def test_qr_chi_pattern():
    B = theorem_B_qr_chi()
    assert B["certified"] is True
    assert B["proved"] is False
    assert B["rows"]["5"]["odd_iff_QR"] is True


def test_chi_pullback_miss():
    C = theorem_C_chi_miss()
    assert C["proved"] is True
    assert C["phi3_divides"] is True
    assert C["eps"] == 0
    assert C["W2_p_law"] is False


def test_open_and_main():
    D = theorem_D_open()
    assert D["W1_p_eq_1"] is False
    assert 3 in D["stay_hit_p5"]
    assert len(D["stay_hit_p29"]) >= 1
    from e1_gmin_m4_prop15620 import main

    out = main()
    assert out["proved"]["sN_killed_p29"] is True
    assert out["proved"]["sN_eps_p_law"] is False
    assert out["proved"]["walsh_general_p"] is False
    assert out["L_status"] == "OPEN"
