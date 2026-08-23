"""Tests for Prop 15.618 — Φ=ε; s_N pullback; 1_M coprime to g."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15618 import (
    theorem_A_phi_eq_eps,
    theorem_B_pullback,
    theorem_C_one_M,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_phi_eq_eps():
    A = theorem_A_phi_eq_eps()
    assert A["proved"] is True
    for rec in A["rows"].values():
        assert rec["gamma_wt_is_2p"] is True
        assert rec["gamma_0"] == 1
        assert rec["scale_1"] is True
        assert rec["zDz_match"] is True
        assert rec["QR_wt_DI"] == 0
        assert rec["QNR_wt_DI"] == 0
    assert A["rows"]["5"]["sN_eps"] == 1
    assert A["rows"]["5"]["sN_QR_odd"] == 1
    assert A["rows"]["5"]["sN_QNR_odd"] == 0


def test_pullback():
    B = theorem_B_pullback()
    assert B["pullback_p_law"] is True
    assert B["orbit_pattern_p_law"] is False
    assert B["orbit_pattern_certified"] is True
    assert B["eps_p_law"] is False
    assert B["proved"] is False
    for rec in B["rows"].values():
        assert rec["mismatches"] == 0
        assert rec["qnr_odd_all_even"] is True
        assert rec["f_dot_qr_odd"] == 1


def test_one_M():
    C = theorem_C_one_M()
    assert C["proved"] is True
    assert C["W2_p_law"] is False
    assert C["gcd_X1_g_is_1"] is True


def test_open():
    D = theorem_D_open()
    assert D["walsh_general_p"] is False
    assert D["W1_p_eq_1"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15618 import main

    out = main()
    assert out["proved"]["phi_eq_eps"] is True
    assert out["proved"]["sN_eps_p_law"] is False
    assert out["proved"]["walsh_general_p"] is False
    assert out["L_status"] == "OPEN"
