"""Tests for Prop 15.604 — QR/QNR in H0; Walsh still open."""
from __future__ import annotations

from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty, multilevel_ND_k_ge_4p_proved
from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15604 import (
    theorem_A_qr_in_h0,
    theorem_B_dilation_fixed,
    theorem_C_restriction_census,
    theorem_D_open,
)


def test_flags_untouched():
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert multilevel_ND_k_ge_4p_proved() is False
    assert type_I_aut_e_3AB_positive_general() is False


def test_qr_qnr_by_p_mod_4():
    A = theorem_A_qr_in_h0()
    assert A["proved"] is True
    for p_s, rec in A["rows"].items():
        p = int(p_s)
        if p % 4 == 1:
            assert rec["qr_in_H0"] is True
            assert rec["qnr_in_H0"] is False
        else:
            assert p % 4 == 3
            assert rec["qr_in_H0"] is False
            assert rec["qnr_in_H0"] is True
        assert rec["chi_Fp_one"] is True
        assert rec["e0_einf_in_H0"] is False
        assert rec["thru0_L_cap_QR"] == p - 1
        assert rec["off0_L_cap_QR"] == (p - 1) // 2


def test_fail_swapped_congruence():
    """The swapped law 1_QR in H0 at p≡3 is false."""
    A = theorem_A_qr_in_h0(primes=(3, 7, 11))
    for rec in A["rows"].values():
        assert rec["qr_in_H0"] is False
        assert rec["qnr_in_H0"] is True


def test_dilation_fixed_dim2():
    B = theorem_B_dilation_fixed()
    assert B["proved"] is True
    for rec in B["rows"].values():
        assert rec["D_N_is_id"] is True
        assert rec["dim_ker_DminusI_cap_H0"] == 2
        assert rec["dim_ker_DminusI_cap_H0"] != 1
        assert rec["dim_ker_DminusI_cap_H0"] != 4
        assert rec["extra_in_H0"] is True
        assert rec["span_card"] == 4  # 2-dim F2-space has 4 vectors


def test_restriction_not_a_theorem():
    C = theorem_C_restriction_census()
    assert C["proved"] is False
    assert C["surjective_in_general"] is False
    assert C["onto_at_listed_primes"] is False
    for rec in C["rows"].values():
        assert rec["ker_to_QR"] >= 2
        assert rec["rank_to_QR"] < rec["n_QR"]


def test_irred_and_walsh_open():
    D = theorem_D_open()
    assert D["proved"] is False
    assert D["walsh_general_p"] is False
    assert D["H0_quotient_irreducible"] is False
    assert D["restriction_QR_surjective"] is False
    assert D["residual_ii_k_eq_4p_empty"] is False


def test_main_does_not_flip():
    from e1_gmin_m4_prop15604 import main

    out = main()
    assert out["proved"]["qr_qnr_in_H0"] is True
    assert out["proved"]["dilation_fixed_dim2"] is True
    assert out["proved"]["restriction_QR_surjective"] is False
    assert out["proved"]["H0_quotient_irreducible"] is False
    assert out["proved"]["walsh_general_p"] is False
    assert out["walsh_15_406_E"] == "OPEN"
    assert out["L_status"] == "OPEN"
