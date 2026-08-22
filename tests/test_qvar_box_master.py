"""Aut-inv master+box cannot prove QVAR: LP min < 0 at p=5, true pairing > 0."""
from __future__ import annotations

import inspect

import e1_gmin_global_qvar as G
import e1_gmin_leftover1_qvar_principal as L
import e1_gmin_qvar_box_master as B
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15275 import type_I_multilevel_bad_case_ND_closed
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general


def test_torb_is_a_true_quotient_at_p5():
    lp = B.aut_box_master_lp(5)
    assert lp["Torb_kappa_err"] < 1e-10
    assert lp["kap_orbit_std"] < 1e-10
    assert lp["n_orbits"] >= 3
    assert lp["ker_dim"] >= 1
    assert lp["lp_ok"]


def test_lp_min_is_negative_true_pairing_positive():
    T = B.theorem_box_master_aut_cannot_prove_qvar()
    assert T["proved"]
    assert T["inequality_proved"] is False
    assert T["claim_lp_min_ge_0"] is False
    assert T["p5_lp_min_is_neg285_over_4"]
    assert T["p5_lp"]["min_pairing"] < 0
    assert T["p5_true"]["pairing"] > 0
    assert T["p5_true"]["clears_QVAR"]
    assert T["p5_true"]["pairing"] > T["p5_lp"]["min_pairing"]


def test_fail_when_claim_box_master_suffices():
    T = B.theorem_box_master_aut_cannot_prove_qvar()
    # Fail-eq: "LP min ≥ 0 would prove QVAR by linear 4-point theory."
    assert T["claim_lp_min_ge_0"] is False
    assert T["p5_lp"]["min_pairing"] < -1.0


def test_does_not_flip_qvar_or_leftovers():
    assert G.global_qvar_proved_general() is False
    assert L.leftover1_qvar_and_principal_proved() is False
    assert phi_F_ge_6_proved_general() is False
    assert residual_ii_k_eq_4p_empty() is False
    assert type_I_multilevel_bad_case_ND_closed() is False
    src = inspect.getsource(G.global_qvar_proved_general)
    assert "return True" not in src
    src_b = inspect.getsource(B.theorem_box_master_aut_cannot_prove_qvar)
    assert "inequality_proved" in src_b
