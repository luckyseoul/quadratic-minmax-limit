"""Leftover 1 QVAR+principal hinge: identities fail-when-wrong; flags are live."""
from __future__ import annotations

import inspect
from fractions import Fraction

import pytest

import e1_gmin_leftover1_qvar_principal as L
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty
from e1_gmin_m4_prop15275 import type_I_multilevel_bad_case_ND_closed
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15589 import (
    lambda_exc_from_quartic_variance,
    quartic_variance_floor_threshold,
    spherical_quartic_variance,
)


def test_reductions_live_and_estimate_flags_are_imported():
    A = L.theorem_A_qvar_iff()
    B = L.theorem_B_spherical_exceeds_qvar()
    C = L.theorem_C_principal_room_reduction()
    D = L.theorem_D_floor_iff_m4_pairing()
    assert A["proved"]
    assert B["proved"]
    assert C["proved"]
    assert D["proved"]
    assert D["pointwise_q2_ge_3By2"] is False
    assert A["qvar_general"] is False  # identity A does not prove the estimate
    assert C["principal_moment_general"] is False  # room formula ≠ moment bound
    assert L.leftover1_reductions_ok()
    leftover1 = L.leftover1_qvar_and_principal_proved()
    assert leftover1 is (
        L.global_qvar_proved_general()
        and L.principal_delta_room_moment_proved()
    )
    src_and = inspect.getsource(L.leftover1_qvar_and_principal_proved)
    assert "global_qvar_proved_general" in src_and
    assert "qvar_k_ge_7_proved_general" not in src_and
    assert phi_F_ge_6_proved_general() is leftover1


def test_phi_F_imports_both_blocks_not_handwritten_true():
    import inspect

    assert phi_F_ge_6_proved_general() is L.leftover1_qvar_and_principal_proved()
    src = inspect.getsource(phi_F_ge_6_proved_general)
    assert "leftover1_qvar_and_principal_proved" in src
    assert "return True" not in src


@pytest.mark.parametrize("p", [5, 7, 11, 13, 17, 19, 23])
def test_FWW_qvar_threshold_drop_16(p):
    thr = quartic_variance_floor_threshold(p)
    assert L.qvar_threshold_wrong_drop16(p) != thr
    assert 6 * L.hs_norm_A_psi(p) == thr


@pytest.mark.parametrize("p", [5, 7])
def test_FWW_lambda_exc_32_not_16(p):
    rec = L.theorem_E_exceptional_quartic_variance()["by_p"][str(p)]
    var = Fraction(rec["E_abs_Zpsi_sq"])
    got = lambda_exc_from_quartic_variance(p, var)
    assert got == Fraction(rec["lambda_exc"])
    assert L.lambda_exc_wrong_16(p, var) != got


@pytest.mark.parametrize("p", [5, 7, 11, 13, 17, 19])
def test_FWW_spherical_drop_q5(p):
    vs = spherical_quartic_variance(p)
    assert L.V_sph_wrong_drop_q5(p) != vs
    assert vs > quartic_variance_floor_threshold(p)


@pytest.mark.parametrize("p", [5, 7, 11, 13, 17, 19])
def test_FWW_principal_room_n14_not_n6(p):
    from e1_gmin_m4_prop15100 import n_of
    from e1_gmin_m4_prop15589 import (
        delta2_room_principal,
        delta2_room_principal_after_exception,
    )

    n = n_of(p)
    new = delta2_room_principal_after_exception(p)
    old = delta2_room_principal(p)
    assert new / old == Fraction(n - 6, n - 14)
    assert new != old
    assert L.crude_Es4_2n3(p) > L.principal_Es4_budget_after_exception(p)


def test_other_leftovers_are_live_units_not_baked_false():
    from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
    from e1_gmin_m4_prop15275 import type_I_aut_e_3AB_positive_general
    from e1_main_chain_status import four_e1_units_closed, run_main_chain

    dump = L.dump_leftover_predicates()
    assert dump["gsum_disj_lb_proved_general"] is gsum_disj_lb_proved_general()
    assert dump["type_I_aut_e_3AB_positive_general"] is (
        type_I_aut_e_3AB_positive_general()
    )
    assert dump["phi_F_ge_6_proved_general"] is phi_F_ge_6_proved_general()
    assert dump["residual_ii_k_eq_4p_empty"] is residual_ii_k_eq_4p_empty()
    assert dump["type_I_multilevel_bad_case_ND_closed"] is (
        type_I_multilevel_bad_case_ND_closed()
    )
    units = four_e1_units_closed()
    out = run_main_chain()
    expect = "CLOSED" if units["closed"] else "OPEN"
    assert dump["L_status"] == expect
    assert out["L_status"] == expect
    # live e1 AND is the old incomplete wiring — not this leftover close
    assert e1_closed_general() is True
