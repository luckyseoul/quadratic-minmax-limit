"""Prop 15.193 — residual (ii) exhaustiveness audit; full residual (ii) OPEN."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15171 import deep_s2_freeness_fail_k_ge_3p_ND_closed
from e1_gmin_m4_prop15179 import residual_ii_dual_twolevel_affine_closed
from e1_gmin_m4_prop15193 import (
    desired_exhaustiveness_lemma_statement,
    freeness_fail_forces_f_e_plus_on_U2,
    freeness_fail_kills_Max_plus_weak_ND_witness,
    freeness_fail_possible_when_k_ge_3p,
    hinge_status_193,
    main,
    multi_level_not_ruled_out,
    non_affine_two_level_not_ruled_out,
    residual_ii_affine_branch_closed,
    residual_ii_exhaustiveness_proved,
    residual_ii_full_closed,
    residual_ii_open_subcases,
)


def test_proved_pieces():
    assert freeness_fail_forces_f_e_plus_on_U2()["proved"] is True
    assert freeness_fail_kills_Max_plus_weak_ND_witness()["proved"] is True
    for p in (5, 7, 11, 13, 17):
        assert freeness_fail_possible_when_k_ge_3p(p)["proved"] is True


def test_exhaustiveness_not_proved():
    assert residual_ii_exhaustiveness_proved() is False
    lemma = desired_exhaustiveness_lemma_statement()
    assert lemma["proved"] is False
    E = multi_level_not_ruled_out()
    assert E["exhaustiveness_proved"] is False
    assert E["multi_level_mass_feasible_examples_exist"] is True
    F = non_affine_two_level_not_ruled_out()
    assert F["exhaustiveness_proved"] is False
    assert F["non_affine_S_H_not_const"] is True


def test_affine_closed_full_closed():
    """Affine (15.179) + (ii-b) 15.236 + (ii-a) 15.237 ⇒ residual (ii) full."""
    assert residual_ii_dual_twolevel_affine_closed() is True
    assert residual_ii_affine_branch_closed() is True
    assert residual_ii_full_closed() is True
    assert deep_s2_freeness_fail_k_ge_3p_ND_closed() is True
    opens = residual_ii_open_subcases()
    assert not any("ii-a" in s for s in opens)
    assert not any("ii-b" in s for s in opens)


def test_e1_still_open():
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_193()
    assert h["bound_proved_general"] is False
    assert h["residual_ii_full"] == "CLOSED"
    out = main()
    assert out["proved"]["residual_ii_full_closed"] is True
    assert out["proved"]["residual_ii_affine_branch_closed"] is True
    assert out["proved"]["residual_ii_exhaustiveness_proved"] is False
    assert out["proved"]["L_closed"] is False
    assert out["proved"]["E1_closed"] is False
