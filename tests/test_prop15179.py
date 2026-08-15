"""Prop 15.179 — residual (ii) freeze-to-tight; residual (i) still open."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15171 import deep_s2_freeness_fail_k_ge_3p_ND_closed
from e1_gmin_m4_prop15179 import (
    dual_twolevel_affine_impossible_k_ge_3p,
    first_moment_forces_k,
    freeze_S_H_equals_3,
    hinge_status_179,
    main,
    residual_ii_dual_twolevel_affine_closed,
    theorem_freeze_all_primes,
)


def test_freeze_and_first_moment():
    assert freeze_S_H_equals_3()["proved"] is True
    for p in (5, 7, 11, 13, 17, 19, 23):
        fm = first_moment_forces_k(p)
        assert fm["k_forced"] == 3 * p - 1
        assert fm["proved"] is True
        imp = dual_twolevel_affine_impossible_k_ge_3p(p)
        assert imp["proved"] is True
        assert imp["no_k_ge_3p_satisfies_first_moment"] is True


def test_theorem_all_primes():
    th = theorem_freeze_all_primes([5, 7, 11, 13, 17, 19, 23, 29, 31])
    assert th["proved"] is True


def test_residual_ii_affine_closed_without_gsum():
    """Affine branch closed by freeze (no Gsum); full residual (ii) still open."""
    assert gsum_disj_lb_proved_general() is False  # residual (i) still open
    assert residual_ii_dual_twolevel_affine_closed() is True
    # Full residual (ii) closed by 15.179+236+237 (ND, not exhaustiveness)
    assert deep_s2_freeness_fail_k_ge_3p_ND_closed() is True
    assert e1_closed_general() is True
    h = hinge_status_179()
    assert h["bound_proved_general"] is False
    out = main()
    assert out["proved"]["L_closed"] is False
    assert out["proved"]["residual_i_closed"] is False
    assert out["proved"]["dual_twolevel_affine_impossible_k_ge_3p"] is True
    assert out["proved"]["residual_ii_dual_twolevel_affine_closed"] is True
