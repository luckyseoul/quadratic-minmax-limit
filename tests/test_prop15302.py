"""Tests for Prop 15.302 — torus decomposition of S_k."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general
from e1_gmin_m4_prop15290 import live_Q_on_T
from e1_gmin_m4_prop15302 import (
    S_from_torus,
    F_star_named,
    coarse_Q_means,
    form_factor,
    live_S,
    main,
    nontrivial_s,
    prove_decomposition,
    prove_floor_rewrite,
    prove_Fstar_named,
    prove_open,
    sigma_k,
)


def test_sigma_values():
    assert sigma_k(5, 0) == 1  # (5-3)/2
    assert sigma_k(5, 4) == 1
    assert sigma_k(5, 2) == -1
    assert sigma_k(7, 0) == 2
    assert sigma_k(7, 12) == 2
    assert sigma_k(7, 8) == -1
    assert sigma_k(7, 2) == -1


def test_decomposition_recovers_S():
    A = prove_decomposition()
    assert A["proved"] is True
    assert A["flip_hit"] >= 1
    assert A["s1_hit"] >= 1
    for p in (5, 7):
        F = _kernel_field(p)
        Q = live_Q_on_T(p)
        rec = S_from_torus(F, Q, 2)
        live = live_S(F, Q, 2)
        assert abs(rec - live) < 1e-8


def test_flip_sigma_breaks_HD():
    F = _kernel_field(7)
    Q = live_Q_on_T(7)
    live = live_S(F, Q, 8)
    bad = S_from_torus(F, Q, 8, flip_sigma=True)
    assert abs(bad - live) > 0.5


def test_s_equals_1_breaks():
    F = _kernel_field(5)
    Q = live_Q_on_T(5)
    live = live_S(F, Q, 2)
    bad = S_from_torus(F, Q, 2, s=1)
    assert abs(bad - live) > 0.5


def test_Fstar_named_in_Q():
    B = prove_Fstar_named()
    assert B["proved"] is True
    assert B["wrong_hit"] >= 1
    F = _kernel_field(7)
    Q = live_Q_on_T(7)
    means = coarse_Q_means(F, Q)
    s = nontrivial_s(F)
    live = form_factor(F, Q, s, 2).real
    named = F_star_named(7, 2, means["++"], means["N*"])
    assert abs(live - named) < 1e-8
    assert abs(named - 2 * (means["++"] - means["N*"])) < 1e-8


def test_floor_is_dQ_circulant():
    C = prove_floor_rewrite()
    assert C["proved"] is True
    assert C["mins"]["5"]["ge6"] is True
    assert C["mins"]["7"]["ge6"] is True
    assert C["phi_F_imported"] is False


def test_floor_still_open():
    D = prove_open()
    assert D["proved"] is False
    assert phi_F_ge_6_proved_general() is False
    assert rhat_ge_budget_proved_general() is False
    assert D["S_k_named_in_p"] is False
    assert D["Q_tau_named_in_p"] is False


def test_flags_untouched():
    assert e1_closed_general() is True
    assert gsum_disj_lb_proved_general() is False
    assert phi_F_ge_6_proved_general() is False


def test_main():
    out = main()
    assert out["proved"]["torus_decomposition"] is True
    assert out["proved"]["Fstar_named_in_Q"] is True
    assert out["proved"]["floor_is_even_circulant_of_dQ"] is True
    assert out["proved"]["phi_F_ge_6_proved_general"] is False
    assert out["proved"]["S_k_named_in_p"] is False
