"""Tests for Prop 15.260 — ∑μκ, ∑μφ, layer counts."""
from fractions import Fraction
from math import comb

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general, n_of
from e1_gmin_m4_prop15260 import (
    hinge_status_260,
    layer_counts,
    residual_i_closed_via_260,
    sum_mu,
    sum_mu_kappa,
    sum_mu_phi,
    theorem_layer_counts,
    theorem_sum_mu_kappa,
    theorem_sum_mu_phi,
)


def test_theorems_proved():
    assert theorem_sum_mu_kappa()["proved"] is True
    assert theorem_sum_mu_phi()["proved"] is True
    assert theorem_layer_counts()["proved"] is True


def test_formulas():
    for p in (5, 7, 11, 13):
        n = n_of(p)
        assert sum_mu_phi(p) == n * sum_mu(p)
        assert sum_mu_kappa(p) == Fraction(n * p * p * (p * p - 1), 8)
        n1, n3 = layer_counts(p)
        B = comb(n, 4)
        sk2 = n * (n - 1) * (n - 2) * (n - 5) // 8
        assert n1 + n3 == B
        assert n1 + 9 * n3 == sk2


def test_predicates_open():
    assert residual_i_closed_via_260() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_260()
    assert h["pi_covariance_general"] is False
    assert h["residual_i_closed"] is False
