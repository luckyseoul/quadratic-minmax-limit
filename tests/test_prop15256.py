"""Tests for Prop 15.256 — global ∑μ and extension (κ,φ) sums."""
from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15256 import (
    E4_balanced,
    f4_extension_sum,
    hinge_status_256,
    residual_i_closed_via_256,
    sum_mu_closed_form,
    theorem_extension_sum_kappa,
    theorem_extension_sum_phi,
    theorem_f4_extension_nonzero,
    theorem_global_sum_mu,
)


def test_theorems_proved():
    assert theorem_global_sum_mu()["proved"] is True
    assert theorem_extension_sum_kappa()["proved"] is True
    assert theorem_extension_sum_phi()["proved"] is True
    assert theorem_f4_extension_nonzero()["proved"] is True


def test_sum_mu_formula():
    for p in (5, 7, 11, 13):
        sm = sum_mu_closed_form(p)
        assert sm == Fraction(-(p * p) * (p * p - 1), 12)
        n = p * p + 1
        E4 = E4_balanced(p)
        assert E4 == p**4 + 6 * p * p + 1
        assert E4 == n + 3 * n * (n - 1) + 24 * sm


def test_f4_extension_nonzero():
    for p in (5, 7, 11, 13, 17):
        s = f4_extension_sum(p)
        assert s != 0
        assert s == Fraction(9 - p * p, p * (p * p + 1))


def test_predicates_open():
    assert residual_i_closed_via_256() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_256()
    assert h["vanishing_sum_general"] is False
    assert h["residual_i_closed"] is False
