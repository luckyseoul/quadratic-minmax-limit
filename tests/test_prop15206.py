"""Prop 15.206: n3=(n-6)/4 Max+-free; residual_i stays open."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15206 import (
    free_e_bound_proved_general,
    mu_bound_proved_general,
    n1_formula_kappa1,
    n3_formula_kappa1,
    residual_i_closed_via_206,
    theorem_n3_degree,
    theorem_paley_divisibility,
    certify_n3_counts,
)


def test_n3_degree_proved():
    r = theorem_n3_degree()
    assert r["proved"] is True


def test_paley_divisibility():
    r = theorem_paley_divisibility()
    assert r["proved"] is True


def test_n3_n1_formulas():
    # p=5: n=26, n3=5, n1=17
    assert n3_formula_kappa1(5) == 5
    assert n1_formula_kappa1(5) == 17
    assert n3_formula_kappa1(7) == 11
    assert n1_formula_kappa1(7) == 35
    assert n3_formula_kappa1(3) == 1
    assert n1_formula_kappa1(3) == 5


def test_n3_census_p3_p5():
    c = certify_n3_counts(primes=[3, 5])
    assert c["certified"] is True
    assert c["by_p"]["3"]["bad"] == 0
    assert c["by_p"]["5"]["bad"] == 0


def test_predicates_remain_open():
    assert free_e_bound_proved_general() is False
    assert mu_bound_proved_general() is False
    assert residual_i_closed_via_206() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
