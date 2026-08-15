"""Prop 15.205: M_cand≤1/(2p); ratio algebra; residual_i stays open."""
from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15205 import (
    M_cand,
    free_e_bound_proved_general,
    free_e_scheme_ratio_candidate,
    mu_bound_proved_general,
    residual_i_closed_via_205,
    theorem_Mcand_le_mu4_thr,
    theorem_ratio_le_three_halves,
    certify_ratio_at_small_p,
    certify_mu_le_Mcand_small,
)


def test_Mcand_le_1_over_2p():
    r = theorem_Mcand_le_mu4_thr()
    assert r["proved"] is True
    assert M_cand(5) == Fraction(3, 65)
    assert M_cand(5) <= Fraction(1, 10)


def test_ratio_le_three_halves():
    r = theorem_ratio_le_three_halves()
    assert r["proved"] is True
    assert free_e_scheme_ratio_candidate(5) == Fraction(41, 28)
    assert free_e_scheme_ratio_candidate(7) == Fraction(163, 113)


def test_ratio_matches_census_p5_p7():
    c = certify_ratio_at_small_p()
    assert c["certified"] is True
    assert c["by_p"]["5"]["match"] is True
    assert c["by_p"]["7"]["match"] is True


def test_mu_le_Mcand_census():
    c = certify_mu_le_Mcand_small()
    assert c["certified"] is True
    assert c["by_p"]["5"]["sharp"] is True


def test_predicates_remain_open():
    assert free_e_bound_proved_general() is False
    assert mu_bound_proved_general() is False
    assert residual_i_closed_via_205() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
