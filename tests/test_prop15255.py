"""Tests for Prop 15.255 — m4± expansion on ∞-sets."""
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15255 import (
    hinge_status_255,
    residual_i_closed_via_255,
    theorem_m4_eq_reduction,
    theorem_m4minus_expansion_infty,
    theorem_m4plus_expansion_infty,
    theorem_reflection_as_sum_bound,
)


def test_expansions_proved():
    assert theorem_m4plus_expansion_infty()["proved"] is True
    assert theorem_m4minus_expansion_infty()["proved"] is True


def test_reduction_criterion():
    t = theorem_m4_eq_reduction()
    assert t["proved"] is True
    assert t["hypothesis_proved_general"] is False


def test_reflection_sum_criterion():
    t = theorem_reflection_as_sum_bound()
    assert t["proved_criterion"] is True
    assert t["hypothesis_proved_general"] is False


def test_predicates_open():
    assert residual_i_closed_via_255() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_255()
    assert h["vanishing_sum_general"] is False
    assert h["residual_i_closed"] is False
