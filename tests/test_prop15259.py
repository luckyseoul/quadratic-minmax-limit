"""Tests for Prop 15.259 — transport lemma and π-covariance reduction."""
from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15259 import (
    hinge_status_259,
    residual_i_closed_via_259,
    theorem_pi_type_transport,
    theorem_reduction_to_pi_covariance,
    theorem_transport_lemma,
)


def test_theorems_proved():
    assert theorem_transport_lemma()["proved"] is True
    assert theorem_reduction_to_pi_covariance()["proved"] is True
    assert theorem_pi_type_transport()["proved"] is True
    assert theorem_reduction_to_pi_covariance()["hypothesis_proved_general"] is False


def test_predicates_open():
    assert residual_i_closed_via_259() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
    h = hinge_status_259()
    assert h["pi_covariance_general"] is False
    assert h["residual_i_closed"] is False
