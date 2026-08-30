"""Prop 15.207: ker=sc reduction; residual_i stays open."""
from __future__ import annotations

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15207 import (
    certify_Gplus_pd_on_subspace,
    certify_ker_Qplus_eq_scheme,
    free_e_bound_proved_general,
    ker_sc_proved_general,
    residual_i_closed_via_207,
    theorem_ker_sc_reduction,
    theorem_scheme_plusplus_formula,
    theorem_wick_identity_zero_diag_plusplus,
)


def test_scheme_plusplus_formula():
    assert theorem_scheme_plusplus_formula()["proved"] is True


def test_wick_identity():
    assert theorem_wick_identity_zero_diag_plusplus()["proved"] is True


def test_ker_sc_reduction():
    assert theorem_ker_sc_reduction()["proved"] is True


def test_Gplus_pd_certified_p3_p5():
    c = certify_Gplus_pd_on_subspace([3, 5])
    assert c["certified"] is True
    assert c["by_p"]["3"]["min_eig_Gplus"] > 7.9
    assert c["by_p"]["5"]["min_eig_Gplus"] > 3.0


def test_ker_Qplus_eq_scheme_p3_p5():
    c = certify_ker_Qplus_eq_scheme([3, 5])
    assert c["certified"] is True
    assert c["by_p"]["5"]["ker_Qplus_dim"] == 25


def test_ker_sc_via_gplus():
    assert ker_sc_proved_general() is True
    assert free_e_bound_proved_general() is False
    assert residual_i_closed_via_207() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is False
