"""Tests for Prop 15.264 — scheme Ext eq; Ext[δ] eigen; ExtΠ=ΠExt."""
from fractions import Fraction

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general, is_prime
from e1_gmin_m4_prop15264 import (
    hinge_status_264,
    mu_part_ab,
    residual_i_closed_via_264,
    theorem_delta_ext_eigen,
    theorem_ext_commutes_with_Pi,
    theorem_mupart_ext_equation,
)


def test_theorems_proved():
    assert theorem_mupart_ext_equation()["proved"] is True
    assert theorem_mupart_ext_equation()["all_primes_5_to_200"] is True
    assert theorem_delta_ext_eigen()["proved"] is True
    assert theorem_ext_commutes_with_Pi()["proved"] is True


def test_coeff_identity():
    for p in (5, 7, 11, 13, 17, 19, 23):
        assert is_prime(p)
        a, b = mu_part_ab(p)
        n = p * p + 1
        lhs = -6 * a + 2 * n * b - b * (p**4 - 1)
        assert lhs == 2


def test_predicates_open():
    assert residual_i_closed_via_264() is False
    assert gsum_disj_lb_proved_general() is False
    assert e1_closed_general() is True
    h = hinge_status_264()
    assert h["pi_odd_kernel_trivial"] is False
    assert h["residual_i_closed"] is False
