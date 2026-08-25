from fractions import Fraction

import pytest

from src.e1_gmin_m4_prop15635 import (
    odd_nonminimum_scaled_floor,
    p11_third_shell_audit,
    third_norm,
    third_pair_harmonic_coefficient,
    third_pair_signed_count,
    third_scaled_norm,
    third_shell_theorem,
)


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23, 29])
def test_uniform_third_norm_and_odd_gap(p):
    assert third_scaled_norm(p) == 2 * (p + 1)
    assert third_norm(p) == Fraction(p + 1, p)
    assert odd_nonminimum_scaled_floor(p) == 3 * p - 6
    assert odd_nonminimum_scaled_floor(p) > third_scaled_norm(p)


@pytest.mark.parametrize("p", [11, 13, 17, 19])
def test_pair_count_and_harmonic_scalar(p):
    assert third_pair_signed_count(p) == p * p * (p * p + 1)
    assert third_pair_harmonic_coefficient(p) == -Fraction(
        p * p + 4 * p - 3, 4 * (p * p + 5)
    )
    assert third_pair_harmonic_coefficient(p) < 0


def test_exact_p11_cumulative_count_exhausts_the_third_shell():
    audit = p11_third_shell_audit()
    assert audit["signed_cumulative_count"] == 31_110
    assert audit["third_shell_signed_count"] == 14_762
    assert audit["predicted_pair_signed_count"] == 14_762
    assert audit["complete_third_shell_is_pair_orbit"]


def test_scope_does_not_soft_close_r1_or_all_prime_shell_classification():
    theorem = third_shell_theorem()
    assert theorem["proved"]
    assert "complete third-shell classification additionally at p=11" in theorem["scope"]
