from fractions import Fraction

import pytest

from src.e1_gmin_m4_prop15636 import (
    complete_third_shell_theorem,
    exceptional_profile_contradiction,
    hasse_forced_power,
)


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23, 29, 31, 37, 41, 43])
def test_uniform_hasse_derivative_contradiction(p):
    m = (p - 1) // 2
    assert hasse_forced_power(p, m - 1) == (1, pow(3, -1, p))
    assert hasse_forced_power(p, m - 2) == (2, pow(5, -1, p))
    cert = exceptional_profile_contradiction(p)
    assert cert["forced_c"] == cert["one_third"]
    assert cert["forced_c_squared"] == cert["one_fifth"]
    assert cert["contradiction_residue"] != 0
    assert cert["equivalent_nonzero_numerator"] == 4
    assert cert["exceptional_profile_impossible"]


def test_complete_third_shell_all_prime_scope():
    theorem = complete_third_shell_theorem()
    assert theorem["proved"]
    assert theorem["scope"] == (
        "complete third dual shell for every odd prime p>=11"
    )
    for p_text, row in theorem["rows"].items():
        p = int(p_text)
        assert row["third_norm"] == str(Fraction(p + 1, p))
        assert row["complete_signed_count"] == p * p * (p * p + 1)
        assert Fraction(row["harmonic_scalar"]) < 0
        assert row["checks"]


def test_small_or_even_inputs_are_outside_stated_scope():
    for p in (2, 3, 5, 7):
        with pytest.raises(ValueError):
            exceptional_profile_contradiction(p)
