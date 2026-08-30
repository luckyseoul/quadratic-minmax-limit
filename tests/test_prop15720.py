"""Prop 15.720: direct degree-congruence bi-tight obstruction."""
from __future__ import annotations

import inspect
from fractions import Fraction

from e1_gmin_m4_prop15720 import (
    bitight_level_obstruction,
    centered_density,
    degree_modulus,
    level_2_arithmetic_obstruction_general,
    level_3_arithmetic_obstruction_general,
    level_4_arithmetic_obstruction_general,
    required_bitight_levels_empty_all_primes,
    scheme_coordinate_from_degree,
    theorem_degree_congruence,
    theorem_required_bitight_levels,
)


def test_centering_and_scheme_coordinate_are_exact():
    for p in (5, 7, 11, 13):
        n = p * p + 1
        for s in (2, 3, 4):
            assert centered_density(p, s) == Fraction(2 * s, n * p)
            assert scheme_coordinate_from_degree(p, s, 0) == -Fraction(
                s, p * (p * p - 1)
            )
            assert degree_modulus(p) == (p * p - 1) // 2


def test_level_2_tail_is_symbolic_not_a_prime_census():
    out = level_2_arithmetic_obstruction_general()
    assert out["proved"] is True
    assert out["base_polynomial_p5"] == 4


def test_level_3_tail_and_p5_exception_are_both_closed():
    out = level_3_arithmetic_obstruction_general()
    assert out["proved"] is True
    assert out["p5"]["empty"] is True
    assert out["p5"]["remainders"] == {0: 6, 1: 4}


def test_level_4_tail_and_small_residues_are_closed():
    out = level_4_arithmetic_obstruction_general()
    assert out["proved"] is True
    assert out["small"][5]["remainders"] == {0: 4, 1: 2}
    assert out["small"][7]["remainders"] == {0: 8, 1: 6}


def test_required_levels_close_for_general_p():
    assert theorem_degree_congruence()["proved"] is True
    assert required_bitight_levels_empty_all_primes() is True
    out = theorem_required_bitight_levels()
    assert out["proved"] is True
    assert out["required_levels"] == [2, 3]
    assert out["proved_bi_tight_levels"] == [2, 3, 4]
    assert out["all_bi_tight_levels_claimed"] is False
    assert out["one_sided_tight_level_4_claimed"] is False
    assert out["spectral_floor_used"] is False
    for p in (5, 7, 11, 101):
        for s in (2, 3, 4):
            assert bitight_level_obstruction(p, s)["bi_tight_empty"] is True


def test_level_4_is_not_mislabeled_as_a_required_or_one_sided_close():
    out = bitight_level_obstruction(5, 4)
    assert out["bi_tight_empty"] is True
    assert out["required_level"] is False
    assert out["level_4_corollary"] is True


def test_new_gate_does_not_import_the_retracted_spectral_predicate():
    src = inspect.getsource(required_bitight_levels_empty_all_primes)
    assert "15167" not in src
    assert "phi_F_ge_6" not in src
    assert "qvar" not in src.lower()
