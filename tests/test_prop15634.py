from fractions import Fraction

import numpy as np
import pytest

from src.e1_gmin_m4_prop15634 import (
    audit_circle_scheme,
    circle_design_parameters,
    circle_evaluation_operator_spectrum,
    circle_graph_spectrum,
    circle_operator_theorem,
    circle_tensor_gram_spectrum,
    common_two_secant_neighbours,
    intersection_valencies,
    n_of,
    projected_tensor_gram_entry,
    second_shadow_negative_definite,
    second_shadow_scalar_offset,
    second_shadow_spectrum,
    square_circle_count,
    z_dimension,
)


@pytest.mark.parametrize("p", [5, 7, 11, 13, 17, 19])
def test_closed_circle_design_and_spectrum_counts(p):
    D = circle_design_parameters(p)
    assert D["blocks"] == square_circle_count(p)
    assert D["blocks"] * D["block_size"] == D["points"] * D["replication"]
    assert (
        D["blocks"] * D["block_size"] * (D["block_size"] - 1)
        == D["points"] * (D["points"] - 1) * D["pair_multiplicity"]
    )
    assert sum(intersection_valencies(p).values()) + 1 == square_circle_count(p)
    assert sum(row["multiplicity"] for row in circle_graph_spectrum(p)) == square_circle_count(p)
    assert sum(row["multiplicity"] for row in circle_tensor_gram_spectrum(p)) == square_circle_count(p)
    assert sum(row["multiplicity"] for row in circle_evaluation_operator_spectrum(p)) == z_dimension(p)


def test_common_neighbour_counts_are_the_matrix_identity_entries():
    for p in (5, 7, 11, 13):
        alpha = Fraction(p * p - 1, 8)
        beta = Fraction((p - 1) ** 2 * (p + 1), 8)
        for j in (0, 1, 2):
            expected = alpha * j + beta - (p if j == 2 else 0)
            assert expected.denominator == 1
            assert common_two_secant_neighbours(p, j) == expected


def test_projected_tensor_gram_entries_match_known_p11_values():
    p = 11
    assert projected_tensor_gram_entry(p, p + 1) == 11858
    assert projected_tensor_gram_entry(p, 2) == Fraction(-605, 3)
    assert projected_tensor_gram_entry(p, 1) == Fraction(-2299, 30)
    assert projected_tensor_gram_entry(p, 0) == Fraction(1452, 5)


def test_circle_gram_nonzero_values_and_multiplicities():
    for p in (5, 7, 11):
        n = n_of(p)
        rows = circle_tensor_gram_spectrum(p)
        assert rows == [
            {"eigenvalue": 0, "multiplicity": n},
            {
                "eigenvalue": p**3 * (p - 1),
                "multiplicity": n * (p - 1) // 4,
            },
            {
                "eigenvalue": p**3 * (p + 1),
                "multiplicity": n * (p - 3) // 4,
            },
        ]


def test_full_second_shadow_closed_forms_and_sign_transition():
    for p in (5, 7, 11, 13, 17, 19):
        a = second_shadow_scalar_offset(p)
        rows = second_shadow_spectrum(p)
        assert rows[0]["eigenvalue"] == -Fraction(
            (p + 2) * (p * p - 4 * p + 1), 4 * p * (p * p + 5)
        )
        assert rows[1]["eigenvalue"] == -Fraction(
            p**3 - 3 * p * p - 19 * p + 9, 8 * p * (p * p + 5)
        )
        assert rows[2]["eigenvalue"] == -Fraction(
            p**3 - 5 * p * p - 19 * p - 1, 8 * p * (p * p + 5)
        )
        assert rows[0]["eigenvalue"] == a
        assert second_shadow_negative_definite(p) == (p >= 11)


@pytest.mark.parametrize("p", [3, 5, 7])
def test_explicit_circle_scheme_exact_small_prime_audit(p):
    record = audit_circle_scheme(p)
    assert record["checks"]
    assert record["quadratic_adjacency_identity_max_error"] == 0


def test_theorem_keeps_r1_open():
    theorem = circle_operator_theorem()
    assert theorem["proved"]
    assert all(theorem["rows"][str(p)]["checks"] for p in (5, 7, 11, 13, 17, 19))
    assert not second_shadow_negative_definite(7)
    assert second_shadow_negative_definite(11)
