import math

import numpy as np
import pytest

from src.e1_gmin_m4_prop15639 import (
    circle_point_signed_count,
    dense_profile_moment_differences,
    dense_profile_moment_factors,
    dependency_certificate,
    direction_count,
    equality_active_profile_counts,
    family_coordinate_signatures,
    first_nonminimal_odd_is_fourth_norm,
    first_nonminimal_odd_norm,
    first_nonminimal_odd_scaled_norm,
    first_nonminimal_odd_shell_classified,
    first_nonminimal_odd_shell_signed_count,
    first_nonminimal_odd_shell_signed_count_closed,
    first_nonminimal_odd_shell_theorem,
    n_of,
    negative_triangle_signed_count,
    p11_first_nonminimal_odd_exact_audit,
    scaled_coordinate_square_sum,
    transformed_t2_scaled_norm,
    unit_coordinate_forced,
)
from src.minmax_quadratic import paley_conference_prime_power


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23, 29, 31, 43])
def test_first_nonminimal_odd_norm_and_unit_coordinate_reduction(p):
    assert first_nonminimal_odd_scaled_norm(p) == 3 * p - 6
    assert first_nonminimal_odd_norm(p).numerator == 3 * p - 6
    assert first_nonminimal_odd_norm(p).denominator == 2 * p
    assert scaled_coordinate_square_sum(p) == 2 * p * (3 * p - 6)
    assert scaled_coordinate_square_sum(p) < 9 * n_of(p)
    assert unit_coordinate_forced(p)


def test_scaled_norm_3p_minus_6_is_fourth_only_at_two_smallest_primes():
    assert first_nonminimal_odd_is_fourth_norm(11)
    assert first_nonminimal_odd_is_fourth_norm(13)
    assert not first_nonminimal_odd_is_fourth_norm(17)
    assert not first_nonminimal_odd_is_fourth_norm(19)


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23, 29, 31])
def test_t1_equality_has_only_two_active_counts(p):
    assert equality_active_profile_counts(p) == (1, direction_count(p) - 2)


def test_dense_profile_moment_factorization():
    for values in ((0, 1, 2), (3, -2, 7), (11, 4, -5), (-8, 9, 13)):
        assert dense_profile_moment_differences(
            *values
        ) == dense_profile_moment_factors(*values)


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23])
def test_dense_equality_maps_to_three_already_classified_even_shells(p):
    assert [transformed_t2_scaled_norm(p, r) for r in range(3)] == [
        2 * (p - 1),
        2 * (p + 1),
        2 * (p + 3),
    ]


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23, 29])
def test_first_nonminimal_odd_shell_family_counts(p):
    n = n_of(p)
    assert negative_triangle_signed_count(p) == math.comb(n, 3)
    assert circle_point_signed_count(p) == p * p * (p - 1) * n
    assert (
        first_nonminimal_odd_shell_signed_count(p)
        == first_nonminimal_odd_shell_signed_count_closed(p)
    )
    assert first_nonminimal_odd_shell_signed_count(p) == (
        p * p * (p - 1) * (p + 7) * n // 6
    )


def test_p11_conference_has_half_negative_triangles():
    p = 11
    C = np.rint(paley_conference_prime_power(p)).astype(np.int64)
    total = math.comb(len(C), 3)
    signed_triangle_sum = int(np.trace(C @ C @ C)) // 6
    negative = (total - signed_triangle_sum) // 2
    assert signed_triangle_sum == 0
    assert 2 * negative == negative_triangle_signed_count(p) == total


def test_p11_exact_qfminim_audit_matches_complete_s27_count():
    audit = p11_first_nonminimal_odd_exact_audit()
    assert audit["scaled_bound"] == 28
    assert audit["signed_cumulative_count"] == 473_970
    assert audit["maximum_scaled_norm"] == 27
    assert audit["scaled_norm_3p_minus_6_shell_signed_count"] == 442_860
    assert audit["predicted_two_family_signed_count"] == 442_860
    assert audit["exact_count_matches_classification"]


@pytest.mark.parametrize("p", [11, 13, 17, 19, 23])
def test_coordinate_signatures_make_families_disjoint(p):
    cert = family_coordinate_signatures(p)
    assert cert["large_magnitude"] == p - 2
    assert cert["negative_triangle_large_coordinates"] == 3
    assert cert["circle_point_large_coordinates"] == 1
    assert cert["other_coordinate_magnitudes_at_most"] == 3
    assert cert["separated"]


@pytest.mark.parametrize("p", [11, 13, 17, 19])
def test_shell_dependencies_and_classification(p):
    assert all(dependency_certificate(p).values())
    assert first_nonminimal_odd_shell_classified(p)


def test_theorem_keeps_harmonic_tail_open():
    theorem = first_nonminimal_odd_shell_theorem()
    assert theorem["proved"]
    assert "complete scaled-norm 3p-6 shell" in theorem["scope"]
    assert all(row["checks"] for row in theorem["rows"].values())
