import itertools

from scripts.p7_size_four_slack_classify import classify_three_odd_fibres
from src.e1_gmin_m4_prop15653 import (
    CERTIFICATE_ARCHIVE_SHA256,
    CERTIFICATE_AUDIT_SHA256,
    ORBIT_SOURCE_SHA256,
    p7_inf_three_orbit_certificate,
    reconstructed_direction_coefficients,
    target_coefficients,
    theorem_p7_infinity_three_exclusion,
)


def test_three_odd_fibre_slack_classification_is_exact():
    result = classify_three_odd_fibres()
    assert result["degree_at_most_two_rank"] == 21
    assert result["left_kernel_dimension"] == 14
    assert result["sparse_correction_candidates"] == 630
    assert result["survivor_count"] == 1
    assert result["unique_formula"] == "A(X)=(|X cap B|-2)^2"
    assert result["proved"] is True


def test_target_polynomials_are_the_forced_slacks():
    for B in ({2}, {0, 2, 5}):
        constant, linear, pairs = target_coefficients(B)
        for chosen in itertools.combinations(range(7), 4):
            z = tuple(1 if s in chosen else -1 for s in range(7))
            value = constant + sum(linear[s] * z[s] for s in range(7)) + sum(
                coefficient * z[s] * z[t]
                for (s, t), coefficient in pairs.items()
            )
            intersection = len(set(chosen) & B)
            slack = intersection if len(B) == 1 else (intersection - 2) ** 2
            assert value == 3 + 2 * slack


def test_coefficient_kernel_reconstructs_target_on_the_middle_slice():
    star_counts = (0, 1, 2, 0, 1, 0, 3)
    for B in ({1}, {0, 3, 6}):
        target_constant, target_linear, target_pairs = target_coefficients(B)
        for kernel_parameter in (-2, 0, 5):
            parallel, linear, cross = reconstructed_direction_coefficients(
                B, star_counts, kernel_parameter
            )
            for chosen in itertools.combinations(range(7), 4):
                z = tuple(1 if s in chosen else -1 for s in range(7))
                reconstructed = parallel + sum(
                    linear[s] * z[s] for s in range(7)
                ) + sum(
                    coefficient * z[s] * z[t]
                    for (s, t), coefficient in cross.items()
                )
                target = target_constant + sum(
                    target_linear[s] * z[s] for s in range(7)
                ) + sum(
                    coefficient * z[s] * z[t]
                    for (s, t), coefficient in target_pairs.items()
                )
                assert reconstructed == target


def test_orbit_certificate_and_theorem_scope():
    coverage = p7_inf_three_orbit_certificate()
    assert coverage["boundary_orbits"] == 416
    assert coverage["orbit_size_sum"] == 18424
    assert coverage["fixed_boundary_infeasible"] == 416
    assert coverage["unknown"] == coverage["feasible"] == 0
    theorem = theorem_p7_infinity_three_exclusion()
    assert theorem["proved"] is True
    assert theorem["p7_infinity_plus_three_finite_both_product_signs"] == "CLOSED"
    assert theorem["p7_four_finite_points"] == "OPEN"
    assert theorem["p5_size_four"] == "OPEN"
    assert theorem["closes_all_p7_size_four"] is False
    assert theorem["closes_residual_ii"] is False
    assert all(
        len(value) == 64
        for value in (
            CERTIFICATE_ARCHIVE_SHA256,
            CERTIFICATE_AUDIT_SHA256,
            ORBIT_SOURCE_SHA256,
        )
    )
