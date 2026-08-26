import math

import numpy as np

from scripts.p7_size_four_slack_classify import (
    classify_four_odd_fibres_phase_one,
    classify_four_odd_fibres_phase_zero,
)
from src.e1_gmin_m4_prop15654 import (
    CERTIFICATE_ARCHIVE_SHA256,
    CERTIFICATE_AUDIT_SHA256,
    ORBIT_SOURCE_SHA256,
    p7_nonsquare_signed_permutation,
    p7_saturated_four_finite_certificate,
    theorem_p7_saturated_four_finite_exclusion,
)
from src.minmax_quadratic import paley_conference_prime_power


def test_complete_saturated_slack_catalogs():
    phase_zero = classify_four_odd_fibres_phase_zero()
    phase_one = classify_four_odd_fibres_phase_one()
    assert phase_zero["proved"] is True
    assert phase_zero["survivor_count"] == 1
    assert phase_zero["unique_formula"] == "A(X)=(|X cap B|-2)^2"
    assert phase_one["proved"] is True
    assert phase_one["survivor_count"] == 36
    assert len({tuple(row["slack_values"]) for row in phase_one["catalog"]}) == 36


def test_nonsquare_signed_permutation_is_exact_anti_isometry():
    result = p7_nonsquare_signed_permutation()
    assert result["fixes_distinguished_edge"] is True
    assert result["signed_conference_anti_isometry"] is True
    assert sorted(result["finite_permutation"]) == list(range(49))
    assert sorted(result["vertex_permutation"]) == list(range(50))


def test_nonsquare_map_flips_product_for_a_29_edge_four_finite_boundary():
    result = p7_nonsquare_signed_permutation()
    permutation = result["vertex_permutation"]
    C = np.rint(paley_conference_prime_power(7)).astype(np.int8)
    edges = {(0, 1), (0, 2), (3, 4)}
    for start in (5, 18):
        cycle = tuple(range(start, start + 13))
        edges.update(
            tuple(sorted((cycle[index], cycle[(index + 1) % len(cycle)])))
            for index in range(len(cycle))
        )
    assert len(edges) == 29
    degrees = [0] * 50
    for a, b in edges:
        degrees[a] += 1
        degrees[b] += 1
    assert tuple(index for index, degree in enumerate(degrees) if degree & 1) == (
        1,
        2,
        3,
        4,
    )
    image_edges = {
        tuple(sorted((permutation[a], permutation[b]))) for a, b in edges
    }
    original_product = math.prod(int(C[a, b]) for a, b in edges)
    image_product = math.prod(int(C[a, b]) for a, b in image_edges)
    assert image_product == -original_product


def test_certificate_counts_and_theorem_scope():
    certificate = p7_saturated_four_finite_certificate()
    assert certificate["saturated_boundaries_per_product_sign"] == 58800
    assert certificate["saturated_orbits_per_product_sign"] == 1225
    assert certificate["fixed_boundary_infeasible"] == 1225
    assert certificate["unknown"] == certificate["feasible"] == 0
    assert certificate["remaining_unsaturated_boundaries_per_product_sign"] == 23520
    assert certificate["remaining_unsaturated_orbits_per_product_sign"] == 518
    theorem = theorem_p7_saturated_four_finite_exclusion()
    assert theorem["proved"] is True
    assert theorem["p7_four_finite_doubly_saturated_both_product_signs"] == "CLOSED"
    assert theorem["p7_four_finite_unsaturated"] == "OPEN"
    assert theorem["p7_all_four_finite_points"] == "OPEN"
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
