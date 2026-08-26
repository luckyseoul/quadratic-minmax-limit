import math

import numpy as np

from scripts.p5_size_four_full_shell_mod7_exception import build_problem_data
from src.e1_gmin_m4_prop15656 import (
    CERTIFICATE_ARCHIVE_SHA256,
    CERTIFICATE_AUDIT_SHA256,
    MOD5_BATCH_SHA256,
    MOD7_EXCEPTION_SHA256,
    ORBIT_SOURCE_SHA256,
    p5_full_shell_four_point_certificate,
    p5_nonsquare_signed_permutation,
    theorem_p5_four_point_exclusion,
)
from src.minmax_quadratic import paley_conference_prime_power


def test_full_shell_exception_matrix_and_parity_are_rebuilt_exactly():
    data = build_problem_data()
    assert data["matrix"].shape == (132, 325)
    assert data["rank"] == 67
    assert data["dependencies"].shape == (65, 132)
    assert np.all(data["dependencies"] @ (data["matrix"] % 7) % 7 == 0)
    assert data["parity_mass"] == 56
    assert data["lift_mass"] == 11


def test_p5_nonsquare_map_is_an_exact_signed_anti_isometry():
    result = p5_nonsquare_signed_permutation()
    assert result["fixes_distinguished_edge"] is True
    assert result["signed_conference_anti_isometry"] is True
    assert sorted(result["finite_permutation"]) == list(range(25))
    assert sorted(result["vertex_permutation"]) == list(range(26))


def test_nonsquare_map_flips_product_for_21_edges_and_four_finite_boundary():
    result = p5_nonsquare_signed_permutation()
    permutation = result["vertex_permutation"]
    C = np.rint(paley_conference_prime_power(5)).astype(np.int8)
    edges = {(0, 1), (0, 2), (3, 4)}
    for start in (5, 14):
        cycle = tuple(range(start, start + 9))
        edges.update(
            tuple(sorted((cycle[index], cycle[(index + 1) % len(cycle)])))
            for index in range(len(cycle))
        )
    assert len(edges) == 21
    degrees = [0] * 26
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
    certificate = p5_full_shell_four_point_certificate()
    assert certificate["floor_surviving_orbit_product_sign_cases"] == 1202
    assert certificate["floor_surviving_boundary_product_sign_cases"] == 26450
    assert certificate["direct_mod5_infeasible_orbits"] == 712
    assert certificate["direct_mod7_infeasible_orbits"] == 1
    assert certificate["transferred_no_infinity_orbits"] == 489
    assert certificate["unknown"] == certificate["feasible"] == 0

    theorem = theorem_p5_four_point_exclusion()
    assert theorem["proved"] is True
    assert theorem["p5_all_size_four_boundaries_with_prop15632"] == "CLOSED"
    assert theorem["closes_all_size_four_for_p_at_least_5"] is True
    assert theorem["first_open_boundary_size_for_p_at_least_5"] == 6
    assert theorem["closes_residual_ii"] is False
    assert theorem["closes_R1"] is False
    assert theorem["L_status"] == "OPEN"
    hashes = (
        CERTIFICATE_ARCHIVE_SHA256,
        CERTIFICATE_AUDIT_SHA256,
        MOD7_EXCEPTION_SHA256,
        *MOD5_BATCH_SHA256.values(),
        *ORBIT_SOURCE_SHA256.values(),
    )
    assert all(len(value) == 64 for value in hashes)
