import itertools

import numpy as np

from scripts.p7_no_infinity_unsaturated_cpsat import direction_target_options
from scripts.p7_unsaturated_mod7_batch import mapped_slack_catalog
from scripts.p7_unsaturated_modular_catalog_filter import (
    equation_matrix,
    left_dependencies,
    scan_catalog,
)
from src.e1_gmin_m4_prop15655 import (
    CERTIFICATE_ARCHIVE_SHA256,
    CERTIFICATE_AUDIT_SHA256,
    CERTIFICATE_BATCH_SHA256,
    ORBIT_SOURCE_SHA256,
    p7_unsaturated_mod7_certificate,
    theorem_p7_unsaturated_four_finite_exclusion,
)


def test_mod7_common_matrix_has_audited_left_nullspace():
    matrix = equation_matrix()
    rank, dependencies = left_dependencies(matrix, 7)
    assert matrix.shape == (282, 1225)
    assert rank == 147
    assert dependencies.shape == (135, 282)
    assert np.all(dependencies @ (matrix % 7) % 7 == 0)


def test_direct_slack_relabelling_matches_interpolated_target_catalog():
    B = {1, 3, 5, 6}
    direct = mapped_slack_catalog(4, 1, 14, B)
    options = direction_target_options(4, 1, B, 32, False)
    pairs = tuple(itertools.combinations(range(7), 2))
    reconstructed = []
    for option in options:
        constant = int(option[1])
        coefficients = tuple(int(value) for value in option[2:])
        values = []
        for X in itertools.combinations(range(7), 4):
            X_set = set(X)
            target = constant + sum(
                coefficients[index]
                * (1 if ((s in X_set) == (t in X_set)) else -1)
                for index, (s, t) in enumerate(pairs)
            )
            assert target >= 3 and target % 2 == 1
            values.append((target - 3) // 2)
        reconstructed.append(tuple(values))
    assert len(direct) == len(reconstructed) == 36
    assert {tuple(row) for row in direct} == set(reconstructed)


def test_hard_orbit_145_catalog_is_completely_rejected_mod7():
    result = scan_catalog(
        -1,
        (1, 9, 24, 34),
        (1,),
        1,
        {},
        (7,),
    )
    assert result["catalog_total"] == 1764
    assert result["moduli"] == [
        {
            "modulus": 7,
            "rank": 147,
            "dependency_dimension": 135,
            "passing_catalog_rows": 0,
            "cumulative_survivors": 0,
        }
    ]
    assert result["surviving_catalog_rows"] == 0


def test_certificate_counts_and_theorem_scope():
    certificate = p7_unsaturated_mod7_certificate()
    assert certificate["unsaturated_boundaries"] == 23520
    assert certificate["unsaturated_orbits"] == 518
    assert certificate["fixed_elevation_cases"] == 2408
    assert certificate["catalog_tuples_excluded"] == 1716742440
    assert certificate["mod7_infeasible_cases"] == 2408
    assert certificate["surviving_cases"] == 0
    assert certificate["independent_audit"] is True

    theorem = theorem_p7_unsaturated_four_finite_exclusion()
    assert theorem["proved"] is True
    assert theorem["p7_four_finite_unsaturated_both_product_signs"] == "CLOSED"
    assert theorem["p7_all_four_finite_points_with_prop15654"] == "CLOSED"
    assert theorem["p7_all_size_four_with_prop15653"] == "CLOSED"
    assert theorem["closes_all_p7_size_four"] is True
    assert theorem["p5_size_four"] == "OPEN"
    assert theorem["closes_residual_ii"] is False
    assert theorem["L_status"] == "OPEN"
    assert all(
        len(value) == 64
        for value in (
            CERTIFICATE_ARCHIVE_SHA256,
            CERTIFICATE_BATCH_SHA256,
            CERTIFICATE_AUDIT_SHA256,
            ORBIT_SOURCE_SHA256,
        )
    )
