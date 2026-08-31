import json
from pathlib import Path

from e1_gmin_m4_prop15742 import (
    EXPECTED_ROW_DIGESTS,
    SIX_DILATE_CUTS,
    aggregate_dependencies_certificate,
    independent_cpsat_energy_exclusion,
    proposition_15742,
    row_catalog_certificate,
)


ROOT = Path(__file__).resolve().parents[1]


def _energy(row):
    return sum(value * value for value in row)


def _assert_row_is_in_relaxation(row, *, total, l1_bound, energy_bound, cut_bound):
    assert sum(row) == total
    assert sum(map(abs, row)) <= l1_bound
    assert _energy(row) <= energy_bound
    assert sum((distance + 1) ** 2 * row[distance] for distance in range(6)) % 13 == 0
    assert all(
        sum(coefficient * value for coefficient, value in zip(cut, row))
        <= cut_bound
        for cut in SIX_DILATE_CUTS
    )


def test_imported_M2_cut_energy_and_parseval_dependencies_are_live():
    row = aggregate_dependencies_certificate()
    assert row["p"] == 13
    assert row["four_exact_stars_force_global_M2_zero"] is True
    assert row["M2_dependency"] == {
        "exact_star_certificate": True,
        "degree_2_four_roots_force_zero": True,
    }
    assert row["signed_total_T_over_h"] == 17
    assert row["row_sums"] == {"elevated": 11, "opposite": -20}
    assert row["l1_bounds"] == {"elevated": 53, "opposite": 56}
    assert "59-P_L" in row["l1_derivation"]
    assert row["cut_bounds"] == {"elevated": 91, "opposite": -130}
    assert len(row["six_interval_dilate_cuts"]) == 6
    assert len({tuple(cut) for cut in row["six_interval_dilate_cuts"]}) == 6
    assert row["six_cut_set_matches_15_741"] is True
    assert row["six_cuts_in_full_74_catalog"] is True
    assert row["cut_image_is_even"] is True
    assert row["cut_column_sums"] == [42] * 6
    assert row["slacks"] == {
        "elevated": {
            "definition": "y=91*1-Cq",
            "parity_and_sign": "odd positive",
            "sum": 84,
        },
        "opposite": {
            "definition": "y=-130*1-Cq",
            "parity_and_sign": "even nonnegative",
            "sum": 60,
        },
    }
    assert row["imported_integer_energy_bounds"] == {
        "elevated": 86,
        "opposite": 106,
    }
    assert row["row_parameters_equal_imported_live_values"] is True
    assert row["imported_nonstar_parseval"] == "707+26*C"
    assert row["collision_parameter_nonnegative"] is True
    assert row["proved"] is True


def test_recursive_elevated_catalog_has_sharp_energy_thirty_one():
    row = row_catalog_certificate("elevated")
    assert row["coordinate_range_from_sum_and_energy"] == [-5, 9]
    assert row["pre_cut_row_count"] == 5844
    assert row["surviving_row_count"] == 30
    assert row["surviving_rows_sha256"] == EXPECTED_ROW_DIGESTS["elevated"]
    assert "superset of realizable" in row["catalog_scope"]
    assert row["prior_energy_cap_redundant_by_independent_cpsat"] is True
    assert row["sharp_energy_maximum"] == 31
    assert row["maximizer_count"] == 6
    assert row["maximizers"] == [
        [0, 3, 1, 4, 1, 2],
        [1, 1, 3, 2, 0, 4],
        [1, 2, 4, 0, 3, 1],
        [2, 0, 1, 3, 4, 1],
        [3, 4, 2, 1, 1, 0],
        [4, 1, 0, 1, 2, 3],
    ]
    assert len(row["surviving_rows"]) == 30
    for candidate in row["surviving_rows"]:
        _assert_row_is_in_relaxation(
            candidate,
            total=11,
            l1_bound=53,
            energy_bound=86,
            cut_bound=91,
        )
    assert max(map(_energy, row["surviving_rows"])) == 31
    assert row["proved"] is True


def test_recursive_opposite_catalog_has_sharp_energy_eighty_two():
    row = row_catalog_certificate("opposite")
    assert row["coordinate_range_from_sum_and_energy"] == [-9, 2]
    assert row["pre_cut_row_count"] == 1704
    assert row["surviving_row_count"] == 24
    assert row["surviving_rows_sha256"] == EXPECTED_ROW_DIGESTS["opposite"]
    assert "superset of realizable" in row["catalog_scope"]
    assert row["prior_energy_cap_redundant_by_independent_cpsat"] is True
    assert row["sharp_energy_maximum"] == 82
    assert row["maximizer_count"] == 6
    assert row["maximizers"] == [
        [-6, -1, -4, -2, -4, -3],
        [-4, -4, -1, -3, -6, -2],
        [-4, -3, -2, -6, -1, -4],
        [-3, -6, -4, -1, -2, -4],
        [-2, -4, -6, -4, -3, -1],
        [-1, -2, -3, -4, -4, -6],
    ]
    assert len(row["surviving_rows"]) == 24
    for candidate in row["surviving_rows"]:
        _assert_row_is_in_relaxation(
            candidate,
            total=-20,
            l1_bound=56,
            energy_bound=106,
            cut_bound=-130,
        )
    assert max(map(_energy, row["surviving_rows"])) == 82
    assert row["proved"] is True


def test_independent_nineteen_variable_models_exclude_the_next_energies():
    elevated = independent_cpsat_energy_exclusion("elevated")
    opposite = independent_cpsat_energy_exclusion("opposite")
    assert elevated["forbidden_energy_floor"] == 32
    assert opposite["forbidden_energy_floor"] == 83
    assert elevated["status"] == opposite["status"] == "INFEASIBLE"
    assert elevated["infeasible"] is opposite["infeasible"] is True
    assert elevated["workers"] == opposite["workers"] == 1
    assert elevated["seed"] == opposite["seed"] == 0
    assert elevated["model_validation"] == opposite["model_validation"] == ""
    assert elevated["variables"] == opposite["variables"] == 19
    assert elevated["constraints"] == opposite["constraints"] == 22
    assert elevated["model_proto_sha256"] == (
        "557e700271596217961bf7f5a6db8107bc32dbdea718e2f62e0cbf4ad8765db3"
    )
    assert opposite["model_proto_sha256"] == (
        "72df1b51c8f369bce8d4133491a74c1a290cb6b878031e014ec7c2b3fc3b0603"
    )
    assert elevated["M2_quotient_domain"] == opposite["M2_quotient_domain"] == [
        -200,
        200,
    ]
    assert elevated["M2_quotient_absolute_bound_from_l1"] == 147
    assert opposite["M2_quotient_absolute_bound_from_l1"] == 156
    assert elevated["prior_energy_upper_constraint_used"] is False
    assert opposite["prior_energy_upper_constraint_used"] is False
    assert elevated["proved"] is opposite["proved"] is True


def test_package_closes_only_the_p13_fourth_shell():
    row = proposition_15742()
    assert row["prop"] == "15.742"
    assert row["result_status"] == "exhaustive finite certificate"
    assert row["remaining_hard_quotient_partition"] == [1, 1, 1, 1, 2, 2, 2]
    assert row["global_energy"] == {
        "three_elevated_upper": 93,
        "seven_opposite_upper": 574,
        "nonstar_upper": 667,
        "exact_parseval": "707+26*C",
        "collision_parameter_lower_bound": 0,
        "parseval_lower": 707,
        "gap": 40,
        "contradiction": True,
    }
    assert row["quartic_code_used"] is False
    assert row["root_quartet_case_split_used"] is False
    assert row["hard_sign_normalization_used"] is False
    assert row["directional_coefficient_matrix_census_used"] is False
    assert row["binary_midpoint_lift_used"] is False
    assert row["common_graph_exists"] is False
    assert row["p13_generic_four_exact_partition_closed"] is True
    assert row["p13_generic_t3_branch_closed"] is True
    assert row["p13_t3_exceptional_u3_closed_by_15_739"] is True
    assert row["p13_k_eq_58_closed"] is True
    assert row["generic_p_ge_17_t3_branch_closed"] is False
    assert row["k_eq_4p_plus_6_shell_closed"] is False
    assert row["residual_ii_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert "p=13 at k>=60" in row["remaining_scope"]
    assert row["proved"] is True
    evidence = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15742.json").read_text()
    )
    assert evidence == row
