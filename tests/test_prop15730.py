from copy import deepcopy

import pytest

from e1_gmin_m4_prop15730 import (
    LINE_TYPE_KEYS,
    bivariate_line_census,
    common_completion_row,
    direction_refinement_row,
    point_signature_rows,
    proposition_15730,
    verify_bivariate_line_census,
    verify_point_signatures,
)


def test_exact_three_by_three_census_for_p31_mixed_blocks():
    census = bivariate_line_census(31, 3)
    assert census == {
        "a0_t0": 452,
        "a0_t1": 26,
        "a0_t2": 42,
        "a1_t0": 42,
        "a1_t1": 200,
        "a1_t2": 0,
        "a2_t0": 224,
        "a2_t1": 4,
        "a2_t2": 3,
    }
    assert tuple(census) == LINE_TYPE_KEYS

    audit = verify_bivariate_line_census(31, 3, census)
    assert audit["actual"] == audit["expected"]
    assert all(audit["checks"].values())
    assert audit["proved"] is True


def test_every_admissible_block_count_in_both_residues_has_exact_moments():
    for p in (31, 41, 43, 47):
        R = (p - 1) // 3
        for y in range(R // 2 + 1):
            row = common_completion_row(p, y)
            census = row["line_type_census"]
            audit = row["line_type_census_audit"]
            assert set(census) == set(LINE_TYPE_KEYS)
            assert all(value >= 0 for value in census.values())
            assert audit["actual"]["total_projective_lines"] == p * p + p + 1
            assert audit["actual"]["A_point_line_incidences"] == (
                row["repair_arc_size_k"] * (p + 1)
            )
            assert audit["actual"]["T_point_line_incidences"] == R * (p + 1)
            assert audit["actual"]["A_T_cross_pair_moment"] == (
                row["repair_arc_size_k"] * R
            )
            assert audit["checks"]["D_occupancy_census"] is True
            assert audit["proved"] is True
            assert row["proved"] is True


def test_census_predicate_fails_when_individual_claims_are_corrupted():
    correct = bivariate_line_census(31, 3)

    wrong_cross = deepcopy(correct)
    wrong_cross["a1_t1"] += 1
    audit = verify_bivariate_line_census(31, 3, wrong_cross)
    assert audit["checks"]["A_T_cross_pair_moment"] is False
    assert audit["proved"] is False

    wrong_empty = deepcopy(correct)
    wrong_empty["a0_t0"] -= 1
    audit = verify_bivariate_line_census(31, 3, wrong_empty)
    assert audit["checks"]["total_projective_lines"] is False
    assert audit["proved"] is False

    missing = deepcopy(correct)
    del missing["a1_t2"]
    audit = verify_bivariate_line_census(31, 3, missing)
    assert audit["exact_keys"] is False
    assert audit["proved"] is False


def test_all_maximum_repairs_and_index_one_complements_are_simultaneous():
    row = common_completion_row(31, 3)
    assert row["R"] == 10
    assert row["repair_arc_size_k"] == 22
    assert row["trisecants_x"] == 4
    assert row["four_secants_y"] == 3
    assert row["rich_block_count"] == 7
    assert row["points_on_rich_blocks"] == 24
    assert row["singleton_points"] == 8

    family = row["repair_family"]
    assert family["maximum_D_subarc_count"] == 3**4 * 6**3 == 17_496
    assert family["maximum_D_subarc_count_formula"] == "3^x 6^y"
    assert family["these_are_all_maximum_D_subarcs"] is True

    complement = row["complement_family"]
    assert complement["complement_size"] == 10
    assert complement["complement_is_an_arc"] is True
    assert complement["every_complement_point_has_A_secant_index"] == 1
    assert complement["unique_secant_blocks_form_matching_on_A"] is True
    assert complement["unique_secant_fibre_counts"] == {1: 4, 2: 3}

    simultaneous = row["simultaneous_unique_trisecants"]
    assert simultaneous["count"] == 10
    assert simultaneous["size"] == 23
    assert simultaneous["all_affine"] is True
    assert simultaneous["each_has_exactly_one_trisecant"] is True


def test_point_signatures_reaggregate_to_every_line_type():
    rows = point_signature_rows(31, 3)
    assert [row["multiplicity"] for row in rows] == [8, 8, 6, 4, 6]
    assert rows[0]["incident_line_types"] == {
        "a2_t0": 21,
        "a1_t1": 10,
        "a1_t0": 1,
    }
    assert rows[1]["incident_line_types"] == {
        "a2_t1": 1,
        "a2_t0": 20,
        "a1_t1": 9,
        "a1_t0": 2,
    }
    assert rows[2]["incident_line_types"] == {
        "a2_t2": 1,
        "a2_t0": 20,
        "a1_t1": 8,
        "a1_t0": 3,
    }
    assert rows[3]["incident_line_types"] == {
        "a2_t1": 1,
        "a1_t1": 20,
        "a0_t2": 9,
        "a0_t1": 2,
    }
    assert rows[4]["incident_line_types"] == {
        "a2_t2": 1,
        "a1_t1": 20,
        "a0_t2": 8,
        "a0_t1": 3,
    }
    audit = verify_point_signatures(31, 3)
    assert audit["A_point_count"] == 22
    assert audit["T_point_count"] == 10
    assert all(audit["checks"].values())
    assert audit["proved"] is True


def test_cotangent_base_counts_are_scoped_to_D_and_one_fixed_repair():
    row = common_completion_row(31, 3)
    bases = row["cotangent_deletion_bases"]
    assert bases["scope"] == (
        "for each fixed maximum repair A; extension counts are within D"
    )
    assert bases["base_size"] == 21
    assert bases["one_within_D_extension_base_count"] == 8
    assert bases["two_within_D_cotangent_extension_base_count"] == 8
    assert bases["three_within_D_cotangent_extension_base_count"] == 6
    assert bases["all_one_point_deletion_bases"] == 22
    assert bases["displayed_D_extensions_individually_preserve_the_arc"] is True
    assert (
        bases["displayed_D_extensions_on_each_tangent_are_pairwise_incompatible"]
        is True
    )
    assert bases["additional_extension_points_outside_D_excluded"] is False


def test_direction_refinement_recovers_nonrich_b2_and_a_rich_profile():
    nonrich = direction_refinement_row(
        31,
        3,
        A_secants=8,
        T_secants=2,
        rich_trisecants=0,
        rich_four_secants=0,
        ordinary_A_T_lines=5,
    )
    assert nonrich["affine_line_type_census"] == {
        "a0_t0": 14,
        "a0_t1": 1,
        "a0_t2": 2,
        "a1_t0": 1,
        "a1_t1": 5,
        "a1_t2": 0,
        "a2_t0": 8,
        "a2_t1": 0,
        "a2_t2": 0,
    }
    assert nonrich["odd_D_fibres_b"] == 2
    assert nonrich["D_fibre_profile"] == {0: 14, 1: 2, 2: 15, 3: 0, 4: 0}
    assert nonrich["line_scope"].startswith("the p affine lines")
    assert nonrich["geometric_realization_claimed"] is False
    assert nonrich["proved"] is True

    rich = direction_refinement_row(
        31,
        3,
        A_secants=8,
        T_secants=2,
        rich_trisecants=1,
        rich_four_secants=1,
        ordinary_A_T_lines=4,
    )
    assert rich["odd_D_fibres_b"] == 4
    assert rich["direction_slack"] == 3
    assert rich["D_fibre_profile"] == {0: 15, 1: 3, 2: 11, 3: 1, 4: 1}
    assert rich["proved"] is True

    impossible = direction_refinement_row(
        31,
        3,
        A_secants=8,
        T_secants=2,
        rich_trisecants=0,
        rich_four_secants=0,
        ordinary_A_T_lines=7,
    )
    assert impossible["checks"]["line_types_nonnegative"] is False
    assert impossible["proved"] is False

    with pytest.raises(ValueError):
        direction_refinement_row(
            31,
            3,
            A_secants=-1,
            T_secants=2,
            rich_trisecants=0,
            rich_four_secants=0,
            ordinary_A_T_lines=5,
        )


def test_proposition_package_is_a_narrowing_not_an_endpoint_close():
    row = proposition_15730()
    assert row["prop"] == "15.730"
    assert row["result_status"] == "proved simultaneous necessary normal form"
    assert row["source_correction"]["proposition_15_729_core_affected"] is False
    assert row["finite_configuration_search_used"] is False
    assert row["endpoint_excluded"] is False
    assert row["p_plus_one_shell_closed"] is False
    assert row["non_walsh_residual_ii_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True
