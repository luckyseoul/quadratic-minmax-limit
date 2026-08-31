import pytest

from e1_gmin_m4_prop15729 import (
    DISTINGUISHED_BLOCK_TYPES,
    cotangent_extension_row,
    endpoint_residue_size_row,
    endpoint_unique_trisecant_construction,
    endpoint_unique_trisecant_rows,
    proposition_15729,
)


def test_both_residue_classes_have_the_claimed_exact_sizes():
    expected = {
        19: (6, 1, 15, 13, "2R+3", "2R+1"),
        31: (10, 1, 23, 21, "2R+3", "2R+1"),
        37: (12, 1, 27, 25, "2R+3", "2R+1"),
        17: (5, 2, 14, 12, "2R+4", "2R+2"),
        23: (7, 2, 18, 16, "2R+4", "2R+2"),
        41: (13, 2, 30, 28, "2R+4", "2R+2"),
    }
    for p, (R, c, unique_size, arc_size, unique_formula, arc_formula) in expected.items():
        row = endpoint_residue_size_row(p)
        assert row["R"] == R
        assert row["c"] == c
        assert row["unique_trisecant_size"] == unique_size == p + 2 - R
        assert row["cotangent_arc_size"] == arc_size == p - R
        assert row["unique_trisecant_size_formula"] == unique_formula
        assert row["cotangent_arc_size_formula"] == arc_formula
        assert row["size_drop_from_D_to_unique_trisecant_set"] == R - 1
        assert row["size_drop_from_unique_trisecant_set_to_arc"] == 2
        assert row["proved"] is True


def test_distinguished_trisecant_deletes_exactly_R_minus_one_points():
    row = endpoint_unique_trisecant_construction(31, 3, "trisecant")
    assert row["trisecants_x"] == 4
    assert row["four_secants_y"] == 3
    assert row["singleton_points"] == 8
    assert row["deleted_on_distinguished_block"] == 0
    assert row["deleted_on_other_trisecants"] == 3
    assert row["deleted_on_other_four_secants"] == 6
    assert row["total_deletions"] == row["target_deletions"] == 9
    assert row["retained_singleton_points"] == 8
    assert row["retained_on_distinguished_block"] == 3
    assert row["other_rich_block_count"] == 6
    assert row["retained_on_every_other_rich_block"] == 2
    assert row["retained_point_count"] == 23
    assert row["all_points_affine"] is True
    assert row["maximum_line_occupancy"] == 3
    assert row["trisecant_count"] == 1
    assert row["finite_configuration_search_used"] is False
    assert row["proved"] is True


def test_distinguished_four_secant_deletes_exactly_R_minus_one_points():
    row = endpoint_unique_trisecant_construction(31, 3, "4-secant")
    assert row["trisecants_x"] == 4
    assert row["four_secants_y"] == 3
    assert row["deleted_on_distinguished_block"] == 1
    assert row["deleted_on_other_trisecants"] == 4
    assert row["deleted_on_other_four_secants"] == 4
    assert row["total_deletions"] == row["target_deletions"] == 9
    assert row["retained_point_count"] == 23
    assert row["maximum_line_occupancy"] == 3
    assert row["trisecant_count"] == 1
    assert row["proved"] is True


def test_every_endpoint_block_row_and_available_distinguished_type_is_covered():
    for p in (17, 19, 23, 29, 31, 37, 41, 43, 47, 53):
        R = (p - 1) // 3
        rows = endpoint_unique_trisecant_rows(p)
        expected_types = []
        for y in range(R // 2 + 1):
            x = R - 2 * y
            if x:
                expected_types.append((y, "trisecant"))
            if y:
                expected_types.append((y, "4-secant"))
        assert [
            (row["four_secants_y"], row["distinguished_block_type"])
            for row in rows
        ] == expected_types
        assert all(row["total_deletions"] == R - 1 for row in rows)
        assert all(row["retained_point_count"] == p + 2 - R for row in rows)
        assert all(row["trisecant_count"] == 1 for row in rows)
        assert all(row["proved"] is True for row in rows)


def test_removing_any_two_triple_points_gives_two_cotangent_extensions():
    for block_type in DISTINGUISHED_BLOCK_TYPES:
        row = cotangent_extension_row(31, 3, block_type)
        assert row["unique_trisecant_set_size"] == 23
        assert row["delete_any_two_points_of_unique_trisecant"] is True
        assert row["remaining_set_is_an_arc"] is True
        assert row["arc_size"] == 21
        assert row["arc_size_formula"] == "p-R"
        assert row["number_of_choices_for_remaining_trisecant_point"] == 3
        assert row["deleted_points_are_distinct_affine_extension_points"] is True
        assert row["adding_either_deleted_point_preserves_the_arc"] is True
        assert row["common_line_meets_arc_only_in_remaining_trisecant_point"] is True
        assert row["extensions_lie_on_one_arc_tangent"] is True
        assert row["proved"] is True


def test_proposition_package_is_a_reduction_not_an_endpoint_close():
    row = proposition_15729()
    assert row["prop"] == "15.729"
    assert row["result_status"] == "proved structural reduction"
    assert row["conclusion"] == {
        "affine_unique_trisecant_set": (
            "a (p+2-R,3)-arc with exactly one trisecant"
        ),
        "affine_cotangent_arc": (
            "a (p-R)-arc with two distinct extension points on one tangent"
        ),
        "residue_c_1_sizes": {"unique_trisecant": "2R+3", "arc": "2R+1"},
        "residue_c_2_sizes": {"unique_trisecant": "2R+4", "arc": "2R+2"},
    }
    assert row["finite_configuration_search_used"] is False
    assert row["new_classification_used"] is False
    assert row["endpoint_excluded"] is False
    assert row["p_plus_one_shell_closed"] is False
    assert row["non_walsh_residual_ii_closed"] is False
    assert row["multi_level_type_I_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True


def test_parameter_validation_rejects_missing_or_unknown_distinguished_blocks():
    with pytest.raises(ValueError):
        endpoint_unique_trisecant_construction(31, 0, "4-secant")
    with pytest.raises(ValueError):
        endpoint_unique_trisecant_construction(37, 6, "trisecant")
    with pytest.raises(ValueError):
        endpoint_unique_trisecant_construction(31, 3, "five-secant")
    with pytest.raises(ValueError):
        endpoint_residue_size_row(25)
