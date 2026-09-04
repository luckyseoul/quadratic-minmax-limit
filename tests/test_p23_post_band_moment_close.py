"""Fail-when-wrong tests for the exceptional p=23 post-band close."""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

from e1_gmin_m4_p23_post_band_moment_close import (
    half_mean_height_three_cube_classification,
    p23_first_post_band_ledger,
    p23_hard_moment_root_certificate,
    p23_k5_moment_sieve,
    p23_local_height_equality,
    p23_post_band_moment_close,
    p23_sharp_hard_family_catalog,
    p23_slice_half_mean_classification,
    write_evidence,
)


ROOT = Path(__file__).resolve().parents[1]


def test_all_dimension_half_mean_equality_is_exactly_f4_or_f5() -> None:
    row = half_mean_height_three_cube_classification()
    assert row["origin_facet_mean_options"] == ["1/2", "3/4"]
    assert row["induction_half_mean_facet"]["F4_extension_affine_pairs_c_u"] == [
        [-2, 1],
        [0, 0],
    ]
    assert row["induction_half_mean_facet"]["F5_extensions"] == ["dummy"]
    assert row["all_three_quarter_facets"]["possible_dimensions"] == [4]
    assert row["all_three_quarter_facets"]["forced_coefficients"] == (
        "a_i=-2 and b_ij=1"
    )
    assert row["classified_forms"]["F4"]["layer_values"] == [3, 1, 0, 0, 1]
    assert row["classified_forms"]["F5"]["layer_values"] == [
        3,
        1,
        0,
        0,
        1,
        3,
    ]
    assert row["classification_exhaustive_in_every_dimension"] is True
    assert row["proved"] is True


def test_local_mass_36_pins_height_three_and_all_half_mean_cubes() -> None:
    row = p23_local_height_equality()
    assert row["scaled_mass_4p_E_C"] == 36
    assert row["mean"] == "9/23"
    assert Fraction(
        row["height_one_exclusion"]["largest_zero_class_complement_bound"]
    ) == Fraction(19404, 2645) < 8
    assert row["height_one_exclusion"]["target_density_absent"] is True
    assert row["height_at_least_two"]["stabilizer_height_upper_bound"] == "9"
    assert row["height_at_least_two"][
        "three_quarter_cube_maximum_upper_bound"
    ] == 6
    assert row["forced_height"] == 3
    assert row["paired_cube_average_at_forced_height"] == "1/2"
    assert row["every_paired_cube_through_a_maximizer_has_mean"] == "1/2"
    assert row["proved"] is True


def test_pairing_compatibility_globalizes_every_slice_equality_form() -> None:
    row = p23_slice_half_mean_classification()
    assert row["all_additive_two_by_two_minors_vanish"] is True
    rectangles = row["additive_binary_rectangle_catalog"]
    assert len(rectangles) == 6
    assert all(
        rectangle["row_constant"] or rectangle["column_constant"]
        for rectangle in rectangles
    )
    assert row["additive_binary_rectangle_classification_proved"] is True
    assert row["binary_additive_matrix_classification"] == [
        "column-only",
        "row-only",
    ]
    assert row["active_column_counts"] == [4, 5]
    assert row["active_row_counts"] == [5]
    coverage = row["paired_cube_coverage_by_swap_count"]
    assert [entry["swap_count"] for entry in coverage] == list(range(12))
    assert all(entry["partial_matching_extends"] for entry in coverage)
    assert row["every_slice_point_lies_in_a_paired_cube_through_X"] is True
    forms = row["global_slice_forms"]
    assert [form["name"] for form in forms] == ["F4", "F5"]
    assert [form["slice_mean"] for form in forms] == ["9/23", "9/23"]
    assert [form["coefficient_offset"] for form in forms] == [1, 5]
    assert row["compatible_hard_opposite_rows"] == [
        {"hard_P": 4, "opposite_Q": 5, "form": "F5", "offset": 5}
    ]
    assert row["all_forced_local_equality_forms_enumerated"] is True
    assert row["proved"] is True


def test_sharp_hard_lift_catalog_has_four_offsets() -> None:
    row = p23_sharp_hard_family_catalog()
    assert row["sharp_scaled_mass"] == 20
    assert row["corrected_transposition_junta_bound"] == "2772/529"
    assert row["selected_table_count"] == 10
    assert row["selected_pair_table_count"] == 6
    assert row["all_equal_triple_table_count"] == 4
    assert row["selected_tables_sha256"] == (
        "a5a39bd5fd61be245c18f2e37b99ba251faaefea4477ec07ed3a96bfff5d61a1"
    )
    assert sorted(
        family["coefficient_offset"] for family in row["hard_families"]
    ) == [2, 3, 4, 5]
    assert row["proved"] is True


def test_p23_endpoint_ledger_leaves_only_p4_q5_f5_before_moments() -> None:
    row = p23_first_post_band_ledger()
    assert row["layer_index_t"] == 9
    assert row["original_k"] == 110
    assert row["H_edge_count"] == 111
    assert row["guaranteed_isolated_vertices"] == 308
    assert row["surviving_residues"] == [9, 11]
    assert all(branch["excluded"] for branch in row["old_branch_exclusions"])
    new_rows = row["new_sharp_branch_family_ledgers"]
    assert sorted(item["forced_hard_parallel_count"] for item in new_rows) == [
        2,
        3,
        4,
        5,
    ]
    assert sorted(item["forced_low_Q"] for item in new_rows) == [4, 5, 6, 7]
    assert all(item["forced_low_scaled_mean"] == 36 for item in new_rows)
    assert all(
        item["directions_at_forced_low_Q_at_least"] == 9 for item in new_rows
    )
    assert row["unique_new_survivor_before_moments"] == {
        "hard_P": 4,
        "opposite_Q": 5,
        "hard_family": "complement-literal plus all-equal-triple",
        "opposite_form": "F5",
        "opposite_Q5_directions_at_least": 9,
    }
    assert row["proved"] is True


def test_twelve_hard_roots_annihilate_both_global_forms() -> None:
    row = p23_hard_moment_root_certificate()
    assert row["full_star_power_sums_degrees_2_4_6_8"] == {
        "2": 0,
        "4": 0,
        "6": 0,
        "8": 0,
    }
    assert row["triangle_pairs_checked"] == 23**2
    assert row["form_degrees"] == {"G4": 4, "G8": 8}
    assert row["distinct_hard_projective_roots"] == 12
    assert row["both_forms_identically_zero"] is True
    assert row["opposite_F5_evaluations"] == {
        "G4": "-2*S4-S2^2",
        "G8": "-24*S8-32*S2*S6+5*S2^4",
    }
    assert row["proved"] is True


def test_exact_k5_sieve_has_disjoint_quartic_and_octic_zero_orbits() -> None:
    row = p23_k5_moment_sieve()
    assert row["five_sets_checked"] == 33649
    assert row["G4_zero_count"] == 1518
    assert row["G8_zero_count"] == 2024
    assert row["simultaneous_zero_count"] == 0
    assert row["G4_zero_sets_sha256"] == (
        "82460f67f3414a1f461b24605c108861d215f970063c0d0af82772de21240c1a"
    )
    assert row["G8_zero_sets_sha256"] == (
        "733bc62c7ad8d0d7083388480d307ad7298d56b4f9e1fcd12562848350c8d6c7"
    )
    assert row["affine_orbit_count"] == 69
    assert row["affine_orbit_size_histogram"] == {"253": 5, "506": 64}
    assert row["affine_orbit_representatives_sha256"] == (
        "34eeb59b625d24907758658f78c0f966291728a72cebb0426a3d4a883fb2022a"
    )
    assert [entry["representative"] for entry in row["G4_zero_orbits"]] == [
        [0, 1, 2, 3, 12],
        [0, 1, 2, 4, 15],
        [0, 1, 2, 7, 17],
    ]
    assert [entry["representative"] for entry in row["G8_zero_orbits"]] == [
        [0, 1, 2, 3, 10],
        [0, 1, 2, 4, 17],
        [0, 1, 2, 4, 18],
        [0, 1, 2, 7, 10],
    ]
    assert row["zero_orbit_representatives_disjoint"] is True
    assert row["proved"] is True


def test_accelerator_replay_script_is_hash_pinned() -> None:
    row = p23_k5_moment_sieve()["independent_accelerator_replay"]
    path = ROOT / row["script"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == row["script_sha256"]
    assert row["each_reported_counts"] == [33649, 1518, 2024, 0]
    assert row["authoritative_certificate"] is False


def test_package_closes_only_the_p23_endpoint_and_replays_evidence(
    tmp_path: Path,
) -> None:
    row = p23_post_band_moment_close()
    assert row["p23_k110_closed"] is True
    assert row["all_boundary_sizes_excluded"] is True
    assert row["finite_graph_or_residual_configuration_census_used"] is False
    assert row["fixed_five_set_coefficient_certificate_used"] is True
    assert row["later_layers_closed"] is False
    assert row["residual_ii_closed_globally"] is False
    assert row["E1_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["proved"] is True

    expected = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_p23_post_band_moment_close.json").read_text()
    )
    assert expected == row
    replay = tmp_path / "p23.json"
    assert write_evidence(replay) == replay
    assert json.loads(replay.read_text()) == row
    assert not list(tmp_path.glob("*.tmp"))
