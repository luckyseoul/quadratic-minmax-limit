import json
from pathlib import Path

from e1_gmin_m4_prop15771 import (
    mean_46_contact_quadratures,
    mean_46_hard_family_catalog,
    mean_46_small_support_equality_catalog,
    middle_boundary_equality_exclusion,
    middle_boundary_swap_cube,
    p23_mass_24_lift_catalog,
    p23_phase_zero_mass_exclusions,
    p23_third_post_band_residue_ledger,
    p23_u11_all_one_common_row_exclusion,
    p23_u11_common_row_exclusion,
    p23_u11_zero_quotient_exclusion,
    p23_u9_two_unit_carry_exclusion,
    proposition_15771,
)


def test_third_post_band_residue_ledger_includes_the_zero_quotient_split():
    row = p23_third_post_band_residue_ledger()
    assert row["proved"]
    assert row["guaranteed_isolated_vertices"] == 300
    assert row["phase_one_quotient_sum"] == "sum k_L=23-u"
    assert row["arithmetic_surviving_residues"] == [9, 10, 11]
    by_u = {item["u"]: item for item in row["rows"]}
    assert by_u[9]["forced_quotient_one_count_if_no_quotient_zero"] == 10
    assert by_u[10]["forced_quotient_one_count_if_no_quotient_zero"] == 11
    assert by_u[11]["forced_quotient_one_count_if_no_quotient_zero"] == 12
    assert by_u[11]["quotient_zero_live_rows"] == [
        {"b": 2, "classification": "exact"},
        {"b": 22, "classification": "exact"},
    ]


def test_all_broad_floor_contact_quadratures_are_exact_and_positive():
    row = mean_46_contact_quadratures()
    assert row["proved"]
    assert row["all_weights_strictly_positive"]
    by_b = {item["b"]: item for item in row["rows"]}
    assert sorted(by_b) == list(range(4, 22, 2))
    assert by_b[4]["nodes"] == [0, 2, 4]
    assert by_b[20]["nodes"] == [10, 12]
    assert by_b[20]["weights"] == ["18/23", "5/23"]
    assert all(item["moments_0_1_2"][0] == "1" for item in row["rows"])


def test_small_support_equalities_include_the_offset_eight_endpoint():
    row = mean_46_small_support_equality_catalog()
    assert row["proved"]
    assert row["b=0"]["coefficient_offset"] == 5
    assert row["b=4"]["labeled_form_count"] == 10
    assert row["b=4"]["orbit_histogram"] == {"2200": 6, "4000": 4}
    assert row["b=4"]["coefficient_offsets"] == [5]
    assert row["b=20"]["labeled_form_count"] == 10
    assert row["b=20"]["orbit_histogram"] == {
        "000;4": 1,
        "200;2": 3,
        "220;0": 3,
        "400;0": 3,
    }
    assert row["b=20"]["coefficient_offsets"] == [4, 6, 8]
    endpoint = next(
        form for form in row["b=20"]["forms"] if form["orbit_type"] == "000;4"
    )
    assert endpoint["triple_value"] == 4
    assert endpoint["signed_linear_coefficients"] == [1, 1, 1]
    assert endpoint["coefficient_offset"] == 8


def test_middle_boundaries_die_by_full_even_half_rank():
    row = middle_boundary_equality_exclusion()
    assert row["proved"]
    assert [item["b"] for item in row["rows"]] == [6, 8, 10, 12, 14, 16, 18]
    assert [item["smaller_parity_side_size"] for item in row["rows"]] == [
        6,
        8,
        10,
        11,
        9,
        7,
        5,
    ]
    assert [item["degree_at_most_two_dimension"] for item in row["rows"]] == [
        22,
        37,
        56,
        67,
        46,
        29,
        16,
    ]
    assert all(
        item["even_half_evaluation_rank_mod_1000003"]
        == item["degree_at_most_two_dimension"]
        for item in row["rows"]
    )
    assert all(item["every_slice_point_has_cross_boundary_swap_cube"] for item in row["rows"])


def test_constructed_swap_cubes_cover_each_intersection_class():
    for boundary in range(6, 20, 2):
        boundary_set = set(range(boundary))
        small = boundary_set if boundary <= 11 else set(range(23)) - boundary_set
        outside = set(range(23)) - small
        dimension = len(small)
        for weight in range(dimension + 1):
            selected = tuple(sorted(sorted(small)[:weight] + sorted(outside)[:12 - weight]))
            cube = middle_boundary_swap_cube(boundary, selected)
            assert cube["dimension"] == dimension
            pairs = cube["pairs"]
            fixed = set(cube["fixed_selected"])
            for mask in range(1 << dimension):
                point = fixed | {pair[0] if mask & (1 << index) else pair[1]
                                 for index, pair in enumerate(pairs)}
                assert len(point) == 12
                intersection = len(point & boundary_set)
                assert intersection == (mask.bit_count() if boundary <= 11 else 12 - mask.bit_count())
                if mask == cube["base_point_mask"]:
                    assert point == set(selected)


def test_mass_24_lift_is_boolean_and_has_exact_three_family_catalog():
    row = p23_mass_24_lift_catalog()
    assert row["proved"]
    assert row["height_at_least_two"]["therefore_H"] == 6
    assert row["height_at_least_two"]["every_paired_cube_through_a_maximizer_has_mean"] == "1/2"
    assert row["height_at_least_two"]["excluded"]
    assert row["therefore_boolean"]
    assert row["corrected_johnson_junta_bound"] == "15708/2645"
    assert row["target_table_count"] == 30
    assert row["target_tables_sha256"] == (
        "b519f1118c7a375b37f974597db6d4539efaee3872b24f03080d2dbd26b60a51"
    )
    assert [
        (
            item["family"],
            item["table_count"],
            item["four_L_offset_increment_after_slice_complement"],
        )
        for item in row["families"]
    ] == [
        ("selected_pair_on_original_slice", 6, 3),
        ("oriented_pair_on_original_slice", 12, 1),
        ("compact_triangle", 12, 1),
    ]


def test_mean_46_catalog_accounts_for_every_even_boundary_and_offset():
    row = mean_46_hard_family_catalog()
    assert row["proved"]
    assert row["all_even_boundary_values_accounted_for"] == list(range(0, 23, 2))
    assert row["excluded_boundary_values"] == [6, 8, 10, 12, 14, 16, 18]
    assert row["possible_coefficient_offsets"] == [4, 5, 6, 7, 8]


def test_u9_two_unit_carry_still_has_enough_roots():
    row = p23_u9_two_unit_carry_exclusion()
    assert row["proved"]
    assert row["low_triangle_minus_star_projective_roots_at_least"] == 10
    assert row["maximum_common_form_degree"] == 8
    assert row["roots_force_G4_and_G8_identically_zero"]
    assert row["unique_survivor_before_moments"] == {
        "hard_P": 4,
        "opposite_Q": 5,
        "opposite_form": "F5",
    }
    assert row["opposite_K5_simultaneous_zero_count"] == 0
    assert all(
        item["directions_at_forced_low_Q_at_least"] == 7
        for item in row["family_ledgers"]
    )


def test_both_u11_quotient_profile_branches_force_five_mass_32_rows():
    zero = p23_u11_zero_quotient_exclusion()
    all_one = p23_u11_all_one_common_row_exclusion()
    package = p23_u11_common_row_exclusion()
    assert zero["proved"] and all_one["proved"] and package["proved"]
    assert [
        (item["quotient_zero_boundary_b"], item["common_base_c"])
        for item in zero["baseline_ledgers"]
    ] == [(2, 4), (22, 3)]
    for item in zero["baseline_ledgers"] + all_one["offset_ledgers"]:
        assert item["forbidden_scaled_mass"] == 8
        assert item["forced_scaled_mass"] == 32
        assert item["opposite_surplus"] == 7
        assert item["directions_at_forced_Q_at_least"] == 5


def test_phase_zero_bridge_checks_all_boundaries_and_nonnegative_lifts():
    row = p23_phase_zero_mass_exclusions()
    assert row["proved"] and row["all_even_boundaries_excluded"]
    assert row["baseline_subtraction_is_nonnegative_by_parity"]
    assert row["pointwise_baselines"] == {
        "b=0": "0", "b=2": "(x_i-x_j)^2", "b=22": "x_j",
    }
    assert row["sharp_positive_lift_mass_floor"] == 20
    masses = {item["scaled_mean"]: item for item in row["rows"]}
    assert masses[8]["floor_surviving_boundaries"] == [0]
    assert masses[32]["floor_surviving_boundaries"] == [0, 2, 22]
    for mass in (8, 32):
        assert [item["b"] for item in masses[mass]["boundary_rows"]] == list(range(0, 23, 2))
        assert all(item["excluded"] for item in masses[mass]["boundary_rows"])
    assert [(item["b"], item["lift_mass"]) for item in masses[32]["boundary_rows"]
            if item["lift_mass"] is not None] == [(0, 32), (2, 8), (22, 8)]
    for package in (p23_u11_zero_quotient_exclusion(), p23_u11_all_one_common_row_exclusion()):
        assert package["all_boundary_opposite_mass_bridge"] == row


def test_proposition_15771_closes_only_the_reviewed_endpoint():
    row = proposition_15771()
    assert row["certificate_checks_passed"]
    assert row["status"].startswith("PROVED")
    assert row["proof_review_complete"]
    assert not row["pending_proof_bridges"]
    assert len(row["completed_proof_bridges"]) == 3
    assert row["proved"]
    assert row["p23_k114_closed"]
    assert row["all_boundary_sizes_excluded"]
    assert not row["new_graph_or_residual_configuration_census_used"]
    assert not row["later_p23_layers_closed"]
    assert not row["residual_ii_closed_globally"]
    assert not row["E1_closed"]
    assert not row["quadratic_minmax_limit_closed"]


def test_reviewed_endpoint_evidence_and_proof_package_are_synchronized():
    root = Path(__file__).resolve().parents[1]
    row = proposition_15771()
    saved = json.loads((root / "evidence/e1_gmin_m4_prop15771.json").read_text())
    assert saved == json.loads(json.dumps(row))
    for relative in row["proof_notes"]:
        assert (root / relative).is_file()
    for name in ("STATUS.md", "README.md", "LONG_HORIZON_GOAL.md", "AGENTS.md"):
        text = (root / name).read_text()
        assert "15.771" in text
        assert "p=23,t>=12" in text
        assert "p=23,t>=11" not in text
