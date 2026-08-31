import pytest

from e1_gmin_m4_prop15727 import (
    CLASSIFIED_ENDPOINT_PRIMES,
    active_first_possible_slack,
    classified_endpoint_exclusions,
    endpoint_block_row,
    endpoint_concavity_ledger,
    endpoint_residue_data,
    first_unexcluded_endpoint,
    linewise_endpoint_equality,
    p17_endpoint_exclusion,
    p19_endpoint_exclusion,
    p23_endpoint_exclusion,
    p29_complete_arc_spectrum,
    p29_complete_twenty_one_arc_certificate,
    p29_endpoint_exclusion,
    p29_klein_twenty_four_arc_certificate,
    proper_deletion_size_contradiction,
    proposition_15727,
    universal_endpoint_rigidity,
)


def test_endpoint_residues_and_current_first_possible_slacks():
    expected = {
        17: (5, 2, 6),
        19: (6, 1, 7),
        23: (7, 2, 8),
        29: (9, 2, 10),
        31: (10, 1, 10),
        37: (12, 1, 12),
        41: (13, 2, 13),
    }
    assert CLASSIFIED_ENDPOINT_PRIMES == (17, 19, 23, 29)
    for p, (R, c, current) in expected.items():
        assert first_unexcluded_endpoint(p) == R
        assert endpoint_residue_data(p) == {"p": p, "R": R, "c": c}
        assert active_first_possible_slack(p) == current


def test_every_proper_deletion_size_is_still_excluded_at_the_endpoint():
    for p in (17, 19, 23, 29, 31, 37, 41, 43, 47, 53):
        R = first_unexcluded_endpoint(p)
        ledger = endpoint_concavity_ledger(p)
        assert ledger["proper_deletion_interval"] == [1, R - 1]
        assert ledger["envelope_size_margin_at_t_R_minus_one"] in (1, 2)
        assert ledger["twice_F_1_minus_R"] > 0
        assert ledger["twice_F_R_minus_one_minus_R"] > 0
        assert ledger["all_proper_deletion_sizes_excluded"] is True
        assert ledger["proved"] is True
        for t in range(1, R):
            row = proper_deletion_size_contradiction(p, t)
            assert row["size_hypothesis_met"] is True
            assert row["twice_incidence_lower_bound"] > 2 * R
            assert row["contradiction"] is True


def test_linewise_equality_has_only_trisecant_and_four_secant_rich_cases():
    rich_cases = []
    for a in range(3):
        for u in range(40):
            row = linewise_endpoint_equality(a, u)
            assert row["classification_correct"] is True
            if row["rich_equality_case"]:
                rich_cases.append((a, u, row["line_occupancy"]))
    assert rich_cases == [(2, 1, 3), (2, 2, 4)]


def test_disjoint_block_equation_and_singleton_counts():
    for p in (17, 19, 23, 29, 31, 37, 41):
        R = first_unexcluded_endpoint(p)
        c = p - 3 * R
        for y in range(R // 2 + 1):
            row = endpoint_block_row(p, y)
            assert row["trisecants_x"] + 2 * row["four_secants_y"] == R
            assert row["singleton_points"] == c + 1 + 2 * y
            assert row["minimum_deletions"] == R
            assert row["maximum_arc_size"] == p + 1 - R
            assert row["maximum_arc_choice_count"] > 0
            assert row["rich_lines_pairwise_D_disjoint"] is True
            assert row["line_moments_match"] is True
            line_counts = row["projective_line_occupancy_counts"]
            assert line_counts[3] == row["trisecants_x"]
            assert line_counts[4] == row["four_secants_y"]
            core = row["regular_trisecant_core"]
            assert core["point_count"] == 3 * (R - y)
            assert core["trisecant_count"] == R - y
            assert core["every_point_on_exactly_one_trisecant"] is True
            assert core["point_count"] + core["tangents_per_point"] == p + 3
            if c == 2:
                assert row["trisecants_x"] >= 1
            assert row["proved"] is True


def test_universal_endpoint_reduction_forces_R_index_one_points_but_not_close():
    for p in (17, 19, 23, 29, 31, 37, 41):
        row = universal_endpoint_rigidity(p)
        R = (p - 1) // 3
        assert row["minimum_cardinality_arc_deletion_size"] == R
        assert row["incidence_sum_I"] == R
        assert row["every_deleted_point_secant_index"] == 1
        assert row["required_c1_of_arc_at_least"] == R
        assert row["block_equation"] == "x+2y=R"
        assert row["rich_lines_pairwise_D_disjoint"] is True
        assert row["finite_search_used_for_universal_reduction"] is False
        assert row["endpoint_excluded_by_universal_reduction_alone"] is False
        assert row["result_status"] == "proved structural reduction"
        assert row["proved"] is True


def test_p17_endpoint_fails_all_complete_and_incomplete_thirteen_arc_branches():
    row = p17_endpoint_exclusion()
    assert row["required_c1"] == 5
    assert row["complete_13_arc_c1_maximum"] == 3
    assert row["complete_14_minus_one_c1_maximum"] == 4
    assert row["conic_contained_c1"] == 0
    assert row["maximum_available_c1"] == 4
    assert row["complete_14_classification"]["complete_14_arc_class_count"] == 1
    assert row["fifteen_arc_classification"]["pgl_class_count_of_15_arcs"] == 1
    assert row["conic_branch"]["off_conic_retained_secants_at_least"] == 3
    assert row["excluded"] is True
    assert row["proved"] is True


def test_p19_endpoint_fails_exhaustive_fourteen_arc_c1_bound():
    row = p19_endpoint_exclusion()
    assert row["required_c1"] == 6
    assert row["maximum_available_c1"] == 4
    assert row["classification"]["projective_fourteen_arc_classes"] == 83
    assert row["classification"]["maximum_c1_over_all_fourteen_arcs"] == 4
    assert row["excluded"] is True
    assert row["proved"] is True


def test_p23_endpoint_fails_complete_seventeen_arc_and_conic_branches():
    row = p23_endpoint_exclusion()
    assert row["required_c1"] == 7
    assert row["complete_17_arc_c1_maximum"] == 1
    assert row["conic_contained_c1"] == 0
    assert row["complete_arc_spectrum"]["no_complete_arc_sizes"] == [
        18,
        19,
        20,
        21,
        22,
        23,
    ]
    assert row["conic_branch"]["off_conic_retained_secants_at_least"] == 4
    assert row["excluded"] is True
    assert row["proved"] is True


def test_p29_complete_twenty_one_representatives_have_no_index_one_points():
    spectrum = p29_complete_arc_spectrum()
    row = p29_complete_twenty_one_arc_certificate()
    assert spectrum["complete_arc_counts_used"] == {21: 2, 24: 1, 30: 1}
    assert spectrum["no_complete_arc_sizes"] == [22, 23, 25, 26, 27, 28, 29]
    assert spectrum["maximum_arc_size"] == 30
    assert row["classified_projective_class_count"] == 2
    assert row["verified_representative_count"] == 2
    assert row["pairwise_distinct_secant_index_histograms"] is True
    assert row["therefore_exhaustive"] is True
    assert row["index_one_point_counts_by_class"] == [0, 0]
    assert [
        representative["outside_secant_index_histogram"]
        for representative in row["representatives"]
    ] == [
        {4: 18, 5: 75, 6: 190, 7: 312, 8: 189, 9: 63, 10: 3},
        {3: 3, 4: 21, 5: 66, 6: 187, 7: 294, 8: 243, 9: 27, 10: 9},
    ]
    assert row["proved"] is True


def test_p29_incomplete_twenty_one_arcs_extend_to_klein_or_conic():
    klein = p29_klein_twenty_four_arc_certificate()
    assert klein["point_count"] == 24
    assert klein["outside_secant_index_histogram"] == {
        6: 28,
        8: 126,
        9: 504,
        10: 84,
        11: 84,
        12: 21,
    }
    assert klein["minimum_outside_secant_index"] == 6
    assert klein["is_arc"] is True
    assert klein["is_complete"] is True
    assert klein["proved"] is True

    row = p29_endpoint_exclusion()
    assert row["required_c1"] == 9
    assert row["complete_21_arc_c1_maximum"] == 0
    assert row["klein_24_minus_three_outside_secant_index_minimum"] == 3
    assert row["conic_branch"]["off_conic_retained_secants_at_least"] == 5
    assert row["excluded"] is True
    assert row["proved"] is True


def test_classified_endpoint_package_and_scope_are_honest():
    rows = classified_endpoint_exclusions()
    assert set(rows) == {17, 19, 23, 29}
    assert all(row["excluded"] and row["proved"] for row in rows.values())

    result = proposition_15727()
    assert result["prop"] == "15.727"
    assert result["first_possible_positive_slack_after"] == {
        "17": 6,
        "19": 7,
        "23": 8,
        "29": 10,
        "31": 10,
        "37": 12,
        "41": 13,
    }
    assert result["first_prime_not_endpoint_excluded_here"] == 31
    assert result["new_long_solver_run_used"] is False
    assert result["external_classification_assisted"] is True
    assert result["p_plus_one_shell_closed"] is False
    assert result["non_walsh_residual_ii_closed"] is False
    assert result["multi_level_type_I_closed"] is False
    assert result["quadratic_minmax_limit_closed"] is False
    assert result["top_level_gates_changed"] is False
    assert "from p=31 onward" in result["remaining_scope"]
    assert result["proved"] is True


def test_parameter_validation_rejects_wrong_ranges():
    for invalid in (2, 3, 15, 16, 21, 25, 27, 35, 49):
        with pytest.raises(ValueError):
            first_unexcluded_endpoint(invalid)
    with pytest.raises(ValueError):
        proper_deletion_size_contradiction(17, 0)
    with pytest.raises(ValueError):
        proper_deletion_size_contradiction(17, 5)
    with pytest.raises(ValueError):
        linewise_endpoint_equality(3, 1)
    with pytest.raises(ValueError):
        linewise_endpoint_equality(2, -1)
    with pytest.raises(ValueError):
        endpoint_block_row(17, 3)
