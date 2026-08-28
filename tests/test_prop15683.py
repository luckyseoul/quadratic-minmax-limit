from e1_gmin_m4_prop15683 import (
    p41_arc_envelope_exclusion,
    p41_endpoint_theorem,
    p41_near_arc_envelope_exclusion,
    p41_residue_zero_profiles,
    tangent_envelope_input,
)


def test_exact_p41_residue_zero_profiles_are_seven_arcs_and_two_near_arcs():
    row = p41_residue_zero_profiles()
    assert row["profile_count"] == 9
    assert row["distinct_global_shape_count"] == 4
    assert row["arc_profile_count"] == 7
    assert row["near_arc_profile_count"] == 2
    assert row["minimum_undetermined_directions"] == 5
    assert {item["pair_slack"] for item in row["profiles"]} == {0, 4}
    for item in row["profiles"]:
        floor_pairs = sum(
            int(secants) * count
            for secants, count in item["global_secant_distribution"].items()
        )
        assert 2 * floor_pairs + item["pair_slack"] == 34 * 33
        assert item["global_secant_distribution"]["17"] == 14
        assert item["global_secant_distribution"]["16"] == 20


def test_ball_lavrauw_tangent_envelope_input_has_the_needed_square_restriction():
    row = tangent_envelope_input()
    assert row["doi"] == "10.1016/j.jcta.2018.06.015"
    assert row["arxiv"] == "1705.10940v4"
    assert row["theorem"] == 11
    assert "degree-2t" in row["odd_order_statement"]
    assert "square" in row["odd_order_statement"]
    assert row["proved"] is True


def test_degree_eighteen_envelope_excludes_all_seven_arc_profiles():
    row = p41_arc_envelope_exclusion()
    assert row["profile_count"] == 7
    assert row["arc_size"] == 34
    assert row["tangents_per_point"] == 9
    assert row["envelope_degree"] == 18
    assert row["minimum_tangents_per_exceptional_direction"] == 28
    assert row["forced_double_direction_component_degree"] == 16
    assert row["residual_curve_degree"] == 2
    assert row["exceptional_secant_edges"] == 3
    assert row["minimum_incident_arc_points"] == 3
    assert row["double_zero_multiplicity_on_point_pencil"] == 4
    assert row["proved"] is True


def test_deleted_triple_degree_twenty_envelope_excludes_both_near_arcs():
    row = p41_near_arc_envelope_exclusion()
    assert row["profile_count"] == 2
    deletion = row["unique_triple_deletion"]
    assert deletion["resulting_arc_size"] == 33
    assert deletion["seven_no_secant_direction_tangents"] == 33
    assert deletion["surviving_exceptional_pair_direction_tangents"] == 31
    assert row["tangents_per_point"] == 10
    assert row["envelope_degree"] == 20
    assert row["forced_double_direction_component_degree"] == 16
    assert row["residual_curve_degree"] == 4
    assert row["surviving_high_pair_endpoints"] == 2
    assert row["low_tangents_at_each_high_pair_endpoint"] == 3
    assert row["residual_after_two_point_pencils"] == 2
    assert row["other_arc_points"] == 31
    assert row["proved"] is True


def test_theorem_closes_only_the_p41_second_endpoint():
    row = p41_endpoint_theorem()
    assert row["positive_residues_excluded_by_prop15681"] is True
    assert row["residue_zero_profile_count"] == 9
    assert row["arc_profiles_excluded"] == 7
    assert row["near_arc_profiles_excluded"] == 2
    assert row["second_all_finite_endpoint_closed"] is True
    assert row["remaining_same_boundary_primes"] == [17, 19, 23]
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True
