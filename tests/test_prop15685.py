from e1_gmin_m4_prop15685 import (
    COMPLETE_17_ARC_REPRESENTATIVES,
    EXPECTED_SECANT_MULTIPLICITY_HISTOGRAMS,
    complete_17_arc_classification_certificate,
    complete_arc_certificate,
    p23_slack_twelve_exclusion,
    projective_points,
    slack_twelve_repair_certificate,
)


def test_projective_plane_has_553_points():
    points = projective_points()
    assert len(points) == 23 * 23 + 23 + 1 == 553
    assert len(set(points)) == 553


def test_all_five_representatives_are_complete_17_arcs():
    assert len(COMPLETE_17_ARC_REPRESENTATIVES) == 5
    rows = [
        complete_arc_certificate(representative)
        for representative in COMPLETE_17_ARC_REPRESENTATIVES
    ]
    assert all(row["point_count"] == 17 for row in rows)
    assert all(row["secant_line_count"] == 136 for row in rows)
    assert all(row["outside_point_count"] == 536 for row in rows)
    assert all(row["secant_outside_incidence_count"] == 2992 for row in rows)
    assert all(row["is_arc"] is True for row in rows)
    assert all(row["is_complete"] is True for row in rows)


def test_secant_multiplicity_histograms_are_exact_and_distinct():
    rows = [
        complete_arc_certificate(representative)
        for representative in COMPLETE_17_ARC_REPRESENTATIVES
    ]
    observed = tuple(row["secant_multiplicity_histogram"] for row in rows)
    assert observed == EXPECTED_SECANT_MULTIPLICITY_HISTOGRAMS
    assert len({tuple(row.items()) for row in observed}) == 5
    assert [row.get(1, 0) for row in observed] == [0, 0, 1, 0, 0]


def test_five_invariant_distinct_representatives_exhaust_classification():
    row = complete_17_arc_classification_certificate()
    assert row["classified_projective_class_count"] == 5
    assert row["verified_representative_count"] == 5
    assert row["pairwise_distinct_invariants"] is True
    assert row["therefore_exhaustive"] is True
    assert row["maximum_one_secant_point_count"] == 1
    assert row["proved"] is True


def test_slack_twelve_forces_three_one_secant_points():
    row = slack_twelve_repair_certificate()
    assert row["unique_profile"]["pair_slack"] == 12
    assert row["repair_deletion_bound"] == 3
    assert row["three_deletion_branch"]["arc_size"] == 17
    assert row["three_deletion_branch"]["therefore_arc_complete"] is True
    assert [item["line_pair_slack"] for item in row["secant_line_slack_rows"]] == [
        4,
        8,
        16,
    ]
    assert row["required_one_secant_points"] == 3
    assert row["proved"] is True


def test_prop15685_excludes_one_profile_without_overclaiming():
    row = p23_slack_twelve_exclusion()
    assert row["slack_twelve_profile_excluded"] is True
    assert row["required_one_secant_points"] == 3
    assert row["maximum_available_in_any_complete_17_arc_class"] == 1
    assert row["p23_profile_count_before"] == 203
    assert row["p23_profiles_excluded_here"] == 1
    assert row["p23_profile_count_after"] == 202
    assert 12 not in row["remaining_pair_slack_histogram"]
    assert sum(row["remaining_pair_slack_histogram"].values()) == 202
    assert row["p23_second_all_finite_endpoint_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert row["proved"] is True
