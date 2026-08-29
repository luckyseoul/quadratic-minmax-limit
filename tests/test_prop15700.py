from e1_gmin_m4_prop15700 import (
    conic_minus_two_affine_profile_census,
    p17_second_boundary_profile_census,
    p17_slack_zero_conic_reduction,
    residue_and_lift_ledger,
)


def test_p17_second_boundary_residue_and_profile_ledger():
    ledger = residue_and_lift_ledger()
    assert ledger["remaining_phase_zero_residues"] == [0, 7, 8]
    assert ledger["remaining_phase_one_residues"] == [0, 8]
    assert all(
        row["excluded"]
        for row in ledger["phase_zero_positive_residues_excluded"]
    )

    census = p17_second_boundary_profile_census()
    assert census["phase_labelled_profile_count"] == 1575
    assert census["slack_zero_profile_count"] == 247
    assert census["slack_zero_residue_pair_histogram"] == {
        (0, 8): 234,
        (7, 8): 4,
        (8, 8): 9,
    }


def test_p17_unique_16_arc_affine_profile_census():
    census = conic_minus_two_affine_profile_census()
    assert census["projective_line_count"] == 307
    assert census["raw_affine_case_count"] == 21267
    assert census["raw_case_count_by_line_intersection"] == {
        0: 20808,
        1: 306,
        2: 153,
    }
    assert census["phase_labelled_profile_count_including_swap"] == 53


def test_prop15700_reduces_slack_zero_to_two_tangent_profiles():
    theorem = p17_slack_zero_conic_reduction()
    assert theorem["profile_count_before"] == 1575
    assert theorem["profiles_excluded_here"] == 245
    assert theorem["profile_count_after"] == 1330
    assert theorem["slack_zero_profile_count_after"] == 2
    assert all(
        row["conic_example"]["line_conic_intersection_size"] == 1
        for row in theorem["surviving_slack_zero_profiles"]
    )
    assert not theorem["p17_second_all_finite_endpoint_closed"]
