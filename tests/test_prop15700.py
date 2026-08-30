from e1_gmin_m4_prop15700 import (
    conic_minus_two_affine_profile_census,
    p17_second_boundary_profile_census,
    p17_slack_zero_conic_reduction,
    residue_and_lift_ledger,
)
from e1_gmin_m4_prop15723 import floor_excess_admissible


def test_p17_second_boundary_residue_and_profile_ledger():
    ledger = residue_and_lift_ledger()
    assert ledger["remaining_phase_zero_residues"] == [0, 7, 8]
    assert ledger["remaining_phase_one_residues"] == [0, 8]
    assert all(
        row["excluded"]
        for row in ledger["phase_zero_positive_residues_excluded"]
    )
    assert ledger["phase_profile_row_counts_at_full_cap"]["0"][0] == 565
    assert ledger["phase_profile_row_counts_at_full_cap"]["0"][6] == 40

    # Complemented reduced-size 3/4 cells and both genuine p=17 equality
    # cells remain admissible; endpoint and generic middle obstructions do not.
    assert floor_excess_admissible(17, 4, 0, 2) is True
    assert floor_excess_admissible(17, 14, 1, 2) is True
    assert floor_excess_admissible(17, 6, 1, 2) is True
    assert floor_excess_admissible(17, 12, 0, 2) is True
    assert floor_excess_admissible(17, 2, 0, 2) is False
    assert floor_excess_admissible(17, 8, 0, 2) is False

    census = p17_second_boundary_profile_census()
    assert census["phase_labelled_profile_count"] == 2503
    assert census["slack_zero_profile_count"] == 286
    assert {
        slack: census["pair_slack_histogram"][slack]
        for slack in (0, 4, 8, 12, 16, 20, 24)
    } == {0: 286, 4: 292, 8: 292, 12: 267, 16: 227, 20: 193, 24: 151}
    assert census["residue_pair_histogram"] == {
        (0, 0): 275,
        (0, 8): 1896,
        (7, 0): 9,
        (7, 8): 9,
        (8, 0): 37,
        (8, 8): 277,
    }
    assert census["slack_zero_residue_pair_histogram"] == {
        (0, 8): 273,
        (7, 8): 4,
        (8, 8): 9,
    }
    assert census["canonical_profile_sha256"] == (
        "48632c09fdf9ed38d4f8608aeb0251bd29af2ac7b5fb81d090657a8ed20793b9"
    )


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
    assert theorem["profile_count_before"] == 2503
    assert theorem["profiles_excluded_here"] == 284
    assert theorem["profile_count_after"] == 2219
    assert theorem["slack_zero_profile_count_after"] == 2
    assert all(
        row["conic_example"]["line_conic_intersection_size"] == 1
        for row in theorem["surviving_slack_zero_profiles"]
    )
    assert len(theorem["remaining_profile_indices"]) == 2219
    assert not theorem["p17_second_all_finite_endpoint_closed"]
