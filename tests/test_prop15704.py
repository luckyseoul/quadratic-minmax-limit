from e1_gmin_m4_prop15704 import (
    one_direction_complete_arc_certificate,
    p17_slack_sixteen_free_direction_exclusion,
    slack_sixteen_profile_ledger,
    slack_sixteen_repair_lemma,
)


def test_slack_sixteen_profile_ledger():
    ledger = slack_sixteen_profile_ledger()
    assert ledger["undetermined_direction_histogram"] == {0: 87, 1: 88, 2: 47, 3: 5}
    assert ledger["profiles_with_at_least_one_undetermined_direction"] == 140
    assert ledger["remaining_zero_direction_profiles"] == 87


def test_slack_sixteen_repair_lemma():
    lemma = slack_sixteen_repair_lemma()
    assert lemma["repair_deletion_bound"] == 4
    assert lemma["minimum_positive_conic_core_slack"] == 20


def test_one_direction_complete_arc_certificate():
    certificate = one_direction_complete_arc_certificate()
    assert certificate["maximum_complete_13_index_one_point_count"] == 3
    assert certificate["raw_candidate_infinity_placement_slack_histogram"] == {
        16: 2,
        28: 16,
        32: 8,
    }
    assert certificate["genuinely_undetermined_infinity_placement_slack_histogram"] == {32: 8}


def test_prop15704_accounting():
    theorem = p17_slack_sixteen_free_direction_exclusion()
    assert theorem["profile_count_before"] == 1368
    assert theorem["profiles_excluded_here"] == 140
    assert theorem["profile_count_after"] == 1228
    assert theorem["remaining_slack_sixteen_profiles"] == 87
    assert len(theorem["remaining_profile_indices"]) == 1228
