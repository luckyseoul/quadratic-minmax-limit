from e1_gmin_m4_prop15703 import (
    complete_fourteen_minus_one_certificate,
    complete_thirteen_arc_certificate,
    p17_complete_thirteen_arc_classification,
    p17_slack_twelve_exclusion,
)


def test_p17_complete_thirteen_arc_classification_fingerprint():
    row = p17_complete_thirteen_arc_classification()
    assert row["complete_13_arc_class_count"] == 8
    assert row["stabilizer_group_orders"] == [1, 2, 2, 2, 2, 3, 4, 6]


def test_all_complete_thirteen_arc_index_one_triples_fail():
    certificate = complete_thirteen_arc_certificate()
    assert certificate["class_count"] == 8
    assert certificate["pairwise_invariant_distinct_secant_index_histograms"]
    assert certificate["index_one_point_counts"] == [0, 0, 0, 0, 0, 0, 2, 3]
    assert certificate["all_candidate_index_one_triple_slacks"] == [16]


def test_complete_fourteen_minus_one_index_one_triples_fail():
    certificate = complete_fourteen_minus_one_certificate()
    assert certificate["index_one_count_histogram"] == {0: 4, 1: 8, 4: 2}
    assert certificate["candidate_index_one_triple_count"] == 8
    assert certificate["candidate_index_one_triple_slack_histogram"] == {20: 8}


def test_prop15703_closes_slack_twelve_and_leaves_753():
    theorem = p17_slack_twelve_exclusion()
    assert theorem["profile_count_before"] == 786
    assert theorem["profiles_excluded_here"] == 33
    assert theorem["profile_count_after"] == 753
    assert 12 not in theorem["remaining_pair_slack_histogram"]
    assert theorem["remaining_pair_slack_histogram"][0] == 2
    assert theorem["remaining_profiles_of_slack_at_least_sixteen"] == 751
    assert not theorem["p17_second_all_finite_endpoint_closed"]
