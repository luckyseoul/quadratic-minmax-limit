from e1_gmin_m4_prop15701 import (
    p17_conic_core_slack_lemma,
    p17_fifteen_arc_classification,
    p17_low_positive_slack_conic_reduction,
    p17_low_positive_slack_profile_ledger,
)


def test_unique_p17_fifteen_arc_is_conic_contained():
    classification = p17_fifteen_arc_classification()
    assert classification["pgl_class_count_of_15_arcs"] == 1
    assert classification["consequence"] == (
        "every arc of size at least 15 in PG(2,17) is conic-contained"
    )


def test_p17_low_positive_slack_profile_ledger():
    ledger = p17_low_positive_slack_profile_ledger()
    assert ledger["excluded_profile_count"] == 475
    observed = {
        row["pair_slack"]: (
            row["profile_count"],
            row["excluded_profile_count"],
            row["remaining_profile_count"],
        )
        for row in ledger["rows"]
    }
    assert observed == {
        4: (292, 292, 0),
        8: (292, 140, 152),
        12: (267, 43, 224),
    }
    assert {
        row["pair_slack"]: row["undetermined_direction_histogram"]
        for row in ledger["rows"]
    } == {
        4: {0: 178, 1: 102, 2: 12},
        8: {0: 152, 1: 116, 2: 24},
        12: {0: 113, 1: 111, 2: 43},
    }


def test_p17_conic_core_forces_slack_at_least_twenty():
    lemma = p17_conic_core_slack_lemma()
    assert [
        lemma["off_conic_count_rows"][h]["pair_slack_floor"]
        for h in range(1, 4)
    ] == [20, 32, 36]
    assert lemma["positive_slack_below_twenty_impossible_after_conic_core"]


def test_prop15701_reduces_p17_remainder_to_1744():
    theorem = p17_low_positive_slack_conic_reduction()
    assert theorem["profile_count_before"] == 2219
    assert theorem["profiles_excluded_here"] == 475
    assert theorem["profile_count_after"] == 1744
    assert theorem["remaining_pair_slack_histogram"][0] == 2
    assert 4 not in theorem["remaining_pair_slack_histogram"]
    assert theorem["remaining_pair_slack_histogram"][8] == 152
    assert theorem["remaining_pair_slack_histogram"][12] == 224
    assert len(theorem["remaining_profile_indices"]) == 1744
    assert not theorem["p17_second_all_finite_endpoint_closed"]
