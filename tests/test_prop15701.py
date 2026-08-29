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
    assert ledger["excluded_profile_count"] == 398
    observed = {
        row["pair_slack"]: (
            row["profile_count"],
            row["excluded_profile_count"],
            row["remaining_profile_count"],
        )
        for row in ledger["rows"]
    }
    assert observed == {
        4: (227, 227, 0),
        8: (195, 128, 67),
        12: (155, 43, 112),
    }


def test_p17_conic_core_forces_slack_at_least_twenty():
    lemma = p17_conic_core_slack_lemma()
    assert [
        lemma["off_conic_count_rows"][h]["pair_slack_floor"]
        for h in range(1, 4)
    ] == [20, 32, 36]
    assert lemma["positive_slack_below_twenty_impossible_after_conic_core"]


def test_prop15701_reduces_p17_remainder_to_932():
    theorem = p17_low_positive_slack_conic_reduction()
    assert theorem["profile_count_before"] == 1330
    assert theorem["profiles_excluded_here"] == 398
    assert theorem["profile_count_after"] == 932
    assert theorem["remaining_pair_slack_histogram"][0] == 2
    assert 4 not in theorem["remaining_pair_slack_histogram"]
    assert theorem["remaining_pair_slack_histogram"][8] == 67
    assert theorem["remaining_pair_slack_histogram"][12] == 112
    assert not theorem["p17_second_all_finite_endpoint_closed"]
