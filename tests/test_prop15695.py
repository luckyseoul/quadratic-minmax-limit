from e1_gmin_m4_prop15695 import (
    p19_b14_floor_equality_ledger,
    p19_b14_layer_rank_certificate,
    p19_slack_twenty_b14_exclusion,
)


def test_b14_phase_one_floors_saturate_the_type_budget():
    row = p19_b14_floor_equality_ledger()
    assert row["phase_profile"] == {2: 9, 14: 1}
    assert row["type_floor_sum"] == row["type_budget"] == 200
    assert row["symmetrized_minimizer"] == "q(t)=1"
    assert row["forced_pointwise_one_layers"] == [6, 8, 10]
    assert row["contradicting_even_parity_layer"] == 5


def test_three_forced_layers_determine_every_quadratic():
    row = p19_b14_layer_rank_certificate()
    assert row["degree_at_most_two_dimension"] == 171
    assert row["witness_layer_histogram"] == {6: 5, 8: 75, 10: 91}
    assert row["finite_field_modulus"] == 101
    assert row["finite_field_rank"] == 171
    assert row["therefore_rational_rank"] == 171


def test_b14_profiles_are_excluded_but_endpoint_stays_open():
    row = p19_slack_twenty_b14_exclusion()
    assert row["slack_twenty_profiles_before"] == 4
    assert row["excluded_b14_profiles"] == 2
    assert row["slack_twenty_profiles_after"] == 2
    assert row["p19_profiles_before"] == 7
    assert row["p19_profiles_after"] == 5
    assert row["remaining_slack_histogram"] == {20: 2, 24: 1, 28: 1, 32: 1}
    assert row["p19_second_all_finite_endpoint_closed"] is False
    assert row["L_status"] == "OPEN"
