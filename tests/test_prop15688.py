from fractions import Fraction

from e1_gmin_m4_prop15688 import (
    p19_second_boundary_reduction,
    p19_residue_zero_profiles,
    sharp_integral_quadratic_lift_floor,
    theorem_sharp_lift_and_p19,
)


def test_sharp_integral_quadratic_lift_floor_and_equality_example():
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 101):
        row = sharp_integral_quadratic_lift_floor(p)
        assert row["sharp_scaled_floor"] == p - 3
        assert row["sharp_mass_floor"] == Fraction(p - 3, 4 * p)
        assert row["equality_example_scaled_mass"] == p - 3
        assert row["H_at_least_two_scaled_floor"] > p - 3
        assert row["proved"] is True


def test_p19_positive_residues_are_removed_and_zero_remains():
    row = p19_second_boundary_reduction()
    assert row["boundary_size"] == 16
    assert row["pair_deficit_budget"] == 240
    assert [r["u0"] for r in row["pair_survivors_before_new_floor"]] == [
        0,
        2,
        3,
        4,
        6,
    ]
    assert {r["u1"] for r in row["pair_survivors_before_new_floor"]} == {9}
    assert [r["zero_quotient_scaled_mean"] for r in row["positive_residue_rows"]] == [
        4,
        6,
        8,
        12,
    ]
    assert all(r["therefore_b_zero"] for r in row["positive_residue_rows"])
    assert all(r["excluded"] for r in row["positive_residue_rows"])
    assert row["residue_zero_minimum_row"]["phase_zero_profile"] == {0: 5, 16: 5}
    assert row["residue_zero_minimum_row"]["phase_one_profile"] == {2: 9, 16: 1}
    assert row["residue_zero_minimum_pair_slack"] == 34
    assert row["residue_zero_minimum_rejected_modulo_four"] is True
    census = row["residue_zero_exact_census"]
    assert census["phase_labelled_profile_count"] == 143
    assert census["global_shape_count"] == 75
    assert census["pair_slack_histogram"] == {
        0: 54,
        4: 37,
        8: 25,
        12: 13,
        16: 7,
        20: 4,
        24: 1,
        28: 1,
        32: 1,
    }
    assert row["p19_second_all_finite_endpoint_closed"] is False
    assert row["top_level_gates_changed"] is False


def test_p19_residue_zero_profile_census_is_complete_and_structured():
    census = p19_residue_zero_profiles()
    assert census["phase_zero_row_count"] == 60
    assert census["phase_one_row_count"] == 9
    assert census["undetermined_direction_histogram_by_slack"][0] == {
        1: 2,
        2: 27,
        3: 25,
    }
    assert census["undetermined_direction_histogram_by_slack"][12] == {
        2: 1,
        3: 6,
        4: 6,
    }
    assert all(int(row["pair_slack"]) % 4 == 0 for row in census["profiles"])


def test_prop15688_does_not_soft_close_any_top_level_gate():
    row = theorem_sharp_lift_and_p19()
    assert row["sharp_floor_all_odd_p_at_least_five"] is True
    assert row["closes_p19_endpoint"] is False
    assert row["closes_residual_ii"] is False
    assert row["closes_R1"] is False
    assert row["closes_type_I"] is False
    assert row["L_status"] == "OPEN"
    assert row["proved"] is True
