from fractions import Fraction

from e1_gmin_m4_prop15680 import (
    _exact_profile_dp,
    p37_mass_ten_exclusion,
    p37_pair_and_lift_ledger,
    paired_cube_boolean_quadratic_floor,
    slice_distance_support_floor,
    theorem_record,
)


def test_paired_cube_transition_and_boolean_floor():
    for p in (5, 7, 11, 17, 37, 41, 101):
        row = paired_cube_boolean_quadratic_floor(p)
        assert row["rho"] == Fraction(1, p + 1)
        assert row["cube_distance_floor"] == Fraction(1, 4)
        assert row["boolean_quadratic_density_floor"] == Fraction(
            p - 3, 4 * p
        )
        assert row["proved"] is True


def test_degree_two_and_four_slice_distance_values():
    assert slice_distance_support_floor(37, 2) == Fraction(171, 2590)
    assert slice_distance_support_floor(37, 4) == Fraction(1938, 441595)


def test_p37_mass_ten_lift_is_impossible():
    row = p37_mass_ten_exclusion()
    assert row["target_mean"] == Fraction(5, 74)
    assert row["maximum_point_value"] == 2
    assert row["value_two_density_upper"] == Fraction(2, 1295)
    assert row["degree_four_gap"] == Fraction(1256, 441595)
    assert row["value_two_excluded"] is True
    assert row["therefore_B_is_boolean"] is True
    assert row["paired_cube_boolean_floor"] == Fraction(17, 74)
    assert row["mass_ten_excluded"] is True


def test_exact_p37_pair_ledger_has_only_four_residues():
    phase_one = _exact_profile_dp(1)
    assert phase_one == [
        {
            "u": 18,
            "quotient_sum": 1,
            "minimum_deficit": 504,
            "profile": {2: 18, 30: 1},
        }
    ]
    ledger = p37_pair_and_lift_ledger()
    assert ledger["surviving_residues"] == [2, 3, 4, 5]
    assert [row["pair_slack"] for row in ledger["pair_rows"] if row["survives_pair_budget"]] == [
        38,
        36,
        8,
        6,
    ]
    assert ledger["old_nonzero_lift_floor"] == 10
    assert all(row["forces_quotient_zero"] for row in ledger["lift_rows"])
    assert all(row["therefore_b_zero"] for row in ledger["lift_rows"])
    assert all(row["excluded"] for row in ledger["lift_rows"])
    assert ledger["endpoint_excluded"] is True


def test_theorem_record_is_honestly_scoped():
    record = theorem_record()
    assert record["proved"] is True
    theorem = record["theorem"]
    assert theorem["remaining_smaller_endpoints"] == [17, 19, 23, 29, 31, 41]
    assert theorem["general_residual_ii"] is False
    assert theorem["R1"] is False
    assert theorem["limit_exists"] is False
