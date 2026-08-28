import json
from pathlib import Path

from e1_gmin_m4_prop15632 import scaled_direction_floor
from e1_gmin_m4_prop15672 import (
    arithmetic_contradiction,
    both_signs_collinear_closed,
    exception_normal_form,
    lift_count_ledger,
    opposite_sign_floor_ledger,
    opposite_sign_near_line_exclusion,
    theorem_record,
)


def test_opposite_sign_floor_ledgers():
    for p in (11, 19, 23, 31):
        row = opposite_sign_floor_ledger(p)
        assert row["c_H"] == -1
        assert row["phase"] == 1
        assert row["line_type_surplus"] == p + 1
        assert row["opposite_type_surplus"] == p + 1
        assert scaled_direction_floor(p, 1, 1) == p - 1
        assert scaled_direction_floor(p, p - 2, 1) == p - 1

    for p in (13, 17, 29, 37):
        row = opposite_sign_floor_ledger(p)
        assert row["c_H"] == 1
        assert row["phase"] == 0
        assert row["line_type_surplus"] == p - 1
        assert row["opposite_type_surplus"] == p + 1
        assert scaled_direction_floor(p, 1, 0) == p + 1
        assert scaled_direction_floor(p, p - 2, 0) == p - 1


def test_four_lifts_exceed_every_type_surplus():
    for p in (11, 13, 17, 19, 23, 31, 101):
        row = lift_count_ledger(p)
        assert row["four_lifts_exceed_surplus"] is True
        assert row["four_lift_scaled_lower_bound"] > p + 1
        assert row["maximum_nonbaseline_directions_per_type"] == 3


def test_exact_two_exception_normal_form():
    for p in (11, 13, 17, 19, 23, 31):
        row = exception_normal_form(p)
        assert row["transverse_baseline_exists_in_each_type"] is True
        assert row["same_type_quantum"] == p + 1
        assert row["exceptions_per_type"] == 1
        assert row["exception_a"] == 2 * p
        assert row["exception_parallel_increment"] == 1


def test_baseline_congruence_contradiction():
    for p in (11, 13, 17, 19, 23, 29, 31, 101):
        row = arithmetic_contradiction(p)
        assert row["maximum_x_plus_y"] == 7
        assert row["substituted_congruences"] == [
            "q divides y+1",
            "q divides x+1",
        ]
        assert row["congruence_lower_bound_on_x_plus_y"] >= 8
        assert row["contradiction"] is True
        assert opposite_sign_near_line_exclusion(p)["excluded"] is True


def test_combined_both_sign_closure_and_open_flags():
    assert both_signs_collinear_closed(11)["both_signs_excluded"] is False
    for p in (13, 17, 19, 23, 29, 31, 101):
        assert both_signs_collinear_closed(p)["both_signs_excluded"] is True
    row = theorem_record()
    assert row["proved"] is True
    assert row["theorem"]["noncollinear_boundary"] == "OPEN"
    assert row["theorem"]["general_residual_ii"] is False
    assert row["theorem"]["R1"] is False


def test_generated_evidence_matches_source():
    root = Path(__file__).resolve().parents[1]
    stored = json.loads(
        (root / "evidence" / "e1_gmin_m4_prop15672.json").read_text()
    )
    expected = json.loads(json.dumps(theorem_record(), default=str))
    assert stored == expected
