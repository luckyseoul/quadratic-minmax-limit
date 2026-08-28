import json
from pathlib import Path

from e1_gmin_m4_prop15671 import (
    coefficient_ledger,
    near_line_floor_ledger,
    rigid_near_line_exclusion,
    theorem_record,
)
from e1_gmin_m4_prop15632 import scaled_direction_floor


def test_rigid_floor_ledgers():
    for p in (19, 23, 31, 43):
        row = near_line_floor_ledger(p)
        assert row["c_H"] == 1
        assert row["phase"] == 0
        assert row["line_type_surplus"] == 0
        assert row["opposite_type_surplus"] == 0
        assert row["all_directional_slacks_forced_to_baseline"] is True
        assert scaled_direction_floor(p, 1, 0) == p + 1
        assert scaled_direction_floor(p, p - 2, 0) == p + 1

    for p in (13, 17, 29, 37):
        row = near_line_floor_ledger(p)
        assert row["c_H"] == -1
        assert row["phase"] == 1
        assert row["line_type_surplus"] == 2
        assert row["opposite_type_surplus"] == 0
        assert row["minimum_nonzero_lift_scaled_cost"] > 2
        assert row["all_directional_slacks_forced_to_baseline"] is True
        assert scaled_direction_floor(p, 1, 1) == p - 1
        assert scaled_direction_floor(p, p - 2, 1) == p + 1


def test_summed_coefficient_congruence():
    assert coefficient_ledger(13)["summed_congruence"] == "I = 2 (mod q)"
    assert coefficient_ledger(19)["summed_congruence"] == "I = 4 (mod q)"


def test_symbolic_exclusions_and_small_boundary():
    for p in (13, 17, 29, 37, 41, 101):
        row = rigid_near_line_exclusion(p)
        assert row["applicable"] is True
        assert row["excluded"] is True
        assert row["arithmetic"]["q_even"] is True

    for p in (19, 23, 31, 43, 47, 59, 103):
        row = rigid_near_line_exclusion(p)
        assert row["applicable"] is True
        assert row["excluded"] is True
        assert row["arithmetic"]["l1_forces_all_a_zero"] is True
        assert row["arithmetic"]["edge_count_forces_k0"] == 8

    assert rigid_near_line_exclusion(7)["excluded"] is False
    assert rigid_near_line_exclusion(11)["excluded"] is False


def test_theorem_record_keeps_top_level_flags_open():
    row = theorem_record()
    assert row["proved"] is True
    assert row["theorem"]["opposite_product_sign"] == "OPEN"
    assert row["theorem"]["general_residual_ii"] is False
    assert row["theorem"]["R1"] is False


def test_generated_evidence_matches_source():
    root = Path(__file__).resolve().parents[1]
    stored = json.loads(
        (root / "evidence" / "e1_gmin_m4_prop15671.json").read_text()
    )
    assert stored == theorem_record()
