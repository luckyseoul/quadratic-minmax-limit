import json
from pathlib import Path

from e1_gmin_m4_prop15674 import (
    full_profile_branch_exclusion,
    full_profile_floor_ledger,
    full_profile_type_pair_ledger,
    residue_classification_ledger,
    theorem_record,
)
from e1_gmin_m4_prop15723 import floor_excess_admissible


def test_full_odd_fibre_floor_bands():
    expected = {
        (19, 0): ((1, 17), ()),
        (17, 1): ((15,), (1,)),
        (19, 1): ((), (1, 17)),
        (17, 0): ((1,), (15,)),
    }
    for key, bands in expected.items():
        row = full_profile_floor_ledger(*key)
        assert (row["floor_P_counts"], row["floor_P_minus_2_counts"]) == bands
        period = row["P"]
        intermediate_floors = {
            row["floors"][b] for b in row["intermediate_counts"]
        }
        assert intermediate_floors <= {2 * period - 8, 2 * period - 2}
        assert all(value > period for value in intermediate_floors)


def test_only_zero_and_p_minus_one_residues_can_survive():
    expected_actual = {
        (19, 0): [0],
        (17, 1): [0, 16],
        (19, 1): [18],
        (17, 0): [0, 16],
    }
    for key, residues in expected_actual.items():
        row = residue_classification_ledger(*key)
        assert row["proved"] is True
        assert row["only_possible_residues_before_branch_floor_availability"] == [
            0,
            key[0] - 1,
        ]
        assert row["surviving_residues_in_this_branch"] == residues
        assert row["intermediate_direction_limit_per_type"] == 1


def test_residue_classification_by_independent_relaxed_dp():
    # Enumerate every floor class and every common residue for the first prime
    # in each congruence class.  We allow every lift of size at least four, so
    # this is an independent relaxation of the symbolic proof.
    for p in (17, 19):
        for phase in (0, 1):
            floor = full_profile_floor_ledger(p, phase)
            period = floor["P"]
            m = floor["m"]
            for u in range(m):
                residue = 2 * u
                target = m - u
                states = {(0, 0, 0, 0)}
                for _ in range(m):
                    next_states = set()
                    for total, low, high, intermediate in states:
                        for b, direction_floor in floor["floors"].items():
                            kind = (
                                "low"
                                if direction_floor == period - 2
                                else "high"
                                if direction_floor == period
                                else "intermediate"
                            )
                            for k in range(target - total + 1):
                                excess = residue + period * k - direction_floor
                                if not floor_excess_admissible(
                                    p, b, phase, excess
                                ):
                                    continue
                                next_states.add(
                                    (
                                        total + k,
                                        low + int(kind == "low"),
                                        high + int(kind == "high"),
                                        intermediate + int(kind == "intermediate"),
                                    )
                                )
                    states = {state for state in next_states if state[0] <= target}
                feasible = {state for state in states if state[0] == target}
                if u not in (0, m - 1):
                    assert feasible == set()
                elif u == 0:
                    assert all(
                        low == 0 and high == m and intermediate == 0
                        for _, low, high, intermediate in feasible
                    )
                else:
                    assert all(
                        low >= m - 1 and intermediate <= 1
                        for _, low, _, intermediate in feasible
                    )


def test_p17_floor_plus_two_exceptions_are_retained_but_add_no_residue():
    for b, phase in ((5, 1), (11, 0)):
        floor = full_profile_floor_ledger(17, phase)
        assert 36 - floor["floors"][b] == 2
        residue = residue_classification_ledger(17, phase)
        assert [b, phase] in residue["floor_plus_two_cells_retained"]
        assert residue["surviving_residues_in_this_branch"] == [0, 16]


def test_geometry_forces_opposite_baseline_types_for_full_profiles():
    expected_offsets = {
        (19, 0): 0,
        (17, 1): 1,
        (19, 1): 2,
        (17, 0): 1,
    }
    for key, offset in expected_offsets.items():
        row = full_profile_type_pair_ledger(*key)
        assert row["proved"] is True
        assert row["two_b1_types"]["contradiction"] is True
        assert row["two_complement_types"]["contradiction"] is True
        assert row["forced_pair"] == [
            "b=1 baseline type",
            "b=p-2 baseline type",
        ]
        assert row["finite_edge_offset"] == offset


def test_same_four_arithmetic_rows_close_every_profile():
    for p in (17, 19, 23, 29, 31, 37, 41, 101):
        for phase in (0, 1):
            row = full_profile_branch_exclusion(p, phase)
            assert row["excluded"] is True
            if p == 17 and phase == 0:
                assert row["p17_l1"]["exact_minimum"] == 75
                assert row["p17_l1"]["transverse_edge_budget"] == 57


def test_theorem_scope_and_generated_evidence():
    row = theorem_record()
    assert row["proved"] is True
    assert row["theorem"]["all_directional_odd_fibre_profiles"] is True
    assert row["theorem"]["general_residual_ii"] is False
    assert row["theorem"]["R1"] is False

    root = Path(__file__).resolve().parents[1]
    stored = json.loads(
        (root / "evidence" / "e1_gmin_m4_prop15674.json").read_text()
    )
    assert stored == row
