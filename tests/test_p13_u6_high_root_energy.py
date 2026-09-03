"""Fail-when-wrong tests for the high-root p=13,u=6 energy close."""
from __future__ import annotations

from p13_u6_high_root_energy import (
    EXPECTED_MODEL_HASHES,
    ROW_SPECS,
    balanced_collision,
    collision_floor_certificate,
    high_root_partition_certificate,
    moment_relation_holds,
    nonexact_parseval_base,
    row_energy_certificate,
    row_satisfies,
    signed_collision_floor,
    translated_cuts,
)


def test_pinned_cut_catalog_and_collision_floor_are_exact() -> None:
    assert len(translated_cuts()) == 74
    assert balanced_collision(4, 6) == 0
    assert balanced_collision(7, 6) == 1
    assert balanced_collision(8, 6) == 2
    assert balanced_collision(9, 6) == 3
    assert [signed_collision_floor(value) for value in (-9, -8, -7, 6, 7, 8)] == [2, 1, 0, 0, 1, 2]
    row = collision_floor_certificate()
    assert row["positive_classes_per_nonzero_bucket"] == 6
    assert row["negative_classes_per_nonzero_bucket"] == 7
    assert row["proved"] is True


def test_every_explicit_maximizer_and_upper_replay_is_live() -> None:
    assert len(ROW_SPECS) == 13
    for spec in ROW_SPECS:
        assert row_satisfies(spec, spec.witness)
        assert sum(value * value for value in spec.witness) == spec.expected_energy
        row = row_energy_certificate(spec.name)
        expected_optimization, expected_replay = EXPECTED_MODEL_HASHES[spec.name]
        assert row["sharp_energy"] == spec.expected_energy
        assert row["optimization"] == {
            "status": "OPTIMAL",
            "workers": 1,
            "model_proto_sha256": expected_optimization,
            "solver_version": row["optimization"]["solver_version"],
        }
        assert row["upper_replay"]["status"] == "INFEASIBLE"
        assert row["upper_replay"]["workers"] == 1
        assert row["upper_replay"]["forbidden_energy_floor"] == spec.expected_energy + 1
        assert row["upper_replay"]["model_proto_sha256"] == expected_replay
        assert row["proved"] is True


def test_quartic_sign_is_fail_when_wrong() -> None:
    hard = next(spec for spec in ROW_SPECS if spec.name == "hard_e1_quartic")
    opposite = next(spec for spec in ROW_SPECS if spec.name == "opposite_q4_quartic")
    assert moment_relation_holds(hard, hard.witness)
    assert moment_relation_holds(opposite, opposite.witness)
    assert not moment_relation_holds(opposite, hard.witness)


def test_all_partition_bases_are_derived_from_common_parseval() -> None:
    expected = {
        (1, 1, 1, 1, 1): 303,
        (2, 1, 1, 1): 298,
        (2, 2, 1): 293,
        (3, 1, 1): 289,
        (3, 2): 284,
        (4, 1): 276,
        (5,): 259,
    }
    assert {partition: nonexact_parseval_base(partition) for partition in expected} == expected


def test_high_root_partition_ledgers_are_strict_and_narrowly_scoped() -> None:
    row = high_root_partition_certificate()
    assert row["closed_partitions"] == [[3, 2], [4, 1], [5]]
    assert [entry["excluded"] for entry in row["partitions"]] == [True, True, True]
    assert row["partitions"][0]["gap"] == 22
    assert [entry["gap"] for entry in row["partitions"][1]["collision_cases"]] == [28, 50]
    assert [entry["gap"] for entry in row["partitions"][2]["collision_cases"]] == [38, 42, 56, 58, 82]
    assert row["p13_t4_u6_fully_closed"] is False
    assert row["graph_or_configuration_census_used"] is False
    assert row["finite_six_bin_aggregate_models_used"] is True
    assert row["proved"] is True
