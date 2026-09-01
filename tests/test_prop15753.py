"""Fail-when-wrong tests for Proposition 15.753."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

import e1_gmin_m4_prop15753 as prop15753
from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
from e1_gmin_m4_prop15753 import (
    ROW_SPECS,
    ROW_SPEC_BY_NAME,
    branch_normalization_certificate,
    branch_ledger_certificate,
    cut_catalog_certificate,
    endpoint_arithmetic,
    hard_residue_ledger,
    moment_relation_holds,
    p17_opposite_sign_regression,
    proposition_15753,
    replay_all_row_certificates,
    replay_row_certificate,
    row_model_sha256,
    validate_witness,
    write_evidence,
)


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "evidence" / "e1_gmin_m4_prop15753.json"


def test_exact_translated_cut_catalogs_are_pinned() -> None:
    expected = {
        17: (698, 72, "a8ac7349cb601db5163ef1526949587c766914d774fe26858fe93eac1d940708"),
        19: (2338, 90, "5f07e9ced107e6dc1551b806043a92147c00d80eb009b70d0cbfd3ce9631c5b7"),
    }
    for p, (count, row_sum, digest) in expected.items():
        row = cut_catalog_certificate(p)
        assert row["distinct_translated_cut_vectors"] == count
        assert row["every_vector_sum"] == row_sum
        assert row["catalog_sha256"] == digest
        assert row["proved"] is True


def test_all_claimed_maximizers_replay_in_exact_integer_arithmetic() -> None:
    finite = 0
    empty = 0
    for spec in ROW_SPECS:
        if spec.expected_energy is None:
            empty += 1
            assert spec.witness is None
            continue
        finite += 1
        assert spec.witness is not None
        row = validate_witness(spec, spec.witness)
        assert row["energy"] == spec.expected_energy
        assert row["maximum_translated_cut"] <= spec.cut_upper
        assert row["proved"] is True
    assert finite == 17
    assert empty == 2


def test_p17_opposite_quartic_uses_the_required_minus_sign() -> None:
    row = p17_opposite_sign_regression()
    spec = ROW_SPEC_BY_NAME["p17_A_opposite_Q4_quartic"]
    wrong_plus_sign = tuple(row["wrong_plus_sign_witness"])
    assert spec.moment_relation == "opposite_quartic"
    assert moment_relation_holds(spec, spec.witness) is True
    assert moment_relation_holds(spec, wrong_plus_sign) is False
    assert row["correct_opposite_energy"] == 11
    assert row["wrong_plus_sign_energy"] == 15
    assert row["wrong_plus_sign_witness_rejected_by_correct_relation"] is True
    with pytest.raises(ArithmeticError):
        validate_witness(spec, wrong_plus_sign)


def test_endpoint_arithmetic_is_exactly_the_fifth_shell() -> None:
    p17 = endpoint_arithmetic(17)
    p19 = endpoint_arithmetic(19)
    assert (p17["layer_t"], p17["original_k"], p17["H_edge_count"]) == (4, 76, 77)
    assert (p19["layer_t"], p19["original_k"], p19["H_edge_count"]) == (4, 84, 85)
    assert p17["maximum_baseline_lift_excess"] < p17["sharp_integral_lift_floor"]
    assert p19["maximum_baseline_lift_excess"] < p19["sharp_integral_lift_floor"]


def test_hard_residue_ledger_derives_exactly_the_four_named_branches() -> None:
    p17 = hard_residue_ledger(17)
    p19 = hard_residue_ledger(19)
    assert p17["possible_branches"] == ["A_XNOR", "B_LITERAL"]
    assert p17["u_0_through_t_rows"][0]["surviving_branch"] == "B_LITERAL"
    assert all(
        row["surviving_branch"] is None for row in p17["u_0_through_t_rows"][1:]
    )
    assert p17["u_equals_m_minus_1_endpoint_b_candidates"] == [2]
    assert all(
        row["k_zero_forbidden"]
        and (
            row["excluded_by_quotient_sum"]
            if row["surviving_branch"] is None
            else row["minimum_k_after_lift_sieve"] == 1
        )
        for row in p17["u_0_through_t_rows"]
    )
    assert p19["possible_branches"] == ["A_XNOR", "C_COMPLEMENT_LITERAL"]
    assert all(row["surviving_branch"] is None for row in p19["u_0_through_t_rows"])
    assert p19["u_equals_m_minus_1_endpoint_b_candidates"] == [2, 18]
    assert p19["equal_mean_endpoint_cells_cannot_mix"] is True
    assert all(
        row["k_zero_forbidden"]
        and row["minimum_k_after_lift_sieve"] == 2
        and row["excluded_by_quotient_sum"]
        for row in p19["u_0_through_t_rows"]
    )
    assert all(
        row["k_zero_forbidden"]
        and row["quotient_sum_less_than_direction_count"]
        for endpoint in (p17, p19)
        for row in endpoint["intermediate_u_rows"]
    )


@pytest.mark.parametrize(
    "p,branch,hT,affine,energy",
    [
        (17, "A_XNOR", 5, "P_L=4+k_L", 1),
        (17, "B_LITERAL", 21, "P_L=4+k_L", 32),
        (19, "A_XNOR", 5, "P_L=4+k_L", 1),
        (19, "C_COMPLEMENT_LITERAL", -15, "P_L=3+k_L", 36),
    ],
)
def test_common_edge_total_forces_each_parallel_normalization(
    p: int, branch: str, hT: int, affine: str, energy: int
) -> None:
    row = branch_normalization_certificate(p, branch)
    assert row["hT_from_edge_split"] == hT
    assert row["forced_affine_normalization"] == affine
    assert row["exact_baseline_row_energy"] == energy
    assert row["exact_baseline_agrees_with_independent_edge_split"] is True
    assert all(
        entry["matching_parallel_counts_in_full_range"]
        == [entry["forced_parallel_count"]]
        and entry["local_sum"] == entry["common_sum"]
        for entry in row["normalization_rows"]
    )


def test_all_four_branch_ledgers_are_strict() -> None:
    row = branch_ledger_certificate()
    p17_A = row["p17"]["branch_A"]
    p17_B = row["p17"]["branch_B"]
    p19_A = row["p19"]["branch_A"]
    p19_C = row["p19"]["branch_C"]
    assert [entry["strict_gap_at_C_zero"] for entry in p17_A["partition_ledgers"]] == [
        342,
        312,
        282,
        212,
        182,
        138,
        162,
    ]
    assert [entry["strict_gap_at_C_zero"] for entry in p17_B["remaining_partition_ledgers"]] == [428, 302]
    assert [entry["strict_gap_at_C_zero"] for entry in p19_A["partition_ledgers"]] == [
        520,
        490,
        460,
        420,
        390,
        312,
        162,
    ]
    assert p19_C["opposite_Q5_row_system_infeasible"] is True
    for branch, nonzero_floor, lift_floor in (
        (p17_A, 16, 14),
        (p19_A, 20, 16),
    ):
        exclusion = branch["opposite_Q3_mean_8_exclusion"]
        assert exclusion == {
            "nonzero_b_floor": nonzero_floor,
            "b_zero_integral_lift_floor": lift_floor,
            "proved": True,
        }
    assert p19_C["opposite_Q4_mean_8_exclusion"] == {
        "nonzero_b_floor": 20,
        "b_zero_integral_lift_floor": 16,
        "proved": True,
    }
    assert row["p17"]["p17_k76_closed"] is True
    assert row["p19"]["p19_k84_closed"] is True


@pytest.mark.parametrize(
    "name",
    tuple(ROW_SPEC_BY_NAME),
)
def test_every_row_model_replays_live_with_one_worker(name: str) -> None:
    row = replay_row_certificate(name)
    assert row["num_search_workers"] == 1
    assert row["status"] == "INFEASIBLE"
    assert row["model_proto_sha256"] == row_model_sha256(ROW_SPEC_BY_NAME[name])
    assert row["proved"] is True


def test_parallel_replay_covers_the_exact_model_set() -> None:
    rows = replay_all_row_certificates()
    assert tuple(rows) == tuple(ROW_SPEC_BY_NAME)
    assert all(
        row["status"] == "INFEASIBLE"
        and row["num_search_workers"] == 1
        and row["proved"] is True
        for row in rows.values()
    )


def test_checked_in_evidence_contains_every_exact_one_worker_replay() -> None:
    row = json.loads(EVIDENCE.read_text())
    assert row["all_row_models_replayed"] is True
    assert row["independent_row_models_run_concurrently"] is True
    assert row["maximum_concurrent_row_models"] == len(ROW_SPECS)
    assert set(row["exact_one_worker_row_replays"]) == set(ROW_SPEC_BY_NAME)
    assert all(
        replay["status"] == "INFEASIBLE"
        and replay["num_search_workers"] == 1
        and replay["proved"] is True
        for replay in row["exact_one_worker_row_replays"].values()
    )
    assert row["finite_prime_aggregate_census_used"] is True
    assert row["graph_or_configuration_census_used"] is False
    assert row["p17_k76_closed"] is True
    assert row["p19_k84_closed"] is True


def test_evidence_model_and_supporting_artifact_hashes_are_live() -> None:
    row = json.loads(EVIDENCE.read_text())
    manifest = row["certificate_manifest"]
    assert manifest["row_model_proto_sha256"] == {
        spec.name: row_model_sha256(spec) for spec in ROW_SPECS
    }
    for relative, expected_hash in manifest["supporting_artifact_sha256"].items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == expected_hash


def test_endpoint_close_does_not_flip_the_global_predicate() -> None:
    row = proposition_15753()
    assert row["p17_k76_closed"] is True
    assert row["p19_k84_closed"] is True
    assert row["finite_prime_aggregate_census_used"] is True
    assert row["graph_or_configuration_census_used"] is False
    assert row["residual_ii_k_ge_4p_ND_closed"] is False
    assert row["E1_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert residual_ii_k_ge_4p_ND_closed() is False


def test_write_evidence_is_atomic(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    payload = {"proved": True, "sentinel": "prop15753"}
    monkeypatch.setattr(prop15753, "proposition_15753", lambda **_: payload)
    path = tmp_path / "prop15753.json"
    assert write_evidence(path) == path
    assert json.loads(path.read_text()) == payload
    assert not list(tmp_path.glob("*.tmp"))
