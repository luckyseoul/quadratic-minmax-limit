"""Fail-when-wrong tests for Proposition 15.754."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

import pytest

import e1_gmin_m4_prop15754 as prop15754
from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
from e1_gmin_m4_prop15754 import (
    ARTIFACT_PATHS,
    EVIDENCE_PATH,
    EXPECTED_ARTIFACT_SHA256,
    EXPECTED_PARTITIONS,
    EXPECTED_ROW_CATALOGS,
    EXPECTED_SUPPORTING_CODE_SHA256,
    complete_partition_ledger,
    necessary_moment_artifact_certificate,
    p13_u6_normalization_certificate,
    proposition_15754,
    translated_cut_coordinate_certificate,
    validate_four_root_payload,
    validate_low_root_payload,
    write_evidence,
)


ROOT = Path(__file__).resolve().parents[1]


def _load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())


def test_normalization_derives_exactly_the_seven_partitions() -> None:
    row = p13_u6_normalization_certificate()
    assert (row["p"], row["t"], row["k"], row["u"]) == (13, 4, 60, 6)
    assert row["prior_remaining_residues"] == [6]
    assert row["hard_mean"] == "a_L=12+14*k_L"
    assert row["hard_quotient_sum"] == 5
    assert row["exact_floor_cell"] == "b=2 XNOR"
    assert row["exact_parallel_candidates"] == [4]
    assert row["hT"] == 5
    assert row["hard_parallel_edge_total"] == 33
    assert row["opposite_parallel_edge_total"] == 28
    assert row["opposite_parallel_counts"] == [4] * 7
    partitions = row["hard_excess_partitions"]
    assert [tuple(entry["hard_excess_partition"]) for entry in partitions] == list(
        EXPECTED_PARTITIONS
    )
    assert [entry["exact_XNOR_root_count"] for entry in partitions] == [
        2,
        3,
        4,
        4,
        5,
        5,
        6,
    ]
    assert [entry["collision_minimum"] for entry in partitions] == [0, 0, 0, 1, 1, 2, 3]
    assert [entry["parseval_nonexact_base"] for entry in partitions] == [
        303,
        298,
        293,
        289,
        284,
        276,
        259,
    ]
    assert row["proved"] is True


def test_translated_cut_duals_make_every_row_box_complete() -> None:
    row = translated_cut_coordinate_certificate()
    assert row["translated_cut_count"] == 74
    assert row["lower_dual_identity"] == ["1", "0", "0", "0", "0", "0"]
    assert row["upper_dual_identity"] == ["-1", "0", "0", "0", "0", "0"]
    assert row["cut_catalog_invariant_under_distance_action"] is True
    assert row["distance_action_transitive"] is True
    assert {
        name: value["integral_bounds"] for name, value in row["row_bounds"].items()
    } == {"H1": [-3, 3], "H2": [-5, 5], "H3": [-7, 7], "O": [-4, 1]}
    assert {
        name: value["enumerated_superset_box"]
        for name, value in row["H3_collision_enumeration_boxes"].items()
    } == {
        "H3_D0": [-7, 6],
        "H3_D1": [-8, 7],
        "H3_D2": [-9, 8],
        "H3_D3": [-10, 9],
    }
    assert row["proved"] is True


def test_low_root_artifact_closes_all_six_collision_strata() -> None:
    row = validate_low_root_payload(_load(ARTIFACT_PATHS["low_root_UG"]))
    assert row["closed_partitions"] == [[1, 1, 1, 1, 1], [2, 1, 1, 1]]
    assert [case["fixed_models_checked"] for case in row["cases"]] == [42, 280]
    assert [case["coefficient_pairs_checked"] for case in row["cases"]] == [
        43_184_232,
        10_221_120,
    ]
    assert [case["maximum_separable_nonexact_energy"] for case in row["cases"]] == [
        293,
        290,
    ]
    assert [case["strict_deficits"] for case in row["cases"]] == [
        [10, 36, 62],
        [8, 34, 60],
    ]
    assert row["graph_or_configuration_census_used"] is False
    assert row["proved"] is True


@pytest.mark.parametrize(
    "name,partition,collisions,maximum,compatible",
    [
        ("four_root_221", (2, 2, 1), [0, 1, 2, 3], 193, 84),
        ("four_root_311", (3, 1, 1), [1, 2, 3, 4], -1, 0),
    ],
)
def test_four_root_artifacts_exhaust_every_common_form(
    name: str,
    partition: tuple[int, ...],
    collisions: list[int],
    maximum: int,
    compatible: int,
) -> None:
    row = validate_four_root_payload(_load(ARTIFACT_PATHS[name]), partition)
    assert row["closed_collision_counts"] == collisions
    assert row["maximum_separable_nonexact_energy"] == maximum
    assert row["hard_locally_compatible_form_assignments_per_collision"] == compatible
    assert row["coverage"] == {
        "UGJ6_coefficient_triples_before_row_pruning": 218_320_284,
        "UGJ6_pairs_after_all_seven_opposite_rows": 336,
        "U_coefficients_passing_exact_root_QR": 7_644,
        "hard_excess_assignment_checks": 1_008,
        "hard_sign_root_sets_checked": 70,
    }
    assert row["sign_safe_opposite_key"] == (
        "(U,G,J6)=(-N2,-N4-N2^2,-N6+N2^3)"
    )
    assert row["W6_regression"] == [1, 12, 1, 1, 12, 12]
    assert row["proved"] is True


def test_artifact_and_supporting_code_hashes_are_pinned_live() -> None:
    row = necessary_moment_artifact_certificate()
    assert row["artifact_sha256"] == EXPECTED_ARTIFACT_SHA256
    assert row["supporting_code_sha256"] == EXPECTED_SUPPORTING_CODE_SHA256
    assert {
        name: hashlib.sha256(path.read_bytes()).hexdigest()
        for name, path in ARTIFACT_PATHS.items()
    } == EXPECTED_ARTIFACT_SHA256
    assert {
        name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        for name in EXPECTED_SUPPORTING_CODE_SHA256
    } == EXPECTED_SUPPORTING_CODE_SHA256


def test_four_root_sign_and_w6_regressions_fail_loudly() -> None:
    original = _load(ARTIFACT_PATHS["four_root_221"])
    wrong_sign = copy.deepcopy(original)
    wrong_sign["common_forms"]["opposite_key"] = (
        "(U,G,J6)=(N2,N4-N2^2,N6-N2^3)"
    )
    with pytest.raises(ArithmeticError):
        validate_four_root_payload(wrong_sign, (2, 2, 1))

    overflowed = copy.deepcopy(original)
    overflowed["W6_regression"][-1] = 3
    with pytest.raises(ArithmeticError):
        validate_four_root_payload(overflowed, (2, 2, 1))


def test_low_root_manifest_corruption_fails_loudly() -> None:
    payload = _load(ARTIFACT_PATHS["low_root_UG"])
    payload["source_checker_sha256"] = "0" * 64
    with pytest.raises(ArithmeticError):
        validate_low_root_payload(payload)


def test_all_seven_partitions_and_only_their_possible_collisions_close() -> None:
    row = complete_partition_ledger()
    ledgers = row["partition_ledgers"]
    assert [tuple(entry["partition"]) for entry in ledgers] == list(
        EXPECTED_PARTITIONS
    )
    assert [entry["raw_independent_energy_upper"] for entry in ledgers[:4]] == [
        357,
        364,
        371,
        417,
    ]
    assert [entry["raw_collision_maximum"] for entry in ledgers[:4]] == [2, 2, 3, 4]
    assert row["independent_top_collision_cross_checks"] == {
        "partition_1^5_C2": True,
        "partition_2_2_1_C3": True,
    }
    assert set(row["row_energy_maxima_from_complete_catalogs"]) == set(
        EXPECTED_ROW_CATALOGS
    )
    assert row["all_seven_partitions_closed"] is True
    assert row["proved"] is True


def test_finite_endpoint_close_does_not_flip_the_global_predicate() -> None:
    row = proposition_15754()
    assert row["result_status"] == (
        "exhaustive finite aggregate/common-form certificate and proved endpoint theorem"
    )
    assert row["p13_t4_u6_closed"] is True
    assert row["p13_k_eq_60_closed"] is True
    assert row["remaining_p13_t4_residues"] == []
    assert row["fifth_shell_k_eq_4p_plus_8_closed_for_every_prime_p_ge_13"] is True
    assert row["finite_prime_aggregate_census_used"] is True
    assert row["finite_common_form_coefficient_certificate_used"] is True
    assert row["graph_or_configuration_census_used"] is False
    assert row["orbit_census_used"] is False
    assert row["coefficient_cell_census_used"] is False
    assert row["residual_ii_k_ge_4p_ND_closed"] is False
    assert row["E1_closed"] is False
    assert row["quadratic_minmax_limit_closed"] is False
    assert row["top_level_gates_changed"] is False
    assert residual_ii_k_ge_4p_ND_closed() is False
    assert row["proved"] is True


def test_checked_in_main_evidence_equals_the_live_certificate() -> None:
    assert _load(EVIDENCE_PATH) == proposition_15754()


def test_write_evidence_is_atomic(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = {"proved": True, "sentinel": "prop15754"}
    monkeypatch.setattr(prop15754, "proposition_15754", lambda: payload)
    target = tmp_path / "prop15754.json"
    assert write_evidence(target) == target
    assert json.loads(target.read_text()) == payload
    assert not list(tmp_path.glob("*.tmp"))


def test_binding_docs_record_15754_without_a_global_overclaim() -> None:
    paths = [
        "AGENTS.md",
        "STATUS.md",
        "HANDOFF.md",
        "LONG_HORIZON_GOAL.md",
        "README.md",
        "solution.md",
        "evidence/PROPOSITION_DEDUP_AUDIT_2026-08-30.md",
        "evidence/share/denseness_path_package.md",
        "evidence/P0_ENGINEERING_GRAPH.md",
    ]
    for relative in paths:
        text = (ROOT / relative).read_text()
        assert "15.754" in text, relative
    status = (ROOT / "STATUS.md").read_text()[:20_000]
    handoff = (ROOT / "HANDOFF.md").read_text()[:25_000]
    assert "p=13,k=60,u=6" in status
    assert "p=13,k=60,u=6" in handoff
    assert "Residual (ii)" in status and "OPEN" in status
    assert "Residual (ii)" in handoff and "OPEN" in handoff
