import json
import math
from pathlib import Path

from e1_gmin_m4_prop15670 import (
    FULL_COUNT,
    NORMALIZED_COUNT,
    PINNED_HISTOGRAM_SHA256,
    TYPE_BUDGET,
    normalization_ledger,
    theorem_record,
)


ROOT = Path(__file__).resolve().parents[1]


def test_ordered_pair_affine_normalization_has_exact_coverage_count():
    row = normalization_ledger()
    assert FULL_COUNT == math.comb(121, 8) == 899_749_479_915
    assert NORMALIZED_COUNT == math.comb(119, 6) == 3_470_108_187
    assert FULL_COUNT * 8 * 7 == NORMALIZED_COUNT * 121 * 120
    assert row["counting_identity"]
    assert row["nonsquare_scalar_swaps_direction_types"]
    assert row["both_c_H_signs_cover_the_swap"]
    similarity = row["similarity_audit"]
    assert similarity["valid"]
    assert similarity["similarity_direction_pairs_checked"] == 120 * 12
    assert similarity["translation_direction_pairs_checked"] == 121 * 12
    assert similarity["phase_transfer_cases_checked"] == 120 * 12 * 2


def test_independent_gpu_histograms_and_representatives_reaudit_exactly():
    row = theorem_record()
    assert row["proved"]
    assert row["pinned_file_sha256_match"]
    assert row["independent_gpu_core_fields_match"]
    assert {audit["c_H"] for audit in row["sign_audits"]} == {-1, 1}
    for audit in row["sign_audits"]:
        assert audit["valid"]
        assert audit["histogram_total"] == NORMALIZED_COUNT
        assert audit["histogram_sha256"] == PINNED_HISTOGRAM_SHA256[
            audit["c_H"]
        ]
        assert audit["floor_survivors"] == 0
        assert audit["minimum_maximum_type_cost"] == 76
        assert audit["minimum_budget_excess"] == 4


def test_prop15670_closes_only_the_finite_p11_size8_branch():
    row = theorem_record()
    theorem = row["theorem"]
    assert theorem["every_finite_p11_size8_boundary"] == "IMPOSSIBLE"
    assert theorem["exact_type_budget"] == TYPE_BUDGET == 72
    assert theorem["contradiction_gap"] == 4
    assert theorem["p11_infinity_plus_9"] == "OPEN"
    assert theorem["larger_boundaries"] == "OPEN"
    assert theorem["general_residual_ii"] is False
    assert theorem["R1"] is False
    assert theorem["global_QVAR"] is False
    assert theorem["type_I"] is False
    assert theorem["limit_exists"] is False


def test_generated_prop15670_record_matches_live_verifier():
    committed = json.loads(
        (ROOT / "evidence" / "e1_gmin_m4_prop15670.json").read_text()
    )
    assert committed == theorem_record()
