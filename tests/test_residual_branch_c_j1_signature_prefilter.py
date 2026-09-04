import hashlib
import json
from pathlib import Path

from scripts.residual_branch_c_j1_signature_prefilter import (
    EXPECTED_OPTION_RECORD_SHA256,
    EXPECTED_TRANSFORMED_DISTANCE_HISTOGRAM,
    HARD_TARGET_ORDER,
    closed_option_catalog,
    transformed_distances,
)


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "evidence" / "resii_p31_j1_signature_prefilter_v1.json"
EVIDENCE_SHA256 = (
    "78aa77fb2554905865f5d0fc32182ff2ce42f4f4b778d6859b83ac103792226c"
)


def test_closed_option_table_has_pinned_counts_and_record_hash() -> None:
    options, digest = closed_option_catalog()
    assert tuple(map(len, options)) == (480, 450, 450, 480) + (450,) * 12
    assert sum(map(len, options)) == 7_260
    assert len(options) == len(HARD_TARGET_ORDER) == 16
    assert digest == EXPECTED_OPTION_RECORD_SHA256


def test_three_weight_six_counterexamples_replay_the_transformed_distance() -> None:
    row = json.loads(EVIDENCE.read_text())
    records = row["weight_six_aggregate_designs"]
    assert [record["component_id"] for record in records] == [5497, 8971, 9077]
    assert [record["aggregate_signature_hex"] for record in records] == [
        "0010a052",
        "18200083",
        "20d40400",
    ]
    for record in records:
        aggregate = int(record["aggregate_signature_hex"], 16)
        distances = transformed_distances(aggregate)
        assert aggregate.bit_count() == 6
        assert min(distances) == 6
        assert record["minimum_after_fixed_and_origin_double"] == 6
        assert record["minimizing_origin_directions"] == [
            index for index, value in enumerate(distances) if value == 6
        ]


def test_bounded_prefilter_evidence_and_scope_are_pinned() -> None:
    raw = EVIDENCE.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == EVIDENCE_SHA256
    row = json.loads(raw)
    assert row["component_catalog_file_sha256"] == (
        "196588c21a37c7788565b64c5b2a7dbfcafaedbd864dadf7e51b8b278895ae5b"
    )
    assert row["option_catalog_file_sha256"] == (
        "6c2c9dd2ca12d007f865c4499dbe038ef53215688eb43bb3759aef7d39daa599"
    )
    assert row["option_catalog_record_sha256"] == EXPECTED_OPTION_RECORD_SHA256
    assert row["minimum_transformed_distance_histogram"] == {
        str(key): value
        for key, value in EXPECTED_TRANSFORMED_DISTANCE_HISTOGRAM.items()
    }
    assert row["checkpoint_design_count"] == 10_000
    assert row["excluded_checkpoint_design_count"] == 10_000
    assert row["passing_checkpoint_design_count"] == 0
    assert row["bounded_checkpoint_f1_d1_excluded"]
    assert not row["component_exhausted"]
    assert not row["full_bfs_component_exhaustively_enumerated"]
    assert not row["global_f1_d1_branch_closed"]
    assert not row["j1_f3_d0_branch_addressed"]
    assert not row["residual_ii_closed"]
