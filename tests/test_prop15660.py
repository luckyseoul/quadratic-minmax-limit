import hashlib
from pathlib import Path

from e1_gmin_m4_prop15660 import (
    CLASS_AUDIT_SHA256,
    GLOBAL_AUDIT_SHA256,
    p5_size_six_class_certificates,
    p5_size_six_global_reduction,
    theorem_p5_size_six_exclusion,
)

ROOT = Path(__file__).resolve().parents[1]


def test_global_reduction_has_exactly_six_classes():
    row = p5_size_six_global_reduction()
    assert row["catalogs_rebuilt_from_definitions"]
    assert row["selection_reduction_proved"]
    assert row["residual_classes_after_symmetry"] == 6
    assert row["no_infinity_sign_transfer_bijective"]


def test_all_six_class_audits_are_recorded():
    row = p5_size_six_class_certificates()
    assert row["proved"]
    assert row["closed_class_count"] == 6
    assert set(row["classes"]) == set(CLASS_AUDIT_SHA256)


def test_theorem_closes_only_p5_size_six():
    row = theorem_p5_size_six_exclusion()
    assert row["proved"]
    assert row["p5_size_six"] == "CLOSED"
    assert row["p7_six_finite"] == "OPEN"
    assert not row["closes_all_size_six"]
    assert not row["closes_residual_ii"]
    assert not row["closes_R1"]
    assert row["L_status"] == "OPEN"


def test_compact_audit_hashes_are_pinned():
    global_path = ROOT / "evidence" / "p5_size6_global_circle_coverage_audit.json"
    assert hashlib.sha256(global_path.read_bytes()).hexdigest() == GLOBAL_AUDIT_SHA256
    for index, expected in CLASS_AUDIT_SHA256.items():
        name = (
            "p5_size6_circle_scip_audit_i0.json"
            if index == 0
            else f"p5_size6_general_circle_scip_audit_i{index}.json"
        )
        assert hashlib.sha256((ROOT / "evidence" / name).read_bytes()).hexdigest() == expected
