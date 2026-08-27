import hashlib
from pathlib import Path

from e1_gmin_m4_prop15662 import (
    GLOBAL_CONIC_AUDIT_SHA256,
    HIGH_MEAN_EXCEPTIONAL_AUDIT_SHA256,
    ORDINARY_EXCEPTIONAL_AUDIT_SHA256,
    p7_size_eight_conic_certificate,
    theorem_p7_size_eight_conic_subbranch,
)

ROOT = Path(__file__).resolve().parents[1]


def test_p7_size_eight_conic_partition_and_closure():
    row = p7_size_eight_conic_certificate()
    assert row["all_size_eight_boundaries_per_sign"] == 450_978_066
    assert row["floor_surviving_boundaries_per_sign"] == 108_754_569
    assert row["minimum_odd_secant_conic_boundaries_per_sign"] == 6_174
    assert row["floor_rejected_conic_boundaries_per_sign"] == 4_851
    assert row["floor_surviving_conic_boundaries_per_sign"] == 1_323
    assert row["saturated_orbits"] + row["exceptional_orbits"] == 32
    assert row["saturated_boundaries"] + row["exceptional_boundaries"] == 1_323
    assert (
        row["saturated_initial_cp_exclusions"]
        + row["saturated_long_cp_exclusions"]
        + row["saturated_catalog_join_exclusions"]
        == row["saturated_mean_allocations"]
    )
    assert (
        row["exceptional_initial_cp_exclusions"]
        + row["exceptional_ordinary_gpu_exclusions"]
        + row["exceptional_high_direction_omission_gpu_exclusions"]
        == row["exceptional_mean_allocations"]
    )
    assert row["remaining_conic_mean_allocations"] == 0
    assert row["nonsquare_sign_transfer"] is True
    assert row["nonconic_floor_survivors_per_sign"] == 108_753_246


def test_prop15662_scope_is_exact():
    row = theorem_p7_size_eight_conic_subbranch()
    assert row["proved"] is True
    assert row["p7_size_eight_minimum_odd_secant_conic_subbranch_both_signs"] == "CLOSED"
    assert row["all_32_floor_surviving_conic_orbits_both_signs"] == "CLOSED"
    assert row["full_p7_size_eight"] == "OPEN"
    assert row["closes_all_p7_size_eight"] is False
    assert row["closes_residual_ii"] is False
    assert row["closes_R1"] is False
    assert row["L_status"] == "OPEN"


def test_prop15662_compact_audit_hashes_are_pinned():
    paths = {
        ROOT / "evidence" / "p7_size8_cminus1" / "global_conic_audit.json": GLOBAL_CONIC_AUDIT_SHA256,
        ROOT / "evidence" / "p7_exceptional_mod7triple_all" / "audit.json": ORDINARY_EXCEPTIONAL_AUDIT_SHA256,
        ROOT / "evidence" / "p7_exceptional_high_direction_omission" / "audit.json": HIGH_MEAN_EXCEPTIONAL_AUDIT_SHA256,
    }
    for path, expected in paths.items():
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
