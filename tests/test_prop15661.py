import hashlib
from pathlib import Path

from e1_gmin_m4_prop15661 import (
    GLOBAL_AUDIT_SHA256,
    NUKA_SUMMARY_SHA256,
    p7_six_finite_certificate,
    theorem_size_six_all_odd_primes,
)

ROOT = Path(__file__).resolve().parents[1]


def test_p7_six_finite_partition_and_closure():
    row = p7_six_finite_certificate()
    assert row["all_boundaries"] == 13_983_816
    assert row["floor_survivors"] == 3_856_300
    assert row["square_semilinear_orbits"] == 80_704
    assert row["ordinary_orbits"] + row["deep_orbits"] == 80_704
    assert row["ordinary_elevation_cases"] == 160_745
    assert row["ordinary_survivors"] == 0
    assert row["deep_initial_infeasible_orbits"] == 92
    assert row["deep_initial_unknown_orbits"] == 93
    assert row["deep_allocation_infeasible_leaves"] == 810
    assert row["deep_low_catalog_join_leaves"] == 120
    assert row["deep_low_catalog_join_survivors"] == 0
    assert row["nonsquare_sign_transfer"] is True
    assert row["proved"] is True


def test_all_size_six_scope_is_exact():
    row = theorem_size_six_all_odd_primes()
    assert row["proved"] is True
    assert row["p7_six_finite_both_product_signs"] == "CLOSED"
    assert row["all_size_six_boundaries_for_odd_p_ge_5"] == "CLOSED"
    assert row["boundaries_size_at_least_eight"] == "OPEN"
    assert row["closes_residual_ii"] is False
    assert row["closes_R1"] is False
    assert row["L_status"] == "OPEN"


def test_compact_evidence_hashes_are_pinned():
    global_path = ROOT / "evidence" / "p7_size6_finite_global_audit.json"
    nuka_path = ROOT / "evidence" / "p7_size6_nuka_independent_summary.json"
    assert hashlib.sha256(global_path.read_bytes()).hexdigest() == GLOBAL_AUDIT_SHA256
    assert hashlib.sha256(nuka_path.read_bytes()).hexdigest() == NUKA_SUMMARY_SHA256
