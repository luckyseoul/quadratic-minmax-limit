#!/usr/bin/env python3
"""Prop. 15.714 -- close positive p7 infinity+7 with z=0.

For z=0 every direction has phase-zero mean eight and a unique complete
Johnson-slice slack catalog. A complete CUDA scan of all C(49,7) finite
boundaries tests the 135 exact mod-seven left dependencies of the common
edge system. All 79,447,032 z=0 boundaries fail; none reaches an edge solve.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15713 import p7_positive_infinity_plus_seven_direction_reduction


ROOT = Path(__file__).resolve().parents[1]
SCAN = ROOT / "evidence" / "p7_infinity7_positive_z0_v100.json"
AUDIT = ROOT / "evidence" / "p7_infinity7_positive_z0_v100_audit.json"


def p7_positive_infinity_plus_seven_z0_exclusion() -> dict[str, object]:
    previous = p7_positive_infinity_plus_seven_direction_reduction()
    scan = json.loads(SCAN.read_text())
    audit = json.loads(AUDIT.read_text())
    expected_histogram = {
        "0": 79_447_032,
        "1": 6_324_528,
        "2": 123_480,
        "3": 5_488,
        "7": 56,
    }
    linear = scan["linear_system"]
    if (
        scan["experiment"] != "p7_infinity7_positive_z0_mod7_gpu"
        or scan["status"] != "complete_exact_mod_seven_z0_boundary_exhaustion"
        or int(scan["checked_boundaries"]) != 85_900_584
        or int(scan["all_boundaries"]) != 85_900_584
        or scan["undetermined_direction_histogram"] != expected_histogram
        or sum(expected_histogram.values()) != 85_900_584
        or int(scan["z0_boundaries"]) != 79_447_032
        or int(scan["mod7_survivors"]) != 0
        or scan["z0_branch_excluded"] is not True
        or int(linear["equations"]) != 282
        or int(linear["edge_variables"]) != 1225
        or int(linear["rank_mod_7"]) != 147
        or int(linear["left_dependency_dimension"]) != 135
        or linear["left_null_audit"] is not True
        or scan["cpu_prefix_verification"]["survivor_ranks"] != []
    ):
        raise ArithmeticError("p7 positive z0 complete scan changed")
    matching_keys = (
        "all_boundaries",
        "checked_boundaries",
        "z0_boundaries",
        "undetermined_direction_histogram",
        "mod7_survivors",
        "base_sha256",
        "tables_sha256",
    )
    if any(audit[key] != scan[key] for key in matching_keys) or audit["blocks"] == scan["blocks"]:
        raise ArithmeticError("independent-grid p7 z0 audit changed")

    projected_before = previous["remaining_undetermined_direction_histogram"]
    if projected_before != {0: 217, 1: 300, 2: 280, 3: 210, 7: 2}:
        raise ArithmeticError("pre-15.714 projected profile envelope changed")
    projected_after = {key: value for key, value in projected_before.items() if key}
    projected_count_after = sum(projected_after.values())
    actual_after = int(scan["all_boundaries"]) - int(scan["z0_boundaries"])
    if projected_count_after != 792 or actual_after != 6_453_552:
        raise ArithmeticError("post-15.714 positive branch count changed")

    return {
        "proposition": "15.714",
        "p": 7,
        "boundary": "infinity plus seven finite points",
        "product_sign": "positive",
        "undetermined_direction_count_excluded": 0,
        "all_finite_boundaries": int(scan["all_boundaries"]),
        "z0_boundaries_excluded": int(scan["z0_boundaries"]),
        "actual_boundary_count_after_z0_exclusion": actual_after,
        "actual_undetermined_direction_histogram": expected_histogram,
        "projected_b_profile_count_before": int(previous["projected_b_profile_count_after"]),
        "projected_b_profiles_excluded_here": 217,
        "projected_b_profile_count_after": projected_count_after,
        "remaining_projected_undetermined_direction_histogram": projected_after,
        "unique_catalog_means": 8,
        "modulus": 7,
        "left_dependency_count": int(linear["left_dependency_dimension"]),
        "mod7_survivors": int(scan["mod7_survivors"]),
        "scan_evidence": str(SCAN.relative_to(ROOT)),
        "independent_grid_audit_evidence": str(AUDIT.relative_to(ROOT)),
        "scan_script": "scripts/p7_infinity7_positive_z0_mod7_gpu.py",
        "positive_z0_branch_closed": True,
        "positive_p7_infinity_plus_seven_closed": False,
        "negative_p7_infinity_plus_seven_changed": False,
        "top_level_gates_changed": False,
        "proved_by_complete_exact_finite_scan": True,
    }


def main() -> None:
    theorem = p7_positive_infinity_plus_seven_z0_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15714.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.714: positive p7 infinity+7 z=0 boundaries "
        f"{theorem['z0_boundaries_excluded']} -> 0"
    )


if __name__ == "__main__":
    main()
