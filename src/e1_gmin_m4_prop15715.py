#!/usr/bin/env python3
"""Prop. 15.715 -- close positive p7 infinity+7 with z=1.

There are exactly four mean allocations for each z=1 boundary.  A complete
CUDA scan projects them onto 23 exact mod-seven dependencies, and the host
checks every projected boundary against the complete catalogs on all 135
dependencies.  None of the 6,324,528 z=1 boundaries survives.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15714 import p7_positive_infinity_plus_seven_z0_exclusion


ROOT = Path(__file__).resolve().parents[1]
SCAN = ROOT / "evidence" / "p7_infinity7_positive_z1_v100.json"
AUDIT = ROOT / "evidence" / "p7_infinity7_positive_z1_v100_audit.json"


def p7_positive_infinity_plus_seven_z1_exclusion() -> dict[str, object]:
    previous = p7_positive_infinity_plus_seven_z0_exclusion()
    scan = json.loads(SCAN.read_text())
    audit = json.loads(AUDIT.read_text())
    expected_digest = "23de1f85d34f641d06279e8cdbc17fe6615fcd98198885e18f695d1812982b4c"
    linear = scan["linear_system"]
    if (
        scan["experiment"] != "p7_infinity7_positive_z1_mod7_gpu"
        or scan["status"] != "complete_projected_then_exact_mod_seven_z1_exhaustion"
        or int(scan["all_boundaries"]) != 85_900_584
        or int(scan["checked_boundaries"]) != 85_900_584
        or int(scan["z1_boundaries"]) != 6_324_528
        or int(scan["mean_allocation_count_per_boundary"]) != 4
        or int(scan["projected_dependency_count"]) != 23
        or int(scan["projected_survivors"]) != 1_326
        or scan["projected_survivor_rank_sha256"] != expected_digest
        or int(scan["all_dependency_survivors"]) != 0
        or scan["all_dependency_survivor_ranks"] != []
        or scan["z1_branch_excluded"] is not True
        or int(linear["equations"]) != 282
        or int(linear["edge_variables"]) != 1_225
        or int(linear["rank_mod_7"]) != 147
        or int(linear["left_dependency_dimension"]) != 135
        or linear["left_null_audit"] is not True
        or scan["catalog_row_histogram_by_direction_mask"]
        != {"1764": 232, "2233": 280}
    ):
        raise ArithmeticError("p7 positive z1 complete scan changed")

    matching_keys = (
        "all_boundaries",
        "checked_boundaries",
        "z1_boundaries",
        "mean_allocation_count_per_boundary",
        "projected_dependency_count",
        "projected_survivors",
        "projected_survivor_rank_sha256",
        "first_projected_survivor_ranks",
        "all_dependency_survivors",
        "all_dependency_survivor_ranks",
        "catalog_row_histogram_by_direction_mask",
    )
    if any(audit[key] != scan[key] for key in matching_keys) or audit["blocks"] == scan["blocks"]:
        raise ArithmeticError("independent-grid p7 z1 audit changed")

    actual_before = int(previous["actual_boundary_count_after_z0_exclusion"])
    actual_after = actual_before - int(scan["z1_boundaries"])
    projected_before = previous["remaining_projected_undetermined_direction_histogram"]
    if projected_before != {1: 300, 2: 280, 3: 210, 7: 2}:
        raise ArithmeticError("pre-15.715 projected profile envelope changed")
    projected_after = {key: value for key, value in projected_before.items() if key != 1}
    projected_count_after = sum(projected_after.values())
    if actual_after != 129_024 or projected_count_after != 492:
        raise ArithmeticError("post-15.715 positive branch count changed")

    return {
        "proposition": "15.715",
        "p": 7,
        "boundary": "infinity plus seven finite points",
        "product_sign": "positive",
        "undetermined_direction_count_excluded": 1,
        "actual_boundary_count_before": actual_before,
        "z1_boundaries_excluded": int(scan["z1_boundaries"]),
        "mean_allocations_per_boundary": int(scan["mean_allocation_count_per_boundary"]),
        "actual_boundary_count_after_z1_exclusion": actual_after,
        "projected_b_profile_count_before": int(previous["projected_b_profile_count_after"]),
        "projected_b_profiles_excluded_here": int(projected_before[1]),
        "projected_b_profile_count_after": projected_count_after,
        "remaining_projected_undetermined_direction_histogram": projected_after,
        "projected_mod7_boundary_candidates": int(scan["projected_survivors"]),
        "full_mod7_survivors": int(scan["all_dependency_survivors"]),
        "modulus": 7,
        "projected_dependency_count": int(scan["projected_dependency_count"]),
        "left_dependency_count": int(linear["left_dependency_dimension"]),
        "scan_evidence": str(SCAN.relative_to(ROOT)),
        "independent_grid_audit_evidence": str(AUDIT.relative_to(ROOT)),
        "scan_script": "scripts/p7_infinity7_positive_z1_mod7_gpu.py",
        "positive_z1_branch_closed": True,
        "positive_p7_infinity_plus_seven_closed": False,
        "negative_p7_infinity_plus_seven_changed": False,
        "top_level_gates_changed": False,
        "proved_by_complete_exact_finite_scan": True,
    }


def main() -> None:
    theorem = p7_positive_infinity_plus_seven_z1_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15715.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.715: positive p7 infinity+7 z=1 boundaries "
        f"{theorem['z1_boundaries_excluded']} -> 0"
    )


if __name__ == "__main__":
    main()
