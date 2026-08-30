#!/usr/bin/env python3
"""Prop. 15.714 -- close positive p7 infinity+7 with z=0.

For z=0 every direction has phase-zero mean eight and a unique complete
Johnson-slice slack catalog. A complete CUDA scan of all C(49,7) finite
boundaries tests the 135 exact mod-seven left dependencies of the common
edge system. All 79,447,032 z=0 boundaries fail; none reaches an edge solve.
"""
from __future__ import annotations

import functools
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15713 import p7_positive_infinity_plus_seven_direction_reduction


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))
SCAN = ROOT / "evidence" / "p7_infinity7_positive_z0_v100.json"
AUDIT = ROOT / "evidence" / "p7_infinity7_positive_z0_v100_audit.json"

SCAN_SHA256 = "ebedf805251bc418f0034e316d37b7d3101a94922d6d4532c06f867cdfa42a76"
AUDIT_SHA256 = "9653275f4f66358fb7a015ffe1d8b22e8fc794a7057b9cbaafd6f1376bfb9ebe"

EXPECTED_INPUT_SHA256 = {
    "equation_matrix": "32b378e8bd6c55deb9b6b546c73ee869b69a5e0d7037f0ed474be5ae882fbc1a",
    "left_dependencies": "0405fad25d2295ed722bd8ee15ebd6592907d8abdf9d97ddfa866500816dbad2",
    "base": "1bd5ca1015bbb735222fe100b1e8b41bb54b13df2aa8667480b11108ad4cea05",
    "floor_tables": "2b66d3b0184b8bb3cc85452ebffa91a819aa26dcf83eac59cae2091c08277656",
    "finite_labels": "96c23cbbee81215029c045ddc326ef7950db23fdc61036d4714f798d9db8e895",
}

EXPECTED_SOURCE_SHA256 = {
    "scripts/p7_infinity7_positive_z0_mod7_gpu.py": (
        "f748005d9be3286094cf7693941544d762f8bef3ae30a27c19b32fb781e7d951"
    ),
    "scripts/p7_size6_positive_infinity_mod7_gpu.py": (
        "d0e54d2749a1fcd2841674134301fef241acbed2abdc5777aedfc5c36e87330a"
    ),
    "scripts/p7_unsaturated_modular_catalog_filter.py": (
        "f9b2781984ab3e2336977d43b657fe337bb09b37baea600b6fdb1f94483d135a"
    ),
    "src/e1_gmin_m4_prop15632.py": (
        "eda17b867dfd9654eb69d5b9f8b6dbb1b9791dd4ef79acc3060ddbf980822e0c"
    ),
}

EXPECTED_SCAN_PREFIX = {"checked": 100_000, "z0": 88_715, "survivor_ranks": []}
EXPECTED_AUDIT_PREFIX = {"checked": 200_000, "z0": 178_533, "survivor_ranks": []}
AUDITED_HIGH_Z_COUNTS = {2: 123_480, 3: 5_488, 7: 56}


def _array_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def _read_pinned_json(path: Path, expected_sha256: str) -> dict[str, object]:
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise ArithmeticError(f"p7 positive z0 evidence bytes changed: {path.name}")
    return json.loads(raw)


@functools.cache
def _current_z0_input_integrity() -> dict[str, object]:
    """Rebuild and pin every cheap mathematical array consumed by the scan."""
    from p7_size6_positive_infinity_mod7_gpu import dependency_tables, finite_labels
    from p7_unsaturated_modular_catalog_filter import equation_matrix, left_dependencies

    matrix = equation_matrix()
    rank, dependencies = left_dependencies(matrix, 7)
    base, floor_tables, linear = dependency_tables()
    labels = finite_labels()
    actual = {
        "equation_matrix": _array_sha256(matrix),
        "left_dependencies": _array_sha256(dependencies),
        "base": _array_sha256(base),
        "floor_tables": _array_sha256(floor_tables),
        "finite_labels": _array_sha256(labels),
    }
    if actual != EXPECTED_INPUT_SHA256:
        raise ArithmeticError("p7 positive z0 mathematical input arrays changed")
    if rank != 147 or linear != {
        "equations": 282,
        "edge_variables": 1_225,
        "rank_mod_7": 147,
        "left_dependency_dimension": 135,
        "left_null_audit": True,
        "valid_odd_fibre_masks": 63,
    }:
        raise ArithmeticError("p7 positive z0 linear-system reconstruction changed")

    source_sha256 = {
        relative: hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        for relative in EXPECTED_SOURCE_SHA256
    }
    if source_sha256 != EXPECTED_SOURCE_SHA256:
        raise ArithmeticError("p7 positive z0 scan source changed")
    return {
        "array_sha256": actual,
        "source_sha256": source_sha256,
        "recomputed_from_current_checkout": True,
    }


def _z01_count_identity() -> dict[str, object]:
    """Algebraically audit N0/N1 from the z>=2 census and direction incidences."""
    all_boundaries = math.comb(49, 7)
    direction_incidences = 8 * 7**7
    z1 = direction_incidences - sum(
        z * count for z, count in AUDITED_HIGH_Z_COUNTS.items()
    )
    z0 = all_boundaries - z1 - sum(AUDITED_HIGH_Z_COUNTS.values())
    if z0 != 79_447_032 or z1 != 6_324_528:
        raise ArithmeticError("p7 z0/z1 direction-incidence identity changed")
    return {
        "all_boundaries": all_boundaries,
        "undetermined_direction_incidences": direction_incidences,
        "audited_high_z_counts": AUDITED_HIGH_Z_COUNTS,
        "derived_z1_boundaries": z1,
        "derived_z0_boundaries": z0,
        "identity": "sum_z z*N_z = 8*7^7",
    }


def p7_positive_infinity_plus_seven_z0_exclusion() -> dict[str, object]:
    previous = p7_positive_infinity_plus_seven_direction_reduction()
    scan = _read_pinned_json(SCAN, SCAN_SHA256)
    audit = _read_pinned_json(AUDIT, AUDIT_SHA256)
    input_integrity = _current_z0_input_integrity()
    count_identity = _z01_count_identity()
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
        or scan["cpu_prefix_verification"] != EXPECTED_SCAN_PREFIX
        or int(scan["blocks"]) != 65_535
        or scan["base_sha256"] != EXPECTED_INPUT_SHA256["base"]
        or scan["tables_sha256"] != EXPECTED_INPUT_SHA256["floor_tables"]
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
        "linear_system",
    )
    if (
        any(audit[key] != scan[key] for key in matching_keys)
        or int(audit["blocks"]) != 32_768
        or audit["cpu_prefix_verification"] != EXPECTED_AUDIT_PREFIX
    ):
        raise ArithmeticError("different-grid p7 z0 rerun changed")
    if (
        int(scan["z0_boundaries"]) != count_identity["derived_z0_boundaries"]
        or int(scan["undetermined_direction_histogram"]["1"])
        != count_identity["derived_z1_boundaries"]
    ):
        raise ArithmeticError("p7 z0/z1 scan counts fail direction-incidence audit")

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
        "different_grid_rerun_evidence": str(AUDIT.relative_to(ROOT)),
        "evidence_sha256": {
            "primary_scan": SCAN_SHA256,
            "different_grid_rerun": AUDIT_SHA256,
        },
        "scan_script": "scripts/p7_infinity7_positive_z0_mod7_gpu.py",
        "input_integrity": input_integrity,
        "count_identity_audit": count_identity,
        "different_grid_rerun_validated": True,
        "independent_implementation_validated": False,
        "validation_scope": (
            "complete exact CUDA scan plus a same-implementation different-grid rerun; "
            "no independent implementation has validated the zero-survivor claim"
        ),
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
