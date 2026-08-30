#!/usr/bin/env python3
"""Prop. 15.717 -- close positive p7 infinity+7 with z=3.

The ten exact affine-semilinear boundary orbits have 400 corrected mean
leaves.  A complete mod-seven catalog exhaustion rejects 398 leaves and
leaves only eight catalog-row tuples across two leaves.  Reusing each same
integer right side modulo three rejects all eight remaining tuples.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15716 import p7_positive_infinity_plus_seven_z2_exclusion


ROOT = Path(__file__).resolve().parents[1]
MOD7 = ROOT / "evidence" / "p7_infinity7_positive_z3_mod7_join_nuka.json"
MULTIMOD = ROOT / "evidence" / "p7_infinity7_positive_z3_multimod_nuka.json"


def p7_positive_infinity_plus_seven_z3_exclusion() -> dict[str, object]:
    previous = p7_positive_infinity_plus_seven_z2_exclusion()
    mod7 = json.loads(MOD7.read_text())
    multimod = json.loads(MULTIMOD.read_text())

    if (
        mod7["experiment"] != "p7_infinity7_positive_z3_mod7_join"
        or mod7["status"] != "complete_rigorous_mod7_necessary_sieve_with_survivors"
        or mod7["full_run"] is not True
        or int(mod7["processed_orbits"]) != 10
        or int(mod7["processed_exact_mean_leaves"]) != 400
        or int(mod7["processed_weighted_boundary_allocation_cases"]) != 225_792
        or int(mod7["rejected_weighted_boundary_allocation_cases"]) != 225_008
        or int(mod7["surviving_cases"]) != 2
        or int(mod7["surviving_weighted_boundary_allocation_cases"]) != 784
        or int(mod7["weighted_exact_mod7_catalog_tuples"]) != 3_136
        or mod7["z3_branch_excluded"] is not False
        or mod7["affine_span_relaxation_used"] is not False
        or mod7["all_case_decisions_sha256"]
        != "c39fb7f530a6380c09d0bf300d6d249df2304370867066e6ce74de62906f275f"
        or mod7["mean_leaf_coverage"]["exact_mean_leaves"] != 400
        or mod7["mean_leaf_coverage"]["residue_histogram"]
        != {"00": 360, "04": 20, "40": 20}
        or mod7["orbit_reduction"]["orbit_count"] != 10
        or mod7["orbit_reduction"]["orbit_size_sum"] != 5_488
    ):
        raise ArithmeticError("p7 positive z3 mod-seven exhaustion changed")

    if (
        multimod["experiment"] != "p7_infinity7_positive_z3_multimod_join"
        or multimod["status"] != "complete_rigorous_multimod_exclusion"
        or int(multimod["input_surviving_mean_cases"]) != 2
        or int(multimod["extracted_mod7_catalog_tuples"]) != 8
        or int(multimod["expected_mod7_catalog_tuples"]) != 8
        or int(multimod["additional_prime_surviving_tuples"]) != 0
        or multimod["all_mod7_tuples_fail_at_least_one_additional_prime"] is not True
        or multimod["same_integer_rhs_tested_at_every_modulus"] is not True
        or multimod["z3_branch_excluded"] is not True
        or multimod["mod7_row_index_certificate_sha256"]
        != "ababbe1e75ad4c913d3262c6df58f5ad44b8755889f5d3614978463ecb826a98"
        or multimod["multimod_decision_certificate_sha256"]
        != "c7e06a92a4581309b1b3e59a2e1e84f1d48e5988c53337864e73521f29debc82"
        or multimod["audited_z3_module"]["sha256"]
        != "4d4e3e52b3baa1838bd7a2d2b0e55fc05e65b4296a86e48374ddfa9be47aba25"
    ):
        raise ArithmeticError("p7 positive z3 same-tuple multimod exhaustion changed")

    linear = {
        int(row["modulus"]): (int(row["rank"]), int(row["dependency_dimension"]))
        for row in multimod["linear_system"]["complete_left_dependency_bases"]
    }
    if linear != {3: (161, 120), 5: (167, 114), 7: (146, 135), 11: (167, 114)}:
        raise ArithmeticError("p7 positive z3 multimodular ranks changed")

    eliminators = [
        tuple(int(value) for value in row["eliminating_additional_primes"])
        for case in multimod["cases"]
        for row in case["tuple_decisions"]
    ]
    if eliminators != [(3,)] * 8:
        raise ArithmeticError("the eight z3 tuple eliminators changed")

    actual_before = int(previous["actual_boundary_count_after_z2_exclusion"])
    actual_after = actual_before - 5_488
    projected_before = previous["remaining_projected_undetermined_direction_histogram"]
    if projected_before != {3: 210, 7: 2}:
        raise ArithmeticError("pre-15.717 projected profile envelope changed")
    projected_after = {7: int(projected_before[7])}
    if actual_after != 56 or sum(projected_after.values()) != 2:
        raise ArithmeticError("post-15.717 positive branch count changed")

    return {
        "proposition": "15.717",
        "p": 7,
        "boundary": "infinity plus seven finite points",
        "product_sign": "positive",
        "undetermined_direction_count_excluded": 3,
        "actual_boundary_count_before": actual_before,
        "z3_boundaries_excluded": 5_488,
        "z3_boundary_orbits": 10,
        "exact_mean_leaves": 400,
        "mean_leaves_rejected_mod7": 398,
        "mod7_surviving_mean_leaves": 2,
        "extracted_exact_mod7_catalog_tuples": 8,
        "same_tuple_mod3_survivors": 0,
        "actual_boundary_count_after_z3_exclusion": actual_after,
        "remaining_actual_undetermined_direction_histogram": {7: 56},
        "remaining_actual_boundary_orbits": {7: 2},
        "projected_b_profile_count_before": int(previous["projected_b_profile_count_after"]),
        "projected_b_profiles_excluded_here": int(projected_before[3]),
        "projected_b_profile_count_after": sum(projected_after.values()),
        "remaining_projected_undetermined_direction_histogram": projected_after,
        "tested_moduli": list(multimod["tested_moduli"]),
        "mod7_evidence": str(MOD7.relative_to(ROOT)),
        "multimod_evidence": str(MULTIMOD.relative_to(ROOT)),
        "mod7_script": "scripts/p7_infinity7_positive_z3_mod7_join.py",
        "multimod_script": "scripts/p7_infinity7_positive_z3_multimod_join.py",
        "positive_z3_branch_closed": True,
        "positive_p7_infinity_plus_seven_closed": False,
        "negative_p7_infinity_plus_seven_changed": False,
        "top_level_gates_changed": False,
        "proved_by_complete_exact_orbit_catalog_and_same_tuple_multimod_exhaustion": True,
    }


def main() -> None:
    theorem = p7_positive_infinity_plus_seven_z3_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15717.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.717: positive p7 infinity+7 z=3 boundaries "
        f"{theorem['z3_boundaries_excluded']} -> 0"
    )


if __name__ == "__main__":
    main()
