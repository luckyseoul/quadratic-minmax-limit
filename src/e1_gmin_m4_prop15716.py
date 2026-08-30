#!/usr/bin/env python3
"""Prop. 15.716 -- close positive p7 infinity+7 with z=2.

Pair transversals reduce 123,480 actual boundaries to 92 exact affine
square-semilinear orbits.  The translation-equivariant 281-row edge system
then rejects all 1,232 exact mean leaves modulo seven, including the
residue-four family.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15715 import p7_positive_infinity_plus_seven_z1_exclusion


ROOT = Path(__file__).resolve().parents[1]
ORBITS = ROOT / "evidence" / "p7_infinity7_positive_zge2_orbits.json"
JOIN = ROOT / "evidence" / "p7_infinity7_positive_z2_mod7_join_nuka.json"


def p7_positive_infinity_plus_seven_z2_exclusion() -> dict[str, object]:
    previous = p7_positive_infinity_plus_seven_z1_exclusion()
    orbits = json.loads(ORBITS.read_text())
    join = json.loads(JOIN.read_text())
    census = orbits["census"]
    linear = join["linear_system"]
    leaf = join["mean_leaf_coverage"]

    if (
        orbits["experiment"] != "p7_infinity7_positive_zge2_orbits"
        or orbits["status"] != "complete_exact_pair_transversal_orbit_census"
        or orbits["all_required_audits_passed"] is not True
        or int(orbits["group_audit"]["group_size"]) != 2_352
        or census["boundary_count_by_z"] != {"2": 123_480, "3": 5_488, "7": 56}
        or census["orbit_count_by_z"] != {"2": 92, "3": 10, "7": 2}
        or int(census["pair_transversal_incidences"]) != 141_120
        or join["experiment"] != "p7_infinity7_positive_z2_mod7_join"
        or join["status"] != "complete_rigorous_mod7_exclusion"
        or join["full_run"] is not True
        or int(join["processed_orbits"]) != 92
        or int(join["processed_exact_mean_leaves"]) != 1_232
        or join["processed_kind_histogram"]
        != {"one_high_q2": 192, "residue4_four_catalog": 48, "two_q1": 992}
        or join["rejected_kind_histogram"] != join["processed_kind_histogram"]
        or int(join["surviving_cases"]) != 0
        or join["z2_branch_excluded"] is not True
        or join["affine_span_relaxation_used"] is not False
        or join["q2_full_direction_block_image_relaxation_used"] is not True
        or join["all_case_decisions_sha256"]
        != "1c39bc61a34cb1f96e7fcafb90a21ccc0e40ad39e27e0110bd87a50224b1f5c6"
        or int(linear["equations"]) != 281
        or int(linear["edge_variables"]) != 1_225
        or int(linear["rank"]) != 146
        or int(linear["left_dependency_dimension"]) != 135
        or int(linear["direction_block_offset"]) != 1
        or int(linear["edge_count_rhs"]) != 29
        or linear["left_null_audit"] is not True
        or int(leaf["exact_mean_leaves"]) != 1_232
        or leaf["obsolete_1184_count_rejected"] is not True
    ):
        raise ArithmeticError("p7 positive z2 orbit/catalog exhaustion changed")

    join_generation = join["boundary_generation"]
    join_orbits = join["orbit_reduction"]
    if (
        join_generation["boundary_histogram_by_z"]
        != {"2": 123_480, "3": 5_488, "7": 56}
        or int(join_generation["pair_transversal_incidences"]) != 141_120
        or int(join_orbits["group_size"]) != 2_352
        or int(join_orbits["orbit_count"]) != 92
        or int(join_orbits["orbit_size_sum"]) != 123_480
        or join_orbits["orbit_size_histogram"]
        != {"588": 18, "1176": 52, "2352": 22}
        or int(join_orbits["same_type_orbits"]) != 48
        or int(join_orbits["split_type_orbits"]) != 44
    ):
        raise ArithmeticError("independent z2 orbit reconstruction changed")

    actual_before = int(previous["actual_boundary_count_after_z1_exclusion"])
    actual_after = actual_before - 123_480
    projected_before = previous["remaining_projected_undetermined_direction_histogram"]
    if projected_before != {2: 280, 3: 210, 7: 2}:
        raise ArithmeticError("pre-15.716 projected profile envelope changed")
    projected_after = {key: value for key, value in projected_before.items() if key != 2}
    if actual_after != 5_544 or sum(projected_after.values()) != 212:
        raise ArithmeticError("post-15.716 positive branch count changed")

    return {
        "proposition": "15.716",
        "p": 7,
        "boundary": "infinity plus seven finite points",
        "product_sign": "positive",
        "undetermined_direction_count_excluded": 2,
        "actual_boundary_count_before": actual_before,
        "z2_boundaries_excluded": 123_480,
        "z2_boundary_orbits": 92,
        "exact_mean_leaves_excluded": int(join["processed_exact_mean_leaves"]),
        "actual_boundary_count_after_z2_exclusion": actual_after,
        "remaining_actual_undetermined_direction_histogram": {3: 5_488, 7: 56},
        "remaining_actual_boundary_orbits": {3: 10, 7: 2},
        "projected_b_profile_count_before": int(previous["projected_b_profile_count_after"]),
        "projected_b_profiles_excluded_here": int(projected_before[2]),
        "projected_b_profile_count_after": sum(projected_after.values()),
        "remaining_projected_undetermined_direction_histogram": projected_after,
        "modulus": 7,
        "translation_equivariant_equations": int(linear["equations"]),
        "left_dependency_count": int(linear["left_dependency_dimension"]),
        "mod7_surviving_mean_leaves": int(join["surviving_cases"]),
        "orbit_evidence": str(ORBITS.relative_to(ROOT)),
        "join_evidence": str(JOIN.relative_to(ROOT)),
        "orbit_script": "scripts/p7_infinity7_positive_zge2_orbits.py",
        "join_script": "scripts/p7_infinity7_positive_z2_mod7_join.py",
        "positive_z2_branch_closed": True,
        "positive_p7_infinity_plus_seven_closed": False,
        "negative_p7_infinity_plus_seven_changed": False,
        "top_level_gates_changed": False,
        "proved_by_complete_exact_orbit_catalog_exhaustion": True,
    }


def main() -> None:
    theorem = p7_positive_infinity_plus_seven_z2_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15716.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.716: positive p7 infinity+7 z=2 boundaries "
        f"{theorem['z2_boundaries_excluded']} -> 0"
    )


if __name__ == "__main__":
    main()
