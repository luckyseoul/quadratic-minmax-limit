#!/usr/bin/env python3
"""Prop. 15.718 -- partial global reduction of positive p7 infinity+7 z=7.

The exact parent affine-hull sieve rejects 3,024 of 4,320 pointed branch
cases and leaves 1,296.  Exact affine symmetry partitions those survivors
into 324 four-case classes.  A same-catalog-row mod-3/mod-7 global join
rigorously rejects 87 representatives; 159 representatives survive only a
necessary relaxation and 78 are skipped at the certified memory cap.

Independently, the Johnson incidence semigroup has a complete 896-row
binary Hilbert basis: 56 generators in grade one, 168 new generators in
grade two, and 672 new generators in grade three.  Its exact layer census
through grade eight closes the high-catalog *structure*, not case
feasibility.  Consequently no actual line boundary is subtracted here:
positive z=7 and the quadratic-minmax-limit theorem remain open.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from e1_gmin_m4_prop15717 import p7_positive_infinity_plus_seven_z3_exclusion


ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "evidence" / "p7_infinity7_positive_z7_global_semigroup_summary.json"
SUMMARY_SHA256 = "8c7a9e388fb4826f4f0872d3718ad6ed1a35018a9ccb32b07d3bfcc01add82d9"


def _integer_keyed(row: dict[str, object]) -> dict[int, int]:
    return {int(key): int(value) for key, value in row.items()}


def p7_positive_infinity_plus_seven_z7_global_semigroup_reduction() -> dict[str, object]:
    """Validate and expose the exact partial z=7 sieve/semigroup certificate."""
    previous = p7_positive_infinity_plus_seven_z3_exclusion()
    raw_summary = SUMMARY.read_bytes()
    if hashlib.sha256(raw_summary).hexdigest() != SUMMARY_SHA256:
        raise ArithmeticError("p7 positive z7 global-semigroup summary changed")
    summary = json.loads(raw_summary)

    if (
        summary["experiment"]
        != "p7_infinity7_positive_z7_global_semigroup_summary"
        or summary["status"]
        != "complete_rigorous_partial_z7_reduction_and_semigroup_certificate"
        or int(summary["p"]) != 7
        or int(summary["z"]) != 7
        or int(summary["phase"]) != 0
        or int(summary["c_H"]) != 1
    ):
        raise ArithmeticError("p7 positive z7 summary identity changed")

    parent = summary["parent_affine_hull_sieve"]
    if (
        int(parent["processed_pointed_branch_cases"]) != 4_320
        or int(parent["rejected_pointed_branch_cases"]) != 3_024
        or int(parent["surviving_pointed_branch_cases"]) != 1_296
        or int(parent["rejected_pointed_branch_cases"])
        + int(parent["surviving_pointed_branch_cases"])
        != int(parent["processed_pointed_branch_cases"])
        or parent["artifact_sha256"]
        != "fa1f1e52d24389a9863274cfa6d2b251a4b06e5e9f2b05e624a0b2587ec65f79"
        or parent["all_case_results_sha256"]
        != "3f65e57cf2f09bc4c674711e3bda3a46503c49169ca54aee4728de5619976aaf"
        or parent["survivor_case_keys_sha256"]
        != "f756e0128e12c78ad7a17f85dd621e4e9ff00f0be80c06ac60a71f74045fc784"
        or parent["script"]
        != "scripts/p7_infinity7_positive_z7_pointed_affine_hull_multimod.py"
        or parent["script_sha256"]
        != "0ae46bf3a0ad64975b7e0ac55aea562c9efcdee43828cb205e4e7cabd196d8b8"
    ):
        raise ArithmeticError("p7 positive z7 parent affine-hull sieve changed")

    symmetry = summary["four_case_symmetry"]
    if (
        symmetry["affine_maps"] != ["u -> 8u", "u -> 32u + 24"]
        or symmetry["direction_permutation"] != [4, 7, 1, 0, 2, 3, 5, 6]
        or symmetry["all_four_augmented_compact_row_spaces_identical"] is not True
        or symmetry["survivor_classes_partition_all_1296_cases"] is not True
        or symmetry["transfer_preserves_exact_same_row_mod3_mod7_catalog_join"]
        is not True
        or int(symmetry["complete_four_case_classes"]) != 324
        or int(symmetry["partial_four_case_classes"]) != 0
        or 4 * int(symmetry["complete_four_case_classes"])
        != int(parent["surviving_pointed_branch_cases"])
        or _integer_keyed(symmetry["row_space_ranks_by_modulus"])
        != {3: 120, 5: 114, 7: 135, 11: 114}
        or symmetry["artifact_sha256"]
        != "0df450d3be1b6897d2bf5e7ae55473ad8c53facffdab86256745d37d5e49cdf3"
        or symmetry["leaf_permutation_sha256_int64"]
        != "cc5af4ca13a2d16713a2b84e98f693c9825e20d7dba5844c87c05207e0642bae"
        or symmetry["script_sha256"]
        != "42ce41934612e0cde03338744975f470ec1c507816b1ee29dea000d19d73950c"
    ):
        raise ArithmeticError("p7 positive z7 four-case symmetry certificate changed")
    expected_rref_hashes = {
        "3": "4896adb1116b115941a784a3ae3bf5310f2dc73125a70b4a65f5e9e422197ba4",
        "5": "f4ec7f92643c521afe188f77d247b9ae2970d7d398e6e07ae1192a46b3c58810",
        "7": "10c0a4c13efbd4573bd75c2f24225a907eddeccbdbed37aead5ebbd3c174d609",
        "11": "5c9ed22d7a49521cb7323e0b707ce70dc319d4fceb49b745564b8f5b82d853f8",
    }
    if symmetry["canonical_rref_sha256_by_modulus"] != expected_rref_hashes:
        raise ArithmeticError("p7 positive z7 canonical row spaces changed")

    join = summary["global_catalog_join"]
    representative_counts = {
        "selected": int(join["selected_representatives"]),
        "processed": int(join["exact_global_join_processed_representatives"]),
        "rejected": int(join["rigorously_rejected_representatives"]),
        "necessary_survivors": int(join["necessary_only_survivors_representatives"]),
        "budget_skips": int(join["budget_skips_representatives"]),
    }
    if (
        representative_counts
        != {
            "selected": 324,
            "processed": 246,
            "rejected": 87,
            "necessary_survivors": 159,
            "budget_skips": 78,
        }
        or representative_counts["processed"]
        != representative_counts["rejected"]
        + representative_counts["necessary_survivors"]
        or representative_counts["selected"]
        != representative_counts["processed"] + representative_counts["budget_skips"]
        or int(join["execution_side_state_cap"]) != 6_000_000
        or join["same_catalog_row_index_used_mod3_mod7"] is not True
        or join["zero_join_is_a_rigorous_rejection"] is not True
        or join["signature_hash_collision_assumption_used"] is not False
        or join["high_catalogs_relaxed_to_exact_affine_hulls"] is not True
        or join["z7_branch_excluded"] is not False
        or join["artifact_sha256"]
        != "5d3bb1a3385fe848932f3405032eb45501741afc8a3aef36d76b126f5859b93a"
        or join["case_results_sha256"]
        != "c34cc913c27910e3876e1b78aed0e9c8c2f42cb2f4368f95054bcd6ead1db7a7"
        or join["rejection_certificates_sha256"]
        != "03ef66e6b4529d05c0351762dd8a190eb9b400e3b050852fafa8f2ea16d0da78"
        or join["survivor_witnesses_sha256"]
        != "b52fbafed16a04514ce7d403d562e39d12fde0ac6ba31845caf3e20b8c8a8dd0"
        or join["skipped_cases_sha256"]
        != "8adc451b2befc8be933f6a40f363b55661b26ead0932a45b65be7d989961b784"
        or join["transferred_case_decisions_sha256"]
        != "7b045433541f00e6b866e443ccd23187906550bfe77c1985ee512ad04a41f210"
        or join["script_sha256"]
        != "86ca9a8055ba20f129f284b6e9001478880a56f20185b701c708000577f34bd8"
    ):
        raise ArithmeticError("p7 positive z7 exact global catalog join changed")

    transferred = {key: int(value) for key, value in join["transferred_counts"].items()}
    if (
        transferred != {"processed": 984, "rejected": 348, "skipped": 312, "surviving": 636}
        or transferred["rejected"] != 4 * representative_counts["rejected"]
        or transferred["surviving"] != 4 * representative_counts["necessary_survivors"]
        or transferred["skipped"] != 4 * representative_counts["budget_skips"]
        or transferred["processed"]
        != transferred["rejected"] + transferred["surviving"]
        or sum(transferred[key] for key in ("rejected", "surviving", "skipped"))
        != int(parent["surviving_pointed_branch_cases"])
    ):
        raise ArithmeticError("p7 positive z7 transferred join accounting changed")

    hostile = summary["hostile_independent_audit"]
    if hostile != {
        "all_324_balanced_partitions_reoptimized": True,
        "exact_sign_convention_rechecked": True,
        "false_negative_path_found": False,
        "randomized_bruteforce_join_checks": 250,
        "verdict": "PASS",
    }:
        raise ArithmeticError("p7 positive z7 hostile audit changed")

    no_go = summary["mod5_mod11_no_go"]
    if (
        int(no_go["candidate_cases_preflighted"]) != 159
        or no_go["all_159_candidate_projected_bases_zero_mod5_mod11"] is not True
        or no_go["all_complete_catalog_contributions_zero_mod5_mod11"] is not True
        or no_go["four_prime_dedup_equals_mod3_mod7_dedup"] is not True
        or int(no_go["new_active_candidate_count"]) != 0
        or int(no_go["new_rigorous_rejections"]) != 0
        or int(no_go["zero_mean_exact_liftable_hull_rank_all_primes"]) != 20
        or no_go["artifact_sha256"]
        != "e59f5609a74ee84ed2e020003ee912ad72b930c36c8d63d98621eeb29e965bcd"
        or no_go["script_sha256"]
        != "a3bbb8b7dcc47d5acf89da888e214a468053b3f77ae61e68f66b88ff9cbabd3d"
    ):
        raise ArithmeticError("p7 positive z7 mod-5/mod-11 no-go audit changed")

    semigroup = summary["johnson_semigroup"]
    generator_histogram = _integer_keyed(semigroup["generator_grade_histogram"])
    layer_counts = _integer_keyed(
        semigroup["complete_semigroup_layer_counts_through_grade_eight"]
    )
    expected_layers = {
        0: 1,
        1: 56,
        2: 1_764,
        3: 37_856,
        4: 575_407,
        5: 6_496_938,
        6: 57_232_105,
        7: 410_200_367,
        8: 2_474_264_653,
    }
    if (
        int(semigroup["hilbert_basis_rows"]) != 896
        or generator_histogram != {1: 56, 2: 168, 3: 672}
        or sum(generator_histogram.values()) != int(semigroup["hilbert_basis_rows"])
        or int(semigroup["maximum_primitive_generator_grade"]) != 3
        or int(semigroup["maximum_generator_coordinate"]) != 1
        or semigroup["S56_equals_complete_grade_one_layer"] is not True
        or int(semigroup["complete_degree_three_rows"]) != 37_856
        or int(semigroup["new_primitive_grade_three_rows"]) != 672
        or semigroup["M1764_decomposition"]
        != {
            "decomposable_S_plus_S_rows": 1_596,
            "primitive_grade_two_rows": 168,
            "union_equals_complete_M1764": True,
        }
        or layer_counts != expected_layers
        or semigroup["required_high_grades"] != [3, 4, 5, 6, 8]
        or semigroup["coordinate_cap_automatic_through_grade_six"] is not True
        or semigroup["grade_eight_requires_explicit_coordinate_cap"] is not True
        or semigroup["input_sha256"]
        != "14b670086a068da3f07c3bc1fb6c8ece339ae9dd78a6b1e412219a286cf45a5f"
        or semigroup["generator_file_sha256"]
        != "3b582d6a0e7c83cb8ed41a421e4950be2645ce2d3aa18dab17432628b787b789"
        or semigroup["hilbert_series_output_sha256"]
        != "92e8a4f648b9cc13395e392164dd0d7e7632e0027cd2954cf7915622a3e31bdf"
        or semigroup["audit_artifact_sha256"]
        != "09c9ffaeace3eab5fa3bcca99a78d2d63496c6a23bf7f1b5f162769395798a24"
        or semigroup["script_sha256"]
        != "465ecf21bfd586379959a292493e717453a81c672e8f40d61e381ab8fb4ed0c5"
    ):
        raise ArithmeticError("p7 positive z7 Johnson semigroup certificate changed")

    attack_order = {
        key: int(value)
        for key, value in summary["remaining_representative_attack_order"].items()
    }
    expected_attack_order = {
        "cap_sensitive_grade_eight": 8,
        "grade_three_only": 51,
        "H0_S0_M7_calibration": 4,
        "maximum_grade_five": 24,
        "maximum_grade_four": 137,
        "maximum_grade_six": 13,
        "total": 237,
    }
    if (
        attack_order != expected_attack_order
        or attack_order["total"]
        != representative_counts["necessary_survivors"]
        + representative_counts["budget_skips"]
        or sum(value for key, value in attack_order.items() if key != "total")
        != attack_order["total"]
    ):
        raise ArithmeticError("p7 positive z7 remaining high-grade census changed")

    semantics = summary["logical_semantics"]
    actual = summary["actual_boundary_status"]
    if semantics != {
        "global_join_survivor_is_binary_edge_feasibility": False,
        "global_join_survivor_is_exact_high_catalog_feasibility": False,
        "no_actual_line_boundary_count_subtracted_here": True,
        "positive_z7_closed": False,
        "semigroup_generation_closes_high_catalog_structure_not_case_feasibility": True,
    }:
        raise ArithmeticError("p7 positive z7 logical semantics changed")
    if actual != {
        "positive_p7_infinity_plus_seven_closed": False,
        "remaining_actual_line_boundaries": 56,
        "remaining_actual_line_orbits": 2,
        "remaining_projected_profiles": 2,
        "top_level_gates_changed": False,
    }:
        raise ArithmeticError("p7 positive z7 actual-boundary status changed")
    if (
        int(previous["actual_boundary_count_after_z3_exclusion"]) != 56
        or previous["remaining_actual_undetermined_direction_histogram"] != {7: 56}
        or previous["remaining_actual_boundary_orbits"] != {7: 2}
        or int(previous["projected_b_profile_count_after"]) != 2
        or previous["remaining_projected_undetermined_direction_histogram"] != {7: 2}
        or previous["positive_p7_infinity_plus_seven_closed"] is not False
    ):
        raise ArithmeticError("pre-15.718 positive z7 remainder changed")

    unresolved_representatives = (
        representative_counts["necessary_survivors"] + representative_counts["budget_skips"]
    )
    return {
        "proposition": "15.718",
        "p": 7,
        "boundary": "infinity plus seven finite points",
        "product_sign": "positive",
        "undetermined_direction_count": 7,
        "summary_evidence": str(SUMMARY.relative_to(ROOT)),
        "summary_evidence_sha256": SUMMARY_SHA256,
        "pointed_branch_cases_before_affine_sieve": 4_320,
        "affine_hull_rigorously_rejected_pointed_branch_cases": 3_024,
        "pointed_branch_cases_after_affine_sieve": 1_296,
        "four_case_symmetry_class_size": 4,
        "four_case_symmetry_representatives": 324,
        "global_join_processed_representatives": 246,
        "global_join_rigorously_rejected_representatives": 87,
        "global_join_necessary_only_survivor_representatives": 159,
        "global_join_budget_skip_representatives": 78,
        "global_join_unresolved_representatives": unresolved_representatives,
        "transferred_pointed_case_counts": transferred,
        "same_catalog_row_mod3_mod7_join": True,
        "mod5_mod11_additional_rejections": 0,
        "johnson_semigroup_hilbert_basis_rows": 896,
        "johnson_semigroup_generator_grade_histogram": generator_histogram,
        "complete_semigroup_layer_counts_through_grade_eight": layer_counts,
        "required_high_grades": list(semigroup["required_high_grades"]),
        "remaining_representative_high_grade_census": attack_order,
        "coordinate_cap_automatic_through_grade_six": True,
        "grade_eight_requires_explicit_coordinate_cap": True,
        "actual_line_boundary_count_before": 56,
        "actual_line_boundaries_excluded_here": 0,
        "actual_line_boundary_count_after": 56,
        "remaining_actual_undetermined_direction_histogram": {7: 56},
        "remaining_actual_boundary_orbits": {7: 2},
        "remaining_projected_b_profile_count": 2,
        "affine_and_global_zero_join_rejections_are_rigorous": True,
        "global_join_survivor_is_feasibility_certificate": False,
        "semigroup_certificate_closes_high_catalog_structure_only": True,
        "positive_z7_branch_closed": False,
        "positive_p7_infinity_plus_seven_closed": False,
        "quadratic_minmax_limit_theorem_closed": False,
        "theorem_remains_open": True,
        "negative_p7_infinity_plus_seven_changed": False,
        "top_level_gates_changed": False,
        "proved_by_exact_affine_sieve_symmetry_global_join_and_semigroup_census": True,
    }


def main() -> None:
    theorem = p7_positive_infinity_plus_seven_z7_global_semigroup_reduction()
    print(
        "Prop. 15.718: positive p7 infinity+7 z=7 pointed cases "
        f"{theorem['pointed_branch_cases_before_affine_sieve']} -> "
        f"{theorem['pointed_branch_cases_after_affine_sieve']}; "
        "z=7 remains open"
    )


if __name__ == "__main__":
    main()
