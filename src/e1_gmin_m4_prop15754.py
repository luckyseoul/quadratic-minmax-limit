#!/usr/bin/env python3
r"""Prop. 15.754 -- the p=13 fifth-shell endpoint is empty.

At ``p=13,t=4,k=60,u=6`` the hard means are ``12+14 k_L`` with
``sum_L k_L=5``.  Exact floor equality is the XNOR row, so at least two
hard directions are exact.  Common signed-edge normalization then gives

    hT=5,  P_L=4+k_L,  Q_L=4.

The seven partitions of five are exhausted by three complementary exact
aggregate mechanisms:

* two or three exact XNOR roots: exhaustive common ``U=h*M2`` and
  ``G=h*M4-M2^2`` separable-energy bounds;
* four exact roots: exhaustive common ``U/G/J6`` coefficient joins, with
  the sign-correct opposite key ``(-N2,-N4-N2^2,-N6+N2^3)``; and
* five or six exact roots: ``G`` vanishes identically, after which all 74
  translated cuts and the global collision budget give strict energy gaps.

This is an exhaustive finite aggregate/common-form certificate and a proved
finite endpoint theorem.  It is not a graph, orbit, coefficient-cell, or
common-realization census, and it does not close residual (ii) globally.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from functools import lru_cache
from math import ceil, floor
from pathlib import Path
from typing import Iterable

from e1_gmin_m4_prop15274 import residual_ii_k_ge_4p_ND_closed
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15734 import (
    BRANCH_B2,
    baseline_coefficient_rules,
    residual_even_floor_table,
)
from e1_gmin_m4_prop15744 import t4_all_residue_sieve
from e1_gmin_m4_prop15749 import proposition_15749
from io_atomic import write_json_atomic
from p13_u6_cut_equalities import (
    translated_cut_vectors,
    u6_energy_ledger_certificate,
)
from p13_u6_high_root_energy import (
    balanced_collision,
    high_root_partition_certificate,
    nonexact_parseval_base,
    signed_collision_floor,
)


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE_PATH = ROOT / "evidence" / "e1_gmin_m4_prop15754.json"
P = 13
M = 7
T = 4
K = 60
U = 6
H_EDGE_COUNT = 61
EXPECTED_REPOSITORY_HEAD = "cf32d2137d9b40ce631f21bcdf8b2cb6e72d0c81"

ARTIFACT_PATHS = {
    "low_root_UG": ROOT / "evidence" / "e1_gmin_m4_prop15754_low_root_ug.json",
    "four_root_221": ROOT / "evidence" / "e1_gmin_m4_prop15754_four_root_221.json",
    "four_root_311": ROOT / "evidence" / "e1_gmin_m4_prop15754_four_root_311.json",
}
EXPECTED_ARTIFACT_SHA256 = {
    "low_root_UG": "9d9a3b75a00410706df94dd02086357393b00bb463b3f9286596a08eb5faa0a4",
    "four_root_221": "83cfbe4f684c4406a2d23216e4049f083c438ed07fda23214bcffbe1c0d2449b",
    "four_root_311": "d924059b403f2ddb33282ff60ef01c32ac7ce3806a276083157eb0b073fec2bc",
}
EXPECTED_SUPPORTING_CODE_SHA256 = {
    "scripts/p13_u6_joint_ug_tables.py": (
        "811eb1833d5551f3fadce62ec6a4302296e301a2cd001bedd399df723afaaebe"
    ),
    "scripts/p13_u6_low_root_ug_bound.py": (
        "213a435e6c8848879e1fb485fa9fdd82832c72b17b863090183c1cd65e62bfaa"
    ),
    "scripts/p13_u6_four_root_ugj.py": (
        "71754b7e9e0f5df38819620f5f90903614152ca4a1516e544fa6fd86f4c0aef1"
    ),
    "src/p13_u6_cut_equalities.py": (
        "376e3a221f35c0e78a131c8fd45a27498a284e40e7c49f027c0b9a45fa649112"
    ),
    "src/p13_u6_high_root_energy.py": (
        "432030353ac92ef4c5000e83ece5a0a3655ef8ae2a29fdd1be5c3eaca5e337d2"
    ),
}
EXPECTED_IMPORTED_HELPER_SHA256 = {
    "src/e1_gmin_m4_prop15740.py": (
        "7a2cfbd12a7057971a0cbaaf523d16f65cff989a47d594a79f641062d02439c3"
    ),
    "scripts/p13_p5_literal_interpolation.py": (
        "31ba186632780d62b3352c214b1112a6beafa48b8a2571882da29e94faf618e5"
    ),
}

EXPECTED_LOW_TABLE_MAXIMUM_SHA256 = {
    "H1": "8267a2d0f1af8a167d27b0a222567c669f7e799f905ff82953faa922dc97b5ff",
    "H2": "14dd0b5429531cfcc52a155a4d7bea5a0d6c653a3e8cf11caa37f59b113cebb3",
    "O": "9c7de941ab58babfa0e78169892fefd47607221208f70e16520a137a3aab925a",
}
EXPECTED_ROW_CATALOGS = {
    "H1": (372, "93ff5fb9059e2cb83242b5716b06d5664e62577baf3d88a5ff4de678d0022a89", 366, "7f758aa8fc9f341e97c1c653cfe5822a9a00f7e073d6716e8c213b67f7e25128", 0, 28),
    "H2": (3078, "a5ddf5eceff6cbeea822a9fb4a1f3a63aa3519a75615f8a50c4fb5fca7270c4d", 2974, "01b06567a085dd53fd05e5953e1621e29f7baab89308dc598c6f7282fcd9033c", 1, 63),
    "H3_D0": (13773, "9f919262fc6ee4d9861e1f9b689b84870b26b77be410d9f5cd7d4fae6776b3be", 12195, "0858ca982e65a5d9ce5e1920c284b80529176c6ff7b93a715bb72c609451b555", 2, 140),
    "H3_D1": (13809, "5a4cb3fb3ee3adfa27856a198e62a61325ad91ce6f121caa998f8ae842b0011b", 12231, "01477cf479321b03abb063cde40bda2b6c0a9b3a7105eac1f782723a59da18e9", 2, 144),
    "H3_D2": (13809, "5a4cb3fb3ee3adfa27856a198e62a61325ad91ce6f121caa998f8ae842b0011b", 12231, "01477cf479321b03abb063cde40bda2b6c0a9b3a7105eac1f782723a59da18e9", 2, 144),
    "H3_D3": (13809, "5a4cb3fb3ee3adfa27856a198e62a61325ad91ce6f121caa998f8ae842b0011b", 12231, "01477cf479321b03abb063cde40bda2b6c0a9b3a7105eac1f782723a59da18e9", 2, 144),
    "O": (242, "7beae84e08778e0231541ba8e815bfb5590086d3b6d5083f07993b1cf60f48ad", 242, "35b0d2d5ddf3b7426c65e79882f81041953a4a035991736db41ec03ecf1dbae0", 15, 31),
}

EXPECTED_PARTITIONS = (
    (1, 1, 1, 1, 1),
    (2, 1, 1, 1),
    (2, 2, 1),
    (3, 1, 1),
    (3, 2),
    (4, 1),
    (5,),
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _integer_partitions(total: int, maximum: int | None = None) -> tuple[tuple[int, ...], ...]:
    if total == 0:
        return ((),)
    upper = min(total, total if maximum is None else maximum)
    rows = []
    for first in range(upper, 0, -1):
        for tail in _integer_partitions(total - first, first):
            rows.append((first,) + tail)
    return tuple(rows)


@lru_cache(maxsize=1)
def p13_u6_normalization_certificate() -> dict[str, object]:
    """Derive the sole residue, parallel normalization, and seven partitions."""
    sieve = t4_all_residue_sieve()
    prior = proposition_15749()
    residue_row = next(row for row in sieve["residue_rows"] if row["u"] == U)
    baseline = baseline_coefficient_rules(P)[BRANCH_B2]
    floors = residual_even_floor_table(P)
    lift = sharp_integral_quadratic_lift_floor(P)

    hard_quotient_sum = M + T - U
    exact_parallel_upper = (H_EDGE_COUNT - hard_quotient_sum) // M
    exact_parallel_candidates = tuple(
        value
        for value in range(exact_parallel_upper + 1)
        if (value - int(baseline["offset"])) % 6 == 0
    )
    exact_parallel = exact_parallel_candidates[0]
    hT = 14 * exact_parallel - 51
    hard_parallel_total = M * exact_parallel + hard_quotient_sum
    opposite_parallel_total = H_EDGE_COUNT - hard_parallel_total

    opposite_mean = lambda q: 14 * q + hT - 39
    q3_mean = opposite_mean(3)
    minimum_opposite_parallel = 4
    opposite_parallel_counts = (minimum_opposite_parallel,) * M

    generated = set(_integer_partitions(hard_quotient_sum))
    partition_rows = []
    for partition in EXPECTED_PARTITIONS:
        exact_count = M - len(partition)
        collision_minimum = sum(
            balanced_collision(exact_parallel + excess, 6)
            for excess in partition
        )
        partition_rows.append(
            {
                "hard_excess_partition": list(partition),
                "hard_quotients": [0] * exact_count + list(partition),
                "hard_parallel_counts": [exact_parallel] * exact_count
                + [exact_parallel + excess for excess in partition],
                "exact_XNOR_root_count": exact_count,
                "collision_minimum": collision_minimum,
                "parseval_nonexact_base": nonexact_parseval_base(partition),
            }
        )

    proved = bool(
        sieve["proved"]
        and sieve["hard_mean_form"] == "a=2u+14k"
        and sieve["hard_quotient_identity"] == "sum k=11-u"
        and residue_row["forced_low_mean"] == 12
        and residue_row["surviving_low_cells"]
        == [{"b": 2, "floor": 12, "excess": 0, "status": "exact baseline", "survives": True}]
        and prior["proved"]
        and prior["remaining_p13_t4_residues"] == [6]
        and hard_quotient_sum == 5
        and baseline["offset"] == 4
        and baseline["congruence"] == "6 divides I+P-4"
        and exact_parallel_upper == 8
        and exact_parallel_candidates == (4,)
        and hT == 5
        and hard_parallel_total == 33
        and opposite_parallel_total == 28
        and q3_mean == 8
        and q3_mean < int(lift["sharp_scaled_floor"]) == 10
        and q3_mean < int(floors["least_nonzero_phase_zero_floor"]) == 12
        and sum(opposite_parallel_counts) == opposite_parallel_total
        and generated == set(EXPECTED_PARTITIONS)
        and [row["exact_XNOR_root_count"] for row in partition_rows]
        == [2, 3, 4, 4, 5, 5, 6]
        and [row["collision_minimum"] for row in partition_rows]
        == [0, 0, 0, 1, 1, 2, 3]
        and [row["parseval_nonexact_base"] for row in partition_rows]
        == [303, 298, 293, 289, 284, 276, 259]
    )
    _require(proved, "the p13 u=6 normalization changed")
    return {
        "p": P,
        "t": T,
        "k": K,
        "u": U,
        "prior_remaining_residues": [6],
        "hard_mean": "a_L=12+14*k_L",
        "hard_quotient_sum": hard_quotient_sum,
        "exact_floor_cell": "b=2 XNOR",
        "exact_XNOR_row_energy": 1,
        "exact_XNOR_rows_at_least": 2,
        "exact_parallel_unspecialized_identity": "hT=14*P0-51",
        "general_hard_parallel_identity": "P_L=P0+k_L",
        "exact_parallel_congruence": baseline["congruence"],
        "exact_parallel_candidates": list(exact_parallel_candidates),
        "exact_parallel_count": exact_parallel,
        "hT": hT,
        "hard_parallel_edge_total": hard_parallel_total,
        "opposite_parallel_edge_total": opposite_parallel_total,
        "opposite_mean_formula": "a(Q)=14*Q+hT-39=14*Q-34",
        "Q3_mean": q3_mean,
        "Q3_excluded_by_phase_zero_floors": True,
        "opposite_parallel_counts": list(opposite_parallel_counts),
        "hard_excess_partitions": partition_rows,
        "proved": proved,
    }


def _distance_permutation(multiplier: int) -> tuple[int, ...]:
    return tuple(
        min(multiplier * value % P, -multiplier * value % P) - 1
        for value in range(1, 7)
    )


@lru_cache(maxsize=1)
def translated_cut_coordinate_certificate() -> dict[str, object]:
    """Prove that every row catalog's finite coordinate box is complete."""
    cuts = translated_cut_vectors()
    cut_set = set(cuts)
    lower_terms = (
        (0, Fraction(-1, 18)),
        (6, Fraction(-1, 6)),
        (34, Fraction(-1, 18)),
    )
    upper_terms = (
        (63, Fraction(-1, 15)),
        (69, Fraction(-1, 30)),
        (71, Fraction(-1, 6)),
        (73, Fraction(-1, 30)),
    )
    lower_identity = tuple(
        Fraction(19, 9)
        + sum(weight * cuts[index][coordinate] for index, weight in lower_terms)
        for coordinate in range(6)
    )
    upper_identity = tuple(
        Fraction(29, 15)
        + sum(weight * cuts[index][coordinate] for index, weight in upper_terms)
        for coordinate in range(6)
    )
    permutations = tuple(_distance_permutation(value) for value in range(1, 7))
    invariant = all(
        tuple(cut[index] for index in permutation) in cut_set
        for permutation in permutations
        for cut in cuts
    )
    transitive = {permutation[0] for permutation in permutations} == set(range(6))

    specifications = {
        "H1": (0, 13, (-3, 3)),
        "H2": (-1, 13, (-5, 5)),
        "H3": (-2, 13, (-7, 7)),
        "O": (-9, -52, (-4, 1)),
    }
    rows = {}
    for name, (total, cut_upper, expected) in specifications.items():
        rational_lower = Fraction(19, 9) * total + sum(
            weight * cut_upper for _index, weight in lower_terms
        )
        lower = ceil(rational_lower)
        minus_coordinate_lower = Fraction(29, 15) * total + sum(
            weight * cut_upper for _index, weight in upper_terms
        )
        rational_upper = -minus_coordinate_lower
        upper = floor(rational_upper)
        _require((lower, upper) == expected, f"the {name} coordinate box changed")
        rows[name] = {
            "sum": total,
            "cut_upper": cut_upper,
            "rational_bounds": [str(rational_lower), str(rational_upper)],
            "integral_bounds": [lower, upper],
        }

    h3_enumeration_boxes = {}
    for extra in range(4):
        collision_lower = min(
            value for value in range(-20, 1) if signed_collision_floor(value) <= extra
        )
        collision_upper = max(
            value for value in range(0, 20) if signed_collision_floor(value) <= extra
        )
        enumeration = (-7 - extra, 6 + extra)
        exact_feasible = (
            max(rows["H3"]["integral_bounds"][0], collision_lower),
            min(rows["H3"]["integral_bounds"][1], collision_upper),
        )
        _require(
            enumeration[0] <= exact_feasible[0]
            and enumeration[1] >= exact_feasible[1],
            f"H3_D{extra} enumeration box is incomplete",
        )
        h3_enumeration_boxes[f"H3_D{extra}"] = {
            "extra_collision_budget": extra,
            "signed_collision_envelope": [collision_lower, collision_upper],
            "exact_feasible_bound": list(exact_feasible),
            "enumerated_superset_box": list(enumeration),
        }

    proved = bool(
        len(cuts) == 74
        and lower_identity == (1, 0, 0, 0, 0, 0)
        and upper_identity == (-1, 0, 0, 0, 0, 0)
        and invariant
        and transitive
        and set(rows) == set(specifications)
    )
    _require(proved, "the p13 u=6 coordinate certificate changed")
    return {
        "translated_cut_count": len(cuts),
        "lower_dual_identity": [str(value) for value in lower_identity],
        "upper_dual_identity": [str(value) for value in upper_identity],
        "cut_catalog_invariant_under_distance_action": invariant,
        "distance_action_transitive": transitive,
        "row_bounds": rows,
        "H3_collision_enumeration_boxes": h3_enumeration_boxes,
        "proved": proved,
    }


def _load_json(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ArithmeticError(f"{path}: expected a JSON object")
    return payload


def validate_low_root_payload(payload: dict[str, object]) -> dict[str, object]:
    """Validate and compact the pinned low-root U/G artifact."""
    imported = {
        row["repository_relative_path"]: row["sha256"]
        for row in payload["imported_helper_files"]
    }
    cases = payload["case_results"]
    expected = (
        ((1, 1, 1, 1, 1), 42, 43_184_232, 12_804_624, 293, (303, 329, 355), "69e06d7f3719f3b18e4f2b372ecdf72ea3031638eabb5a95ad853a6c719b5438"),
        ((2, 1, 1, 1), 280, 10_221_120, 3_057_264, 290, (298, 324, 350), "34e82cb53f67cd3f5afe0174e0d62e9cbfc5e55e5f1bf5f3c3bd4cf3d993bf75"),
    )
    summaries = []
    checks = [
        payload["repository_HEAD"] == EXPECTED_REPOSITORY_HEAD,
        payload["source_checker_path"] == "scripts/p13_u6_joint_ug_tables.py",
        payload["source_checker_sha256"]
        == EXPECTED_SUPPORTING_CODE_SHA256["scripts/p13_u6_joint_ug_tables.py"],
        payload["script_path"] == "scripts/p13_u6_low_root_ug_bound.py",
        payload["script_sha256"]
        == EXPECTED_SUPPORTING_CODE_SHA256["scripts/p13_u6_low_root_ug_bound.py"],
        payload["output_path"]
        == "evidence/e1_gmin_m4_prop15754_low_root_ug.json",
        imported == EXPECTED_IMPORTED_HELPER_SHA256,
        payload["table_maximum_sha256"] == EXPECTED_LOW_TABLE_MAXIMUM_SHA256,
        payload["direction_signs"].count(-1) == 7,
        payload["direction_signs"].count(1) == 7,
        payload["all_targets_excluded"] is True,
        payload["proved"] is True,
        len(cases) == len(expected),
    ]
    for case, row in zip(cases, expected):
        partition, fixed, pairs, compatible, maximum, targets, summary_hash = row
        target_rows = case["targets"]
        checks.extend(
            [
                tuple(case["partition"]) == partition,
                case["fixed_models_checked"] == fixed,
                case["expected_fixed_models"] == fixed,
                len(case["fixed_models"]) == fixed,
                case["coefficient_pairs_checked"] == pairs,
                case["locally_compatible_coefficient_pairs"] == compatible,
                case["global_maximum_separable_nonexact_energy"] == maximum,
                case["fixed_model_summary_sha256"] == summary_hash,
                tuple(target["target_nonexact_energy"] for target in target_rows)
                == targets,
                tuple(target["collision"] for target in target_rows) == (0, 1, 2),
                all(target["all_fixed_models_excluded"] for target in target_rows),
                all(target["maximum_separable_nonexact_energy"] == maximum for target in target_rows),
                case["original_CP_table_model_attainment_validation"]["energy_sum"]
                == maximum,
            ]
        )
        summaries.append(
            {
                "partition": list(partition),
                "fixed_models_checked": fixed,
                "coefficient_pairs_checked": pairs,
                "locally_compatible_coefficient_pairs": compatible,
                "maximum_separable_nonexact_energy": maximum,
                "targets": list(targets),
                "strict_deficits": [target - maximum for target in targets],
                "maximum_attained_in_original_CP_table_model": True,
                "fixed_model_summary_sha256": summary_hash,
            }
        )
    proved = all(checks)
    _require(proved, "the low-root U/G artifact changed")
    return {
        "method": "exhaustive common U/G separable-energy upper bound",
        "cases": summaries,
        "closed_partitions": [[1, 1, 1, 1, 1], [2, 1, 1, 1]],
        "graph_or_configuration_census_used": False,
        "proved": proved,
    }


def _row_catalog_manifest(payload: dict[str, object]) -> dict[str, tuple[object, ...]]:
    return {
        name: (
            row["row_count"],
            row["row_sha256"],
            row["table_tuple_count"],
            row["table_sha256"],
            row["minimum_energy"],
            row["maximum_energy"],
        )
        for name, row in payload["row_catalogs"].items()
    }


def validate_four_root_payload(
    payload: dict[str, object], partition: tuple[int, ...]
) -> dict[str, object]:
    """Validate and compact one pinned four-root U/G/J6 artifact."""
    expected_by_partition = {
        (2, 2, 1): ((0, 1, 2, 3), (293, 319, 345, 371), 84, 193, "12e941d504b56f0fe86903710295fef4174dd7174f449137d553e27aa75dc084"),
        (3, 1, 1): ((1, 2, 3, 4), (315, 341, 367, 393), 0, -1, "4f53cda18c2baa0c0354bb5f9a3ecbe5ed12ab4d8e11ba873c2f11161202b945"),
    }
    collisions, targets, compatible, maximum, reachable_hash = expected_by_partition[
        partition
    ]
    rows = payload["collision_rows"]
    coverage = payload["algorithm"]["coverage"]
    checks = [
        payload["repository_HEAD"] == EXPECTED_REPOSITORY_HEAD,
        payload["script_path"] == "scripts/p13_u6_four_root_ugj.py",
        payload["script_sha256"]
        == EXPECTED_SUPPORTING_CODE_SHA256["scripts/p13_u6_four_root_ugj.py"],
        tuple(payload["partition"]) == partition,
        payload["exact_XNOR_root_count"] == 4,
        payload["imported_helper_sha256"] == EXPECTED_IMPORTED_HELPER_SHA256,
        payload["moment_weights"]["W2"] == [1, 4, 9, 3, 12, 10],
        payload["moment_weights"]["W4"] == [1, 3, 3, 9, 1, 9],
        payload["moment_weights"]["W6"] == [1, 12, 1, 1, 12, 12],
        payload["W6_regression"] == [1, 12, 1, 1, 12, 12],
        payload["common_forms"]["opposite_key"]
        == "(U,G,J6)=(-N2,-N4-N2^2,-N6+N2^3)",
        _row_catalog_manifest(payload) == EXPECTED_ROW_CATALOGS,
        coverage
        == {
            "UGJ6_coefficient_triples_before_row_pruning": 218_320_284,
            "UGJ6_pairs_after_all_seven_opposite_rows": 336,
            "U_coefficients_passing_exact_root_QR": 7_644,
            "hard_excess_assignment_checks": 1_008,
            "hard_sign_root_sets_checked": 70,
        },
        tuple(payload["closed_collision_strata"]) == collisions,
        payload["remaining_collision_strata"] == [],
        payload["all_targets_excluded"] is True,
        payload["proved"] is True,
        tuple(row["collision"] for row in rows) == collisions,
        tuple(row["target_nonexact_energy"] for row in rows) == targets,
        all(row["hard_locally_compatible_form_assignments"] == compatible for row in rows),
        all(row["maximum_separable_nonexact_energy"] == maximum for row in rows),
        all(row["reachable_energy_values_sha256"] == reachable_hash for row in rows),
        all(row["joint_UGJ6_compatible"] is False for row in rows),
        all(row["all_form_coefficients_exhausted"] is True for row in rows),
    ]
    proved = all(checks)
    _require(proved, f"the four-root {partition} U/G/J6 artifact changed")
    return {
        "partition": list(partition),
        "method": "exhaustive common U/G/J6 coefficient join and energy bitsets",
        "closed_collision_counts": list(collisions),
        "targets": list(targets),
        "hard_locally_compatible_form_assignments_per_collision": compatible,
        "maximum_separable_nonexact_energy": maximum,
        "coverage": coverage,
        "sign_safe_opposite_key": payload["common_forms"]["opposite_key"],
        "W6_regression": payload["W6_regression"],
        "proved": proved,
    }


@lru_cache(maxsize=1)
def necessary_moment_artifact_certificate() -> dict[str, object]:
    supporting_code = {
        name: file_sha256(ROOT / name) for name in EXPECTED_SUPPORTING_CODE_SHA256
    }
    _require(
        supporting_code == EXPECTED_SUPPORTING_CODE_SHA256,
        "a p13 u=6 supporting verifier changed",
    )
    artifact_hashes = {name: file_sha256(path) for name, path in ARTIFACT_PATHS.items()}
    _require(
        artifact_hashes == EXPECTED_ARTIFACT_SHA256,
        "a p13 u=6 necessary-moment artifact changed",
    )
    low = validate_low_root_payload(_load_json(ARTIFACT_PATHS["low_root_UG"]))
    four_221 = validate_four_root_payload(
        _load_json(ARTIFACT_PATHS["four_root_221"]), (2, 2, 1)
    )
    four_311 = validate_four_root_payload(
        _load_json(ARTIFACT_PATHS["four_root_311"]), (3, 1, 1)
    )
    proved = bool(low["proved"] and four_221["proved"] and four_311["proved"])
    _require(proved, "a low/four-root p13 u=6 partition survived")
    return {
        "low_root_UG": low,
        "four_root_UGJ6": [four_221, four_311],
        "artifact_sha256": artifact_hashes,
        "supporting_code_sha256": supporting_code,
        "closed_partitions": [
            [1, 1, 1, 1, 1],
            [2, 1, 1, 1],
            [2, 2, 1],
            [3, 1, 1],
        ],
        "proved": proved,
    }


@lru_cache(maxsize=1)
def complete_partition_ledger() -> dict[str, object]:
    """Show that every collision stratum in every partition is excluded."""
    normalization = p13_u6_normalization_certificate()
    artifacts = necessary_moment_artifact_certificate()
    high = high_root_partition_certificate()
    equalities = u6_energy_ledger_certificate()
    row_maxima = {
        name: values[-1] for name, values in EXPECTED_ROW_CATALOGS.items()
    }
    low_rows = [
        {
            "partition": [1, 1, 1, 1, 1],
            "collision_minimum": 0,
            "raw_independent_energy_upper": 5 * row_maxima["H1"] + 7 * row_maxima["O"],
            "raw_collision_maximum": 2,
            "closed_collision_counts": [0, 1, 2],
            "method": "common U/G separable-energy bound",
        },
        {
            "partition": [2, 1, 1, 1],
            "collision_minimum": 0,
            "raw_independent_energy_upper": row_maxima["H2"] + 3 * row_maxima["H1"] + 7 * row_maxima["O"],
            "raw_collision_maximum": 2,
            "closed_collision_counts": [0, 1, 2],
            "method": "common U/G separable-energy bound",
        },
        {
            "partition": [2, 2, 1],
            "collision_minimum": 0,
            "raw_independent_energy_upper": 2 * row_maxima["H2"] + row_maxima["H1"] + 7 * row_maxima["O"],
            "raw_collision_maximum": 3,
            "closed_collision_counts": [0, 1, 2, 3],
            "method": "common U/G/J6 coefficient join",
        },
        {
            "partition": [3, 1, 1],
            "collision_minimum": 1,
            "raw_independent_energy_upper": row_maxima["H3_D3"] + 2 * row_maxima["H1"] + 7 * row_maxima["O"],
            "raw_collision_maximum": 4,
            "closed_collision_counts": [1, 2, 3, 4],
            "method": "common U/G/J6 coefficient join",
        },
    ]
    high_rows = high["partitions"]
    ledgers = low_rows + high_rows
    bases = {
        tuple(row["hard_excess_partition"]): row["parseval_nonexact_base"]
        for row in normalization["hard_excess_partitions"]
    }
    proved = bool(
        artifacts["proved"]
        and high["proved"]
        and equalities["proved"]
        and [row["raw_independent_energy_upper"] for row in low_rows]
        == [357, 364, 371, 417]
        and [row["raw_collision_maximum"] for row in low_rows] == [2, 2, 3, 4]
        and bases
        == {
            (1, 1, 1, 1, 1): 303,
            (2, 1, 1, 1): 298,
            (2, 2, 1): 293,
            (3, 1, 1): 289,
            (3, 2): 284,
            (4, 1): 276,
            (5,): 259,
        }
        and [row["partition"] for row in ledgers]
        == [list(partition) for partition in EXPECTED_PARTITIONS]
        and all(row.get("excluded", True) for row in high_rows)
    )
    _require(proved, "a p13 u=6 partition or collision stratum remained")
    return {
        "parseval_identity": (
            "required nonexact energy = partition base + 26*C"
        ),
        "row_energy_maxima_from_complete_catalogs": row_maxima,
        "partition_ledgers": ledgers,
        "independent_top_collision_cross_checks": {
            "partition_1^5_C2": equalities["partition_1^5"]["C_equals_2_excluded"],
            "partition_2_2_1_C3": equalities["partition_2_2_1"]["C_equals_3_excluded"],
        },
        "all_seven_partitions_closed": proved,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def proposition_15754() -> dict[str, object]:
    """Package the finite p13 endpoint close without flipping a global gate."""
    normalization = p13_u6_normalization_certificate()
    coordinates = translated_cut_coordinate_certificate()
    artifacts = necessary_moment_artifact_certificate()
    partitions = complete_partition_ledger()
    residual_open = residual_ii_k_ge_4p_ND_closed()
    proved = bool(
        normalization["proved"]
        and coordinates["proved"]
        and artifacts["proved"]
        and partitions["proved"]
        and not residual_open
    )
    _require(proved, "Proposition 15.754 certificate failed")
    return {
        "proposition": "15.754",
        "title": "Common-form energy closes the p13 fifth-shell endpoint",
        "result_status": (
            "exhaustive finite aggregate/common-form certificate and proved endpoint theorem"
        ),
        "statement": "the residual-(ii) branch p=13,t=4,k=60,u=6 is empty",
        "normalization_and_partition_exhaustiveness": normalization,
        "translated_cut_coordinate_completeness": coordinates,
        "necessary_moment_artifacts": artifacts,
        "complete_partition_ledger": partitions,
        "certificate_manifest": {
            "artifact_sha256": artifacts["artifact_sha256"],
            "supporting_code_sha256": artifacts["supporting_code_sha256"],
            "repository_head_recorded_by_artifacts": EXPECTED_REPOSITORY_HEAD,
        },
        "p13_t4_u6_closed": True,
        "p13_k_eq_60_closed": True,
        "remaining_p13_t4_residues": [],
        "fifth_shell_k_eq_4p_plus_8_closed_for_every_prime_p_ge_13": True,
        "finite_prime_aggregate_census_used": True,
        "finite_common_form_coefficient_certificate_used": True,
        "graph_or_configuration_census_used": False,
        "orbit_census_used": False,
        "coefficient_cell_census_used": False,
        "residual_ii_k_ge_4p_ND_closed": residual_open,
        "E1_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "remaining_scope": (
            "critical p=5,7; p=11 at k>=50; later p13 layers; layers beyond "
            "Proposition 15.752's band; residual (ii), E1, and the limit remain open"
        ),
        "proved": proved,
    }


def write_evidence(path: Path = EVIDENCE_PATH) -> Path:
    """Write the deterministic Proposition 15.754 evidence atomically."""
    write_json_atomic(path, proposition_15754())
    return path


def main() -> None:
    result = proposition_15754()
    path = write_evidence()
    print(
        json.dumps(
            {
                "proposition": result["proposition"],
                "result_status": result["result_status"],
                "p13_t4_u6_closed": result["p13_t4_u6_closed"],
                "p13_k_eq_60_closed": result["p13_k_eq_60_closed"],
                "residual_ii_k_ge_4p_ND_closed": result[
                    "residual_ii_k_ge_4p_ND_closed"
                ],
                "output": str(path),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()


__all__ = [
    "ARTIFACT_PATHS",
    "EVIDENCE_PATH",
    "EXPECTED_ARTIFACT_SHA256",
    "EXPECTED_PARTITIONS",
    "EXPECTED_ROW_CATALOGS",
    "EXPECTED_SUPPORTING_CODE_SHA256",
    "complete_partition_ledger",
    "file_sha256",
    "necessary_moment_artifact_certificate",
    "p13_u6_normalization_certificate",
    "proposition_15754",
    "translated_cut_coordinate_certificate",
    "validate_four_root_payload",
    "validate_low_root_payload",
    "write_evidence",
]
