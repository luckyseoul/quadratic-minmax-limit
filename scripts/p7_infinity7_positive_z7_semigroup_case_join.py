#!/usr/bin/env python3
"""Bounded exact semigroup case joins for the positive p=7, z=7 grade-three cases.

This is the case-level bridge between the audited Johnson-semigroup Hilbert
basis and the 324 orbit-0/branch-A representatives from the existing global
catalog join.  It targets the 51 representatives whose only high catalogs
have excess grade three:

* ``H1_S3_M4`` with high grades ``[3]`` (four cases);
* ``H2_S2_M3`` with high grades ``[3,3]`` (23 cases left after one prior
  rigorous rejection);
* ``H3_S5_M0`` with high grades ``[3,3,3]`` (24 cases).

Every selected projection retains all six characteristic-three torsion
coordinates and a deterministic subset of the 21 characteristic-seven
coordinates.  For every direction, the complete anchor-relative support is
constructed jointly modulo 3 and 7 from the same Hilbert-generator rows.
Grades zero through three are box-exact.  The eight direction supports are
then convolved globally and the exact projected target is tested.

A missing target in a completed projection is a rigorous rejection.  Target
presence, including exact full-group saturation, is necessary only.  Any
state-cap hit discards the partial support and is an explicit skip.  This
first implementation intentionally permits smoke runs only and makes no z=7
closure claim.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import p7_infinity7_positive_z7_global_catalog_join as old_join  # noqa: E402
import p7_infinity7_positive_z7_high_semigroup_support as semigroup  # noqa: E402
import p7_infinity7_positive_z7_pointed_affine_hull_multimod as affine  # noqa: E402
import p7_infinity7_positive_z7_torsion_support_projection as torsion  # noqa: E402


P = 7
AMBIENT_DIMENSION = 35
MODULI = (3, 7)
EXPECTED_REPRESENTATIVES = 324
EXPECTED_CURRENT_COUNTS = {"rejected": 87, "surviving": 159, "skipped": 78}
EXPECTED_CURRENT_CASE_RESULTS_SHA256 = (
    "c34cc913c27910e3876e1b78aed0e9c8c2f42cb2f4368f95054bcd6ead1db7a7"
)
TARGET_PATTERNS = {
    "H1_S3_M4": (3,),
    "H2_S2_M3": (3, 3),
    "H3_S5_M0": (3, 3, 3),
}
EXPECTED_TARGET_PATTERN_COUNTS = {
    "H1_S3_M4": 4,
    "H2_S2_M3": 23,
    "H3_S5_M0": 24,
}
EXPECTED_GRADE_CATALOG_ROWS = {0: 1, 1: 56, 2: 1_764, 3: 37_856}

ARCHIVE_ROOT = Path(
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/2026-08-30-p7-z7"
)
DEFAULT_PARENT_INPUT = ARCHIVE_ROOT / "p7_inf7_z7_affine_hull_multimod_full.json"
DEFAULT_CURRENT_JOIN = (
    ARCHIVE_ROOT / "global-join/p7_z7_global_join_full_explicit.json"
)
DEFAULT_HILBERT_BASIS = ARCHIVE_ROOT / "semigroup/p7_johnson_semigroup.gen"

DEFAULT_MOD7_COORDINATE_COUNT = 3
DEFAULT_PROJECTION_COUNT = 2
DEFAULT_SMOKE_CASES = 3
DEFAULT_STATE_CAP = 300_000
DEFAULT_PAIR_CHUNK_CAP = 100_000
MAX_SMOKE_MOD7_COORDINATES = 5
MAX_SMOKE_PROJECTIONS = 8
MAX_SMOKE_CASES = 51
MAX_SMOKE_STATE_CAP = 2_000_000
MAX_AUTOMATIC_PROJECTION_CANDIDATES = 128


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def array_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def json_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def pattern_name(pattern: tuple[int, int, int]) -> str:
    return f"H{pattern[0]}_S{pattern[1]}_M{pattern[2]}"


def load_current_problem(
    parent_path: Path, current_join_path: Path
) -> tuple[dict, list[dict], dict, list[dict], dict]:
    """Rebuild and bind all 324 representatives to the current decisions."""
    parent_payload, survivors, parent_audit = old_join.load_parent_input(parent_path)
    representatives = old_join.representative_survivors(survivors)
    require(len(representatives) == EXPECTED_REPRESENTATIVES, "representative census changed")

    raw = current_join_path.read_bytes()
    current = json.loads(raw)
    require(
        current.get("experiment") == old_join.EXPERIMENT,
        "current-decision input is not the global catalog join",
    )
    require(current.get("smoke_test") is False, "current decisions came from a smoke run")
    require(
        current.get("complete_selected_universe_run") is True,
        "current decisions do not cover their selected universe",
    )
    require(
        current.get("representative_selection", {}).get("enabled") is True,
        "current decisions are not the 324 representatives",
    )
    rows = current["case_results"]
    require(len(rows) == EXPECTED_REPRESENTATIVES, "current decision census changed")
    require(
        old_join.canonical_case_digest(rows) == current["case_results_sha256"],
        "current case-result certificate failed",
    )
    require(
        current["case_results_sha256"] == EXPECTED_CURRENT_CASE_RESULTS_SHA256,
        "current global-join decision set is not the audited 2026-08-30 set",
    )
    require(
        current["input_provenance"]["file_sha256"] == parent_audit["file_sha256"],
        "current decisions were computed from another parent artifact",
    )
    representative_keys = [str(row["case_key"]) for row in representatives]
    current_keys = [str(row["case_key"]) for row in rows]
    require(current_keys == representative_keys, "current decisions and representatives differ")

    counts = {
        "rejected": sum(bool(row["rigorously_rejected"]) for row in rows),
        "surviving": sum(bool(row["necessary_only_survivor"]) for row in rows),
        "skipped": sum(bool(row["skipped"]) for row in rows),
    }
    require(counts == EXPECTED_CURRENT_COUNTS, "current global-join counts changed")
    require(
        all(
            sum(
                bool(row[field])
                for field in ("rigorously_rejected", "necessary_only_survivor", "skipped")
            )
            == 1
            for row in rows
        ),
        "a current decision has ambiguous status",
    )

    representative_by_key = {str(row["case_key"]): row for row in representatives}
    targets = [
        representative_by_key[str(row["case_key"])]
        for row in rows
        if str(row["catalog_pattern"]) in TARGET_PATTERNS
        and not bool(row["rigorously_rejected"])
    ]
    target_keys = {str(row["case_key"]) for row in targets}
    target_counts = Counter(
        str(row["catalog_pattern"])
        for row in rows
        if str(row["case_key"]) in target_keys
    )
    require(dict(target_counts) == EXPECTED_TARGET_PATTERN_COUNTS, "grade-three target census changed")
    require(len(targets) == 51, "expected exactly 51 remaining grade-three representatives")

    current_by_key = {str(row["case_key"]): row for row in rows}
    return parent_payload, representatives, current_by_key, targets, {
        "parent_input": parent_audit,
        "current_global_join": {
            "path": str(current_join_path.resolve()),
            "file_bytes": len(raw),
            "file_sha256": hashlib.sha256(raw).hexdigest(),
            "case_results_sha256": current["case_results_sha256"],
            "representative_case_count": len(rows),
            "decision_counts": counts,
            "all_324_current_decisions_revalidated": True,
        },
        "target_pattern_counts": dict(sorted(target_counts.items())),
        "target_case_count": len(targets),
        "target_case_keys_sha256": json_sha256([row["case_key"] for row in targets]),
        "one_previously_rejected_H2_S2_M3_case_excluded_from_target": True,
    }


def reconstruct_and_validate(
    parent_payload: dict, representatives: list[dict]
) -> tuple[dict, dict]:
    rebuilt = old_join.reconstruction()
    parent_match = old_join.validate_reconstruction_against_parent(parent_payload, rebuilt)
    validated_keys = []
    for row in representatives:
        old_join.validate_parent_survivor(row, rebuilt)
        validated_keys.append(str(row["case_key"]))
    require(len(validated_keys) == EXPECTED_REPRESENTATIVES, "reconstruction missed a representative")
    return rebuilt, {
        "parent_reconstruction": parent_match,
        "orbit0_branchA_representatives_reconstructed": len(validated_keys),
        "representative_case_keys_sha256": json_sha256(validated_keys),
        "all_representative_leaf_metadata_and_pointed_systems_revalidated": True,
    }


def leaf_grade(orbit: dict, leaf: dict, direction: int) -> int:
    grade, _floor_mean, _floor = semigroup.excess_grade(
        int(orbit["masks"][direction]), int(leaf["scaled_means"][direction])
    )
    return int(grade)


def audit_target_grades(
    targets: list[dict], current_by_key: dict[str, dict], rebuilt: dict
) -> dict:
    class_grade = {"U": 0, "S": 1, "M": 2, "H": 3}
    rows = []
    grade_vectors = []
    for target in targets:
        orbit, leaf, _system, _factory = old_join.validate_parent_survivor(target, rebuilt)
        grades = tuple(leaf_grade(orbit, leaf, direction) for direction in range(P + 1))
        classes = tuple(str(value) for value in leaf["catalog_classes"])
        require(
            all(grades[direction] == class_grade[classes[direction]] for direction in range(P + 1)),
            "target catalog class/excess-grade correspondence changed",
        )
        high_grades = tuple(grades[direction] for direction in leaf["high_directions"])
        current_pattern = str(current_by_key[str(target["case_key"])]["catalog_pattern"])
        require(high_grades == TARGET_PATTERNS[current_pattern], "target high-grade profile changed")
        require(sum(grades) == 14, "target total directional grade changed")
        grade_vectors.append(list(grades))
        rows.append(
            {
                "case_key": target["case_key"],
                "catalog_pattern": current_pattern,
                "directional_excess_grades": list(grades),
                "high_directions": list(leaf["high_directions"]),
                "high_excess_grades": list(high_grades),
                "total_directional_excess_grade": sum(grades),
                "prior_global_join_decision": current_by_key[str(target["case_key"])][
                    "decision_status"
                ],
            }
        )
    return {
        "case_count": len(rows),
        "all_catalog_classes_equal_grades_U0_S1_M2_H3": True,
        "all_high_catalogs_have_box_exact_grade_three": True,
        "all_cases_have_total_directional_grade_fourteen": True,
        "directional_grade_vectors_sha256": json_sha256(grade_vectors),
        "case_grade_records_sha256": json_sha256(rows),
        "records": rows,
    }


def project_equation_vector(
    vector: np.ndarray,
    quotient_data: dict[int, dict],
    mod7_coordinates: tuple[int, ...],
) -> np.ndarray:
    source = np.ascontiguousarray(vector, dtype=np.int64)
    components = []
    for modulus in MODULI:
        complement = np.asarray(quotient_data[modulus]["complement"], dtype=np.int64)
        selected = tuple(range(len(complement))) if modulus == 3 else mod7_coordinates
        component = complement[list(selected)] @ (source % modulus) % modulus
        components.append(np.ascontiguousarray(component, dtype=np.uint8))
    return np.ascontiguousarray(np.concatenate(components), dtype=np.uint8)


def full_generator_projections(
    basis: np.ndarray,
    orbit: dict,
    quotient_data: dict[int, dict],
) -> tuple[dict[int, np.ndarray], dict]:
    rows = {}
    audits = []
    all_mod7 = tuple(range(len(quotient_data[7]["complement"])))
    for direction in range(P + 1):
        mask = int(orbit["masks"][direction])
        mapped = np.ascontiguousarray(
            np.stack([affine.map_canonical_row(mask, row) for row in basis]),
            dtype=np.int64,
        )
        variation = np.ascontiguousarray(-2 * mapped, dtype=np.int64)
        joint = semigroup.project_rows(variation, direction, quotient_data, all_mod7)
        require(joint.shape == (len(basis), 27), "full torsion generator projection shape changed")
        rows[direction] = joint
        audits.append(
            {
                "direction": direction,
                "mask": mask,
                "mapped_basis_sha256_int64": array_sha256(mapped),
                "full_F3_6_F7_21_generator_projection_sha256_uint8": array_sha256(joint),
            }
        )
    return rows, {
        "directions": len(rows),
        "joint_full_projection_shape_per_direction": [len(basis), 27],
        "all_six_mod3_and_all_twenty_one_mod7_coordinates_derived_before_selection": True,
        "same_Hilbert_generator_row_supplies_both_prime_components": True,
        "records_sha256": json_sha256(audits),
        "records": audits,
    }


def projection_intrinsic_score(
    coordinates: tuple[int, ...],
    full_rows: dict[int, np.ndarray],
    degrees: np.ndarray,
) -> tuple[int, ...]:
    columns = tuple(range(6)) + tuple(6 + value for value in coordinates)
    by_degree: dict[int, list[int]] = {1: [], 2: [], 3: []}
    for direction in range(P + 1):
        selected = np.ascontiguousarray(full_rows[direction][:, columns], dtype=np.uint8)
        for degree in (1, 2, 3):
            signatures = selected[degrees == degree]
            keys = signatures.view(np.dtype((np.void, signatures.shape[1]))).reshape(-1)
            by_degree[degree].append(len(np.unique(keys)))
    return (
        min(by_degree[3]),
        sum(by_degree[3]),
        min(by_degree[2]),
        sum(by_degree[2]),
        min(by_degree[1]),
        sum(by_degree[1]),
    )


def evenly_sample_combinations(
    coordinate_count: int, maximum_candidates: int
) -> tuple[list[tuple[int, ...]], int, bool]:
    all_rows = list(itertools.combinations(range(21), coordinate_count))
    total = len(all_rows)
    if total <= maximum_candidates:
        return all_rows, total, False
    indices = {
        round(index * (total - 1) / (maximum_candidates - 1))
        for index in range(maximum_candidates)
    }
    sampled = [all_rows[index] for index in sorted(indices)]
    require(len(sampled) == maximum_candidates, "deterministic combination sampling collided")
    return sampled, total, True


def parse_projection_spec(value: str) -> tuple[tuple[int, ...], ...]:
    projections = []
    for group in value.split(";"):
        group = group.strip()
        if not group:
            continue
        try:
            coordinates = tuple(sorted(int(item.strip()) for item in group.split(",")))
        except ValueError as error:
            raise AssertionError("projection coordinates must be integers") from error
        require(len(coordinates) == len(set(coordinates)), "projection repeats a coordinate")
        require(coordinates and all(0 <= value < 21 for value in coordinates), "bad F7 coordinate")
        projections.append(coordinates)
    require(projections, "explicit projection specification is empty")
    require(len(projections) == len(set(projections)), "explicit projection repeated")
    return tuple(sorted(projections))


def select_projection_subsets(
    *,
    full_rows: dict[int, np.ndarray],
    degrees: np.ndarray,
    coordinate_count: int,
    projection_count: int,
    explicit: tuple[tuple[int, ...], ...] | None,
) -> tuple[tuple[tuple[int, ...], ...], dict]:
    if explicit is not None:
        require(len(explicit) <= MAX_SMOKE_PROJECTIONS, "too many explicit smoke projections")
        widths = {len(row) for row in explicit}
        require(len(widths) == 1, "explicit projections have different widths")
        selected = explicit
        scored = [
            {
                "mod7_coordinates": list(row),
                "intrinsic_score": list(projection_intrinsic_score(row, full_rows, degrees)),
            }
            for row in selected
        ]
        return selected, {
            "selection_mode": "explicit_sorted_coordinate_subsets",
            "candidate_subsets_considered": len(selected),
            "total_coordinate_subsets": math.comb(21, len(selected[0])),
            "candidate_pool_deterministically_sampled": False,
            "score_order": [
                "minimum_grade3_generator_signatures",
                "sum_grade3_generator_signatures",
                "minimum_grade2_generator_signatures",
                "sum_grade2_generator_signatures",
                "minimum_grade1_generator_signatures",
                "sum_grade1_generator_signatures",
            ],
            "selected": scored,
        }

    require(1 <= coordinate_count <= MAX_SMOKE_MOD7_COORDINATES, "automatic F7 width out of range")
    require(1 <= projection_count <= MAX_SMOKE_PROJECTIONS, "automatic projection count out of range")
    candidates, total, sampled = evenly_sample_combinations(
        coordinate_count, MAX_AUTOMATIC_PROJECTION_CANDIDATES
    )
    scored_rows = [
        (projection_intrinsic_score(row, full_rows, degrees), row) for row in candidates
    ]
    scored_rows.sort(key=lambda item: tuple(-value for value in item[0]) + item[1])
    require(projection_count <= len(scored_rows), "projection count exceeds candidate pool")
    selected = tuple(row for _score, row in scored_rows[:projection_count])
    repeated = tuple(row for _score, row in sorted(
        scored_rows, key=lambda item: tuple(-value for value in item[0]) + item[1]
    )[:projection_count])
    require(selected == repeated, "deterministic projection ranking did not reproduce")
    return selected, {
        "selection_mode": "deterministic_generator_signature_ranking",
        "mod7_coordinates_per_projection": coordinate_count,
        "requested_projection_count": projection_count,
        "total_coordinate_subsets": total,
        "candidate_subsets_considered": len(candidates),
        "candidate_pool_deterministically_sampled": sampled,
        "candidate_sampling_rule": (
            "all lexicographic combinations"
            if not sampled
            else f"{MAX_AUTOMATIC_PROJECTION_CANDIDATES} evenly spaced indices in lexicographic combination order"
        ),
        "score_order": [
            "minimum_grade3_generator_signatures",
            "sum_grade3_generator_signatures",
            "minimum_grade2_generator_signatures",
            "sum_grade2_generator_signatures",
            "minimum_grade1_generator_signatures",
            "sum_grade1_generator_signatures",
        ],
        "tie_break": "lexicographically smallest F7 coordinate tuple",
        "selected": [
            {"mod7_coordinates": list(row), "intrinsic_score": list(score)}
            for score, row in scored_rows[:projection_count]
        ],
        "selection_recomputed_identically": True,
    }


def projected_direction_supports(
    *,
    projection: tuple[int, ...],
    full_generator_rows: dict[int, np.ndarray],
    basis: np.ndarray,
    degrees: np.ndarray,
    orbit: dict,
    anchors: affine.AnchorFactory,
    quotient_data: dict[int, dict],
    state_cap: int,
    pair_chunk_cap: int,
) -> tuple[dict[int, dict[int, np.ndarray]], dict]:
    moduli = (3,) * 6 + (7,) * len(projection)
    codec = semigroup.MixedRadixCodec(moduli)
    columns = tuple(range(6)) + tuple(6 + value for value in projection)
    support_table: dict[int, dict[int, np.ndarray]] = {}
    direction_audits = []
    for direction in range(P + 1):
        joint_generators = np.ascontiguousarray(
            full_generator_rows[direction][:, columns], dtype=np.uint8
        )
        generator_codes = {
            degree: codec.unique_codes(joint_generators[degrees == degree])
            for degree in (1, 2, 3)
        }
        raw_supports, recurrence_records = semigroup.support_recurrence(
            generator_codes,
            codec,
            state_cap,
            pair_chunk_cap,
            max_grade=3,
        )
        mask = int(orbit["masks"][direction])
        translated: dict[int, np.ndarray] = {}
        grade_audits = []
        for grade in range(4):
            recurrence_record = recurrence_records[grade]
            if grade not in raw_supports:
                require(not recurrence_record["completed"], "missing support was labeled complete")
                grade_audits.append(
                    {
                        "excess_grade": grade,
                        "decision_status": recurrence_record["decision_status"],
                        "completed": False,
                        "skipped": True,
                        "partial_support_used": False,
                        "rigorous_rejection_allowed": False,
                    }
                )
                continue

            _grade, floor_mean, floor = semigroup.excess_grade(
                mask, 0 if mask.bit_count() == 7 else 8
            )
            mean = floor_mean + 4 * grade
            anchor = anchors.get(mask, mean)
            offset = semigroup.project_rows(
                (anchor - floor)[None, :], direction, quotient_data, projection
            )[0]
            observed = codec.translate(raw_supports[grade], offset)

            # The repository's complete direct catalogs independently calibrate
            # every grade used in this first case-level attack, including H3.
            catalog = affine.mapped_catalog(mask, mean).astype(np.int64)
            require(
                len(catalog) == EXPECTED_GRADE_CATALOG_ROWS[grade],
                f"grade-{grade} complete catalog census changed",
            )
            delta = np.ascontiguousarray(anchor[None, :] - catalog, dtype=np.int64)
            require(not np.any(delta.sum(axis=1)), "direct catalog delta changed mean")
            require(not np.any(anchors.kernel_rows @ delta.T), "direct catalog left degree two")
            require(not np.any(delta % 2), "direct catalog delta changed parity")
            direct = codec.unique_codes(
                semigroup.project_rows(delta, direction, quotient_data, projection)
            )
            require(
                np.array_equal(observed, direct),
                f"direction {direction} grade-{grade} semigroup/direct support mismatch",
            )
            translated[grade] = observed
            grade_audits.append(
                {
                    "excess_grade": grade,
                    "scaled_mean": mean,
                    "complete_catalog_rows": len(catalog),
                    "projected_unique_states": len(observed),
                    "support_sha256_uint64": array_sha256(observed.astype("<u8", copy=False)),
                    "full_projected_group_saturated": len(observed) == codec.group_size,
                    "semigroup_recurrence_equals_complete_direct_catalog_projection": True,
                    "box_exact": True,
                    "same_source_row_used_mod3_mod7": True,
                }
            )
        support_table[direction] = translated
        direction_audits.append(
            {
                "direction": direction,
                "mask": mask,
                "projected_generator_unique_states_by_degree": {
                    str(degree): len(generator_codes[degree]) for degree in (1, 2, 3)
                },
                "projected_joint_generator_sha256_uint8": array_sha256(joint_generators),
                "recurrence_records_sha256": json_sha256(recurrence_records),
                "grades": grade_audits,
                "all_grades_zero_through_three_completed": set(translated) == set(range(4)),
            }
        )
    return support_table, {
        "mod7_coordinates": list(projection),
        "projected_group": f"F3^6 x F7^{len(projection)}",
        "projected_group_size": codec.group_size,
        "state_cap": state_cap,
        "pair_chunk_cap": pair_chunk_cap,
        "directions": direction_audits,
        "all_direction_grades_zero_through_three_completed": all(
            row["all_grades_zero_through_three_completed"] for row in direction_audits
        ),
        "all_completed_grade_supports_match_complete_direct_catalogs": True,
    }


def convolve_support_sequence(
    factors: tuple[tuple[int, np.ndarray], ...],
    codec: semigroup.MixedRadixCodec,
    state_cap: int,
    pair_chunk_cap: int,
) -> tuple[np.ndarray | None, dict]:
    current = np.asarray([0], dtype=np.uint64)
    records = []
    ordered = tuple(sorted(factors, key=lambda row: (len(row[1]), row[0])))
    for position, (direction, support) in enumerate(ordered):
        require(len(support) > 0, "direction support is empty")
        if len(current) == codec.group_size:
            records.append(
                {
                    "direction": direction,
                    "direction_states": len(support),
                    "status": "exact_full_group_translation_shortcut",
                    "input_states": len(current),
                    "output_states": len(current),
                }
            )
            continue
        output, audit = semigroup.bounded_minkowski_union(
            current,
            support,
            np.empty(0, dtype=np.uint64),
            codec,
            state_cap,
            pair_chunk_cap,
        )
        record = {
            "direction": direction,
            "direction_states": len(support),
            "input_states": len(current),
            "output_states_lower_bound_or_exact": len(output),
            **audit,
        }
        records.append(record)
        if audit["status"] == "state_cap_exceeded":
            return None, {
                "status": "skipped_state_cap",
                "completed": False,
                "skipped": True,
                "state_cap": state_cap,
                "partial_support_used": False,
                "first_incomplete_direction": direction,
                "unprocessed_directions": [row[0] for row in ordered[position + 1 :]],
                "steps": records,
            }
        current = output
    require(len(current) <= state_cap or len(current) == codec.group_size, "completed side escaped cap")
    saturated = len(current) == codec.group_size
    return current, {
        "status": "exact_full_group_saturation" if saturated else "complete_exact_support",
        "completed": True,
        "skipped": False,
        "state_count": len(current),
        "support_sha256_uint64": array_sha256(current.astype("<u8", copy=False)),
        "full_projected_group_saturated": saturated,
        "convolution_order": [row[0] for row in ordered],
        "steps": records,
    }


def meet_target(
    left: np.ndarray,
    right: np.ndarray,
    target_digits: np.ndarray,
    codec: semigroup.MixedRadixCodec,
    chunk_states: int,
) -> dict:
    require(len(left) and len(right), "case join side is empty")
    require(np.array_equal(left, np.unique(left)), "left support is not sorted unique")
    target = np.ascontiguousarray(target_digits, dtype=np.uint8)
    require(target.shape == (len(codec.moduli),), "target width changed")
    matches = 0
    first: tuple[int, int] | None = None
    certificate = hashlib.sha256()
    for start in range(0, len(right), chunk_states):
        block = right[start : start + chunk_states]
        rows = codec.decode(block).astype(np.int16)
        needed_digits = (
            target[None, :].astype(np.int16) - rows
        ) % codec.moduli_array[None, :].astype(np.int16)
        needed = codec.encode(needed_digits.astype(np.uint8))
        positions = np.searchsorted(left, needed)
        candidates = np.flatnonzero(positions < len(left))
        if not len(candidates):
            continue
        hits = candidates[left[positions[candidates]] == needed[candidates]]
        if not len(hits):
            continue
        left_indices = positions[hits].astype(np.uint64)
        right_indices = (start + hits).astype(np.uint64)
        matches += len(hits)
        certificate.update(left_indices.astype("<u8", copy=False).tobytes())
        certificate.update(right_indices.astype("<u8", copy=False).tobytes())
        if first is None:
            first = (int(left_indices[0]), int(right_indices[0]))
    return {
        "matching_projected_support_pairs": matches,
        "first_matching_side_indices": list(first) if first is not None else None,
        "matching_pair_index_certificate_sha256": certificate.hexdigest(),
        "exact_sorted_uint64_intersection": True,
        "hash_collision_assumption_used": False,
    }


def evaluate_projected_case(
    *,
    target_row: dict,
    current_row: dict,
    rebuilt: dict,
    common: np.ndarray,
    quotient_data: dict[int, dict],
    projection: tuple[int, ...],
    supports: dict[int, dict[int, np.ndarray]],
    state_cap: int,
    pair_chunk_cap: int,
) -> dict:
    orbit, leaf, system, _factory = old_join.validate_parent_survivor(target_row, rebuilt)
    codec = semigroup.MixedRadixCodec((3,) * 6 + (7,) * len(projection))
    grades = tuple(leaf_grade(orbit, leaf, direction) for direction in range(P + 1))
    factors = []
    missing = []
    for direction, grade in enumerate(grades):
        support = supports[direction].get(grade)
        if support is None:
            missing.append({"direction": direction, "excess_grade": grade})
        else:
            factors.append((direction, support))
    base_record = {
        "mod7_coordinates": list(projection),
        "projected_group": f"F3^6 x F7^{len(projection)}",
        "projected_group_size": codec.group_size,
        "directional_excess_grades": list(grades),
        "direction_support_sizes": [
            len(supports[direction][grade]) if grade in supports[direction] else None
            for direction, grade in enumerate(grades)
        ],
    }
    if missing:
        return {
            **base_record,
            "decision_status": "skipped_incomplete_direction_support",
            "rigorously_rejected": False,
            "necessary_only": False,
            "skipped": True,
            "missing_direction_supports": missing,
            "partial_support_used": False,
        }

    anchor_rhs, _raw_syndromes = affine.anchor_rhs_and_raw_syndromes(
        orbit, leaf, system, rebuilt["anchors"]
    )
    base_digits = project_equation_vector(anchor_rhs, quotient_data, projection)
    target_digits = np.ascontiguousarray(
        (-base_digits.astype(np.int16))
        % np.asarray(codec.moduli, dtype=np.int16),
        dtype=np.uint8,
    )

    # Common exact dependencies are not discarded by assumption: exact grade,
    # parity, and degree-two catalogs force them, and one direct catalog row in
    # every direction checks the resulting syndrome for this exact leaf.
    exact_rhs = anchor_rhs.copy()
    for direction, grade in enumerate(grades):
        mask = int(orbit["masks"][direction])
        mean = int(leaf["scaled_means"][direction])
        catalog_row = affine.mapped_catalog(mask, mean)[0].astype(np.int64)
        anchor = rebuilt["anchors"].get(mask, mean)
        block = slice(1 + AMBIENT_DIMENSION * direction, 1 + AMBIENT_DIMENSION * (direction + 1))
        exact_rhs[block] += anchor - catalog_row
    require(not np.any(common @ exact_rhs), "common exact dependency syndrome did not vanish")

    sizes = tuple(len(row[1]) for row in factors)
    directions = tuple(row[0] for row in factors)
    partition = old_join.balanced_partition(directions, sizes)
    by_direction = {direction: support for direction, support in factors}
    left_factors = tuple(
        (direction, by_direction[direction]) for direction in partition["left_directions"]
    )
    right_factors = tuple(
        (direction, by_direction[direction]) for direction in partition["right_directions"]
    )
    left, left_audit = convolve_support_sequence(
        left_factors, codec, state_cap, pair_chunk_cap
    )
    right, right_audit = convolve_support_sequence(
        right_factors, codec, state_cap, pair_chunk_cap
    )

    # One full side plus any nonempty opposite product is the full group.
    saturated_side = (
        left is not None
        and len(left) == codec.group_size
        or right is not None
        and len(right) == codec.group_size
    )
    common_fields = {
        **base_record,
        "base_digits_sha256_uint8": array_sha256(base_digits),
        "target_digits": target_digits.tolist(),
        "target_digits_sha256_uint8": array_sha256(target_digits),
        "balanced_partition": partition,
        "left_support": left_audit,
        "right_support": right_audit,
        "common_exact_dependency_syndrome_checked_with_direct_catalog_rows": True,
        "same_direction_support_codes_couple_mod3_mod7_before_global_convolution": True,
    }
    if saturated_side:
        return {
            **common_fields,
            "decision_status": "necessary_only_exact_full_group_saturation",
            "rigorously_rejected": False,
            "necessary_only": True,
            "skipped": False,
            "exact_projection_completed": True,
            "target_present_by_exact_full_group_saturation": True,
            "saturation_proves_feasibility": False,
        }
    if left is None or right is None:
        return {
            **common_fields,
            "decision_status": "skipped_case_convolution_state_cap",
            "rigorously_rejected": False,
            "necessary_only": False,
            "skipped": True,
            "exact_projection_completed": False,
            "partial_support_used": False,
        }

    joined = meet_target(left, right, target_digits, codec, pair_chunk_cap)
    rejected = joined["matching_projected_support_pairs"] == 0
    return {
        **common_fields,
        "decision_status": (
            "rigorous_exact_semigroup_projection_rejection"
            if rejected
            else "necessary_only_exact_semigroup_projection_survivor"
        ),
        "rigorously_rejected": rejected,
        "necessary_only": not rejected,
        "skipped": False,
        "exact_projection_completed": True,
        "join": joined,
        "missing_target_is_rigorous_rejection": rejected,
        "target_presence_proves_feasibility": False,
        "prior_global_join_decision": current_row["decision_status"],
    }


def select_smoke_cases(
    targets: list[dict], current_by_key: dict[str, dict], count: int
) -> tuple[list[dict], dict]:
    require(1 <= count <= min(MAX_SMOKE_CASES, len(targets)), "smoke case count out of range")
    buckets = {
        pattern: [
            row
            for row in targets
            if current_by_key[str(row["case_key"])]["catalog_pattern"] == pattern
        ]
        for pattern in TARGET_PATTERNS
    }
    selected = []
    depth = 0
    while len(selected) < count:
        added = False
        for pattern in TARGET_PATTERNS:
            if depth < len(buckets[pattern]) and len(selected) < count:
                selected.append(buckets[pattern][depth])
                added = True
        require(added, "smoke round-robin selection stalled")
        depth += 1
    return selected, {
        "selection_rule": "round-robin in declared target-pattern order, preserving current decision order",
        "selected_case_count": len(selected),
        "selected_pattern_counts": dict(
            sorted(
                Counter(
                    current_by_key[str(row["case_key"])]["catalog_pattern"] for row in selected
                ).items()
            )
        ),
        "selected_case_keys": [row["case_key"] for row in selected],
        "selected_case_keys_sha256": json_sha256([row["case_key"] for row in selected]),
        "full_51_case_coverage_claimed": False,
    }


def manufactured_case_join_audit() -> dict:
    codec = semigroup.MixedRadixCodec((3, 7))

    # A positive and a negative target are compared with direct Cartesian
    # enumeration using joint rows, so independent-prime mixing cannot hide.
    factors = (
        (0, codec.unique_codes(np.asarray([[0, 1], [1, 0]], dtype=np.uint8))),
        (1, codec.unique_codes(np.asarray([[0, 6], [2, 0]], dtype=np.uint8))),
        (2, codec.unique_codes(np.asarray([[1, 1], [2, 2]], dtype=np.uint8))),
    )
    left, left_audit = convolve_support_sequence(factors[:2], codec, 21, 21)
    right, right_audit = convolve_support_sequence(factors[2:], codec, 21, 21)
    require(left is not None and right is not None, "manufactured support unexpectedly capped")
    brute = {
        tuple(
            (
                np.asarray(a, dtype=np.int16)
                + np.asarray(b, dtype=np.int16)
                + np.asarray(c, dtype=np.int16)
            )
            % np.asarray((3, 7), dtype=np.int16)
        )
        for a in codec.decode(factors[0][1]).tolist()
        for b in codec.decode(factors[1][1]).tolist()
        for c in codec.decode(factors[2][1]).tolist()
    }
    positive_target = np.asarray(sorted(brute)[0], dtype=np.uint8)
    all_group = set(itertools.product(range(3), range(7)))
    absent_rows = sorted(all_group - brute)
    require(absent_rows, "manufactured support unexpectedly filled its group")
    negative_target = np.asarray(absent_rows[0], dtype=np.uint8)
    positive = meet_target(left, right, positive_target, codec, 21)
    negative = meet_target(left, right, negative_target, codec, 21)
    require(positive["matching_projected_support_pairs"] > 0, "manufactured witness was lost")
    require(negative["matching_projected_support_pairs"] == 0, "manufactured missing target passed")

    trap_support = codec.unique_codes(np.asarray([[0, 1], [1, 0]], dtype=np.uint8))
    zero = np.asarray([0, 0], dtype=np.uint8)
    trap = meet_target(trap_support, np.asarray([0], dtype=np.uint64), zero, codec, 21)
    require(trap["matching_projected_support_pairs"] == 0, "same-row prime trap failed")
    require(
        np.any(codec.decode(trap_support)[:, 0] == 0)
        and np.any(codec.decode(trap_support)[:, 1] == 0),
        "same-row trap lacks independent-prime marginal hits",
    )

    full = np.arange(codec.group_size, dtype=np.uint64)
    saturated, saturated_audit = convolve_support_sequence(
        ((0, full), (1, trap_support)), codec, 21, 21
    )
    require(
        saturated is not None and len(saturated) == codec.group_size,
        "manufactured saturation was not exact",
    )

    cap_codec = semigroup.MixedRadixCodec((7, 7))
    axis_a = cap_codec.unique_codes(
        np.asarray([[value, 0] for value in range(5)], dtype=np.uint8)
    )
    axis_b = cap_codec.unique_codes(
        np.asarray([[0, value] for value in range(5)], dtype=np.uint8)
    )
    capped, capped_audit = convolve_support_sequence(
        ((0, axis_a), (1, axis_b)), cap_codec, 10, 25
    )
    require(capped is None and capped_audit["status"] == "skipped_state_cap", "cap trap failed")
    require(capped_audit["partial_support_used"] is False, "cap trap retained a partial support")
    return {
        "passed": True,
        "brute_force_joint_target_count": len(brute),
        "positive_target": positive_target.tolist(),
        "negative_target": negative_target.tolist(),
        "positive_target_found": True,
        "negative_target_rigorously_absent": True,
        "same_row_cross_prime_false_positive_trap_rejected": True,
        "full_group_saturation_is_exact_and_only_implies_target_presence": True,
        "state_cap_is_explicit_skip_and_partial_support_is_discarded": True,
        "left_support": left_audit,
        "right_support": right_audit,
        "positive_join": positive,
        "negative_join": negative,
        "same_row_trap_join": trap,
        "saturation_audit": saturated_audit,
        "cap_audit": capped_audit,
    }


def run(
    *,
    smoke: bool,
    parent_path: Path,
    current_join_path: Path,
    hilbert_basis_path: Path,
    output_path: Path,
    smoke_cases: int,
    coordinate_count: int,
    projection_count: int,
    explicit_projections: tuple[tuple[int, ...], ...] | None,
    state_cap: int,
    pair_chunk_cap: int,
) -> dict:
    started = time.time()
    require(smoke, "this first case-level implementation deliberately supports --smoke only")
    require(1 <= state_cap <= MAX_SMOKE_STATE_CAP, "smoke state cap is outside its audited bound")
    require(1 <= pair_chunk_cap <= MAX_SMOKE_STATE_CAP, "pair chunk cap is outside its audited bound")

    parent_payload, representatives, current_by_key, targets, input_audit = load_current_problem(
        parent_path, current_join_path
    )
    rebuilt, reconstruction_audit = reconstruct_and_validate(parent_payload, representatives)
    target_grade_audit = audit_target_grades(targets, current_by_key, rebuilt)
    selected_cases, smoke_selection = select_smoke_cases(targets, current_by_key, smoke_cases)

    system = rebuilt["systems"][0]["A"]
    matrix = np.ascontiguousarray(system["matrix"], dtype=np.int64)
    kernel_rows = np.ascontiguousarray(rebuilt["kernel_rows"], dtype=np.int64)
    basis, degrees, basis_audit = semigroup.load_hilbert_basis(
        hilbert_basis_path, kernel_rows
    )
    common, common_audit = torsion.exact_common_dependency_basis(matrix, kernel_rows)
    quotient_data, quotient_audit = semigroup.derive_torsion_quotients(matrix, common)
    require(
        len(quotient_data[3]["complement"]) == 6
        and len(quotient_data[7]["complement"]) == 21,
        "effective torsion dimensions changed",
    )

    orbit = rebuilt["orbits"][0]
    full_rows, full_projection_audit = full_generator_projections(
        basis, orbit, quotient_data
    )
    projections, projection_selection = select_projection_subsets(
        full_rows=full_rows,
        degrees=degrees,
        coordinate_count=coordinate_count,
        projection_count=projection_count,
        explicit=explicit_projections,
    )
    require(
        all(len(row) <= MAX_SMOKE_MOD7_COORDINATES for row in projections),
        "selected projection is too wide for smoke mode",
    )

    support_tables = {}
    support_audits = []
    for projection in projections:
        table, audit = projected_direction_supports(
            projection=projection,
            full_generator_rows=full_rows,
            basis=basis,
            degrees=degrees,
            orbit=orbit,
            anchors=rebuilt["anchors"],
            quotient_data=quotient_data,
            state_cap=state_cap,
            pair_chunk_cap=pair_chunk_cap,
        )
        support_tables[projection] = table
        support_audits.append(audit)

    case_results = []
    for target_row in selected_cases:
        key = str(target_row["case_key"])
        projection_results = [
            evaluate_projected_case(
                target_row=target_row,
                current_row=current_by_key[key],
                rebuilt=rebuilt,
                common=common,
                quotient_data=quotient_data,
                projection=projection,
                supports=support_tables[projection],
                state_cap=state_cap,
                pair_chunk_cap=pair_chunk_cap,
            )
            for projection in projections
        ]
        rejected = any(row["rigorously_rejected"] for row in projection_results)
        skipped = not rejected and any(row["skipped"] for row in projection_results)
        necessary = not rejected and not skipped
        require(sum((rejected, skipped, necessary)) == 1, "aggregate case decision is ambiguous")
        row = {
            "case_key": key,
            "catalog_pattern": current_by_key[key]["catalog_pattern"],
            "prior_global_join_decision": current_by_key[key]["decision_status"],
            "projection_results": projection_results,
            "decision_status": (
                "rigorous_semigroup_case_projection_rejection"
                if rejected
                else "skipped_at_least_one_selected_projection"
                if skipped
                else "necessary_only_survivor_of_all_selected_projections"
            ),
            "rigorously_rejected": rejected,
            "necessary_only_survivor": necessary,
            "skipped": skipped,
            "all_selected_projections_completed": not any(
                result["skipped"] for result in projection_results
            ),
            "selected_projection_count": len(projection_results),
        }
        row["decision_certificate_sha256"] = json_sha256(row)
        case_results.append(row)

    counts = {
        "selected": len(case_results),
        "rejected": sum(row["rigorously_rejected"] for row in case_results),
        "surviving": sum(row["necessary_only_survivor"] for row in case_results),
        "skipped": sum(row["skipped"] for row in case_results),
    }
    require(
        counts["rejected"] + counts["surviving"] + counts["skipped"] == counts["selected"],
        "result count partition failed",
    )
    script_path = Path(__file__).resolve()
    result = {
        "experiment": "p7_infinity7_positive_z7_semigroup_case_join",
        "status": "bounded_smoke_semigroup_case_join_only",
        "p": 7,
        "c_H": 1,
        "z": 7,
        "phase": 0,
        "smoke_test": True,
        "full_run": False,
        "source_provenance": {
            "this_script_path": str(script_path),
            "this_script_sha256": file_sha256(script_path),
            "high_semigroup_support_script_sha256": file_sha256(Path(semigroup.__file__)),
            "torsion_projection_script_sha256": file_sha256(Path(torsion.__file__)),
            "global_catalog_join_script_sha256": file_sha256(Path(old_join.__file__)),
        },
        "configuration": {
            "smoke_cases": smoke_cases,
            "automatic_mod7_coordinate_count": coordinate_count,
            "automatic_projection_count": projection_count,
            "explicit_projections": (
                [list(row) for row in explicit_projections]
                if explicit_projections is not None
                else None
            ),
            "state_cap": state_cap,
            "pair_chunk_cap": pair_chunk_cap,
        },
        "input_and_current_decision_audit": input_audit,
        "representative_reconstruction_audit": reconstruction_audit,
        "target_grade_audit": target_grade_audit,
        "smoke_case_selection": smoke_selection,
        "normaliz_Hilbert_basis": basis_audit,
        "exact_common_rational_dependency_audit": common_audit,
        "varying_torsion_quotient_audit": quotient_audit,
        "full_generator_projection_audit": full_projection_audit,
        "projection_subset_selection": projection_selection,
        "projected_direction_support_audits": support_audits,
        "manufactured_semigroup_support_audit": semigroup.manufactured_self_audit(),
        "manufactured_case_join_audit": manufactured_case_join_audit(),
        "result_counts": counts,
        "case_results_sha256": old_join.canonical_case_digest(case_results),
        "case_results": case_results,
        "logical_semantics": {
            "all_six_mod3_torsion_coordinates_retained_in_every_projection": True,
            "selected_mod7_coordinates_are_from_the_derived_twenty_one_dimensional_quotient": True,
            "same_Hilbert_or_catalog_row_supplies_mod3_and_mod7_before_deduplication": True,
            "grades_zero_through_three_are_exact_box_bounded_catalog_supports": True,
            "all_eight_direction_supports_are_convolved_in_each_completed_case_projection": True,
            "missing_target_in_completed_projection_is_rigorous_rejection": True,
            "target_presence_is_only_necessary": True,
            "exact_full_group_saturation_is_only_target_presence_not_feasibility": True,
            "state_cap_hit_is_explicit_skip": True,
            "partial_support_after_cap_is_discarded": True,
            "different_projection_witnesses_are_not_assumed_compatible": True,
            "smoke_run_claims_full_51_case_coverage": False,
            "smoke_run_claims_positive_z7_closure": False,
            "binary_edge_feasibility_claimed": False,
        },
        "all_51_target_cases_processed": False,
        "positive_z7_excluded": False,
        "output_path": str(output_path.resolve()),
        "elapsed_seconds": time.time() - started,
    }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--parent-input", type=Path, default=DEFAULT_PARENT_INPUT)
    parser.add_argument("--current-join", type=Path, default=DEFAULT_CURRENT_JOIN)
    parser.add_argument("--hilbert-basis", type=Path, default=DEFAULT_HILBERT_BASIS)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smoke-cases", type=int, default=DEFAULT_SMOKE_CASES)
    parser.add_argument(
        "--mod7-coordinate-count", type=int, default=DEFAULT_MOD7_COORDINATE_COUNT
    )
    parser.add_argument("--projection-count", type=int, default=DEFAULT_PROJECTION_COUNT)
    parser.add_argument(
        "--mod7-projections",
        help="explicit semicolon-separated coordinate subsets, e.g. '0,1,2;3,4,5'",
    )
    parser.add_argument("--state-cap", type=int, default=DEFAULT_STATE_CAP)
    parser.add_argument("--pair-chunk-cap", type=int, default=DEFAULT_PAIR_CHUNK_CAP)
    args = parser.parse_args()
    explicit = parse_projection_spec(args.mod7_projections) if args.mod7_projections else None
    result = run(
        smoke=args.smoke,
        parent_path=args.parent_input,
        current_join_path=args.current_join,
        hilbert_basis_path=args.hilbert_basis,
        output_path=args.output,
        smoke_cases=args.smoke_cases,
        coordinate_count=args.mod7_coordinate_count,
        projection_count=args.projection_count,
        explicit_projections=explicit,
        state_cap=args.state_cap,
        pair_chunk_cap=args.pair_chunk_cap,
    )
    torsion.pointed.atomic_write(args.output, result)
    print(
        json.dumps(
            {
                "status": result["status"],
                "output": str(args.output),
                "target_case_count": result["input_and_current_decision_audit"][
                    "target_case_count"
                ],
                "smoke_selected_cases": result["result_counts"]["selected"],
                "rigorously_rejected": result["result_counts"]["rejected"],
                "necessary_only_survivors": result["result_counts"]["surviving"],
                "skipped": result["result_counts"]["skipped"],
                "selected_projections": [
                    row["mod7_coordinates"]
                    for row in result["projection_subset_selection"]["selected"]
                ],
                "elapsed_seconds": result["elapsed_seconds"],
            },
            indent=2,
            sort_keys=True,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
