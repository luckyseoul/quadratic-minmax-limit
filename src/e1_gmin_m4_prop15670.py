#!/usr/bin/env python3
"""Prop. 15.670 -- every finite p=11 size-eight boundary is impossible.

An ordered pair in any finite eight-set can be sent to field points 0 and 1
by a unique affine map.  A square normalizing scalar preserves quadratic
direction types; a nonsquare scalar swaps the types, equivalently swapping
the two possible values of c_H in the parity-floor test.  It therefore
suffices to test both signs on the C(119,6) sets containing 0 and 1.

Independent CUDA (V100) and HIP (RX 9070 XT) censuses test all 3,470,108,187
normalized sets.  Both complete cost-pair histograms agree exactly, both
signs have zero floor survivors, and the minimum larger type cost is 76,
four above the exact type budget 72.  Independent CPU itertools code agrees
on every histogram entry in a 100,000-set prefix.

This closes the first p=11 floor-plus-pair survivor from Proposition 15.669.
Infinity plus nine, larger boundaries, residual (ii), R1, QVAR, Type I, and
the limit remain open.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15632 import (
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)


ROOT = Path(__file__).resolve().parents[1]
P = 11
Q = P * P
NORMALIZED_COUNT = math.comb(Q - 2, 6)
FULL_COUNT = math.comb(Q, 8)
TYPE_BUDGET = (P + 1) ** 2 // 2
COST_STRIDE = 133
V100_PATH = ROOT / "evidence" / "p11_size8_normalized_floor_v100.json"
RX9070XT_PATH = (
    ROOT / "evidence" / "p11_size8_normalized_floor_rx9070xt.json"
)
PINNED_FILE_SHA256 = {
    "v100": "3ba6bd0574b6dd48087647be9380fb5775c63b65380da3bf4ef3c7368b716a5a",
    "rx9070xt": "9b2cef9f9558cef9f0612383fa6673aee56bd96e0f5303322f4e4e9a28fca16a",
}
PINNED_HISTOGRAM_SHA256 = {
    -1: "106de27687d372ce083a3fcfb1adb3fcf50830f14cd314842eebacd9488fbf01",
    1: "8df9472aa3d7bb1c24db7240a967d4a861b6e63343e07ea48adf07a473958e39",
}


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_record(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())


def histogram_array(sign_row: dict[str, object]) -> np.ndarray:
    histogram = np.zeros((COST_STRIDE, COST_STRIDE), dtype="<u8")
    for row in sign_row["cost_pair_histogram"]:
        histogram[int(row["cost_minus"]), int(row["cost_plus"])] = int(
            row["count"]
        )
    return histogram


def direction_tables() -> tuple[np.ndarray, np.ndarray]:
    labels = []
    epsilons = []
    for direction in projective_directions(P):
        epsilon, row = field_direction_data(P, direction)
        epsilons.append(int(epsilon))
        labels.append([int(value) for value in row])
    return np.asarray(labels, dtype=np.int8), np.asarray(epsilons, dtype=np.int8)


def audit_boundary(boundary: list[int], c_h: int) -> dict[str, object]:
    labels, epsilons = direction_tables()
    totals = {-1: 0, 1: 0}
    b_values = []
    for direction in range(P + 1):
        mask = 0
        for point in boundary:
            mask ^= 1 << int(labels[direction, point])
        b = mask.bit_count()
        phase = int(int(epsilons[direction]) == c_h)
        cost = scaled_direction_floor(P, b, phase)
        totals[int(epsilons[direction])] += cost
        b_values.append(b)
    return {
        "type_costs": {str(key): value for key, value in totals.items()},
        "b_by_direction": b_values,
        "maximum_type_cost": max(totals.values()),
    }


def sign_semantic_audit(sign_row: dict[str, object]) -> dict[str, object]:
    c_h = int(sign_row["c_H"])
    histogram = histogram_array(sign_row)
    semantic_hash = hashlib.sha256(histogram.tobytes()).hexdigest()
    support = np.argwhere(histogram != 0)
    minimum_maximum_cost = min(max(int(a), int(b)) for a, b in support)
    floor_survivors = sum(
        int(histogram[a, b])
        for a in range(TYPE_BUDGET + 1)
        for b in range(TYPE_BUDGET + 1)
    )
    representative = audit_boundary(
        [int(value) for value in sign_row["first_minimum_boundary"]], c_h
    )
    expected_costs = {
        str(key): int(value)
        for key, value in sign_row["first_minimum_type_costs"].items()
    }
    valid = bool(
        int(histogram.sum()) == NORMALIZED_COUNT
        and semantic_hash == str(sign_row["histogram_sha256"])
        and semantic_hash == PINNED_HISTOGRAM_SHA256[c_h]
        and floor_survivors == int(sign_row["floor_survivors"]) == 0
        and minimum_maximum_cost
        == int(sign_row["minimum_maximum_type_cost"])
        == 76
        and representative["type_costs"] == expected_costs
        and representative["b_by_direction"]
        == [int(value) for value in sign_row["first_minimum_b_by_direction"]]
        and representative["maximum_type_cost"] == 76
    )
    return {
        "c_H": c_h,
        "histogram_total": int(histogram.sum()),
        "histogram_sha256": semantic_hash,
        "floor_survivors": floor_survivors,
        "minimum_maximum_type_cost": minimum_maximum_cost,
        "minimum_budget_excess": minimum_maximum_cost - TYPE_BUDGET,
        "representative_recomputed": representative,
        "valid": valid,
    }


def _partition_signature(labels: list[int]) -> tuple[tuple[int, ...], ...]:
    fibres = [[] for _ in range(P)]
    for point, label in enumerate(labels):
        fibres[int(label)].append(point)
    return tuple(sorted(tuple(fibre) for fibre in fibres))


def affine_similarity_audit() -> dict[str, object]:
    """Audit the normalization group on every point and direction.

    Multiplication by ``alpha`` permutes the parallel classes and changes
    every kernel type by ``chi(alpha)``.  If the latter is negative, replacing
    ``c_H`` by ``-c_H`` preserves the phase test ``eps_d == c_H``.  Translation
    only relabels the fibres inside each fixed direction.
    """
    from e1_gmin_m4_prop15598 import field_ctx

    q, multiply, add, character, *_rest = field_ctx(P)
    directions = projective_directions(P)
    direction_rows = [field_direction_data(P, direction) for direction in directions]
    signature_to_direction = {
        _partition_signature(labels): index
        for index, (_epsilon, labels) in enumerate(direction_rows)
    }
    if len(signature_to_direction) != P + 1:
        raise AssertionError("projective direction partitions are not distinct")

    similarity_pairs = 0
    phase_transfer_cases = 0
    for alpha in range(1, q):
        alpha_type = int(character(alpha))
        for epsilon, labels in direction_rows:
            transformed = [0] * q
            for point, label in enumerate(labels):
                transformed[multiply(alpha, point)] = int(label)
            target_index = signature_to_direction[_partition_signature(transformed)]
            target_epsilon = int(direction_rows[target_index][0])
            if target_epsilon != alpha_type * int(epsilon):
                raise AssertionError("quadratic direction type did not transfer")
            for c_h in (-1, 1):
                target_c_h = alpha_type * c_h
                if (int(epsilon) == c_h) != (target_epsilon == target_c_h):
                    raise AssertionError("finite-boundary phase did not transfer")
                phase_transfer_cases += 1
            similarity_pairs += 1

    translation_pairs = 0
    for beta in range(q):
        for _epsilon, labels in direction_rows:
            transformed = [0] * q
            for point, label in enumerate(labels):
                transformed[add(point, beta)] = int(label)
            if _partition_signature(transformed) != _partition_signature(labels):
                raise AssertionError("translation changed a parallel class")
            translation_pairs += 1

    return {
        "all_nonzero_scalars": q - 1,
        "all_translations": q,
        "all_projective_directions": P + 1,
        "similarity_direction_pairs_checked": similarity_pairs,
        "translation_direction_pairs_checked": translation_pairs,
        "phase_transfer_cases_checked": phase_transfer_cases,
        "direction_type_multiplies_by_scalar_character": True,
        "finite_boundary_phase_is_preserved_after_c_H_transfer": True,
        "valid": True,
    }


def normalization_ledger() -> dict[str, object]:
    ordered_pairs_per_set = 8 * 7
    affine_maps = Q * (Q - 1)
    left = FULL_COUNT * ordered_pairs_per_set
    right = NORMALIZED_COUNT * affine_maps
    similarity = affine_similarity_audit()
    return {
        "all_finite_eight_sets": FULL_COUNT,
        "ordered_pairs_per_set": ordered_pairs_per_set,
        "normalized_sets_containing_0_and_1": NORMALIZED_COUNT,
        "affine_maps_x_to_a_x_plus_b": affine_maps,
        "pointed_set_count_before_normalization": left,
        "pointed_set_count_after_normalization": right,
        "counting_identity": left == right,
        "unique_map_for_ordered_pair": "z -> (z-x)/(y-x)",
        "similarity_audit": similarity,
        "nonsquare_scalar_swaps_direction_types": similarity[
            "direction_type_multiplies_by_scalar_character"
        ],
        "both_c_H_signs_cover_the_swap": similarity[
            "finite_boundary_phase_is_preserved_after_c_H_transfer"
        ],
    }


def theorem_record() -> dict[str, object]:
    records = {
        "v100": load_record(V100_PATH),
        "rx9070xt": load_record(RX9070XT_PATH),
    }
    file_hashes = {
        "v100": file_sha256(V100_PATH),
        "rx9070xt": file_sha256(RX9070XT_PATH),
    }
    hashes_pinned = file_hashes == PINNED_FILE_SHA256
    core_fields_match = all(
        records["v100"][key] == records["rx9070xt"][key]
        for key in (
            "normalized_boundaries",
            "checked_boundaries",
            "type_budget",
            "signs",
            "prefix_verification",
            "all_p11_eight_finite_boundaries_excluded",
        )
    )
    sign_audits = [sign_semantic_audit(row) for row in records["v100"]["signs"]]
    normalization = normalization_ledger()
    devices = {
        name: {
            "device": record["device"],
            "elapsed_seconds": record["elapsed_seconds"],
            "boundaries_per_second": record["boundaries_per_second"],
        }
        for name, record in records.items()
    }
    proved = bool(
        hashes_pinned
        and core_fields_match
        and normalization["counting_identity"]
        and normalization["similarity_audit"]["valid"]
        and all(row["valid"] for row in sign_audits)
        and all(
            int(record["checked_boundaries"]) == NORMALIZED_COUNT
            and record["prefix_verification"][
                "independent_itertools_all_histograms_match"
            ]
            is True
            and int(record["prefix_verification"]["checked"]) == 100_000
            for record in records.values()
        )
    )
    return {
        "prop": "15.670",
        "title": "Finite p=11 size-eight boundary exclusion",
        "proved": proved,
        "normalization": normalization,
        "file_sha256": file_hashes,
        "pinned_file_sha256_match": hashes_pinned,
        "independent_gpu_core_fields_match": core_fields_match,
        "sign_audits": sign_audits,
        "benchmarks_not_used_as_proof": devices,
        "theorem": {
            "every_finite_p11_size8_boundary": "IMPOSSIBLE",
            "exact_minimum_larger_type_cost": 76,
            "exact_type_budget": TYPE_BUDGET,
            "contradiction_gap": 4,
            "p11_infinity_plus_9": "OPEN",
            "larger_boundaries": "OPEN",
            "general_residual_ii": False,
            "R1": False,
            "global_QVAR": False,
            "type_I": False,
            "limit_exists": False,
        },
        "L_status": "OPEN",
    }


def main() -> dict[str, object]:
    record = theorem_record()
    if record["proved"] is not True:
        raise ArithmeticError("Proposition 15.670 audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15670.json"
    destination.write_text(json.dumps(record, indent=2) + "\n")
    print("Prop 15.670 finite p=11 size-eight boundary: proved")
    print(f"  wrote {destination}")
    return record


if __name__ == "__main__":
    main()
