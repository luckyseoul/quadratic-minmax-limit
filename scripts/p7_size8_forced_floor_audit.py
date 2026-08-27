#!/usr/bin/env python3
"""Independent audit of the p=7 size-eight forced-floor CUDA exclusion.

The CUDA result selects the finite boundaries whose two quadratic direction
types both spend their exact budget 32 at the parity floor.  This audit does
not import the CUDA implementation.  It rebuilds the affine score matrix and
its complete mod-seven left kernel, reconstructs every projected survivor's
exact Johnson-slice catalogs, and verifies that none satisfies all 135
dependencies.  It also checks, directly on the eight affine direction
partitions, that the nonsquare Paley anti-isometry bijects the stratum between
the two product signs.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
    scaled_direction_floor,
)
from e1_gmin_m4_prop15654 import p7_nonsquare_signed_permutation  # noqa: E402
from p7_fixed_boundary_catalog_join import mapped_catalog  # noqa: E402
from p7_unsaturated_modular_catalog_filter import (  # noqa: E402
    equation_matrix,
    left_dependencies,
)


ALL_BOUNDARIES = math.comb(49, 8)
EXPECTED_FORCED = 83_770_008
EXPECTED_PROFILES = 2_016
EXPECTED_ODD_HISTOGRAM = {
    16: 254_016,
    20: 8_396_640,
    24: 32_673_984,
    28: 30_465_456,
    32: 10_459_344,
    36: 1_467_648,
    40: 52_920,
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"expected JSON object: {path}")
    return payload


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def array_sha256(value: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(value).tobytes()).hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def modular_rank(matrix: np.ndarray, modulus: int = 7) -> int:
    value = np.asarray(matrix, dtype=np.int64).copy() % modulus
    rank = 0
    for column in range(value.shape[1]):
        candidates = np.flatnonzero(value[rank:, column])
        if not len(candidates):
            continue
        pivot = rank + int(candidates[0])
        if pivot != rank:
            value[[rank, pivot]] = value[[pivot, rank]]
        value[rank] = value[rank] * pow(int(value[rank, column]), -1, modulus) % modulus
        factors = value[:, column].copy()
        factors[rank] = 0
        active = np.flatnonzero(factors)
        if len(active):
            value[active] = (
                value[active] - factors[active, None] * value[rank, None]
            ) % modulus
        rank += 1
        if rank == value.shape[0]:
            break
    return rank


def direction_tables() -> tuple[np.ndarray, tuple[int, ...]]:
    epsilons = []
    labels = []
    for direction in projective_directions(7):
        eps, row = field_direction_data(7, direction)
        epsilons.append(int(eps))
        labels.append(tuple(int(value) for value in row))
    return np.asarray(labels, dtype=np.int8), tuple(epsilons)


def floor_for(c_h: int, eps: int, b: int) -> int:
    return int(scaled_direction_floor(7, b, int(eps == c_h)))


def type_costs(
    profile: tuple[int, ...], epsilons: tuple[int, ...], c_h: int
) -> tuple[int, int]:
    return tuple(
        sum(
            floor_for(c_h, epsilons[index], b)
            for index, b in enumerate(profile)
            if epsilons[index] == eps
        )
        for eps in (-1, 1)
    )


def source_scope(source: dict, c_h: int, epsilons: tuple[int, ...]) -> dict:
    require(
        source.get("experiment") == "p7_size8_floor_profile_gpu"
        and source.get("status") == "complete_exact_floor_profile_census"
        and int(source.get("p", 0)) == 7
        and int(source.get("c_H", 0)) == c_h
        and int(source.get("checked_boundaries", 0)) == ALL_BOUNDARIES,
        f"invalid complete floor source for c_H={c_h}",
    )
    profiles = 0
    boundaries = 0
    histogram: Counter[int] = Counter()
    for row in source["survivor_ordered_profiles"]:
        profile = tuple(int(value) for value in row["b_by_direction"])
        count = int(row["count"])
        if type_costs(profile, epsilons, c_h) == (32, 32):
            profiles += 1
            boundaries += count
            histogram[sum(profile)] += count
    out = {
        "ordered_profile_count": profiles,
        "boundary_count": boundaries,
        "odd_secant_histogram": dict(sorted(histogram.items())),
    }
    require(
        out
        == {
            "ordered_profile_count": EXPECTED_PROFILES,
            "boundary_count": EXPECTED_FORCED,
            "odd_secant_histogram": EXPECTED_ODD_HISTOGRAM,
        },
        f"forced-floor source scope changed for c_H={c_h}",
    )
    return out


def unrank_lex(rank: int) -> tuple[int, ...]:
    require(0 <= rank < ALL_BOUNDARIES, "rank outside C(49,8)")
    out = []
    next_value = 0
    for position in range(8):
        remaining = 7 - position
        for candidate in range(next_value, 49 - remaining):
            ways = math.comb(48 - candidate, remaining)
            if rank < ways:
                out.append(candidate)
                next_value = candidate + 1
                break
            rank -= ways
        else:
            raise AssertionError("lexicographic unranking failed")
    return tuple(out)


def boundary_masks(boundary: tuple[int, ...], labels: np.ndarray) -> tuple[int, ...]:
    masks = []
    for row in labels:
        mask = 0
        for point in boundary:
            mask ^= 1 << int(row[point])
        masks.append(mask)
    return tuple(masks)


def contribution(
    dependency: np.ndarray,
    direction_index: int,
    mask: int,
    eps: int,
    c_h: int,
) -> np.ndarray:
    b = mask.bit_count()
    phase = int(eps == c_h)
    mean = floor_for(c_h, eps, b)
    odd_fibres = {value for value in range(7) if mask & (1 << value)}
    values = mapped_catalog(b, phase, mean, odd_fibres, None).astype(np.int64)
    bad = 13 - values
    block = dependency[
        :, 2 + 35 * direction_index : 2 + 35 * (direction_index + 1)
    ].astype(np.int64)
    return (block @ (bad.T % 7) % 7).astype(np.uint8)


def dependency_matches(
    boundary: tuple[int, ...],
    labels: np.ndarray,
    epsilons: tuple[int, ...],
    c_h: int,
    dependency: np.ndarray,
) -> tuple[int, dict]:
    masks = boundary_masks(boundary, labels)
    profile = tuple(mask.bit_count() for mask in masks)
    require(type_costs(profile, epsilons, c_h) == (32, 32), "candidate left forced stratum")
    syndrome = (
        dependency[:, :2].astype(np.int64)
        @ np.asarray([29, 1], dtype=np.int64)
        % 7
    )
    variable = None
    variable_metadata = None
    catalog_sizes = []
    for direction_index, (eps, mask) in enumerate(zip(epsilons, masks)):
        values = contribution(dependency, direction_index, mask, eps, c_h).astype(np.int64)
        catalog_sizes.append(int(values.shape[1]))
        if values.shape[1] == 1:
            syndrome = (syndrome + values[:, 0]) % 7
        else:
            require(variable is None, "forced stratum has more than one variable catalog")
            variable = values
            variable_metadata = {
                "direction_index": direction_index,
                "odd_fibre_mask": mask,
                "catalog_rows": int(values.shape[1]),
            }
    if variable is None:
        matches = int(not np.any(syndrome % 7))
    else:
        matches = int(
            np.count_nonzero(
                np.all((variable + syndrome[:, None]) % 7 == 0, axis=0)
            )
        )
    return matches, {
        "odd_fibre_masks": list(masks),
        "odd_fibre_profile": list(profile),
        "catalog_sizes": catalog_sizes,
        "variable_catalog": variable_metadata,
    }


def audit_direction_transfer(
    labels: np.ndarray, epsilons: tuple[int, ...], finite_permutation: tuple[int, ...]
) -> dict:
    direction_map = []
    fibre_maps = []
    for old_direction in range(8):
        matches = []
        for new_direction in range(8):
            fibre_map: dict[int, int] = {}
            valid = True
            for point in range(49):
                old_label = int(labels[old_direction, point])
                new_label = int(labels[new_direction, finite_permutation[point]])
                previous = fibre_map.setdefault(old_label, new_label)
                if previous != new_label:
                    valid = False
                    break
            if valid and len(fibre_map) == 7 and len(set(fibre_map.values())) == 7:
                matches.append((new_direction, tuple(fibre_map[value] for value in range(7))))
        require(len(matches) == 1, "nonsquare multiplication did not uniquely map a direction")
        new_direction, fibre_map = matches[0]
        require(
            epsilons[new_direction] == -epsilons[old_direction],
            "nonsquare direction map did not reverse quadratic type",
        )
        direction_map.append(new_direction)
        fibre_maps.append(fibre_map)
    require(sorted(direction_map) == list(range(8)), "direction transfer is not bijective")
    return {
        "direction_map": direction_map,
        "fibre_label_permutations": [list(row) for row in fibre_maps],
        "direction_map_is_bijective": True,
        "quadratic_types_reverse": True,
        "odd_fibre_counts_are_preserved": True,
        "phase_int_eps_equals_cH_is_preserved_when_cH_flips": True,
        "type_floor_sums_swap": True,
        "forced_floor_32_32_stratum_is_bijective": True,
    }


def run(args: argparse.Namespace) -> dict:
    started = time.time()
    result = load_json(args.gpu_result)
    source_minus = load_json(args.floor_minus)
    source_plus = load_json(args.floor_plus)
    labels, epsilons = direction_tables()
    scope_minus = source_scope(source_minus, -1, epsilons)
    scope_plus = source_scope(source_plus, 1, epsilons)
    require(scope_minus == scope_plus, "the two forced-floor source scopes differ")

    require(
        result.get("experiment") == "p7_size8_forced_floor_gpu"
        and result.get("status") == "complete_exact_doubly_saturated_boundary_exhaustion"
        and int(result.get("p", 0)) == 7
        and int(result.get("c_H", 0)) == -1
        and int(result.get("checked_boundaries", 0)) == ALL_BOUNDARIES
        and int(result.get("forced_floor_boundaries", 0)) == EXPECTED_FORCED,
        "invalid complete CUDA result",
    )
    require(result.get("source_sha256") == sha256(args.floor_minus), "CUDA floor-source hash changed")
    require(result.get("source_scope") == {
        "ordered_profile_count": EXPECTED_PROFILES,
        "boundary_count": EXPECTED_FORCED,
        "odd_secant_histogram": {str(k): v for k, v in EXPECTED_ODD_HISTOGRAM.items()},
    }, "CUDA source scope changed")
    require(
        {int(k): int(v) for k, v in result["forced_floor_odd_secant_histogram"].items()}
        == EXPECTED_ODD_HISTOGRAM,
        "CUDA odd-secant histogram changed",
    )

    matrix = equation_matrix()
    rank, dependency = left_dependencies(matrix, 7)
    dependency = dependency.astype(np.uint8)
    require(rank == 147 and dependency.shape == (135, 282), "mod-seven system changed")
    require(
        not np.any(dependency.astype(np.int64) @ (matrix.astype(np.int64) % 7) % 7),
        "reconstructed dependency basis is not left-null",
    )
    require(
        array_sha256(dependency) == result["full_dependency_sha256"],
        "reconstructed dependency basis hash changed",
    )
    coefficients = np.asarray(result["projection"]["coefficients_mod7"], dtype=np.uint8)
    projected = coefficients.astype(np.int64) @ dependency.astype(np.int64) % 7
    require(
        coefficients.shape == (8, 135)
        and modular_rank(coefficients) == 8
        and modular_rank(projected) == 8
        and array_sha256(projected.astype(np.uint8))
        == result["projection"]["projected_dependency_sha256"],
        "projected dependency reconstruction changed",
    )

    ranks = [int(value) for value in result["projected_survivor_ranks"]]
    require(
        ranks == sorted(set(ranks))
        and len(ranks) == int(result["projected_dependency_survivors"]),
        "projected survivor ranks are incomplete or duplicated",
    )
    evidence_rows = {
        int(row["rank"]): row for row in result["full_135_dependency_recheck_rows"]
    }
    require(set(evidence_rows) == set(ranks), "full recheck rows do not cover projected ranks")

    exact_survivors = []
    candidate_digest = hashlib.sha256()
    catalog_histogram: Counter[int] = Counter()
    for candidate_rank in ranks:
        boundary = unrank_lex(candidate_rank)
        projected_matches, _ = dependency_matches(
            boundary, labels, epsilons, -1, projected.astype(np.uint8)
        )
        require(projected_matches > 0, "recorded rank does not pass its projected dependencies")
        full_matches, metadata = dependency_matches(
            boundary, labels, epsilons, -1, dependency
        )
        if full_matches:
            exact_survivors.append(candidate_rank)
        row = evidence_rows[candidate_rank]
        require(
            row["boundary_finite_field"] == list(boundary)
            and row["boundary_vertices"] == [value + 1 for value in boundary]
            and row["odd_fibre_masks"] == metadata["odd_fibre_masks"]
            and row["variable_catalog"] == metadata["variable_catalog"]
            and int(row["full_dependency_matches"]) == full_matches,
            f"candidate recheck row changed at rank {candidate_rank}",
        )
        catalog_histogram.update(metadata["catalog_sizes"])
        candidate_digest.update(np.asarray([candidate_rank], dtype=np.uint64).tobytes())
        candidate_digest.update(np.asarray(metadata["odd_fibre_masks"], dtype=np.uint8).tobytes())
        candidate_digest.update(np.asarray([full_matches], dtype=np.uint16).tobytes())
    require(not exact_survivors, "a full mod-seven dependency survivor remains")
    require(
        int(result.get("full_dependency_survivors", -1)) == 0
        and result.get("full_dependency_survivor_ranks") == []
        and result.get("all_forced_floor_boundaries_mod7_infeasible") is True,
        "CUDA result's exact-closure flags changed",
    )

    symmetry = p7_nonsquare_signed_permutation()
    require(
        symmetry["fixes_distinguished_edge"] is True
        and symmetry["signed_conference_anti_isometry"] is True,
        "Paley nonsquare anti-isometry failed",
    )
    finite_permutation = tuple(int(value) for value in symmetry["finite_permutation"])
    transfer = audit_direction_transfer(labels, epsilons, finite_permutation)

    return {
        "experiment": "p7_size8_forced_floor_audit",
        "status": "passed_independent_complete_forced_floor_audit",
        "p": 7,
        "finite_boundary_size": 8,
        "gpu_result": str(args.gpu_result),
        "gpu_result_sha256": sha256(args.gpu_result),
        "floor_sources": {
            "c_H=-1": str(args.floor_minus),
            "c_H=+1": str(args.floor_plus),
        },
        "floor_source_sha256": {
            "c_H=-1": sha256(args.floor_minus),
            "c_H=+1": sha256(args.floor_plus),
        },
        "source_scope_each_sign": scope_minus,
        "all_boundaries_checked_by_cuda": ALL_BOUNDARIES,
        "forced_floor_boundaries_each_sign": EXPECTED_FORCED,
        "full_linear_system_shape": list(matrix.shape),
        "full_linear_system_rank_mod7": rank,
        "full_dependency_shape": list(dependency.shape),
        "full_dependency_sha256": array_sha256(dependency),
        "projected_dependency_rank_mod7": modular_rank(projected),
        "projected_survivors_reconstructed": len(ranks),
        "projected_candidate_audit_sha256": candidate_digest.hexdigest(),
        "candidate_catalog_size_histogram": dict(sorted(catalog_histogram.items())),
        "full_135_dependency_survivors": len(exact_survivors),
        "all_cminus_forced_floor_boundaries_excluded": not exact_survivors,
        "nonsquare_anti_isometry": {
            "nonsquare_multiplier": int(symmetry["nonsquare_multiplier"]),
            "fixes_distinguished_edge": True,
            "signed_conference_anti_isometry": True,
        },
        "sign_transfer": transfer,
        "all_cplus_forced_floor_boundaries_excluded": not exact_survivors,
        "all_forced_floor_boundaries_both_signs_excluded": not exact_survivors,
        "closes_all_nonconic_size_eight": False,
        "closes_all_p7_size_eight": False,
        "closes_residual_ii": False,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gpu-result", type=Path, required=True)
    parser.add_argument("--floor-minus", type=Path, required=True)
    parser.add_argument("--floor-plus", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = run(args)
    atomic_json(args.output, out)
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
