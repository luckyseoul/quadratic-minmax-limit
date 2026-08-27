#!/usr/bin/env python3
"""Independent audit and mod-3 closure of the p=7 four-allocation stratum.

This script does not import the CUDA scanner.  It rebuilds the common score
matrix and its full left kernels over F_3 and F_7, replays every projected
rank/direction candidate, intersects exact catalog-row solutions across the
two characteristics, identifies the complete mod-seven survivor family as
an affine line plus one off-line point, and verifies the nonsquare
anti-isometry transfer to the opposite product sign.
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
from p7_fixed_boundary_mean_allocation_batch import POINTS, allocations  # noqa: E402
from p7_unsaturated_modular_catalog_filter import (  # noqa: E402
    equation_matrix,
    left_dependencies,
)


ALL_BOUNDARIES = math.comb(49, 8)
EXPECTED_ALLOCATION_BOUNDARIES = {
    4: 23_563_806,
    11: 154_056,
    16: 1_194_816,
    24: 1_176,
    44: 69_384,
}
EXPECTED_ALLOCATION_PROFILES = {4: 2_245, 11: 248, 16: 516, 24: 8, 44: 110}
EXPECTED_FOUR_ODD_HISTOGRAM = {
    16: 691_488,
    20: 5_603_640,
    24: 9_190_146,
    28: 5_990_544,
    32: 1_846_908,
    36: 232_848,
    40: 5_880,
    44: 2_352,
}
EXPECTED_PROJECTED_LEAVES = 1_191
EXPECTED_PROJECTED_BOUNDARIES = 1_177
EXPECTED_MOD7_SURVIVOR_LEAVES = 1_176


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
    labels = []
    epsilons = []
    for direction in projective_directions(7):
        eps, row = field_direction_data(7, direction)
        labels.append(tuple(int(value) for value in row))
        epsilons.append(int(eps))
    label_array = np.asarray(labels, dtype=np.int8)
    require(label_array.shape == (8, 49), "direction-label table changed")
    require(sorted(epsilons) == [-1] * 4 + [1] * 4, "direction types changed")
    return label_array, tuple(epsilons)


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


def abstract_rows(
    profile: tuple[int, ...], epsilons: tuple[int, ...], c_h: int
) -> list[dict]:
    rows = []
    for b, eps in zip(profile, epsilons):
        phase = int(eps == c_h)
        floor = floor_for(c_h, eps, b)
        odd_fibres = set(range(b))
        parity_mass = sum(
            (sum(value in odd_fibres for value in point) + phase) & 1
            for point in POINTS
        )
        allowed = tuple(
            mean
            for mean in range(floor, 33, 2)
            if 5 * mean >= 2 * parity_mass
            and (5 * mean - 2 * parity_mass) % 4 == 0
        )
        rows.append(
            {
                "eps": int(eps),
                "floor": floor,
                "allowed_scaled_means": allowed,
            }
        )
    return rows


def source_scope(source: dict, c_h: int, epsilons: tuple[int, ...]) -> dict:
    require(
        source.get("experiment") == "p7_size8_floor_profile_gpu"
        and source.get("status") == "complete_exact_floor_profile_census"
        and int(source.get("p", 0)) == 7
        and int(source.get("c_H", 0)) == c_h
        and int(source.get("checked_boundaries", 0)) == ALL_BOUNDARIES,
        f"invalid complete floor source for c_H={c_h}",
    )
    allocation_boundaries: Counter[int] = Counter()
    allocation_profiles: Counter[int] = Counter()
    four_floor_pairs: Counter[tuple[int, int]] = Counter()
    four_odd: Counter[int] = Counter()
    for source_row in source["survivor_ordered_profiles"]:
        profile = tuple(int(value) for value in source_row["b_by_direction"])
        count = int(source_row["count"])
        floors = type_costs(profile, epsilons, c_h)
        if sum(profile) == 8 or floors == (32, 32):
            continue
        leaves = allocations(abstract_rows(profile, epsilons, c_h))
        allocation_boundaries[len(leaves)] += count
        allocation_profiles[len(leaves)] += 1
        if len(leaves) == 4:
            deficient_eps = -1 if floors == (24, 32) else 1
            floor_vector = tuple(floor_for(c_h, eps, b) for eps, b in zip(epsilons, profile))
            elevated = []
            for leaf in leaves:
                difference = tuple(value - floor for value, floor in zip(leaf, floor_vector))
                support = [index for index, value in enumerate(difference) if value]
                require(
                    len(support) == 1
                    and difference[support[0]] == 8
                    and epsilons[support[0]] == deficient_eps,
                    "four-allocation leaf shape changed",
                )
                elevated.append(support[0])
            require(
                sorted(elevated)
                == sorted(index for index, eps in enumerate(epsilons) if eps == deficient_eps),
                "four-allocation elevated directions changed",
            )
            four_floor_pairs[floors] += count
            four_odd[sum(profile)] += count
    require(dict(allocation_boundaries) == EXPECTED_ALLOCATION_BOUNDARIES, "boundary census changed")
    require(dict(allocation_profiles) == EXPECTED_ALLOCATION_PROFILES, "profile census changed")
    require(dict(four_odd) == EXPECTED_FOUR_ODD_HISTOGRAM, "odd-secant census changed")
    expected_floor_pairs = (
        {(24, 32): 17_298_078, (32, 24): 6_265_728}
        if c_h == -1
        else {(24, 32): 6_265_728, (32, 24): 17_298_078}
    )
    require(dict(four_floor_pairs) == expected_floor_pairs, "floor-pair census changed")
    return {
        "allocation_count_boundary_histogram": dict(sorted(allocation_boundaries.items())),
        "allocation_count_ordered_profile_histogram": dict(sorted(allocation_profiles.items())),
        "four_allocation_boundaries": allocation_boundaries[4],
        "four_allocation_leaves": 4 * allocation_boundaries[4],
        "four_allocation_floor_pair_histogram": [
            {"type_floor_sums": list(key), "boundaries": value}
            for key, value in sorted(four_floor_pairs.items())
        ],
        "four_allocation_odd_secant_histogram": dict(sorted(four_odd.items())),
    }


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


def rank_lex(boundary: tuple[int, ...]) -> int:
    require(len(boundary) == 8 and tuple(sorted(set(boundary))) == boundary, "bad boundary")
    rank = 0
    next_value = 0
    for position, selected in enumerate(boundary):
        remaining = 7 - position
        for candidate in range(next_value, selected):
            rank += math.comb(48 - candidate, remaining)
        next_value = selected + 1
    require(unrank_lex(rank) == boundary, "rank/unrank audit failed")
    return rank


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
    modulus: int,
    direction: int,
    mask: int,
    eps: int,
    c_h: int,
    mean: int,
) -> np.ndarray:
    odd_fibres = {value for value in range(7) if mask & (1 << value)}
    values = mapped_catalog(
        mask.bit_count(), int(eps == c_h), mean, odd_fibres, None
    ).astype(np.int64)
    bad = 13 - values
    block = dependency[:, 2 + 35 * direction : 2 + 35 * (direction + 1)].astype(
        np.int64
    )
    return (block @ (bad.T % modulus) % modulus).astype(np.uint8)


def catalog_matches(
    syndrome: np.ndarray,
    variables: list[np.ndarray],
    modulus: int,
) -> set[tuple[int, ...]]:
    if not variables:
        return {()} if not np.any(syndrome % modulus) else set()
    if len(variables) == 1:
        hits = np.flatnonzero(
            np.all(
                (variables[0].astype(np.int64) + syndrome[:, None]) % modulus == 0,
                axis=0,
            )
        )
        return {(int(index),) for index in hits}
    require(len(variables) == 2, "more than two variable catalogs")
    lookup: dict[bytes, list[int]] = {}
    for second_index in range(variables[1].shape[1]):
        key = np.ascontiguousarray(variables[1][:, second_index]).tobytes()
        lookup.setdefault(key, []).append(second_index)
    matches = set()
    for first_index in range(variables[0].shape[1]):
        target = (
            -syndrome.astype(np.int64)
            - variables[0][:, first_index].astype(np.int64)
        ) % modulus
        key = np.ascontiguousarray(target.astype(np.uint8)).tobytes()
        for second_index in lookup.get(key, []):
            matches.add((first_index, second_index))
    return matches


def dependency_matches(
    boundary: tuple[int, ...],
    omitted: int,
    labels: np.ndarray,
    epsilons: tuple[int, ...],
    c_h: int,
    dependency: np.ndarray,
    modulus: int,
    cache: dict[tuple[int, int, int, int], np.ndarray],
) -> tuple[set[tuple[int, ...]], dict]:
    masks = boundary_masks(boundary, labels)
    profile = tuple(mask.bit_count() for mask in masks)
    floors = type_costs(profile, epsilons, c_h)
    require(floors in ((24, 32), (32, 24)), "candidate left four-allocation stratum")
    deficient_eps = -1 if floors == (24, 32) else 1
    require(epsilons[omitted] == deficient_eps, "candidate raises saturated type")
    means = [floor_for(c_h, eps, b) for eps, b in zip(epsilons, profile)]
    means[omitted] += 8
    syndrome = (
        dependency[:, :2].astype(np.int64)
        @ np.asarray([29, 1], dtype=np.int64)
        % modulus
    )
    variables = []
    variable_metadata = []
    catalog_sizes = []
    for direction, (mask, eps, mean) in enumerate(zip(masks, epsilons, means)):
        key = (direction, mask, mean, modulus)
        if key not in cache:
            cache[key] = contribution(
                dependency, modulus, direction, mask, eps, c_h, mean
            )
        values = cache[key].astype(np.int64)
        catalog_sizes.append(int(values.shape[1]))
        if values.shape[1] == 1:
            syndrome = (syndrome + values[:, 0]) % modulus
        else:
            variables.append(values)
            variable_metadata.append(
                {
                    "direction_index": direction,
                    "odd_fibre_mask": mask,
                    "scaled_mean": mean,
                    "catalog_rows": int(values.shape[1]),
                }
            )
    matches = catalog_matches(syndrome, variables, modulus)
    return matches, {
        "odd_fibre_masks": list(masks),
        "odd_fibre_profile": list(profile),
        "scaled_means_direction_order": means,
        "catalog_sizes": catalog_sizes,
        "variable_catalogs": variable_metadata,
    }


def projected_passes(
    boundary: tuple[int, ...],
    omitted: int,
    labels: np.ndarray,
    epsilons: tuple[int, ...],
    projected: np.ndarray,
    cache: dict[tuple[int, int, int, int], np.ndarray],
) -> bool:
    masks = boundary_masks(boundary, labels)
    syndrome = (
        projected[:, :2].astype(np.int64)
        @ np.asarray([29, 1], dtype=np.int64)
        % 7
    )
    variable = None
    for direction, (mask, eps) in enumerate(zip(masks, epsilons)):
        if direction == omitted:
            continue
        mean = floor_for(-1, eps, mask.bit_count())
        key = (omitted, direction, mask, mean)
        if key not in cache:
            cache[key] = contribution(projected, 7, direction, mask, eps, -1, mean)
        values = cache[key].astype(np.int64)
        if values.shape[1] == 1:
            syndrome = (syndrome + values[:, 0]) % 7
        else:
            require(variable is None, "projected leaf has two floor variables")
            variable = values
    if variable is None:
        return not np.any(syndrome % 7)
    return bool(np.any(np.all((variable + syndrome[:, None]) % 7 == 0, axis=0)))


def line_plus_point_family(
    labels: np.ndarray,
    epsilons: tuple[int, ...],
    line_type: int,
) -> set[tuple[int, int]]:
    family = set()
    for direction, eps in enumerate(epsilons):
        if eps != line_type:
            continue
        for fibre in range(7):
            line = tuple(point for point in range(49) if int(labels[direction, point]) == fibre)
            require(len(line) == 7, "affine direction fibre is not a line")
            line_set = set(line)
            for point in range(49):
                if point in line_set:
                    continue
                boundary = tuple(sorted(line + (point,)))
                masks = boundary_masks(boundary, labels)
                require(
                    masks[direction].bit_count() == 2
                    and all(
                        mask.bit_count() == 6
                        for index, mask in enumerate(masks)
                        if index != direction
                    ),
                    "line-plus-point odd-fibre profile changed",
                )
                family.add((rank_lex(boundary), direction))
    require(len(family) == 4 * 7 * 42, "line-plus-point family count changed")
    return family


def audit_direction_transfer(
    labels: np.ndarray,
    epsilons: tuple[int, ...],
    finite_permutation: tuple[int, ...],
) -> tuple[dict, tuple[int, ...]]:
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
        require(len(matches) == 1, "nonsquare multiplication did not map a unique direction")
        new_direction, fibre_map = matches[0]
        require(
            epsilons[new_direction] == -epsilons[old_direction],
            "nonsquare direction map did not reverse type",
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
        "four_allocation_stratum_is_bijective": True,
    }, tuple(direction_map)


def run(args: argparse.Namespace) -> dict:
    started = time.time()
    result = load_json(args.gpu_result)
    source_minus = load_json(args.floor_minus)
    source_plus = load_json(args.floor_plus)
    table_summary = load_json(args.table_summary)
    labels, epsilons = direction_tables()
    scope_minus = source_scope(source_minus, -1, epsilons)
    scope_plus = source_scope(source_plus, 1, epsilons)
    require(
        scope_minus["four_allocation_boundaries"]
        == scope_plus["four_allocation_boundaries"]
        == EXPECTED_ALLOCATION_BOUNDARIES[4],
        "two-sign four-allocation scopes differ",
    )

    require(
        result.get("experiment") == "p7_size8_one_elevation_gpu"
        and result.get("status") == "complete_exact_four_allocation_boundary_exhaustion"
        and int(result.get("p", 0)) == 7
        and int(result.get("c_H", 0)) == -1
        and int(result.get("checked_boundaries", 0)) == ALL_BOUNDARIES
        and int(result.get("four_allocation_boundaries", 0))
        == EXPECTED_ALLOCATION_BOUNDARIES[4]
        and int(result.get("four_allocation_leaves", 0))
        == 4 * EXPECTED_ALLOCATION_BOUNDARIES[4],
        "invalid complete CUDA result",
    )
    require(result.get("source_sha256") == sha256(args.floor_minus), "CUDA source hash changed")
    require(
        int(result.get("projected_dependency_survivor_leaves", -1))
        == EXPECTED_PROJECTED_LEAVES
        and int(result.get("projected_dependency_survivor_boundaries", -1))
        == EXPECTED_PROJECTED_BOUNDARIES
        and int(result.get("full_dependency_survivor_leaves", -1))
        == EXPECTED_MOD7_SURVIVOR_LEAVES,
        "CUDA survivor census changed",
    )

    matrix = equation_matrix()
    ranks = {}
    dependencies = {}
    for modulus, expected_rank, expected_dimension in ((3, 162, 120), (7, 147, 135)):
        rank, dependency = left_dependencies(matrix, modulus)
        dependency = dependency.astype(np.uint8)
        require(
            rank == expected_rank and dependency.shape == (expected_dimension, 282),
            f"mod-{modulus} score system changed",
        )
        require(
            not np.any(
                dependency.astype(np.int64)
                @ (matrix.astype(np.int64) % modulus)
                % modulus
            ),
            f"mod-{modulus} dependency basis is not left-null",
        )
        ranks[modulus] = rank
        dependencies[modulus] = dependency
    require(
        array_sha256(dependencies[7]) == result["full_dependency_sha256"],
        "independent mod-seven dependency hash changed",
    )

    require(
        table_summary.get("experiment") == "p7_size8_one_elevation_tables"
        and table_summary.get("status")
        == "complete_exact_elevated_direction_omission_tables"
        and table_summary.get("output_sha256") == sha256(args.tables),
        "invalid conditioned table cache",
    )
    with np.load(args.tables, allow_pickle=False) as handle:
        arrays = {key: handle[key] for key in handle.files}
    for key, value in arrays.items():
        require(
            table_summary["array_sha256"].get(key) == array_sha256(value),
            f"conditioned cache hash changed for {key}",
        )
    require(np.array_equal(arrays["labels"], labels), "cached labels changed")
    require(
        np.array_equal(arrays["epsilons"], np.asarray(epsilons, dtype=np.int8)),
        "cached direction signs changed",
    )
    require(
        np.array_equal(arrays["dependency"], dependencies[7]),
        "cached full dependency basis changed",
    )
    coefficients = arrays["selected_coefficients"].astype(np.int64)
    projected_family = arrays["projected_dependencies"].astype(np.uint8)
    require(
        np.array_equal(
            coefficients @ dependencies[7].astype(np.int64) % 7,
            projected_family,
        ),
        "conditioned dependency reconstruction changed",
    )
    for omitted in range(8):
        columns = slice(2 + 35 * omitted, 2 + 35 * (omitted + 1))
        require(not np.any(projected_family[omitted, :, columns]), "omission block is nonzero")
        require(modular_rank(projected_family[omitted]) == 22, "projection lost rank")

    projected_pairs = [
        (int(row[0]), int(row[1]))
        for row in result["projected_survivor_rank_direction_pairs"]
    ]
    require(
        projected_pairs == sorted(set(projected_pairs))
        and len(projected_pairs) == EXPECTED_PROJECTED_LEAVES,
        "projected pair list is incomplete or duplicated",
    )
    evidence_rows = {
        (int(row["rank"]), int(row["elevated_direction_index"])): row
        for row in result["full_135_dependency_recheck_rows"]
    }
    require(set(evidence_rows) == set(projected_pairs), "full recheck rows miss projected pairs")

    projected_cache: dict[tuple[int, int, int, int], np.ndarray] = {}
    full_caches = {3: {}, 7: {}}
    mod7_survivors = set()
    joint_survivors = set()
    mod7_match_histogram: Counter[int] = Counter()
    mod3_match_histogram: Counter[int] = Counter()
    intersection_histogram: Counter[int] = Counter()
    catalog_shape_histogram: Counter[tuple[int, ...]] = Counter()
    candidate_digest = hashlib.sha256()
    for pair in projected_pairs:
        candidate_rank, omitted = pair
        boundary = unrank_lex(candidate_rank)
        require(
            projected_passes(
                boundary,
                omitted,
                labels,
                epsilons,
                projected_family[omitted],
                projected_cache,
            ),
            f"recorded candidate fails conditioned projection: {pair}",
        )
        matches = {}
        metadata = None
        for modulus in (3, 7):
            matches[modulus], rebuilt_metadata = dependency_matches(
                boundary,
                omitted,
                labels,
                epsilons,
                -1,
                dependencies[modulus],
                modulus,
                full_caches[modulus],
            )
            if metadata is None:
                metadata = rebuilt_metadata
            else:
                require(
                    metadata["odd_fibre_masks"] == rebuilt_metadata["odd_fibre_masks"]
                    and metadata["scaled_means_direction_order"]
                    == rebuilt_metadata["scaled_means_direction_order"]
                    and metadata["variable_catalogs"] == rebuilt_metadata["variable_catalogs"],
                    "catalog metadata depends on modulus",
                )
        intersection = matches[3] & matches[7]
        if matches[7]:
            mod7_survivors.add(pair)
        if intersection:
            joint_survivors.add(pair)
        mod3_match_histogram[len(matches[3])] += 1
        mod7_match_histogram[len(matches[7])] += 1
        intersection_histogram[len(intersection)] += 1
        catalog_shape_histogram[tuple(row["catalog_rows"] for row in metadata["variable_catalogs"])] += 1

        evidence = evidence_rows[pair]
        require(
            evidence["boundary_finite_field"] == list(boundary)
            and evidence["boundary_vertices"] == [value + 1 for value in boundary]
            and evidence["odd_fibre_masks"] == metadata["odd_fibre_masks"]
            and evidence["scaled_means_direction_order"]
            == metadata["scaled_means_direction_order"]
            and evidence["variable_catalogs"] == metadata["variable_catalogs"]
            and int(evidence["full_dependency_match_count"]) == len(matches[7])
            and {tuple(int(value) for value in row) for row in evidence["matching_catalog_rows"]}
            == matches[7],
            f"CUDA full recheck row changed: {pair}",
        )
        candidate_digest.update(np.asarray(pair, dtype=np.uint64).tobytes())
        candidate_digest.update(np.asarray(metadata["odd_fibre_masks"], dtype=np.uint8).tobytes())
        candidate_digest.update(np.asarray([len(matches[3]), len(matches[7]), len(intersection)], dtype=np.uint16).tobytes())

    recorded_mod7 = {
        (int(row[0]), int(row[1]))
        for row in result["full_dependency_survivor_rank_direction_pairs"]
    }
    require(mod7_survivors == recorded_mod7, "independent mod-seven survivor set changed")
    require(len(mod7_survivors) == EXPECTED_MOD7_SURVIVOR_LEAVES, "mod-seven survivor count changed")
    require(not joint_survivors, "a catalog row survives both mod three and mod seven")

    geometric_minus = line_plus_point_family(labels, epsilons, -1)
    require(
        mod7_survivors == geometric_minus,
        "mod-seven survivor set is not exactly the negative-type line-plus-point family",
    )
    require(
        mod7_match_histogram == Counter({0: 15, 2: 1_176}),
        "mod-seven match histogram changed",
    )
    require(intersection_histogram == Counter({0: 1_191}), "joint match histogram changed")

    symmetry = p7_nonsquare_signed_permutation()
    require(
        symmetry["fixes_distinguished_edge"] is True
        and symmetry["signed_conference_anti_isometry"] is True,
        "Paley nonsquare anti-isometry failed",
    )
    finite_permutation = tuple(int(value) for value in symmetry["finite_permutation"])
    transfer, direction_map = audit_direction_transfer(labels, epsilons, finite_permutation)
    geometric_plus = line_plus_point_family(labels, epsilons, 1)
    mapped_family = set()
    for candidate_rank, old_direction in geometric_minus:
        boundary = unrank_lex(candidate_rank)
        mapped_boundary = tuple(sorted(finite_permutation[point] for point in boundary))
        mapped_family.add((rank_lex(mapped_boundary), direction_map[old_direction]))
    require(mapped_family == geometric_plus, "line-plus-point survivor family did not transfer")

    return {
        "experiment": "p7_size8_one_elevation_audit",
        "status": "passed_independent_complete_four_allocation_exclusion_audit",
        "p": 7,
        "finite_boundary_size": 8,
        "gpu_result": str(args.gpu_result),
        "gpu_result_sha256": sha256(args.gpu_result),
        "conditioned_tables": str(args.tables),
        "conditioned_tables_sha256": sha256(args.tables),
        "conditioned_table_summary": str(args.table_summary),
        "conditioned_table_summary_sha256": sha256(args.table_summary),
        "floor_sources": {
            "c_H=-1": str(args.floor_minus),
            "c_H=+1": str(args.floor_plus),
        },
        "floor_source_sha256": {
            "c_H=-1": sha256(args.floor_minus),
            "c_H=+1": sha256(args.floor_plus),
        },
        "source_scope_by_sign": {"c_H=-1": scope_minus, "c_H=+1": scope_plus},
        "all_boundaries_checked_by_cuda": ALL_BOUNDARIES,
        "four_allocation_boundaries_each_sign": EXPECTED_ALLOCATION_BOUNDARIES[4],
        "four_allocation_leaves_each_sign": 4 * EXPECTED_ALLOCATION_BOUNDARIES[4],
        "full_linear_system_shape": list(matrix.shape),
        "full_linear_system_ranks": {str(key): value for key, value in ranks.items()},
        "full_dependency_shapes": {
            str(key): list(value.shape) for key, value in dependencies.items()
        },
        "full_dependency_sha256": {
            str(key): array_sha256(value) for key, value in dependencies.items()
        },
        "projected_survivor_leaves_reconstructed": len(projected_pairs),
        "projected_survivor_boundaries_reconstructed": len(
            {rank for rank, _direction in projected_pairs}
        ),
        "candidate_audit_sha256": candidate_digest.hexdigest(),
        "variable_catalog_shape_histogram": [
            {"catalog_rows": list(key), "candidates": value}
            for key, value in sorted(catalog_shape_histogram.items())
        ],
        "mod3_match_count_histogram": dict(sorted(mod3_match_histogram.items())),
        "mod7_match_count_histogram": dict(sorted(mod7_match_histogram.items())),
        "joint_mod3_mod7_match_count_histogram": dict(sorted(intersection_histogram.items())),
        "mod7_survivor_leaves": len(mod7_survivors),
        "mod7_survivor_geometry": {
            "description": "one affine line plus one point off that line",
            "line_direction_type_for_cminus1": -1,
            "direction_families": 4,
            "parallel_lines_per_direction": 7,
            "off_line_points_per_line": 42,
            "count": len(geometric_minus),
            "all_have_odd_fibre_profile": "2 in the line-normal direction and 6 in each other direction",
            "all_have_two_mod7_catalog_rows": True,
        },
        "joint_mod3_mod7_survivor_leaves": len(joint_survivors),
        "all_cminus_four_allocation_leaves_excluded": not joint_survivors,
        "nonsquare_anti_isometry": {
            "nonsquare_multiplier": int(symmetry["nonsquare_multiplier"]),
            "fixes_distinguished_edge": True,
            "signed_conference_anti_isometry": True,
        },
        "sign_transfer": transfer,
        "line_plus_point_family_transfers_bijectively": True,
        "all_cplus_four_allocation_leaves_excluded": not joint_survivors,
        "all_four_allocation_boundaries_both_signs_excluded": not joint_survivors,
        "remaining_nonconic_floor_survivors_each_sign": 1_419_432,
        "closes_all_nonconic_size_eight": False,
        "closes_all_p7_size_eight": False,
        "closes_residual_ii": False,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gpu-result", type=Path, required=True)
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--table-summary", type=Path, required=True)
    parser.add_argument("--floor-minus", type=Path, required=True)
    parser.add_argument("--floor-plus", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = run(args)
    atomic_json(args.output, out)
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
