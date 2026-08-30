#!/usr/bin/env python3
"""Exact mod-seven catalog exhaustion for positive p=7 infinity+7, z=2.

The distinguished-edge equation is deliberately omitted.  This leaves the
translation-equivariant 281-by-1225 edge system, of rank 146 over F_7 and
with 135 left dependencies.  Consequently the full 2,352-element affine
square-semilinear group may be used for boundary orbit reduction.

All exact mean leaves are covered:

* 192 residue-zero leaves have one q=2 direction.  The full 112-dimensional
  subspace of dependencies annihilating that entire direction block rejects
  them under a necessary block-image relaxation larger than the catalog.
* 992 residue-zero leaves have two complete q=1 catalogs.  Their rows are
  joined exactly in all 135 dependency coordinates.
* 48 same-type leaves have residue four: the two b=7 directions have mean 4,
  the other two directions of that type have mean 12, and the opposite type
  remains at mean 8.  Their four complete catalogs are joined exactly by a
  2+2 meet-in-the-middle join in all 135 dependency coordinates.

Passing a modular test is only a necessary condition for an edge lift.  A
zero survivor count, however, is a rigorous exclusion.  In particular, this
program never reports an affine-span relaxation as an exact catalog result.
"""
from __future__ import annotations

import argparse
import functools
import hashlib
import itertools
import json
import math
import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
)
from p7_unsaturated_modular_catalog_filter import (  # noqa: E402
    equation_matrix,
    left_dependencies,
)
from p7_unsaturated_slack_catalog import exact_slack_catalog_values  # noqa: E402


P = 7
Q = P * P
MODULUS = 7
EDGE_COUNT = 4 * P + 1
EXPECTED_PAIR_INCIDENCES = math.comb(P + 1, 2) * math.factorial(P)
EXPECTED_Z_HISTOGRAM = {2: 123_480, 3: 5_488, 7: 56}
EXPECTED_Z2_BOUNDARIES = 123_480
EXPECTED_GROUP_SIZE = Q * (Q - 1)
EXPECTED_ORBIT_COUNT = 92
EXPECTED_ORBIT_SIZE_HISTOGRAM = {588: 18, 1176: 52, 2352: 22}
EXPECTED_SAME_TYPE_ORBITS = 48
EXPECTED_SPLIT_TYPE_ORBITS = 44
EXPECTED_HIGH_LEAVES = 192
EXPECTED_PAIR_LEAVES = 992
EXPECTED_RESIDUE4_LEAVES = 48
EXPECTED_TOTAL_LEAVES = 1232

POINTS = tuple(itertools.combinations(range(P), 4))
POINT_INDEX = {point: index for index, point in enumerate(POINTS)}
DIRECTIONS = projective_directions(P)
DIRECTION_DATA = tuple(field_direction_data(P, direction) for direction in DIRECTIONS)
DIRECTION_TYPES = tuple(int(row[0]) for row in DIRECTION_DATA)
LABELS = tuple(tuple(int(value) for value in row[1]) for row in DIRECTION_DATA)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def atomic_write(path: Path, payload: dict) -> None:
    """Write JSON by same-directory replace, leaving no partial result."""
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    try:
        temporary.write_text(json.dumps(payload, indent=2) + "\n")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def matrix_sha256(matrix: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(matrix).tobytes()).hexdigest()


def modular_rank(matrix: np.ndarray, modulus: int = MODULUS) -> int:
    work = np.asarray(matrix, dtype=np.int64).copy() % modulus
    require(work.ndim == 2, "modular rank input must be a matrix")
    row = 0
    for column in range(work.shape[1]):
        candidates = np.flatnonzero(work[row:, column])
        if not len(candidates):
            continue
        pivot = row + int(candidates[0])
        work[[row, pivot]] = work[[pivot, row]]
        work[row] = work[row] * pow(int(work[row, column]), -1, modulus) % modulus
        for target in range(row + 1, work.shape[0]):
            if work[target, column]:
                work[target] = (
                    work[target] - work[target, column] * work[row]
                ) % modulus
        row += 1
        if row == work.shape[0]:
            break
    return row


def modular_right_nullspace(
    matrix: np.ndarray, modulus: int = MODULUS
) -> tuple[np.ndarray, int]:
    """Return a row basis of {x: matrix @ x = 0} over F_modulus."""
    work = np.asarray(matrix, dtype=np.int64).copy() % modulus
    require(work.ndim == 2, "nullspace input must be a matrix")
    rows, columns = work.shape
    pivots: list[int] = []
    row = 0
    for column in range(columns):
        candidates = np.flatnonzero(work[row:, column])
        if not len(candidates):
            continue
        pivot = row + int(candidates[0])
        work[[row, pivot]] = work[[pivot, row]]
        work[row] = work[row] * pow(int(work[row, column]), -1, modulus) % modulus
        for target in range(rows):
            if target != row and work[target, column]:
                work[target] = (
                    work[target] - work[target, column] * work[row]
                ) % modulus
        pivots.append(column)
        row += 1
        if row == rows:
            break
    pivot_set = set(pivots)
    basis = []
    for free_column in (column for column in range(columns) if column not in pivot_set):
        vector = np.zeros(columns, dtype=np.int64)
        vector[free_column] = 1
        for pivot_row, pivot_column in enumerate(pivots):
            vector[pivot_column] = -work[pivot_row, free_column] % modulus
        basis.append(vector)
    return np.asarray(basis, dtype=np.int64), row


def translation_equivariant_system() -> tuple[np.ndarray, np.ndarray, dict]:
    """Delete only the distinguished-edge row from the common edge system."""
    full = equation_matrix()
    require(full.shape == (282, 1225), "fixed-edge source system changed shape")
    require(np.all(full[0] == 1), "row zero is no longer the edge-count equation")
    require(
        int(np.count_nonzero(full[1])) == 1 and int(np.sum(full[1])) == 1,
        "row one is no longer the distinguished-edge equation",
    )
    matrix = np.concatenate((full[:1], full[2:]), axis=0)
    rank, dependencies = left_dependencies(matrix, MODULUS)
    require(matrix.shape == (281, 1225), "translation-equivariant shape changed")
    require(rank == 146, f"translation-equivariant rank changed: {rank}")
    require(
        dependencies.shape == (135, 281),
        f"left dependency shape changed: {dependencies.shape}",
    )
    require(
        not np.any(dependencies @ (matrix.astype(np.int64) % MODULUS) % MODULUS),
        "left-null dependency audit failed",
    )
    calibration = (17 * np.arange(matrix.shape[1], dtype=np.int64) + 3) % 2
    manufactured_rhs = matrix.astype(np.int64) @ calibration % MODULUS
    require(
        not np.any(dependencies @ manufactured_rhs % MODULUS),
        "manufactured consistent right side was rejected",
    )
    return matrix, dependencies.astype(np.int64), {
        "source_equations": 282,
        "omitted_source_row": 1,
        "omitted_equation": "distinguished edge",
        "equations": 281,
        "edge_variables": 1225,
        "edge_count_rhs": EDGE_COUNT,
        "direction_block_offset": 1,
        "direction_block_width": 35,
        "modulus": MODULUS,
        "rank": rank,
        "left_dependency_dimension": len(dependencies),
        "left_null_audit": True,
        "manufactured_rhs_calibration": True,
        "matrix_sha256": matrix_sha256(matrix),
        "dependency_sha256": matrix_sha256(dependencies.astype(np.uint8)),
    }


def pair_transversal_boundaries() -> tuple[set[tuple[int, ...]], dict]:
    """Generate every boundary having at least two b=7 directions."""
    permutations = tuple(itertools.permutations(range(P)))
    incidence: Counter[tuple[int, ...]] = Counter()
    for first, second in itertools.combinations(range(P + 1), 2):
        grid = [[-1] * P for _ in range(P)]
        for point in range(Q):
            a = LABELS[first][point]
            b = LABELS[second][point]
            require(grid[a][b] == -1, "two directions failed unique intersection")
            grid[a][b] = point
        require(
            all(value >= 0 for row in grid for value in row),
            "two-direction intersection grid is incomplete",
        )
        for permutation in permutations:
            boundary = tuple(sorted(grid[index][permutation[index]] for index in range(P)))
            incidence[boundary] += 1

    require(
        sum(incidence.values()) == EXPECTED_PAIR_INCIDENCES,
        "pair-transversal incidence count changed",
    )
    z_histogram: Counter[int] = Counter()
    multiplicity_histogram: Counter[int] = Counter()
    for boundary, multiplicity in incidence.items():
        z = 0
        for labels in LABELS:
            mask = 0
            for point in boundary:
                mask ^= 1 << labels[point]
            z += int(mask.bit_count() == P)
        require(z >= 2, "generated boundary has fewer than two transversals")
        require(
            multiplicity == math.comb(z, 2),
            "pair-incidence multiplicity disagrees with direct z count",
        )
        z_histogram[z] += 1
        multiplicity_histogram[multiplicity] += 1
    require(dict(sorted(z_histogram.items())) == EXPECTED_Z_HISTOGRAM, "z census changed")
    z2 = {boundary for boundary, multiplicity in incidence.items() if multiplicity == 1}
    require(len(z2) == EXPECTED_Z2_BOUNDARIES, "z=2 boundary count changed")
    return z2, {
        "direction_pairs": math.comb(P + 1, 2),
        "permutations_per_pair": math.factorial(P),
        "pair_transversal_incidences": sum(incidence.values()),
        "distinct_boundaries_with_z_at_least_2": len(incidence),
        "boundary_histogram_by_z": {
            str(key): value for key, value in sorted(z_histogram.items())
        },
        "incidence_multiplicity_histogram": {
            str(key): value for key, value in sorted(multiplicity_histogram.items())
        },
        "direct_z_and_pair_multiplicity_agree": True,
    }


def affine_square_semilinear_group() -> tuple[tuple[int, ...], ...]:
    q, multiply, add, character, frobenius, _norm, _ia, _ib = field_ctx(P)
    permutations = {
        tuple(
            add(
                multiply(alpha, frobenius(point) if use_frobenius else point),
                translate,
            )
            for point in range(q)
        )
        for translate in range(q)
        for alpha in range(1, q)
        if character(alpha) == 1
        for use_frobenius in (False, True)
    }
    require(len(permutations) == EXPECTED_GROUP_SIZE, "affine group size changed")
    require(
        all(tuple(sorted(permutation)) == tuple(range(q)) for permutation in permutations),
        "affine group contains a non-permutation",
    )
    return tuple(sorted(permutations))


def boundary_masks(boundary: tuple[int, ...]) -> tuple[int, ...]:
    masks = []
    for labels in LABELS:
        mask = 0
        for point in boundary:
            mask ^= 1 << labels[point]
        masks.append(mask)
    return tuple(masks)


def boundary_orbits(z2: set[tuple[int, ...]]) -> tuple[list[dict], dict]:
    group = affine_square_semilinear_group()
    remaining = set(z2)
    orbits = []
    while remaining:
        representative = min(remaining)
        orbit = {
            tuple(sorted(permutation[point] for point in representative))
            for permutation in group
        }
        require(orbit <= z2, "z=2 set is not invariant under the affine group")
        masks = boundary_masks(representative)
        undetermined = tuple(
            direction for direction, mask in enumerate(masks) if mask.bit_count() == P
        )
        require(len(undetermined) == 2, "orbit representative is not z=2")
        same_type = DIRECTION_TYPES[undetermined[0]] == DIRECTION_TYPES[undetermined[1]]
        orbits.append(
            {
                "representative": representative,
                "size": len(orbit),
                "masks": masks,
                "b_values": tuple(mask.bit_count() for mask in masks),
                "undetermined": undetermined,
                "same_type": same_type,
            }
        )
        remaining -= orbit

    size_histogram = Counter(int(orbit["size"]) for orbit in orbits)
    same_type = sum(bool(orbit["same_type"]) for orbit in orbits)
    require(len(orbits) == EXPECTED_ORBIT_COUNT, "z=2 orbit count changed")
    require(sum(int(orbit["size"]) for orbit in orbits) == len(z2), "orbit weights changed")
    require(
        dict(sorted(size_histogram.items())) == EXPECTED_ORBIT_SIZE_HISTOGRAM,
        "orbit-size histogram changed",
    )
    require(same_type == EXPECTED_SAME_TYPE_ORBITS, "same-type orbit count changed")
    require(
        len(orbits) - same_type == EXPECTED_SPLIT_TYPE_ORBITS,
        "split-type orbit count changed",
    )
    return orbits, {
        "group": "full affine square-semilinear group",
        "group_size": len(group),
        "orbit_count": len(orbits),
        "orbit_size_sum": sum(int(orbit["size"]) for orbit in orbits),
        "orbit_size_histogram": {
            str(key): value for key, value in sorted(size_histogram.items())
        },
        "same_type_orbits": same_type,
        "split_type_orbits": len(orbits) - same_type,
        "orbit_invariance_audit": True,
    }


@functools.lru_cache(maxsize=None)
def weak_compositions(total: int, parts: int) -> tuple[tuple[int, ...], ...]:
    if parts == 1:
        return ((total,),)
    return tuple(
        (first, *tail)
        for first in range(total + 1)
        for tail in weak_compositions(total - first, parts - 1)
    )


def exact_mean_leaves(orbits: list[dict]) -> tuple[list[list[dict]], dict]:
    """Return the 1,232 exact leaves, including the residue-four family."""
    by_orbit: list[list[dict]] = []
    kind_histogram: Counter[str] = Counter()
    pair_pattern_histogram: Counter[tuple[int, int]] = Counter()
    residue_pattern_histogram: Counter[tuple[int, ...]] = Counter()
    for orbit_index, orbit in enumerate(orbits):
        masks = orbit["masks"]
        b_values = orbit["b_values"]
        type_options = []
        for quadratic_type in (-1, 1):
            directions = tuple(
                index for index, value in enumerate(DIRECTION_TYPES) if value == quadratic_type
            )
            require(len(directions) == 4, "quadratic direction type is not four-by-four")
            missing_units = sum(b_values[index] == P for index in directions)
            type_options.append(
                tuple((directions, row) for row in weak_compositions(missing_units, 4))
            )

        leaves = []
        for left, right in itertools.product(*type_options):
            q_values = [0] * (P + 1)
            for directions, values in (left, right):
                for direction, value in zip(directions, values):
                    q_values[direction] = value
            nonzero = [value for value in q_values if value]
            means = tuple(
                (0 if b_values[direction] == P else 8) + 8 * q_values[direction]
                for direction in range(P + 1)
            )
            if sorted(nonzero) == [2]:
                kind = "one_high_q2"
                variable_directions = ()
                omitted_direction = q_values.index(2)
            elif sorted(nonzero) == [1, 1]:
                kind = "two_q1"
                variable_directions = tuple(
                    direction for direction, value in enumerate(q_values) if value == 1
                )
                omitted_direction = None
                sizes = tuple(
                    sorted(2233 if b_values[direction] == 3 else 1764 for direction in variable_directions)
                )
                pair_pattern_histogram[sizes] += 1
            else:
                raise AssertionError(f"unexpected residue-zero allocation {q_values}")
            leaves.append(
                {
                    "orbit_index": orbit_index,
                    "kind": kind,
                    "q_values": tuple(q_values),
                    "means": means,
                    "variable_directions": variable_directions,
                    "omitted_direction": omitted_direction,
                }
            )
            kind_histogram[kind] += 1

        if orbit["same_type"]:
            undetermined = orbit["undetermined"]
            active_type = DIRECTION_TYPES[undetermined[0]]
            active = tuple(
                direction
                for direction, value in enumerate(DIRECTION_TYPES)
                if value == active_type
            )
            means = tuple(
                4
                if direction in undetermined
                else 12
                if direction in active
                else 8
                for direction in range(P + 1)
            )
            sizes = tuple(
                sorted(62 if b_values[direction] == 3 else 56 for direction in active)
            )
            residue_pattern_histogram[sizes] += 1
            leaves.append(
                {
                    "orbit_index": orbit_index,
                    "kind": "residue4_four_catalog",
                    "q_values": None,
                    "means": means,
                    "variable_directions": active,
                    "omitted_direction": None,
                }
            )
            kind_histogram["residue4_four_catalog"] += 1
        by_orbit.append(leaves)

    require(kind_histogram["one_high_q2"] == EXPECTED_HIGH_LEAVES, "q=2 coverage changed")
    require(kind_histogram["two_q1"] == EXPECTED_PAIR_LEAVES, "two-q=1 coverage changed")
    require(
        kind_histogram["residue4_four_catalog"] == EXPECTED_RESIDUE4_LEAVES,
        "residue-four coverage changed",
    )
    require(sum(kind_histogram.values()) == EXPECTED_TOTAL_LEAVES, "total leaf coverage changed")
    require(
        pair_pattern_histogram
        == Counter({(1764, 1764): 282, (1764, 2233): 508, (2233, 2233): 202}),
        "two-catalog size-pattern census changed",
    )
    require(
        residue_pattern_histogram
        == Counter(
            {
                (56, 56, 56, 56): 6,
                (56, 56, 56, 62): 20,
                (56, 56, 62, 62): 22,
            }
        ),
        "residue-four catalog-size census changed",
    )
    require(
        all(len(leaves) == (11 if orbit["same_type"] else 16) for orbit, leaves in zip(orbits, by_orbit)),
        "per-orbit exact leaf count changed",
    )
    return by_orbit, {
        "exact_mean_leaves": sum(kind_histogram.values()),
        "kind_histogram": dict(sorted(kind_histogram.items())),
        "two_q1_catalog_size_patterns": {
            "x".join(map(str, key)): value
            for key, value in sorted(pair_pattern_histogram.items())
        },
        "residue4_catalog_size_patterns": {
            "x".join(map(str, key)): value
            for key, value in sorted(residue_pattern_histogram.items())
        },
        "same_type_leaves_per_orbit": 11,
        "split_type_leaves_per_orbit": 16,
        "obsolete_1184_count_rejected": sum(kind_histogram.values()) != 1184,
    }


def parity_row(b: int) -> np.ndarray:
    B = set(range(b))
    return np.asarray(
        [sum(value in B for value in point) & 1 for point in POINTS],
        dtype=np.int16,
    )


@functools.lru_cache(maxsize=None)
def canonical_catalog(b: int, mean: int) -> np.ndarray:
    """Complete phase-zero catalog in canonical fibre coordinates."""
    if b == P:
        # On J(7,4), intersection with all seven fibres has even parity.
        rows = canonical_catalog(0, mean)
    elif b in (1, 5) and mean >= 8:
        # The unique parity floor plus an arbitrary even degree-two excess.
        rows = canonical_catalog(0, mean - 8) + parity_row(b)[None, :]
    else:
        rows = np.asarray(exact_slack_catalog_values(b, 0, mean), dtype=np.int16)
    rows = np.ascontiguousarray(rows, dtype=np.int16)
    require(rows.ndim == 2 and rows.shape[1] == 35, "catalog shape changed")
    require(len({row.tobytes() for row in rows}) == len(rows), "canonical catalog duplicates")
    require(
        np.all(2 * rows.sum(axis=1, dtype=np.int64) == 5 * mean),
        "canonical catalog mean changed",
    )
    canonical_parity = parity_row(b)
    require(
        np.all((rows % 2) == canonical_parity[None, :]),
        "canonical catalog parity changed",
    )
    require(np.all((0 <= rows) & (rows <= 13)), "catalog leaves affine bad-count range")
    return rows


@functools.lru_cache(maxsize=None)
def mapped_catalog(mask: int, mean: int) -> np.ndarray:
    B = {index for index in range(P) if mask & (1 << index)}
    b = len(B)
    canonical = canonical_catalog(b, mean)
    actual_B = sorted(B)
    actual_complement = sorted(set(range(P)) - B)
    permutation = dict(zip(range(b), actual_B)) | dict(
        zip(range(b, P), actual_complement)
    )
    inverse = {target: source for source, target in permutation.items()}
    source_columns = [
        POINT_INDEX[tuple(sorted(inverse[value] for value in point))]
        for point in POINTS
    ]
    mapped = np.ascontiguousarray(canonical[:, source_columns], dtype=np.int16)
    require(len({row.tobytes() for row in mapped}) == len(mapped), "mapped catalog duplicates")
    actual_parity = np.asarray(
        [sum(value in B for value in point) & 1 for point in POINTS],
        dtype=np.int16,
    )
    require(np.all((mapped % 2) == actual_parity[None, :]), "mapped catalog parity changed")
    return mapped


def expected_catalog_audit() -> dict:
    expected = {
        (1, 8): 1,
        (3, 8): 1,
        (5, 8): 1,
        (7, 0): 1,
        (7, 8): 1764,
        (1, 16): 1764,
        (3, 16): 2233,
        (5, 16): 1764,
        (7, 4): 56,
        (1, 12): 56,
        (3, 12): 62,
        (5, 12): 56,
    }
    observed = {(b, mean): len(canonical_catalog(b, mean)) for b, mean in expected}
    require(observed == expected, f"exact catalog sizes changed: {observed}")
    return {
        f"b{b}_mean{mean}": count for (b, mean), count in sorted(observed.items())
    }


def row_keys(rows: np.ndarray) -> np.ndarray:
    rows = np.ascontiguousarray(rows, dtype=np.uint8)
    require(rows.ndim == 2 and rows.shape[1] > 0, "signature rows have bad shape")
    return rows.view(np.dtype((np.void, rows.shape[1]))).reshape(-1)


def summed_rows(contributions: tuple[np.ndarray, ...]) -> np.ndarray:
    """Enumerate exact syndrome sums; matrices are coordinate-by-catalog."""
    require(bool(contributions), "cannot enumerate an empty catalog side")
    rows = np.zeros((1, contributions[0].shape[0]), dtype=np.uint8)
    for matrix in contributions:
        require(matrix.shape[0] == rows.shape[1], "join dimensions disagree")
        rows = (
            rows[:, None, :].astype(np.int16)
            + matrix.T[None, :, :].astype(np.int16)
        ) % MODULUS
        rows = np.ascontiguousarray(rows.reshape(-1, rows.shape[-1]), dtype=np.uint8)
    return rows


def count_matching_rows(first: np.ndarray, needed: np.ndarray) -> int:
    """Count matches using injective raw-byte keys, including multiplicity."""
    first_keys, counts = np.unique(row_keys(first), return_counts=True)
    needed_keys = row_keys(needed)
    positions = np.searchsorted(first_keys, needed_keys)
    valid = positions < len(first_keys)
    if not np.any(valid):
        return 0
    valid_indices = np.flatnonzero(valid)
    equal = first_keys[positions[valid_indices]] == needed_keys[valid_indices]
    hits = valid_indices[equal]
    return int(np.sum(counts[positions[hits]], dtype=np.int64))


def choose_two_plus_two(sizes: tuple[int, int, int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    best = None
    for first in itertools.combinations(range(4), 2):
        if 0 not in first:
            continue
        second = tuple(index for index in range(4) if index not in first)
        products = (math.prod(sizes[index] for index in first), math.prod(sizes[index] for index in second))
        score = (max(products), min(products), first)
        if best is None or score < best[0]:
            best = (score, first, second)
    require(best is not None, "2+2 join partition failed")
    return best[1], best[2]


def exact_join_count(
    base: np.ndarray, contributions: tuple[np.ndarray, ...]
) -> tuple[int, dict]:
    """Count complete-catalog tuples summing with base to zero in F_7^135."""
    require(base.shape == (135,), "join base must have 135 coordinates")
    require(len(contributions) in (2, 4), "only exact 1+1 and 2+2 joins are expected")
    sizes = tuple(int(matrix.shape[1]) for matrix in contributions)
    if len(contributions) == 2:
        first_indices = (0,)
        second_indices = (1,)
    else:
        first_indices, second_indices = choose_two_plus_two(sizes)  # type: ignore[arg-type]
    first = summed_rows(tuple(contributions[index] for index in first_indices))
    second = summed_rows(tuple(contributions[index] for index in second_indices))
    needed = np.ascontiguousarray(
        (-base[None, :].astype(np.int16) - second.astype(np.int16)) % MODULUS,
        dtype=np.uint8,
    )
    count = count_matching_rows(first, needed)
    return count, {
        "partition": [list(first_indices), list(second_indices)],
        "catalog_sizes": list(sizes),
        "enumerated_signatures": [len(first), len(second)],
        "cartesian_catalog_tuples": math.prod(sizes),
        "exact_injective_signature_encoding": "135 raw base-seven bytes",
    }


def join_implementation_self_audit() -> dict:
    """Compare both join shapes with positive-hit duplicate-aware brute force."""
    audits = []
    cases = (
        (
            np.asarray([0], dtype=np.uint8),
            (
                np.asarray([[0, 0, 1]], dtype=np.uint8),
                np.asarray([[0, 6, 6, 1]], dtype=np.uint8),
            ),
        ),
        (
            np.asarray([0], dtype=np.uint8),
            (
                np.asarray([[0, 0]], dtype=np.uint8),
                np.asarray([[0, 1, 1]], dtype=np.uint8),
                np.asarray([[0, 6]], dtype=np.uint8),
                np.asarray([[0, 0, 6]], dtype=np.uint8),
            ),
        ),
    )
    for base, matrices in cases:
        sizes = tuple(matrix.shape[1] for matrix in matrices)
        dimension = len(base)
        # The production helper requires 135 coordinates.  Pad with zeros so
        # this small independently brute-forced instance exercises that path.
        padded_base = np.pad(base, (0, 135 - dimension))
        padded = tuple(np.pad(matrix, ((0, 135 - dimension), (0, 0))) for matrix in matrices)
        joined, metadata = exact_join_count(padded_base, padded)
        brute = 0
        for indices in itertools.product(*(range(size) for size in sizes)):
            syndrome = base.astype(np.int16).copy()
            for matrix, index in zip(matrices, indices):
                syndrome += matrix[:, index]
            brute += int(not np.any(syndrome % MODULUS))
        require(joined == brute, "meet-in-the-middle join failed brute-force audit")
        require(joined > 0, "join self-audit did not exercise a successful match")
        require(
            any(len({bytes(matrix[:, index]) for index in range(matrix.shape[1])}) < matrix.shape[1] for matrix in matrices),
            "join self-audit did not exercise duplicate multiplicity",
        )
        audits.append({"catalog_sizes": list(sizes), "join_count": joined, **metadata})
    return {
        "passed": True,
        "positive_matches_exercised": True,
        "duplicate_multiplicity_exercised": True,
        "cases": audits,
    }


def direction_annihilators(dependencies: np.ndarray) -> tuple[tuple[np.ndarray, ...], dict]:
    coefficients = []
    rows = []
    for direction in range(P + 1):
        columns = slice(1 + 35 * direction, 1 + 35 * (direction + 1))
        block = dependencies[:, columns]
        basis, block_rank = modular_right_nullspace(block.T)
        conditioned = basis @ dependencies % MODULUS
        require(block_rank == 23, f"direction {direction} block rank changed")
        require(basis.shape == (112, 135), "annihilator dimension changed")
        require(modular_rank(basis) == 112, "annihilator coefficient basis lost rank")
        require(modular_rank(conditioned) == 112, "conditioned dependencies lost rank")
        require(not np.any(conditioned[:, columns]), "conditioned block is not zero")
        coefficients.append(basis.astype(np.int64))
        rows.append(
            {
                "direction": direction,
                "direction_block_rank": block_rank,
                "annihilator_dimension": len(basis),
                "conditioned_rank": modular_rank(conditioned),
                "conditioned_block_zero": True,
                "coefficient_sha256": matrix_sha256(basis.astype(np.uint8)),
            }
        )
    return tuple(coefficients), {
        "method": "full right nullspace of each 35-column direction block",
        "affine_span_relaxation_used": False,
        "directions": rows,
    }


class ContributionFactory:
    def __init__(self, dependencies: np.ndarray):
        self.dependencies = dependencies
        self.cache: dict[tuple[int, int, int], np.ndarray] = {}

    def get(self, direction: int, mask: int, mean: int) -> np.ndarray:
        key = (direction, mask, mean)
        if key not in self.cache:
            values = mapped_catalog(mask, mean).astype(np.int64)
            bad_counts = 13 - values
            require(np.all(bad_counts >= 0), "catalog score exceeds affine edge range")
            block = self.dependencies[
                :, 1 + 35 * direction : 1 + 35 * (direction + 1)
            ]
            contribution = block @ bad_counts.T % MODULUS
            self.cache[key] = np.ascontiguousarray(contribution, dtype=np.uint8)
        return self.cache[key]


def render_case(
    orbit_index: int,
    orbit: dict,
    leaf_index: int,
    leaf: dict,
    test: str,
    exact_mod7_catalog_tuples: int | None,
    join: dict | None,
) -> dict:
    row = {
        "orbit_index": orbit_index,
        "orbit_leaf_index": leaf_index,
        "orbit_size": int(orbit["size"]),
        "representative_finite_field": list(orbit["representative"]),
        "direction_types": list(DIRECTION_TYPES),
        "b_values": list(orbit["b_values"]),
        "undetermined_directions": list(orbit["undetermined"]),
        "same_quadratic_type": bool(orbit["same_type"]),
        "leaf_kind": leaf["kind"],
        "q_values": None if leaf["q_values"] is None else list(leaf["q_values"]),
        "scaled_means": list(leaf["means"]),
        "variable_directions": list(leaf["variable_directions"]),
        "omitted_direction": leaf["omitted_direction"],
        "mod7_test": test,
        "exact_mod7_catalog_tuples": exact_mod7_catalog_tuples,
        "weighted_boundary_allocation_cases": int(orbit["size"]),
    }
    if exact_mod7_catalog_tuples is not None:
        row["weighted_exact_mod7_catalog_tuples"] = (
            int(orbit["size"]) * exact_mod7_catalog_tuples
        )
    if join is not None:
        row["join"] = join
    return row


def run(orbit_limit: int = 0) -> dict:
    started = time.time()
    require(Counter(DIRECTION_TYPES) == Counter({-1: 4, 1: 4}), "direction types changed")
    matrix, dependencies, linear = translation_equivariant_system()
    z2, generation = pair_transversal_boundaries()
    orbits, orbit_audit = boundary_orbits(z2)
    leaves_by_orbit, leaf_audit = exact_mean_leaves(orbits)
    catalogs = expected_catalog_audit()
    join_audit = join_implementation_self_audit()
    annihilators, annihilator_audit = direction_annihilators(dependencies)

    if orbit_limit < 0:
        raise ValueError("--orbit-limit must be nonnegative")
    process_count = len(orbits) if orbit_limit == 0 else min(orbit_limit, len(orbits))
    full_run = process_count == len(orbits)
    base_edge_syndrome = (dependencies[:, 0] * EDGE_COUNT % MODULUS).astype(np.uint8)
    factory = ContributionFactory(dependencies)
    survivors = []
    decision_digest = hashlib.sha256()
    processed_kind_histogram: Counter[str] = Counter()
    rejected_kind_histogram: Counter[str] = Counter()
    processed_weighted = 0
    rejected_weighted = 0
    surviving_weighted = 0
    weighted_exact_mod7_tuples = 0
    pair_cartesian_tuples = 0
    residue_cartesian_tuples = 0
    per_orbit = []

    for orbit_index, (orbit, leaves) in enumerate(zip(orbits[:process_count], leaves_by_orbit[:process_count])):
        orbit_survivors = 0
        for leaf_index, leaf in enumerate(leaves):
            kind = leaf["kind"]
            processed_kind_histogram[kind] += 1
            processed_weighted += int(orbit["size"])
            fixed = base_edge_syndrome.astype(np.int16).copy()
            variable = set(leaf["variable_directions"])
            omitted = leaf["omitted_direction"]
            for direction, (mask, mean) in enumerate(zip(orbit["masks"], leaf["means"])):
                if direction in variable or direction == omitted:
                    continue
                contribution = factory.get(direction, mask, mean)
                require(contribution.shape[1] == 1, "fixed direction catalog is not unique")
                fixed += contribution[:, 0]
            fixed %= MODULUS

            exact_tuple_count: int | None
            join_metadata: dict | None
            if kind == "one_high_q2":
                require(omitted is not None and not variable, "bad q=2 leaf structure")
                passing = not np.any(annihilators[omitted] @ fixed % MODULUS)
                exact_tuple_count = None
                join_metadata = None
                test = "full_112_dependency_direction_annihilation"
            else:
                contributions = tuple(
                    factory.get(direction, orbit["masks"][direction], leaf["means"][direction])
                    for direction in leaf["variable_directions"]
                )
                expected_sizes = (
                    tuple(2233 if orbit["b_values"][direction] == 3 else 1764 for direction in leaf["variable_directions"])
                    if kind == "two_q1"
                    else tuple(62 if orbit["b_values"][direction] == 3 else 56 for direction in leaf["variable_directions"])
                )
                require(
                    tuple(matrix.shape[1] for matrix in contributions) == expected_sizes,
                    "variable complete-catalog size changed",
                )
                exact_tuple_count, join_metadata = exact_join_count(
                    fixed.astype(np.uint8), contributions
                )
                passing = exact_tuple_count > 0
                test = "exact_full_135_two_catalog_join" if kind == "two_q1" else "exact_full_135_four_catalog_2plus2_join"
                if kind == "two_q1":
                    pair_cartesian_tuples += int(join_metadata["cartesian_catalog_tuples"])
                else:
                    residue_cartesian_tuples += int(join_metadata["cartesian_catalog_tuples"])

            digest_row = {
                "orbit": orbit_index,
                "leaf": leaf_index,
                "kind": kind,
                "means": leaf["means"],
                "passing": passing,
                "tuples": exact_tuple_count,
            }
            decision_digest.update(json.dumps(digest_row, sort_keys=True, separators=(",", ":")).encode())
            decision_digest.update(b"\n")
            if passing:
                orbit_survivors += 1
                surviving_weighted += int(orbit["size"])
                if exact_tuple_count is not None:
                    weighted_exact_mod7_tuples += int(orbit["size"]) * exact_tuple_count
                survivors.append(
                    render_case(
                        orbit_index,
                        orbit,
                        leaf_index,
                        leaf,
                        test,
                        exact_tuple_count,
                        join_metadata,
                    )
                )
            else:
                rejected_kind_histogram[kind] += 1
                rejected_weighted += int(orbit["size"])

        per_orbit.append(
            {
                "orbit_index": orbit_index,
                "orbit_size": int(orbit["size"]),
                "same_quadratic_type": bool(orbit["same_type"]),
                "exact_mean_leaves": len(leaves),
                "surviving_leaves": orbit_survivors,
            }
        )

    processed_leaves = sum(processed_kind_histogram.values())
    require(processed_leaves == sum(len(rows) for rows in leaves_by_orbit[:process_count]), "processed leaf count mismatch")
    require(processed_weighted == rejected_weighted + surviving_weighted, "weighted decision census mismatch")
    if full_run:
        require(processed_leaves == EXPECTED_TOTAL_LEAVES, "full run missed exact mean leaves")
        require(
            processed_kind_histogram
            == Counter(
                {
                    "one_high_q2": EXPECTED_HIGH_LEAVES,
                    "two_q1": EXPECTED_PAIR_LEAVES,
                    "residue4_four_catalog": EXPECTED_RESIDUE4_LEAVES,
                }
            ),
            "full processed kind census changed",
        )

    excluded = full_run and not survivors
    return {
        "experiment": "p7_infinity7_positive_z2_mod7_join",
        "status": (
            "complete_rigorous_mod7_exclusion"
            if excluded
            else "complete_rigorous_mod7_necessary_sieve_with_survivors"
            if full_run
            else "complete_prefix_rigorous_mod7_necessary_sieve"
        ),
        "p": P,
        "c_H": 1,
        "infinity_in_boundary": True,
        "finite_boundary_points": 7,
        "z": 2,
        "phase": 0,
        "linear_system": linear,
        "boundary_generation": generation,
        "orbit_reduction": orbit_audit,
        "mean_leaf_coverage": leaf_audit,
        "catalog_row_counts": catalogs,
        "catalog_source": "complete exact Johnson-slice catalogs",
        "join_self_audit": join_audit,
        "q2_annihilator_audit": annihilator_audit,
        "affine_span_relaxation_used": False,
        "q2_full_direction_block_image_relaxation_used": True,
        "modular_passing_is_edge_feasibility": False,
        "processed_orbits": process_count,
        "full_orbit_count": len(orbits),
        "full_run": full_run,
        "processed_exact_mean_leaves": processed_leaves,
        "processed_kind_histogram": dict(sorted(processed_kind_histogram.items())),
        "rejected_kind_histogram": dict(sorted(rejected_kind_histogram.items())),
        "processed_weighted_boundary_allocation_cases": processed_weighted,
        "rejected_weighted_boundary_allocation_cases": rejected_weighted,
        "surviving_cases": len(survivors),
        "surviving_weighted_boundary_allocation_cases": surviving_weighted,
        "weighted_exact_mod7_catalog_tuples": weighted_exact_mod7_tuples,
        "two_q1_cartesian_catalog_tuples_joined_on_representatives": pair_cartesian_tuples,
        "residue4_cartesian_catalog_tuples_joined_on_representatives": residue_cartesian_tuples,
        "z2_branch_excluded": excluded,
        "all_case_decisions_sha256": decision_digest.hexdigest(),
        "syndrome_contribution_cache_entries": len(factory.cache),
        "per_orbit_summary": per_orbit,
        "survivor_cases": survivors,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--orbit-limit",
        type=int,
        default=0,
        help="process only the first N of 92 orbits; zero means the full run",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    out = run(args.orbit_limit)
    atomic_write(args.output, out)
    if not args.quiet:
        summary = {key: value for key, value in out.items() if key not in {"per_orbit_summary", "survivor_cases"}}
        print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
