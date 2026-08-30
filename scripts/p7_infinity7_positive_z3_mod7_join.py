#!/usr/bin/env python3
"""Exact mod-seven catalog exhaustion for positive p=7 infinity+7, z=3.

The ten boundary-orbit representatives are read from
``evidence/p7_infinity7_positive_zge2_orbits.json``.  The distinguished-edge
equation is omitted, leaving the translation-equivariant 281-by-1225 system
of rank 146 over F_7 and its 135 left dependencies.

All 400 corrected exact mean leaves are covered: 360 residue-00 leaves and
20 leaves in each of residue classes 04 and 40.  Leaves without a high
catalog are joined exactly in all 135 dependency coordinates, using 1+2 or
a cardinality-balanced 2+3 meet in the middle.  A leaf with one high catalog
uses all 112 dependencies annihilating that direction block and then tests
the remaining zero, one, or three complete catalogs exactly in those
conditioned coordinates.

Every passing modular test is only necessary for an integral edge lift.  A
zero survivor count, however, rigorously closes the z=3 branch.
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
MODULUS = 7
EDGE_COUNT = 4 * P + 1
EVIDENCE_PATH = ROOT / "evidence" / "p7_infinity7_positive_zge2_orbits.json"

EXPECTED_ORBIT_COUNT = 10
EXPECTED_BOUNDARIES = 5_488
EXPECTED_ORBIT_SIZE_HISTOGRAM = {392: 6, 784: 4}
EXPECTED_RESIDUE_HISTOGRAM = {"00": 360, "04": 20, "40": 20}
EXPECTED_TOTAL_LEAVES = 400
EXPECTED_PATTERN_HISTOGRAM = {
    "M1764^3": 76,
    "M1764^2*M2233": 104,
    "M1764*M2233^2": 20,
    "H*M1764": 120,
    "H*M2233": 32,
    "H": 8,
    "S56^2*S62^2*M1764": 8,
    "S56^3*S62*M1764": 12,
    "S56^3*S62*M2233": 4,
    "S56^4*M1764": 4,
    "S56^4*M2233": 4,
    "H*S56^3": 8,
}
EXPECTED_KIND_HISTOGRAM = {
    "residue00_three_catalog": 200,
    "residue00_high_plus_one_catalog": 152,
    "residue00_high_only": 8,
    "residue4_five_catalog": 32,
    "residue4_high_plus_three_catalog": 8,
}
MAX_JOIN_CHUNK_ROWS = 100_000

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
    """Write JSON by same-directory replacement, leaving no partial result."""
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
    for free_column in (
        column for column in range(columns) if column not in pivot_set
    ):
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
    dependencies = dependencies.astype(np.int64)
    return matrix, dependencies, {
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


def boundary_masks(boundary: tuple[int, ...]) -> tuple[int, ...]:
    masks = []
    for labels in LABELS:
        mask = 0
        for point in boundary:
            mask ^= 1 << labels[point]
        masks.append(mask)
    return tuple(masks)


def evidence_z3_orbits() -> tuple[list[dict], dict]:
    """Load and independently audit the ten recorded z=3 representatives."""
    encoded = EVIDENCE_PATH.read_bytes()
    payload = json.loads(encoded)
    require(payload.get("p") == P, "orbit evidence has the wrong prime")
    require(
        payload.get("all_required_audits_passed") is True,
        "orbit evidence is not marked fully audited",
    )
    source_rows = [row for row in payload["orbits"] if int(row["z"]) == 3]
    require(len(source_rows) == EXPECTED_ORBIT_COUNT, "z=3 orbit count changed")
    orbits = []
    for branch_index, row in enumerate(source_rows):
        require(
            int(row["branch_orbit_index"]) == branch_index,
            "z=3 branch orbit indices are not contiguous",
        )
        representative = tuple(int(value) for value in row["representative_finite_field"])
        require(
            len(representative) == P and len(set(representative)) == P,
            "z=3 representative is not a seven-set",
        )
        masks = boundary_masks(representative)
        recorded_masks = tuple(int(value) for value in row["direction_masks"])
        require(masks == recorded_masks, "recorded direction masks changed")
        b_values = tuple(mask.bit_count() for mask in masks)
        undetermined = tuple(
            direction for direction, b in enumerate(b_values) if b == P
        )
        require(len(undetermined) == 3, "evidence representative is not z=3")
        require(
            undetermined == tuple(int(value) for value in row["undetermined_directions"]),
            "recorded undetermined directions changed",
        )
        require(
            tuple(int(value) for value in row["direction_types"]) == DIRECTION_TYPES,
            "recorded quadratic direction types changed",
        )
        type_counts = {
            quadratic_type: sum(
                DIRECTION_TYPES[direction] == quadratic_type
                for direction in undetermined
            )
            for quadratic_type in (-1, 1)
        }
        recorded_counts = {
            int(key): int(value)
            for key, value in row["undetermined_type_counts"].items()
        }
        require(type_counts == recorded_counts, "undetermined type census changed")
        majority_types = [value for value in (-1, 1) if type_counts[value] >= 2]
        require(len(majority_types) == 1, "z=3 orbit has no unique majority type")
        orbits.append(
            {
                "source_orbit_index": int(row["orbit_index"]),
                "branch_orbit_index": branch_index,
                "representative": representative,
                "size": int(row["size"]),
                "stabilizer_size": int(row["stabilizer_size"]),
                "masks": masks,
                "b_values": b_values,
                "undetermined": undetermined,
                "undetermined_type_counts": type_counts,
                "majority_type": majority_types[0],
            }
        )

    size_histogram = Counter(int(orbit["size"]) for orbit in orbits)
    majority_histogram = Counter(int(orbit["majority_type"]) for orbit in orbits)
    require(
        dict(sorted(size_histogram.items())) == EXPECTED_ORBIT_SIZE_HISTOGRAM,
        "z=3 orbit-size histogram changed",
    )
    require(sum(size * count for size, count in size_histogram.items()) == EXPECTED_BOUNDARIES,
            "z=3 orbit weights changed")
    require(majority_histogram == Counter({-1: 5, 1: 5}),
            "z=3 majority-type balance changed")
    return orbits, {
        "source": str(EVIDENCE_PATH.relative_to(ROOT)),
        "source_sha256": hashlib.sha256(encoded).hexdigest(),
        "selection": "all source orbits with z=3",
        "orbit_count": len(orbits),
        "orbit_size_sum": sum(int(orbit["size"]) for orbit in orbits),
        "orbit_size_histogram": {
            str(key): value for key, value in sorted(size_histogram.items())
        },
        "majority_type_histogram": {
            str(key): value for key, value in sorted(majority_histogram.items())
        },
        "representative_mask_recomputation_audit": True,
        "direction_type_audit": True,
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


def catalog_size(role: str, b: int) -> int:
    if role == "M":
        return 2233 if b == 3 else 1764
    if role == "S":
        return 62 if b == 3 else 56
    raise AssertionError(f"no complete catalog size for role {role}")


def catalog_pattern(roles: tuple[str, ...], b_values: tuple[int, ...]) -> str:
    labels = []
    for role, b in zip(roles, b_values):
        if role == "H":
            labels.append("H")
        elif role in ("M", "S"):
            labels.append(f"{role}{catalog_size(role, b)}")
    counts = Counter(labels)
    ordered = []
    for label in ("H", "S56", "S62", "M1764", "M2233"):
        count = counts.pop(label, 0)
        if count:
            ordered.append(label if count == 1 else f"{label}^{count}")
    require(not counts, f"unrendered catalog pattern entries: {counts}")
    return "*".join(ordered)


def make_leaf(
    orbit_index: int,
    residue: str,
    kind: str,
    means: tuple[int, ...],
    roles: tuple[str, ...],
    q_values: tuple[int, ...] | None,
    residue_increments: tuple[int, ...] | None,
) -> dict:
    variable = tuple(index for index, role in enumerate(roles) if role in ("M", "S"))
    high = tuple(index for index, role in enumerate(roles) if role == "H")
    require(len(high) <= 1, "leaf has more than one high direction")
    return {
        "orbit_index": orbit_index,
        "residue": residue,
        "kind": kind,
        "means": means,
        "roles": roles,
        "q_values": q_values,
        "residue_increments": residue_increments,
        "variable_directions": variable,
        "high_direction": high[0] if high else None,
    }


def exact_mean_leaves(orbits: list[dict]) -> tuple[list[list[dict]], dict]:
    """Construct and audit all 400 corrected z=3 exact mean leaves."""
    directions_by_type = {
        quadratic_type: tuple(
            index
            for index, value in enumerate(DIRECTION_TYPES)
            if value == quadratic_type
        )
        for quadratic_type in (-1, 1)
    }
    require(
        all(len(directions) == 4 for directions in directions_by_type.values()),
        "quadratic directions are not four-by-four",
    )
    by_orbit: list[list[dict]] = []
    residue_histogram: Counter[str] = Counter()
    pattern_histogram: Counter[str] = Counter()
    kind_histogram: Counter[str] = Counter()

    for orbit_index, orbit in enumerate(orbits):
        b_values = orbit["b_values"]
        z_by_type = {
            quadratic_type: sum(
                b_values[direction] == P
                for direction in directions_by_type[quadratic_type]
            )
            for quadratic_type in (-1, 1)
        }
        require(sum(z_by_type.values()) == 3, "orbit z census changed")
        leaves = []

        standard_options = []
        for quadratic_type in (-1, 1):
            directions = directions_by_type[quadratic_type]
            standard_options.append(
                tuple(
                    (directions, values)
                    for values in weak_compositions(z_by_type[quadratic_type], 4)
                )
            )
        for left, right in itertools.product(*standard_options):
            q_values = [0] * (P + 1)
            for directions, values in (left, right):
                for direction, value in zip(directions, values):
                    q_values[direction] = value
            means = tuple(
                (0 if b_values[direction] == P else 8) + 8 * q_values[direction]
                for direction in range(P + 1)
            )
            roles = tuple(
                "fixed" if value == 0 else "M" if value == 1 else "H"
                for value in q_values
            )
            positive = sorted(value for value in q_values if value)
            if positive == [1, 1, 1]:
                kind = "residue00_three_catalog"
            elif positive == [1, 2]:
                kind = "residue00_high_plus_one_catalog"
            elif positive == [3]:
                kind = "residue00_high_only"
            else:
                raise AssertionError(f"unexpected residue-00 allocation {q_values}")
            leaves.append(
                make_leaf(
                    orbit_index,
                    "00",
                    kind,
                    means,
                    roles,
                    tuple(q_values),
                    None,
                )
            )

        active_type = int(orbit["majority_type"])
        other_type = -active_type
        require(z_by_type[active_type] in (2, 3), "bad residue-four active z")
        active_directions = directions_by_type[active_type]
        other_directions = directions_by_type[other_type]
        active_options = weak_compositions(z_by_type[active_type] - 2, 4)
        other_options = weak_compositions(z_by_type[other_type], 4)
        residue = "40" if active_type == -1 else "04"
        for active_values, other_values in itertools.product(
            active_options, other_options
        ):
            means = [0] * (P + 1)
            roles = ["fixed"] * (P + 1)
            increments = [0] * (P + 1)
            q_values = [0] * (P + 1)
            for direction, value in zip(active_directions, active_values):
                increments[direction] = value
                means[direction] = (
                    (4 if b_values[direction] == P else 12) + 8 * value
                )
                roles[direction] = "S" if value == 0 else "H"
            for direction, value in zip(other_directions, other_values):
                q_values[direction] = value
                means[direction] = (
                    (0 if b_values[direction] == P else 8) + 8 * value
                )
                roles[direction] = "fixed" if value == 0 else "M"
            high_count = roles.count("H")
            kind = (
                "residue4_high_plus_three_catalog"
                if high_count
                else "residue4_five_catalog"
            )
            leaves.append(
                make_leaf(
                    orbit_index,
                    residue,
                    kind,
                    tuple(means),
                    tuple(roles),
                    tuple(q_values),
                    tuple(increments),
                )
            )

        expected_per_orbit = 24 if 3 in z_by_type.values() else 44
        require(len(leaves) == expected_per_orbit, "per-orbit leaf count changed")
        for leaf in leaves:
            residue_histogram[leaf["residue"]] += 1
            kind_histogram[leaf["kind"]] += 1
            pattern = catalog_pattern(leaf["roles"], b_values)
            leaf["catalog_pattern"] = pattern
            pattern_histogram[pattern] += 1
        by_orbit.append(leaves)

    require(sum(map(len, by_orbit)) == EXPECTED_TOTAL_LEAVES,
            "total corrected leaf coverage changed")
    require(dict(sorted(residue_histogram.items())) == EXPECTED_RESIDUE_HISTOGRAM,
            "residue leaf census changed")
    require(dict(sorted(kind_histogram.items())) == dict(sorted(EXPECTED_KIND_HISTOGRAM.items())),
            "leaf kind census changed")
    require(dict(pattern_histogram) == EXPECTED_PATTERN_HISTOGRAM,
            f"catalog pattern census changed: {dict(pattern_histogram)}")
    return by_orbit, {
        "exact_mean_leaves": sum(map(len, by_orbit)),
        "residue_histogram": dict(sorted(residue_histogram.items())),
        "kind_histogram": dict(sorted(kind_histogram.items())),
        "catalog_pattern_histogram": {
            key: pattern_histogram[key] for key in EXPECTED_PATTERN_HISTOGRAM
        },
        "split_type_leaves_per_orbit": 44,
        "single_type_leaves_per_orbit": 24,
        "all_400_corrected_leaves_covered": True,
    }


def parity_row(b: int) -> np.ndarray:
    odd_fibres = set(range(b))
    return np.asarray(
        [sum(value in odd_fibres for value in point) & 1 for point in POINTS],
        dtype=np.int16,
    )


@functools.lru_cache(maxsize=None)
def canonical_catalog(b: int, mean: int) -> np.ndarray:
    """Complete phase-zero catalog in canonical fibre coordinates."""
    if b == P:
        rows = canonical_catalog(0, mean)
    elif b in (1, 5) and mean >= 8:
        rows = canonical_catalog(0, mean - 8) + parity_row(b)[None, :]
    else:
        rows = np.asarray(exact_slack_catalog_values(b, 0, mean), dtype=np.int16)
    rows = np.ascontiguousarray(rows, dtype=np.int16)
    require(rows.ndim == 2 and rows.shape[1] == 35, "catalog shape changed")
    require(len({row.tobytes() for row in rows}) == len(rows),
            "canonical catalog duplicates")
    require(np.all(2 * rows.sum(axis=1, dtype=np.int64) == 5 * mean),
            "canonical catalog mean changed")
    canonical_parity = parity_row(b)
    require(np.all((rows % 2) == canonical_parity[None, :]),
            "canonical catalog parity changed")
    require(np.all((0 <= rows) & (rows <= 13)),
            "catalog leaves affine bad-count range")
    return rows


@functools.lru_cache(maxsize=None)
def mapped_catalog(mask: int, mean: int) -> np.ndarray:
    odd_fibres = {index for index in range(P) if mask & (1 << index)}
    b = len(odd_fibres)
    canonical = canonical_catalog(b, mean)
    actual_odd = sorted(odd_fibres)
    actual_even = sorted(set(range(P)) - odd_fibres)
    permutation = dict(zip(range(b), actual_odd)) | dict(
        zip(range(b, P), actual_even)
    )
    inverse = {target: source for source, target in permutation.items()}
    source_columns = [
        POINT_INDEX[tuple(sorted(inverse[value] for value in point))]
        for point in POINTS
    ]
    mapped = np.ascontiguousarray(canonical[:, source_columns], dtype=np.int16)
    require(len({row.tobytes() for row in mapped}) == len(mapped),
            "mapped catalog duplicates")
    actual_parity = np.asarray(
        [sum(value in odd_fibres for value in point) & 1 for point in POINTS],
        dtype=np.int16,
    )
    require(np.all((mapped % 2) == actual_parity[None, :]),
            "mapped catalog parity changed")
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
    require(rows.ndim == 2 and rows.shape[1] > 0,
            "signature rows have bad shape")
    return rows.view(np.dtype((np.void, rows.shape[1]))).reshape(-1)


def summed_rows(contributions: tuple[np.ndarray, ...]) -> np.ndarray:
    """Materialize exact syndrome sums; matrices are coordinate-by-catalog."""
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


def iter_summed_rows(
    contributions: tuple[np.ndarray, ...],
) -> itertools.chain | itertools.product:
    """Yield bounded chunks containing every exact sum with multiplicity."""
    require(bool(contributions), "cannot stream an empty catalog side")
    dimension = contributions[0].shape[0]
    require(all(matrix.shape[0] == dimension for matrix in contributions),
            "streamed join dimensions disagree")
    ordered = tuple(
        sorted(contributions, key=lambda matrix: matrix.shape[1], reverse=True)
    )
    total = math.prod(matrix.shape[1] for matrix in ordered)
    if total <= MAX_JOIN_CHUNK_ROWS:
        yield summed_rows(ordered)
        return

    first = ordered[0]
    tail_matrices = ordered[1:]
    if not tail_matrices:
        for start in range(0, first.shape[1], MAX_JOIN_CHUNK_ROWS):
            yield np.ascontiguousarray(
                first[:, start : start + MAX_JOIN_CHUNK_ROWS].T,
                dtype=np.uint8,
            )
        return
    for tail in iter_summed_rows(tail_matrices):
        first_width = max(1, MAX_JOIN_CHUNK_ROWS // len(tail))
        for start in range(0, first.shape[1], first_width):
            first_rows = first[:, start : start + first_width].T
            rows = (
                first_rows[:, None, :].astype(np.int16)
                + tail[None, :, :].astype(np.int16)
            ) % MODULUS
            yield np.ascontiguousarray(
                rows.reshape(-1, dimension), dtype=np.uint8
            )


def matching_key_count(
    available_keys: np.ndarray,
    available_counts: np.ndarray,
    needed: np.ndarray,
) -> int:
    needed_keys = row_keys(needed)
    positions = np.searchsorted(available_keys, needed_keys)
    valid = positions < len(available_keys)
    if not np.any(valid):
        return 0
    valid_indices = np.flatnonzero(valid)
    equal = available_keys[positions[valid_indices]] == needed_keys[valid_indices]
    hits = valid_indices[equal]
    return int(np.sum(available_counts[positions[hits]], dtype=np.int64))


def choose_partition(
    sizes: tuple[int, ...], first_count: int
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    best = None
    for first in itertools.combinations(range(len(sizes)), first_count):
        second = tuple(index for index in range(len(sizes)) if index not in first)
        products = (
            math.prod(sizes[index] for index in first),
            math.prod(sizes[index] for index in second),
        )
        score = (max(products), min(products), first)
        if best is None or score < best[0]:
            best = (score, first, second)
    require(best is not None, "meet-in-the-middle partition failed")
    return best[1], best[2]


def exact_join_count(
    base: np.ndarray, contributions: tuple[np.ndarray, ...]
) -> tuple[int, dict]:
    """Count exact complete-catalog tuples summing with base to zero."""
    base = np.ascontiguousarray(base, dtype=np.uint8)
    require(base.ndim == 1 and len(base) > 0, "join base has bad shape")
    require(all(matrix.shape[0] == len(base) for matrix in contributions),
            "catalog contribution dimension changed")
    sizes = tuple(int(matrix.shape[1]) for matrix in contributions)
    count_catalogs = len(contributions)

    if count_catalogs == 0:
        count = int(not np.any(base))
        return count, {
            "partition": [[], []],
            "catalog_sizes": [],
            "enumerated_signatures": [1, 1],
            "cartesian_catalog_tuples": 1,
            "streamed_side_peak_rows": 1,
            "exact_injective_signature_encoding": f"{len(base)} raw base-seven bytes",
        }

    if count_catalogs == 1:
        available_keys, available_counts = np.unique(
            row_keys(contributions[0].T), return_counts=True
        )
        needed = np.ascontiguousarray(
            (-base[None, :].astype(np.int16)) % MODULUS, dtype=np.uint8
        )
        count = matching_key_count(available_keys, available_counts, needed)
        return count, {
            "partition": [[0], []],
            "catalog_sizes": list(sizes),
            "enumerated_signatures": [sizes[0], 1],
            "cartesian_catalog_tuples": sizes[0],
            "streamed_side_peak_rows": 1,
            "exact_injective_signature_encoding": f"{len(base)} raw base-seven bytes",
        }

    if count_catalogs == 3:
        first_indices, second_indices = choose_partition(sizes, 1)
    elif count_catalogs == 5:
        first_indices, second_indices = choose_partition(sizes, 2)
    else:
        raise AssertionError(
            f"expected zero, one, three, or five catalogs, got {count_catalogs}"
        )

    first_product = math.prod(sizes[index] for index in first_indices)
    second_product = math.prod(sizes[index] for index in second_indices)
    if first_product > second_product:
        first_indices, second_indices = second_indices, first_indices
        first_product, second_product = second_product, first_product
    first_chunks = tuple(
        iter_summed_rows(tuple(contributions[index] for index in first_indices))
    )
    first = np.ascontiguousarray(np.concatenate(first_chunks, axis=0), dtype=np.uint8)
    require(len(first) == first_product, "materialized join-side size changed")
    available_keys, available_counts = np.unique(row_keys(first), return_counts=True)

    count = 0
    peak = 0
    seen = 0
    for second in iter_summed_rows(
        tuple(contributions[index] for index in second_indices)
    ):
        peak = max(peak, len(second))
        seen += len(second)
        needed = np.ascontiguousarray(
            (-base[None, :].astype(np.int16) - second.astype(np.int16)) % MODULUS,
            dtype=np.uint8,
        )
        count += matching_key_count(available_keys, available_counts, needed)
    require(seen == second_product, "streamed join-side size changed")
    return count, {
        "partition": [list(first_indices), list(second_indices)],
        "partition_shape": [len(first_indices), len(second_indices)],
        "catalog_sizes": list(sizes),
        "enumerated_signatures": [first_product, second_product],
        "cartesian_catalog_tuples": math.prod(sizes),
        "streamed_side_peak_rows": peak,
        "exact_injective_signature_encoding": f"{len(base)} raw base-seven bytes",
    }


def join_implementation_self_audit() -> dict:
    """Compare all production join shapes with direct brute force."""
    zero = np.asarray([0, 0], dtype=np.uint8)
    cases: tuple[tuple[np.ndarray, tuple[np.ndarray, ...]], ...] = (
        (zero, ()),
        (zero, (np.asarray([[0, 1, 0], [0, 0, 0]], dtype=np.uint8),)),
        (
            zero,
            (
                np.asarray([[0, 0], [0, 0]], dtype=np.uint8),
                np.asarray([[0, 1, 1], [0, 0, 0]], dtype=np.uint8),
                np.asarray([[0, 6], [0, 0]], dtype=np.uint8),
            ),
        ),
        (
            zero,
            (
                np.asarray([[0, 0], [0, 0]], dtype=np.uint8),
                np.asarray([[0, 1], [0, 0]], dtype=np.uint8),
                np.asarray([[0, 6], [0, 0]], dtype=np.uint8),
                np.asarray([[0, 0], [0, 1]], dtype=np.uint8),
                np.asarray([[0, 0], [0, 6]], dtype=np.uint8),
            ),
        ),
    )
    rows = []
    for base, matrices in cases:
        joined, metadata = exact_join_count(base, matrices)
        brute = 0
        sizes = tuple(matrix.shape[1] for matrix in matrices)
        assignments = itertools.product(*(range(size) for size in sizes))
        for indices in assignments:
            syndrome = base.astype(np.int16).copy()
            for matrix, index in zip(matrices, indices):
                syndrome += matrix[:, index]
            brute += int(not np.any(syndrome % MODULUS))
        require(joined == brute, "meet-in-the-middle join failed brute-force audit")
        require(joined > 0, "join self-audit did not exercise a successful match")
        rows.append({"catalog_count": len(matrices), "join_count": joined, **metadata})
    return {
        "passed": True,
        "catalog_counts_exercised": [0, 1, 3, 5],
        "positive_matches_exercised": True,
        "duplicate_multiplicity_exercised": True,
        "cases": rows,
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
        require(modular_rank(conditioned) == 112,
                "conditioned dependencies lost rank")
        require(not np.any(conditioned[:, columns]),
                "conditioned direction block is not zero")
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
        "conditioned_dependency_dimension": 112,
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
            require(np.all(bad_counts >= 0),
                    "catalog score exceeds affine edge range")
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
    matching_tuples: int,
    join: dict,
) -> dict:
    high_direction = leaf["high_direction"]
    row = {
        "orbit_index": orbit_index,
        "source_orbit_index": int(orbit["source_orbit_index"]),
        "orbit_leaf_index": leaf_index,
        "orbit_size": int(orbit["size"]),
        "representative_finite_field": list(orbit["representative"]),
        "direction_types": list(DIRECTION_TYPES),
        "b_values": list(orbit["b_values"]),
        "undetermined_directions": list(orbit["undetermined"]),
        "undetermined_type_counts": {
            str(key): value
            for key, value in sorted(orbit["undetermined_type_counts"].items())
        },
        "residue": leaf["residue"],
        "leaf_kind": leaf["kind"],
        "catalog_pattern": leaf["catalog_pattern"],
        "q_values": None if leaf["q_values"] is None else list(leaf["q_values"]),
        "residue_increments": (
            None
            if leaf["residue_increments"] is None
            else list(leaf["residue_increments"])
        ),
        "scaled_means": list(leaf["means"]),
        "catalog_roles": list(leaf["roles"]),
        "variable_directions": list(leaf["variable_directions"]),
        "high_direction": high_direction,
        "high_direction_block_relaxed": high_direction is not None,
        "mod7_test": test,
        "matching_remaining_catalog_tuples": matching_tuples,
        "weighted_boundary_allocation_cases": int(orbit["size"]),
        "join": join,
    }
    if high_direction is None:
        row["exact_mod7_catalog_tuples"] = matching_tuples
        row["weighted_exact_mod7_catalog_tuples"] = (
            int(orbit["size"]) * matching_tuples
        )
    return row


def run(orbit_limit: int = 0) -> dict:
    started = time.time()
    require(Counter(DIRECTION_TYPES) == Counter({-1: 4, 1: 4}),
            "direction types changed")
    _matrix, dependencies, linear = translation_equivariant_system()
    orbits, orbit_audit = evidence_z3_orbits()
    leaves_by_orbit, leaf_audit = exact_mean_leaves(orbits)
    catalogs = expected_catalog_audit()
    join_audit = join_implementation_self_audit()
    annihilators, annihilator_audit = direction_annihilators(dependencies)

    if orbit_limit < 0:
        raise ValueError("--orbit-limit must be nonnegative")
    process_count = len(orbits) if orbit_limit == 0 else min(orbit_limit, len(orbits))
    full_run = process_count == len(orbits)
    base_edge_syndrome = (
        dependencies[:, 0] * EDGE_COUNT % MODULUS
    ).astype(np.uint8)
    factory = ContributionFactory(dependencies)
    survivors = []
    decision_digest = hashlib.sha256()
    processed_kind_histogram: Counter[str] = Counter()
    rejected_kind_histogram: Counter[str] = Counter()
    processed_residue_histogram: Counter[str] = Counter()
    processed_pattern_histogram: Counter[str] = Counter()
    processed_weighted = 0
    rejected_weighted = 0
    surviving_weighted = 0
    weighted_exact_mod7_tuples = 0
    exact_cartesian_tuples = 0
    conditioned_cartesian_tuples = 0
    per_orbit = []

    for orbit_index, (orbit, leaves) in enumerate(
        zip(orbits[:process_count], leaves_by_orbit[:process_count])
    ):
        orbit_survivors = 0
        for leaf_index, leaf in enumerate(leaves):
            kind = leaf["kind"]
            roles = leaf["roles"]
            processed_kind_histogram[kind] += 1
            processed_residue_histogram[leaf["residue"]] += 1
            processed_pattern_histogram[leaf["catalog_pattern"]] += 1
            processed_weighted += int(orbit["size"])

            fixed = base_edge_syndrome.astype(np.int16).copy()
            for direction, (mask, mean, role) in enumerate(
                zip(orbit["masks"], leaf["means"], roles)
            ):
                if role != "fixed":
                    continue
                contribution = factory.get(direction, mask, mean)
                require(contribution.shape[1] == 1,
                        "fixed direction catalog is not unique")
                fixed += contribution[:, 0]
            fixed = np.ascontiguousarray(fixed % MODULUS, dtype=np.uint8)

            contributions = tuple(
                factory.get(direction, orbit["masks"][direction], leaf["means"][direction])
                for direction in leaf["variable_directions"]
            )
            expected_sizes = tuple(
                catalog_size(roles[direction], orbit["b_values"][direction])
                for direction in leaf["variable_directions"]
            )
            require(tuple(matrix.shape[1] for matrix in contributions) == expected_sizes,
                    "variable complete-catalog size changed")

            high_direction = leaf["high_direction"]
            if high_direction is None:
                require(len(contributions) in (3, 5),
                        "no-high leaf has the wrong catalog count")
                matching_tuples, join_metadata = exact_join_count(fixed, contributions)
                passing = matching_tuples > 0
                test = (
                    "exact_full_135_three_catalog_1plus2_join"
                    if len(contributions) == 3
                    else "exact_full_135_five_catalog_balanced_2plus3_join"
                )
                exact_cartesian_tuples += int(join_metadata["cartesian_catalog_tuples"])
            else:
                require(roles[high_direction] == "H",
                        "high direction role mismatch")
                require(len(contributions) in (0, 1, 3),
                        "high leaf has the wrong remaining catalog count")
                annihilator = annihilators[high_direction]
                conditioned_base = np.ascontiguousarray(
                    annihilator @ fixed.astype(np.int64) % MODULUS,
                    dtype=np.uint8,
                )
                conditioned_contributions = tuple(
                    np.ascontiguousarray(
                        annihilator @ contribution.astype(np.int64) % MODULUS,
                        dtype=np.uint8,
                    )
                    for contribution in contributions
                )
                matching_tuples, join_metadata = exact_join_count(
                    conditioned_base, conditioned_contributions
                )
                passing = matching_tuples > 0
                test = (
                    f"full_112_dependency_high_direction_annihilation_then_exact_"
                    f"{len(contributions)}_catalog_test"
                )
                conditioned_cartesian_tuples += int(
                    join_metadata["cartesian_catalog_tuples"]
                )

            digest_row = {
                "orbit": orbit_index,
                "leaf": leaf_index,
                "residue": leaf["residue"],
                "pattern": leaf["catalog_pattern"],
                "means": leaf["means"],
                "passing": passing,
                "matching_tuples": matching_tuples,
            }
            decision_digest.update(
                json.dumps(digest_row, sort_keys=True, separators=(",", ":")).encode()
            )
            decision_digest.update(b"\n")
            if passing:
                orbit_survivors += 1
                surviving_weighted += int(orbit["size"])
                if high_direction is None:
                    weighted_exact_mod7_tuples += (
                        int(orbit["size"]) * matching_tuples
                    )
                survivors.append(
                    render_case(
                        orbit_index,
                        orbit,
                        leaf_index,
                        leaf,
                        test,
                        matching_tuples,
                        join_metadata,
                    )
                )
            else:
                rejected_kind_histogram[kind] += 1
                rejected_weighted += int(orbit["size"])

        per_orbit.append(
            {
                "orbit_index": orbit_index,
                "source_orbit_index": int(orbit["source_orbit_index"]),
                "orbit_size": int(orbit["size"]),
                "majority_quadratic_type": int(orbit["majority_type"]),
                "exact_mean_leaves": len(leaves),
                "surviving_leaves": orbit_survivors,
            }
        )

    processed_leaves = sum(processed_kind_histogram.values())
    require(
        processed_leaves == sum(len(rows) for rows in leaves_by_orbit[:process_count]),
        "processed leaf count mismatch",
    )
    require(processed_weighted == rejected_weighted + surviving_weighted,
            "weighted decision census mismatch")
    if full_run:
        require(processed_leaves == EXPECTED_TOTAL_LEAVES,
                "full run missed corrected exact mean leaves")
        require(dict(processed_kind_histogram) == EXPECTED_KIND_HISTOGRAM,
                "full processed kind census changed")
        require(dict(processed_residue_histogram) == EXPECTED_RESIDUE_HISTOGRAM,
                "full processed residue census changed")
        require(dict(processed_pattern_histogram) == EXPECTED_PATTERN_HISTOGRAM,
                "full processed pattern census changed")

    excluded = full_run and not survivors
    return {
        "experiment": "p7_infinity7_positive_z3_mod7_join",
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
        "z": 3,
        "phase": 0,
        "linear_system": linear,
        "orbit_reduction": orbit_audit,
        "mean_leaf_coverage": leaf_audit,
        "catalog_row_counts": catalogs,
        "catalog_source": "complete exact Johnson-slice catalogs",
        "join_self_audit": join_audit,
        "high_direction_annihilator_audit": annihilator_audit,
        "affine_span_relaxation_used": False,
        "high_full_direction_block_image_relaxation_used": True,
        "modular_passing_is_edge_feasibility": False,
        "processed_orbits": process_count,
        "full_orbit_count": len(orbits),
        "full_run": full_run,
        "processed_exact_mean_leaves": processed_leaves,
        "processed_kind_histogram": dict(sorted(processed_kind_histogram.items())),
        "processed_residue_histogram": dict(sorted(processed_residue_histogram.items())),
        "processed_catalog_pattern_histogram": {
            key: processed_pattern_histogram[key]
            for key in EXPECTED_PATTERN_HISTOGRAM
            if processed_pattern_histogram[key]
        },
        "rejected_kind_histogram": dict(sorted(rejected_kind_histogram.items())),
        "processed_weighted_boundary_allocation_cases": processed_weighted,
        "rejected_weighted_boundary_allocation_cases": rejected_weighted,
        "surviving_cases": len(survivors),
        "surviving_weighted_boundary_allocation_cases": surviving_weighted,
        "weighted_exact_mod7_catalog_tuples": weighted_exact_mod7_tuples,
        "exact_full_coordinate_cartesian_catalog_tuples_joined_on_representatives": exact_cartesian_tuples,
        "conditioned_cartesian_catalog_tuples_tested_on_representatives": conditioned_cartesian_tuples,
        "z3_branch_excluded": excluded,
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
        help="process only the first N of 10 orbits; zero means the full run",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    out = run(args.orbit_limit)
    atomic_write(args.output, out)
    if not args.quiet:
        summary = {
            key: value
            for key, value in out.items()
            if key not in {"per_orbit_summary", "survivor_cases"}
        }
        print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
