#!/usr/bin/env python3
"""Affine-span mod-seven filter for p=7 negative infinity-plus-five cases.

Each quadratic direction type has exactly one scaled-mean-14 direction and
three scaled-mean-6 directions.  The two nontrivial mean-14 catalogs are
relaxed to their affine spans over F_7.  Inconsistency of this larger linear
set is a rigorous exclusion of every exact catalog tuple.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import sys
import time
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
)
from p7_unsaturated_mod7_batch import contribution_matrix  # noqa: E402
from p7_unsaturated_modular_catalog_filter import (  # noqa: E402
    equation_matrix,
    left_dependencies,
)


MODULUS = 7


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def add_to_basis(basis: dict[int, np.ndarray], row: np.ndarray) -> bool:
    row = np.asarray(row, dtype=np.int16).copy() % MODULUS
    for pivot in sorted(basis):
        if row[pivot]:
            row = (row - int(row[pivot]) * basis[pivot]) % MODULUS
    nonzero = np.flatnonzero(row)
    if not len(nonzero):
        return False
    pivot = int(nonzero[0])
    row = row * pow(int(row[pivot]), -1, MODULUS) % MODULUS
    for old_pivot, old_row in list(basis.items()):
        if old_row[pivot]:
            basis[old_pivot] = (
                old_row - int(old_row[pivot]) * row
            ) % MODULUS
    basis[pivot] = row
    return True


def reduced_basis(rows) -> tuple[tuple[int, ...], np.ndarray]:
    basis: dict[int, np.ndarray] = {}
    for row in rows:
        add_to_basis(basis, row)
    pivots = tuple(sorted(basis))
    matrix = np.stack([basis[pivot] for pivot in pivots]) if pivots else np.empty((0, 135), dtype=np.int16)
    return pivots, matrix.astype(np.int16)


def in_span(target: np.ndarray, pivots: tuple[int, ...], basis: np.ndarray) -> bool:
    row = np.asarray(target, dtype=np.int16).copy() % MODULUS
    for pivot, vector in zip(pivots, basis):
        if row[pivot]:
            row = (row - int(row[pivot]) * vector) % MODULUS
    return not bool(np.any(row))


def run(limit: int | None = None, source: Path | None = None) -> dict:
    started = time.time()
    matrix = equation_matrix()
    rank, dependencies = left_dependencies(matrix, MODULUS)
    if rank != 147 or dependencies.shape != (135, 282):
        raise AssertionError("unexpected dependency dimensions")
    labels = []
    epsilons = []
    for direction in projective_directions(7):
        eps, row = field_direction_data(7, direction)
        epsilons.append(int(eps))
        labels.append(tuple(int(value) for value in row))
    type_directions = {
        eps: tuple(index for index, value in enumerate(epsilons) if value == eps)
        for eps in (-1, 1)
    }
    edge_base = (
        dependencies[:, :2] @ np.asarray([29, 1], dtype=np.int64) % MODULUS
    ).astype(np.int16)

    minimum_cache: dict[tuple[int, int], np.ndarray] = {}
    variable_cache: dict[tuple[int, int], tuple[np.ndarray, tuple[int, ...], np.ndarray, int]] = {}
    pair_cache: dict[tuple[tuple[int, int], tuple[int, int]], tuple[tuple[int, ...], np.ndarray]] = {}

    def B_from_mask(mask: int) -> set[int]:
        return {value for value in range(7) if (mask >> value) & 1}

    def minimum(direction: int, mask: int) -> np.ndarray:
        key = (direction, mask)
        if key not in minimum_cache:
            b = mask.bit_count()
            if b not in (1, 5):
                raise ValueError("only b=1,5 have a mean-six minimum")
            minimum_cache[key] = contribution_matrix(
                dependencies, direction, b, 1, 6, B_from_mask(mask)
            )[:, 0].astype(np.int16)
        return minimum_cache[key]

    def variable(direction: int, mask: int):
        key = (direction, mask)
        if key not in variable_cache:
            b = mask.bit_count()
            catalog = contribution_matrix(
                dependencies, direction, b, 1, 14, B_from_mask(mask)
            ).astype(np.int16)
            representative = catalog[:, 0]
            pivots, basis = reduced_basis(
                (catalog[:, index] - representative) % MODULUS
                for index in range(1, catalog.shape[1])
            )
            variable_cache[key] = (
                representative,
                pivots,
                basis,
                int(catalog.shape[1]),
            )
        return variable_cache[key]

    def pair_basis(key1: tuple[int, int], key2: tuple[int, int]):
        pair_key = tuple(sorted((key1, key2)))
        if pair_key not in pair_cache:
            _r1, _p1, b1, _n1 = variable(*pair_key[0])
            _r2, _p2, b2, _n2 = variable(*pair_key[1])
            pair_cache[pair_key] = reduced_basis(itertools.chain(b1, b2))
        return pair_cache[pair_key]

    checked_boundaries = 0
    floor_rejected = 0
    elevation_cases = 0
    affine_span_surviving_cases = 0
    surviving_boundaries = []
    surviving_cases = []
    catalog_pattern_counts: dict[str, int] = {}
    if source is None:
        boundary_rows = (
            (finite, 1) for finite in itertools.combinations(range(49), 5)
        )
        source_orbits = None
    else:
        payload = json.loads(source.read_text())
        if (
            int(payload["p"]) != 7
            or int(payload["c_H"]) != -1
            or int(payload["boundary_size"]) != 6
            or int(payload["infinity_value"]) != 1
        ):
            raise ValueError("orbit source has the wrong branch scope")
        boundary_rows = (
            (
                tuple(int(value) for value in orbit["representative_finite_field"]),
                int(orbit["size"]),
            )
            for orbit in payload["orbits"]
        )
        source_orbits = len(payload["orbits"])
    checked_boundary_weight = 0
    surviving_boundary_weight = 0
    for finite, orbit_size in boundary_rows:
        if limit is not None and checked_boundaries >= limit:
            break
        masks = []
        for direction in range(8):
            mask = 0
            for point in finite:
                mask ^= 1 << labels[direction][point]
            if mask.bit_count() not in (1, 3, 5):
                raise AssertionError("five points produced an even fibre mask")
            masks.append(mask)
        options = {}
        invalid = False
        for eps in (-1, 1):
            b3 = [d for d in type_directions[eps] if masks[d].bit_count() == 3]
            if len(b3) > 1:
                invalid = True
                break
            options[eps] = tuple(b3 or type_directions[eps])
        checked_boundaries += 1
        checked_boundary_weight += orbit_size
        if invalid:
            floor_rejected += 1
            continue

        base = edge_base.copy()
        for direction, mask in enumerate(masks):
            if mask.bit_count() in (1, 5):
                base = (base + minimum(direction, mask)) % MODULUS
        boundary_survives = False
        for d_minus in options[-1]:
            for d_plus in options[1]:
                elevation_cases += 1
                keys = ((d_minus, masks[d_minus]), (d_plus, masks[d_plus]))
                case = base.copy()
                pattern = []
                for direction, mask in keys:
                    representative, _pivots, _basis, count = variable(direction, mask)
                    if mask.bit_count() in (1, 5):
                        case = (case - minimum(direction, mask)) % MODULUS
                    case = (case + representative) % MODULUS
                    pattern.append(count)
                pattern_name = "x".join(map(str, sorted(pattern, reverse=True)))
                catalog_pattern_counts[pattern_name] = catalog_pattern_counts.get(pattern_name, 0) + 1
                pivots, basis = pair_basis(*keys)
                if in_span((-case) % MODULUS, pivots, basis):
                    affine_span_surviving_cases += 1
                    boundary_survives = True
                    surviving_cases.append(
                        {
                            "finite_boundary": list(finite),
                            "orbit_size": orbit_size,
                            "elevated_directions": [d_minus, d_plus],
                            "direction_masks": list(masks),
                            "catalog_sizes": pattern,
                            "catalog_pattern": pattern_name,
                        }
                    )
        if boundary_survives:
            surviving_boundaries.append(list(finite))
            surviving_boundary_weight += orbit_size
        if checked_boundaries % 1000 == 0:
            print(
                f"checked={checked_boundaries} floor_rejected={floor_rejected} "
                f"span_survivors={affine_span_surviving_cases}",
                flush=True,
            )

    return {
        "experiment": "p7_size6_negative_infinity_mod7_affine",
        "status": "exact_mod_seven_affine_catalog_span_filter",
        "p": 7,
        "c_H": -1,
        "infinity_value": 1,
        "checked_boundaries": checked_boundaries,
        "checked_boundary_weight": checked_boundary_weight,
        "requested_limit": limit,
        "source": str(source) if source is not None else None,
        "source_sha256": sha256(source) if source is not None else None,
        "source_orbits": source_orbits,
        "complete_boundary_enumeration": limit is None and source is None,
        "complete_orbit_enumeration": limit is None and source is not None,
        "floor_rejected_boundaries": floor_rejected,
        "elevation_cases": elevation_cases,
        "catalog_pattern_counts": catalog_pattern_counts,
        "variable_catalog_keys": len(variable_cache),
        "pair_span_keys": len(pair_cache),
        "affine_span_surviving_cases": affine_span_surviving_cases,
        "affine_span_surviving_boundaries": len(surviving_boundaries),
        "affine_span_surviving_boundary_weight": surviving_boundary_weight,
        "surviving_boundaries_finite_field": surviving_boundaries,
        "surviving_cases": surviving_cases,
        "all_checked_cases_affine_span_infeasible": affine_span_surviving_cases == 0,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = run(args.limit, args.source)
    if args.output is not None:
        atomic_write(args.output, out)
    print(json.dumps({key: value for key, value in out.items() if key != "surviving_boundaries_finite_field"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
