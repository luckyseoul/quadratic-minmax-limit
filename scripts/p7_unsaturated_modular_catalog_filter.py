#!/usr/bin/env python3
"""Odd-prime linear-consistency filter for p=7 slack catalogs.

For a fully fixed directional target, each of the 280 affine scores gives
an exact linear equation ``sum(bad edges) = (29-score)/2``.  Edge count and
the distinguished edge give two more equations.  Reducing this common
282-by-1225 system modulo odd primes exposes any catalog rows whose right
hand sides violate a linear dependency.  Passing is only a necessary
condition; rejection is a rigorous finite exclusion of that catalog row.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import field_direction_data, projective_directions  # noqa: E402
from p7_no_infinity_unsaturated_cpsat import atomic_write, direction_target_options  # noqa: E402
from p7_unsaturated_gf2_catalog_filter import direction_scope  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


POINTS = tuple(itertools.combinations(range(7), 4))
PAIRS = tuple(itertools.combinations(range(7), 2))


def equation_matrix() -> np.ndarray:
    """Return edge count, fixed edge, and all 280 bad-edge equations."""
    data = geometry(7, "affine")
    C = data["C"]
    edges = data["edges"]
    rows: list[np.ndarray] = []

    rows.append(np.ones(len(edges), dtype=np.int16))
    fixed = np.zeros(len(edges), dtype=np.int16)
    fixed[edges.index((0, 1))] = 1
    rows.append(fixed)
    for direction in projective_directions(7):
        eps, labels = field_direction_data(7, direction)
        for X in POINTS:
            X_set = set(X)
            bad = np.zeros(len(edges), dtype=np.int16)
            for edge_index, (a, endpoint) in enumerate(edges):
                y_a = eps if a == 0 else (1 if labels[a - 1] in X_set else -1)
                y_b = 1 if labels[endpoint - 1] in X_set else -1
                if eps * y_a * y_b * int(C[a, endpoint]) < 0:
                    bad[edge_index] = 1
            rows.append(bad)
    matrix = np.stack(rows)
    if matrix.shape != (282, 1225):
        raise AssertionError(f"unexpected equation shape {matrix.shape}")
    return matrix


def left_dependencies(matrix: np.ndarray, modulus: int) -> tuple[int, np.ndarray]:
    """Row-reduce over F_modulus and return rank and left-null witnesses."""
    if modulus < 3:
        raise ValueError("this filter requires an odd prime modulus")
    A = matrix.astype(np.int64) % modulus
    m, n = A.shape
    transform = np.eye(m, dtype=np.int64)
    rank = 0
    for column in range(n):
        candidates = np.flatnonzero(A[rank:, column])
        if not len(candidates):
            continue
        pivot = rank + int(candidates[0])
        if pivot != rank:
            A[[rank, pivot]] = A[[pivot, rank]]
            transform[[rank, pivot]] = transform[[pivot, rank]]
        inverse = pow(int(A[rank, column]), -1, modulus)
        A[rank] = (A[rank] * inverse) % modulus
        transform[rank] = (transform[rank] * inverse) % modulus
        factors = A[:, column].copy()
        factors[rank] = 0
        active = np.flatnonzero(factors)
        if len(active):
            A[active] = (
                A[active] - factors[active, None] * A[rank, None]
            ) % modulus
            transform[active] = (
                transform[active]
                - factors[active, None] * transform[rank, None]
            ) % modulus
        rank += 1
        if rank == m:
            break
    if np.any(A[rank:]):
        raise AssertionError("row reduction did not leave zero dependency rows")
    dependencies = transform[rank:] % modulus
    if len(dependencies) != m - rank:
        raise AssertionError("dependency dimension mismatch")
    return rank, dependencies


def target_bad_counts(option: tuple[int, ...]) -> np.ndarray:
    constant = int(option[1])
    coefficients = tuple(int(value) for value in option[2:])
    out = []
    for X in POINTS:
        X_set = set(X)
        score = constant + sum(
            coefficients[pair_index]
            * (1 if ((s in X_set) == (t in X_set)) else -1)
            for pair_index, (s, t) in enumerate(PAIRS)
        )
        if score < 3 or score > 29 or score % 2 == 0:
            raise AssertionError(f"invalid normalized catalog score {score}")
        out.append((29 - score) // 2)
    return np.asarray(out, dtype=np.int64)


def scan_catalog(
    c_h: int,
    boundary: tuple[int, ...],
    elevated: tuple[int, ...],
    scan_direction: int,
    fixed_indices: dict[int, int],
    moduli: tuple[int, ...],
) -> dict:
    started = time.time()
    direction_rows, type_floors = direction_scope(c_h, boundary, elevated)
    if scan_direction not in elevated:
        raise ValueError("scan direction must be elevated")
    elevated_set = set(elevated)
    option_tables = []
    for direction_index, row in enumerate(direction_rows):
        options = direction_target_options(
            len(row["B"]),
            int(row["phase"]),
            set(row["B"]),
            type_floors[int(row["eps"])],
            direction_index in elevated_set,
        )
        if direction_index != scan_direction and len(options) > 1:
            if direction_index not in fixed_indices:
                raise ValueError(
                    f"direction {direction_index} has {len(options)} choices; "
                    "provide --fixed-index"
                )
            fixed_index = int(fixed_indices[direction_index])
            if not 0 <= fixed_index < len(options):
                raise ValueError(f"fixed index outside direction {direction_index}")
        option_tables.append(options)

    scan_options = option_tables[scan_direction]
    right_sides = np.empty((282, len(scan_options)), dtype=np.int64)
    right_sides[0] = 29
    right_sides[1] = 1
    for direction_index, options in enumerate(option_tables):
        row_slice = slice(2 + 35 * direction_index, 2 + 35 * (direction_index + 1))
        if direction_index == scan_direction:
            right_sides[row_slice] = np.stack(
                [target_bad_counts(option) for option in scan_options], axis=1
            )
        else:
            option_index = int(fixed_indices.get(direction_index, 0))
            right_sides[row_slice] = target_bad_counts(options[option_index])[:, None]

    matrix = equation_matrix()
    survivors = np.ones(len(scan_options), dtype=bool)
    rows = []
    for modulus in moduli:
        rank, dependencies = left_dependencies(matrix, modulus)
        if len(dependencies):
            syndromes = dependencies @ (right_sides % modulus) % modulus
            passes = np.all(syndromes == 0, axis=0)
        else:
            passes = np.ones(len(scan_options), dtype=bool)
        survivors &= passes
        rows.append(
            {
                "modulus": modulus,
                "rank": rank,
                "dependency_dimension": int(len(dependencies)),
                "passing_catalog_rows": int(np.count_nonzero(passes)),
                "cumulative_survivors": int(np.count_nonzero(survivors)),
            }
        )

    survivor_indices = np.flatnonzero(survivors).tolist()
    return {
        "experiment": "p7_unsaturated_modular_catalog_filter",
        "status": "complete_exact_odd_prime_linear_consistency_scan",
        "p": 7,
        "c_H": c_h,
        "fixed_boundary": list(boundary),
        "fixed_elevated_directions": list(elevated),
        "scan_direction": scan_direction,
        "catalog_total": len(scan_options),
        "equations": int(matrix.shape[0]),
        "edge_variables": int(matrix.shape[1]),
        "moduli": rows,
        "surviving_catalog_rows": len(survivor_indices),
        "rejected_catalog_rows": len(scan_options) - len(survivor_indices),
        "surviving_catalog_indices": (
            "ALL" if len(survivor_indices) == len(scan_options) else survivor_indices
        ),
        "fixed_indices": {
            str(direction): index for direction, index in sorted(fixed_indices.items())
        },
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c-h", type=int, choices=(-1, 1), required=True)
    parser.add_argument("--fixed-boundary", type=int, nargs=4, required=True)
    parser.add_argument("--elevated-directions", type=int, nargs="+", required=True)
    parser.add_argument("--scan-direction", type=int, required=True)
    parser.add_argument("--moduli", type=int, nargs="+", default=(3, 5, 7, 11, 13))
    parser.add_argument(
        "--fixed-index",
        type=int,
        nargs=2,
        action="append",
        metavar=("DIRECTION", "INDEX"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    fixed_indices = {}
    for direction, index in args.fixed_index or []:
        if direction in fixed_indices:
            raise ValueError(f"duplicate fixed index for direction {direction}")
        fixed_indices[direction] = index
    out = scan_catalog(
        args.c_h,
        tuple(sorted(args.fixed_boundary)),
        tuple(sorted(args.elevated_directions)),
        args.scan_direction,
        fixed_indices,
        tuple(args.moduli),
    )
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        atomic_write(args.output, out)


if __name__ == "__main__":
    main()
