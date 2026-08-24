#!/usr/bin/env python3
"""Search low-degree affine relations in normalized d=21 boundary scalars.

For six scalar coordinates over F_8, solve

    1 + sum_m c_m m(x_1,...,x_6) = 0

on every normalized census row, using all nonconstant monomials through a
requested total degree.  A solution has value one at the simultaneous zero
triple and would therefore exclude that obstruction.  Searches with fewer
features than rows are non-interpolating; any candidate is still evidence
until derived symbolically and validated on a disjoint range.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT / "src"))

from io_atomic import write_json_atomic  # noqa: E402


MODULUS = 0x57


def multiply(left: int, right: int) -> int:
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left >> 6:
            left ^= MODULUS
    return result


def power(value: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply(result, value)
        value = multiply(value, value)
        exponent >>= 1
    return result


def exponent_vectors(n_variables: int, max_degree: int):
    for degree in range(1, max_degree + 1):
        for indices in itertools.combinations_with_replacement(
            range(n_variables), degree
        ):
            exponents = [0] * n_variables
            for index in indices:
                exponents[index] += 1
            yield tuple(exponents)


def solve_affine(
    features: np.ndarray, multiplication: np.ndarray, inverses: np.ndarray
) -> tuple[list[int] | None, int, int]:
    n_rows, n_columns = features.shape
    matrix = np.concatenate(
        [features.copy(), np.ones((n_rows, 1), dtype=np.uint8)], axis=1
    )
    pivot_columns = []
    pivot_row = 0
    for column in range(n_columns):
        candidates = np.flatnonzero(matrix[pivot_row:, column])
        if not len(candidates):
            continue
        found = pivot_row + int(candidates[0])
        if found != pivot_row:
            matrix[[pivot_row, found]] = matrix[[found, pivot_row]]
        inverse = inverses[matrix[pivot_row, column]]
        matrix[pivot_row] = multiplication[inverse, matrix[pivot_row]]
        active = np.flatnonzero(matrix[:, column])
        active = active[active != pivot_row]
        if len(active):
            scales = matrix[active, column]
            matrix[active] ^= multiplication[
                scales[:, None], matrix[pivot_row][None, :]
            ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == n_rows:
            break
    inconsistent = any(
        not np.any(row[:n_columns]) and row[n_columns]
        for row in matrix[pivot_row:]
    )
    rank = len(pivot_columns)
    augmented_rank = rank + int(inconsistent)
    if inconsistent:
        return None, rank, augmented_rank
    solution = [0] * n_columns
    for row, column in enumerate(pivot_columns):
        solution[column] = int(matrix[row, n_columns])
    return solution, rank, augmented_rank


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--max-degree", type=int, default=5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.max_degree < 1:
        raise ValueError("max-degree must be positive")
    data = json.loads(args.input.read_text())
    values = np.asarray(
        [
            left + right
            for row in data["rows"]
            for left, right in [row["scalar_triples_common_f8"]]
        ],
        dtype=np.uint8,
    )
    multiplication = np.asarray(
        [[multiply(left, right) for right in range(64)] for left in range(64)],
        dtype=np.uint8,
    )
    inverses = np.asarray(
        [0] + [power(value, 62) for value in range(1, 64)], dtype=np.uint8
    )
    powers = np.ones((len(values), 6, args.max_degree + 1), dtype=np.uint8)
    for exponent in range(1, args.max_degree + 1):
        powers[:, :, exponent] = multiplication[
            powers[:, :, exponent - 1], values
        ]

    results = []
    for degree in range(1, args.max_degree + 1):
        exponents = list(exponent_vectors(6, degree))
        columns = []
        for vector in exponents:
            column = np.ones(len(values), dtype=np.uint8)
            for variable, exponent in enumerate(vector):
                if exponent:
                    column = multiplication[
                        column, powers[:, variable, exponent]
                    ]
            columns.append(column)
        features = np.column_stack(columns)
        solution, rank, augmented_rank = solve_affine(
            features, multiplication, inverses
        )
        candidate = solution is not None
        row = {
            "degree": degree,
            "rows": len(values),
            "features": len(exponents),
            "rank": rank,
            "augmented_rank": augmented_rank,
            "excludes_simultaneous_zero_candidate": candidate,
        }
        if candidate:
            row["terms"] = [
                {"exponents": list(exponent), "coefficient": coefficient}
                for exponent, coefficient in zip(exponents, solution)
                if coefficient
            ]
        results.append(row)
        print(json.dumps(row, separators=(",", ":")), flush=True)
        if candidate:
            break
    result = {
        "input": str(args.input),
        "field_modulus": hex(MODULUS),
        "results": results,
        "warning": "candidate relations require disjoint validation and proof",
    }
    if args.output:
        write_json_atomic(args.output, result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
