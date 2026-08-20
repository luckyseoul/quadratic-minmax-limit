#!/usr/bin/env python3
"""Exact coefficient sieve for k=4 at p=19,23,31.

The one-profile energy minimum leaves only a few four-profile energy
partitions.  Degree-two coefficients lie in a one-dimensional homogeneous
kernel and degree-one coefficients in a two-dimensional kernel.  This script
checks the remaining constant reconstruction congruence for every direction
subset, top-kernel scalar, degree-one kernel vector, and admissible energy
type.  No candidate survives, before Boolean endpoint reconstruction.

The top-kernel scalar cannot be zero: four active degree-one profiles already
use 4/3 of the conserved energy, while an active degree-zero endpoint profile
uses the whole total by itself.
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from e1_gmin_m4_prop15588 import directions  # noqa: E402
from e1_gmin_m4_prop15589 import quadratic_profile_min_b  # noqa: E402


def kernel_modp(matrix: np.ndarray, p: int) -> list[np.ndarray]:
    a = np.asarray(matrix, dtype=np.int64).copy() % p
    row = 0
    pivots = []
    for column in range(a.shape[1]):
        pivot = next((r for r in range(row, a.shape[0]) if a[r, column]), None)
        if pivot is None:
            continue
        a[[row, pivot]] = a[[pivot, row]]
        a[row] = a[row] * pow(int(a[row, column]), p - 2, p) % p
        for r in range(a.shape[0]):
            if r != row and a[r, column]:
                a[r] = (a[r] - a[r, column] * a[row]) % p
        pivots.append(column)
        row += 1
        if row == a.shape[0]:
            break
    free = [column for column in range(a.shape[1]) if column not in pivots]
    basis = []
    for free_column in free:
        vector = np.zeros(a.shape[1], dtype=np.int64)
        vector[free_column] = 1
        for r, column in reversed(list(enumerate(pivots))):
            vector[column] = -(a[r] @ vector) % p
        basis.append(vector)
    return basis


def quadratic_types(p: int) -> dict[int, list[dict]]:
    midpoint = (p - 1) // 2
    centered = [value if value <= midpoint else value - p for value in range(p)]
    chi = [0] + [
        1 if pow(value, (p - 1) // 2, p) == 1 else -1
        for value in range(1, p)
    ]
    out = {-1: [], 1: []}
    for leading_character in (-1, 1):
        base = [
            1 if value == 0 else 2 if chi[value] == leading_character else 0
            for value in range(p)
        ]
        for completed_constant in range(p):
            counts = [
                base[(value - completed_constant) % p] for value in range(p)
            ]
            standard_sum = sum(
                counts[value] * centered[value] for value in range(p)
            )
            if standard_sum % p:
                raise RuntimeError("lift sum is not divisible by p")
            endpoint_flips = standard_sum // p
            if not 0 <= endpoint_flips <= counts[midpoint]:
                continue
            energy = sum(
                counts[value] * centered[value] ** 2 for value in range(p)
            ) + endpoint_flips * p
            if energy % (2 * p):
                raise RuntimeError("energy is not divisible by 2p")
            out[leading_character].append(
                {
                    "b": energy // (2 * p),
                    "completed_constant": completed_constant,
                    "endpoint_flips": endpoint_flips,
                    "endpoint_multiplicity": counts[midpoint],
                }
            )
    return out


def energy_partitions(types: dict[int, list[dict]], total: int) -> list[list[int]]:
    values = sorted({record["b"] for records in types.values() for record in records})
    return [
        list(parts)
        for parts in itertools.combinations_with_replacement(values, 4)
        if sum(parts) == total
    ]


def homogeneous_matrix(forms, degree: int, p: int) -> np.ndarray:
    # The row order is immaterial; this matches the profile reconstruction
    # convention used by the exact enumerator.
    return np.asarray(
        [
            [
                pow(forms[j][0], exponent, p)
                * pow(forms[j][1], degree - exponent, p)
                % p
                for j in range(4)
            ]
            for exponent in range(degree + 1)
        ],
        dtype=np.int64,
    )


def scan_prime(p: int) -> dict:
    midpoint = (p - 1) // 2
    normalized_total = (p * p - 1) // 8
    types = quadratic_types(p)
    minimum = quadratic_profile_min_b(p)
    maximum_relevant = normalized_total - 3 * minimum
    relevant = {
        sign: [record for record in records if record["b"] <= maximum_relevant]
        for sign, records in types.items()
    }
    partitions = energy_partitions(relevant, normalized_total)
    square, _ = directions(p)
    forms = [form for _t_of, form in square]
    subsets = list(itertools.combinations(range(len(forms)), 4))
    degree_one_parameters = np.asarray(
        list(itertools.product(range(p), repeat=2)), dtype=np.int64
    )
    chi = np.asarray(
        [0]
        + [
            1 if pow(value, (p - 1) // 2, p) == 1 else -1
            for value in range(1, p)
        ],
        dtype=np.int8,
    )
    required_constant_sum = (-(4 + 1) * pow(2, p - 2, p)) % p
    candidate_histogram = Counter()
    total_candidates = 0

    for subset in subsets:
        selected_forms = [forms[j] for j in subset]
        top_kernel = kernel_modp(homogeneous_matrix(selected_forms, 2, p), p)
        linear_kernel = kernel_modp(homogeneous_matrix(selected_forms, 1, p), p)
        if len(top_kernel) != 1 or len(linear_kernel) != 2:
            raise RuntimeError("unexpected coefficient-kernel dimensions")
        top_vector = top_kernel[0]
        if np.any(top_vector == 0):
            raise RuntimeError("a four-direction top kernel coordinate vanished")
        linear_vectors = degree_one_parameters @ np.asarray(linear_kernel) % p
        subset_candidates = 0

        for scalar in range(1, p):
            leading = scalar * top_vector % p
            inverse_4a = np.asarray(
                [pow(int(4 * value % p), p - 2, p) for value in leading]
            )
            square_completion = linear_vectors * linear_vectors * inverse_4a % p
            choices = [relevant[int(chi[value])] for value in leading]
            for combination in itertools.product(*choices):
                if sum(record["b"] for record in combination) != normalized_total:
                    continue
                completed = np.asarray(
                    [record["completed_constant"] for record in combination]
                )
                finite_constants = (square_completion + completed) % p
                rho_constants = (finite_constants + midpoint) % p
                survives = np.sum(rho_constants, axis=1) % p == required_constant_sum
                subset_candidates += int(np.count_nonzero(survives))

        candidate_histogram[subset_candidates] += 1
        total_candidates += subset_candidates

    return {
        "p": p,
        "n_square_directions": len(forms),
        "n_direction_subsets": len(subsets),
        "normalized_total_T": normalized_total,
        "minimum_quadratic_b": minimum,
        "relevant_energy_partitions": partitions,
        "coefficient_candidate_histogram": {
            str(count): frequency
            for count, frequency in sorted(candidate_histogram.items())
        },
        "total_coefficient_candidates": total_candidates,
        "boolean_endpoint_search_needed": total_candidates != 0,
        "k4_empty": total_candidates == 0,
    }


def main() -> dict:
    report = {str(p): scan_prime(p) for p in (19, 23, 31)}
    expected = {
        "19": (210, [[10, 10, 10, 15]]),
        "23": (495, [[16, 16, 16, 18], [16, 16, 17, 17]]),
        "31": (1820, [[30, 30, 30, 30]]),
    }
    if not all(
        report[key]["n_direction_subsets"] == n_subsets
        and report[key]["relevant_energy_partitions"] == partitions
        and report[key]["coefficient_candidate_histogram"] == {"0": n_subsets}
        and report[key]["k4_empty"]
        for key, (n_subsets, partitions) in expected.items()
    ):
        raise RuntimeError("k=4 coefficient sieve audit failed")

    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
