#!/usr/bin/env python3
"""Fast exact coefficient and Boolean sieve for residual k=6 strata.

For six distinct projective directions, every three columns of the quadratic
coefficient matrix and every two columns of the linear coefficient matrix are
invertible.  For each energy assignment we therefore enumerate records only
in the three smallest profile groups, solve for the other three quadratic
coefficients, enumerate one further linear coefficient, and solve for the
last two.  This is exactly equivalent to the six-way Cartesian product in
``k6_p41_coefficient_sieve.py`` but avoids materializing almost all tuples.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np
from numba import njit

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "evidence"))

from k5_p29_coefficient_sieve import (  # noqa: E402
    cyclic_direction_permutation,
    homogeneous_matrix,
    kernel_modp,
    square_directions,
    subset_orbits,
)
from k5_p23_coefficient_sieve import quartic_kernel  # noqa: E402
from k6_p43_coefficient_sieve import relevant_quartic_types  # noqa: E402


def inverse_modp(matrix: np.ndarray, p: int) -> np.ndarray:
    """Return the inverse of a small square matrix over F_p."""
    n = matrix.shape[0]
    work = np.concatenate(
        [np.asarray(matrix, dtype=np.int64) % p, np.eye(n, dtype=np.int64)],
        axis=1,
    )
    for column in range(n):
        pivot = next(row for row in range(column, n) if work[row, column] % p)
        work[[column, pivot]] = work[[pivot, column]]
        work[column] *= pow(int(work[column, column]), p - 2, p)
        work[column] %= p
        for row in range(n):
            if row != column and work[row, column]:
                work[row] -= work[row, column] * work[column]
                work[row] %= p
    return work[:, n:] % p


def coefficient_solvers(
    quadratic_matrix: np.ndarray, linear_matrix: np.ndarray, p: int
) -> tuple[np.ndarray, np.ndarray]:
    """Precompute ordered 3-column and 2-column inverse matrices."""
    quadratic = np.zeros((6, 6, 6, 3, 3), dtype=np.int64)
    linear = np.zeros((6, 6, 2, 2), dtype=np.int64)
    for ordered in itertools.permutations(range(6), 3):
        quadratic[ordered] = inverse_modp(quadratic_matrix[:, ordered], p)
    for ordered in itertools.permutations(range(6), 2):
        linear[ordered] = inverse_modp(linear_matrix[:, ordered], p)
    return quadratic, linear


def dense_profile_groups(
    p: int, energy_values: list[int], types: list[tuple[int, int, int, int, int]]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Pack profile records and exact q/(q,l) lookup ranges."""
    energy_index = {energy: index for index, energy in enumerate(energy_values)}
    choices: dict[tuple[int, int], list[tuple[int, int, int]]] = defaultdict(list)
    for energy, leading, quadratic, linear, constant in types:
        choices[(leading, energy_index[energy])].append(
            (quadratic, linear, constant)
        )
    maximum = max(map(len, choices.values()), default=0)
    n_energy = len(energy_values)
    records = np.zeros((p, n_energy, maximum, 3), dtype=np.int64)
    counts = np.zeros((p, n_energy), dtype=np.int64)
    q_start = np.full((p, n_energy, p), -1, dtype=np.int64)
    q_end = np.full((p, n_energy, p), -1, dtype=np.int64)
    ql_start = np.full((p, n_energy, p, p), -1, dtype=np.int64)
    ql_end = np.full((p, n_energy, p, p), -1, dtype=np.int64)
    for (leading, energy), group in choices.items():
        group = sorted(group)
        counts[leading, energy] = len(group)
        records[leading, energy, : len(group)] = group
        for index, (quadratic, linear, _constant) in enumerate(group):
            if q_start[leading, energy, quadratic] < 0:
                q_start[leading, energy, quadratic] = index
            q_end[leading, energy, quadratic] = index + 1
            if ql_start[leading, energy, quadratic, linear] < 0:
                ql_start[leading, energy, quadratic, linear] = index
            ql_end[leading, energy, quadratic, linear] = index + 1
    return records, counts, q_start, q_end, ql_start, ql_end


@njit(cache=True)
def solve_coefficient_assignments(
    records,
    counts,
    q_start,
    q_end,
    ql_start,
    ql_end,
    assignments,
    leading,
    quadratic_matrix,
    linear_matrix,
    quadratic_inverse,
    linear_inverse,
    p,
    target_constant,
    output,
):
    """Enumerate all coefficient-compatible tuples for one top scalar."""
    total_candidates = 0
    stored_candidates = 0
    raw_tuples = 0
    group_counts = np.empty(6, dtype=np.int64)
    order = np.empty(6, dtype=np.int64)
    chosen_indices = np.empty(6, dtype=np.int64)
    q_required = np.empty(3, dtype=np.int64)
    linear_rhs = np.empty(2, dtype=np.int64)

    for assignment_index in range(assignments.shape[0]):
        valid = True
        raw = 1
        for direction in range(6):
            count = counts[leading[direction], assignments[assignment_index, direction]]
            group_counts[direction] = count
            if count == 0:
                valid = False
            raw *= count
            order[direction] = direction
        if not valid:
            continue
        raw_tuples += raw

        # Put the three smallest groups first.  They are the free quadratic
        # variables; the other three are then uniquely determined.
        for i in range(1, 6):
            value = order[i]
            j = i - 1
            while j >= 0 and group_counts[order[j]] > group_counts[value]:
                order[j + 1] = order[j]
                j -= 1
            order[j + 1] = value
        a0, a1, a2 = order[0], order[1], order[2]
        b0, b1, b2 = order[3], order[4], order[5]
        e0 = assignments[assignment_index, a0]
        e1 = assignments[assignment_index, a1]
        e2 = assignments[assignment_index, a2]
        eb0 = assignments[assignment_index, b0]
        eb1 = assignments[assignment_index, b1]
        eb2 = assignments[assignment_index, b2]
        l0, l1, l2 = leading[a0], leading[a1], leading[a2]
        lb0, lb1, lb2 = leading[b0], leading[b1], leading[b2]

        for i0 in range(group_counts[a0]):
            q0 = records[l0, e0, i0, 0]
            for i1 in range(group_counts[a1]):
                q1 = records[l1, e1, i1, 0]
                for i2 in range(group_counts[a2]):
                    q2 = records[l2, e2, i2, 0]
                    for row in range(3):
                        rhs = -(
                            quadratic_matrix[row, a0] * q0
                            + quadratic_matrix[row, a1] * q1
                            + quadratic_matrix[row, a2] * q2
                        )
                        q_required[row] = rhs % p
                    rhs0, rhs1, rhs2 = q_required[0], q_required[1], q_required[2]
                    q_b0 = (
                        quadratic_inverse[b0, b1, b2, 0, 0] * rhs0
                        + quadratic_inverse[b0, b1, b2, 0, 1] * rhs1
                        + quadratic_inverse[b0, b1, b2, 0, 2] * rhs2
                    ) % p
                    q_b1 = (
                        quadratic_inverse[b0, b1, b2, 1, 0] * rhs0
                        + quadratic_inverse[b0, b1, b2, 1, 1] * rhs1
                        + quadratic_inverse[b0, b1, b2, 1, 2] * rhs2
                    ) % p
                    q_b2 = (
                        quadratic_inverse[b0, b1, b2, 2, 0] * rhs0
                        + quadratic_inverse[b0, b1, b2, 2, 1] * rhs1
                        + quadratic_inverse[b0, b1, b2, 2, 2] * rhs2
                    ) % p

                    start0 = q_start[lb0, eb0, q_b0]
                    if start0 < 0:
                        continue
                    end0 = q_end[lb0, eb0, q_b0]
                    for ib0 in range(start0, end0):
                        for row in range(2):
                            linear_rhs[row] = -(
                                linear_matrix[row, a0] * records[l0, e0, i0, 1]
                                + linear_matrix[row, a1] * records[l1, e1, i1, 1]
                                + linear_matrix[row, a2] * records[l2, e2, i2, 1]
                                + linear_matrix[row, b0]
                                * records[lb0, eb0, ib0, 1]
                            ) % p
                        l_b1 = (
                            linear_inverse[b1, b2, 0, 0] * linear_rhs[0]
                            + linear_inverse[b1, b2, 0, 1] * linear_rhs[1]
                        ) % p
                        l_b2 = (
                            linear_inverse[b1, b2, 1, 0] * linear_rhs[0]
                            + linear_inverse[b1, b2, 1, 1] * linear_rhs[1]
                        ) % p
                        start1 = ql_start[lb1, eb1, q_b1, l_b1]
                        start2 = ql_start[lb2, eb2, q_b2, l_b2]
                        if start1 < 0 or start2 < 0:
                            continue
                        end1 = ql_end[lb1, eb1, q_b1, l_b1]
                        end2 = ql_end[lb2, eb2, q_b2, l_b2]
                        constant_first = (
                            records[l0, e0, i0, 2]
                            + records[l1, e1, i1, 2]
                            + records[l2, e2, i2, 2]
                            + records[lb0, eb0, ib0, 2]
                        )
                        for ib1 in range(start1, end1):
                            partial = constant_first + records[lb1, eb1, ib1, 2]
                            for ib2 in range(start2, end2):
                                if (
                                    partial + records[lb2, eb2, ib2, 2]
                                ) % p != target_constant:
                                    continue
                                chosen_indices[a0] = i0
                                chosen_indices[a1] = i1
                                chosen_indices[a2] = i2
                                chosen_indices[b0] = ib0
                                chosen_indices[b1] = ib1
                                chosen_indices[b2] = ib2
                                if stored_candidates < output.shape[0]:
                                    for direction in range(6):
                                        energy = assignments[
                                            assignment_index, direction
                                        ]
                                        output[stored_candidates, direction] = energy
                                        index = chosen_indices[direction]
                                        for coefficient in range(3):
                                            output[
                                                stored_candidates,
                                                6 + 3 * direction + coefficient,
                                            ] = records[
                                                leading[direction],
                                                energy,
                                                index,
                                                coefficient,
                                            ]
                                    stored_candidates += 1
                                total_candidates += 1
    return total_candidates, stored_candidates, raw_tuples


def energy_assignments(energy_values: list[int], total: int) -> tuple[list, np.ndarray]:
    partitions = [
        partition
        for partition in itertools.combinations_with_replacement(energy_values, 6)
        if sum(partition) == total
    ]
    assignments = sorted(
        {
            assignment
            for partition in partitions
            for assignment in itertools.permutations(partition)
        }
    )
    index = {energy: i for i, energy in enumerate(energy_values)}
    return partitions, np.asarray(
        [[index[energy] for energy in assignment] for assignment in assignments],
        dtype=np.int64,
    )


def scan_prime_fast(
    p: int,
    minimum: int,
    shard_index: int = 0,
    shard_count: int = 1,
    candidate_capacity: int = 200_000,
) -> dict:
    midpoint = (p - 1) // 2
    total = (p * p - 1) // 8
    cutoff = total - 5 * minimum
    target = (
        (-(6 + 1) * pow(2, p - 2, p)) % p - 6 * midpoint
    ) % p
    types = relevant_quartic_types(p, cutoff)
    histogram = Counter(record[0] for record in types)
    energies = sorted(histogram)
    partitions, assignments = energy_assignments(energies, total)
    packed = dense_profile_groups(p, energies, types)
    records, counts, q_start, q_end, ql_start, ql_end = packed

    square = square_directions(p)
    coordinates = [coordinate for coordinate, _form in square]
    forms = [form for _coordinate, form in square]
    permutation = cyclic_direction_permutation(forms, p)
    all_orbits = subset_orbits(len(forms), 6, permutation)
    kernel_real, kernel_imag = quartic_kernel(p)

    candidate_histogram = Counter()
    quartic_abs_sq_histogram = Counter()
    raw_tuples = total_candidates = endpoint_branches = 0
    boolean_representatives = processed_subsets = 0
    processed_orbits = 0
    output = np.empty((candidate_capacity, 24), dtype=np.int64)

    for orbit_index, (subset, orbit_size) in enumerate(all_orbits):
        if orbit_index % shard_count != shard_index:
            continue
        processed_orbits += 1
        processed_subsets += orbit_size
        selected = [forms[index] for index in subset]
        selected_coordinates = [coordinates[index] for index in subset]
        top_kernel = kernel_modp(homogeneous_matrix(selected, 4, p), p)
        if len(top_kernel) != 1 or np.any(top_kernel[0] == 0):
            raise RuntimeError("unexpected quartic top kernel")
        top = top_kernel[0]
        quadratic_matrix = homogeneous_matrix(selected, 2, p)
        linear_matrix = homogeneous_matrix(selected, 1, p)
        quadratic_inverse, linear_inverse = coefficient_solvers(
            quadratic_matrix, linear_matrix, p
        )
        subset_candidates = 0

        for scalar in range(1, p):
            leading = scalar * top % p
            count, stored, raw = solve_coefficient_assignments(
                records,
                counts,
                q_start,
                q_end,
                ql_start,
                ql_end,
                assignments,
                leading,
                quadratic_matrix,
                linear_matrix,
                quadratic_inverse,
                linear_inverse,
                p,
                target,
                output,
            )
            if stored != count:
                raise RuntimeError(
                    f"candidate capacity {candidate_capacity} too small for {count}"
                )
            raw_tuples += orbit_size * int(raw)
            subset_candidates += int(count)

            s = np.arange(p, dtype=np.int64)
            for candidate in output[:stored]:
                profile_records = candidate[6:].reshape(6, 3)
                polynomial = (
                    leading[:, None] * s[None, :] ** 4
                    + profile_records[:, 0, None] * s[None, :] ** 2
                    + profile_records[:, 1, None] * s[None, :]
                    + profile_records[:, 2, None]
                ) % p
                centered = np.where(
                    polynomial <= midpoint, polynomial, polynomial - p
                )
                replacements = np.sum(centered, axis=1) // p
                endpoint_choices = [
                    list(
                        itertools.combinations(
                            np.where(polynomial[j] == midpoint)[0],
                            int(replacements[j]),
                        )
                    )
                    for j in range(6)
                ]
                for endpoint_sets in itertools.product(*endpoint_choices):
                    endpoint_branches += orbit_size
                    profiles = centered.copy()
                    for j, endpoint_set in enumerate(endpoint_sets):
                        if endpoint_set:
                            profiles[j, list(endpoint_set)] -= p
                    point_sum = sum(
                        profiles[j][selected_coordinates[j]] for j in range(6)
                    )
                    if np.all(
                        (point_sum == midpoint) | (point_sum == -midpoint - 1)
                    ):
                        boolean_representatives += orbit_size
                        negative = (point_sum == -midpoint - 1).astype(np.int64)
                        real = int(negative @ kernel_real @ negative)
                        imag = int(negative @ kernel_imag @ negative)
                        quartic_abs_sq_histogram[real * real + imag * imag] += (
                            orbit_size
                        )

        candidate_histogram[subset_candidates] += orbit_size
        total_candidates += orbit_size * subset_candidates
        if processed_orbits % 10 == 0:
            print(
                f"shard {shard_index}/{shard_count} orbit "
                f"{processed_orbits} candidates={total_candidates}",
                flush=True,
            )

    threshold = Fraction(3 * p * p * (p * p - 1), 16)
    moment = (
        Fraction(
            sum(value * count for value, count in quartic_abs_sq_histogram.items()),
            boolean_representatives,
        )
        if boolean_representatives
        else None
    )
    return {
        "p": p,
        "algorithm": "quadratic-3 plus linear-2 exact elimination",
        "shard_index": shard_index,
        "shard_count": shard_count,
        "processed_cyclic_subset_orbits": processed_orbits,
        "processed_direction_subsets": processed_subsets,
        "n_square_directions": len(forms),
        "n_direction_subsets": sum(size for _subset, size in all_orbits),
        "n_cyclic_subset_orbits": len(all_orbits),
        "normalized_total_T": total,
        "minimum_quartic_b": minimum,
        "energy_partitions": [list(partition) for partition in partitions],
        "relevant_type_histogram": {
            str(energy): count for energy, count in sorted(histogram.items())
        },
        "total_type_tuples_before_coefficient_sieve": raw_tuples,
        "coefficient_candidate_histogram": {
            str(count): frequency
            for count, frequency in sorted(candidate_histogram.items())
        },
        "total_coefficient_candidates": total_candidates,
        "endpoint_branches": endpoint_branches,
        "boolean_representatives_mod_translation": boolean_representatives,
        "eps_plus_count_including_translations": boolean_representatives * p * p,
        "abs_Zpsi_sq_histogram": {
            str(value): count
            for value, count in sorted(quartic_abs_sq_histogram.items())
        },
        "E_abs_Zpsi_sq": str(moment) if moment is not None else None,
        "QVAR_threshold": str(threshold),
        "clears_QVAR": bool(moment >= threshold if moment is not None else True),
        "k6_empty": boolean_representatives == 0,
    }


def main() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("p", type=int)
    parser.add_argument("minimum", type=int)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--candidate-capacity", type=int, default=200_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = scan_prime_fast(
        args.p,
        args.minimum,
        args.shard_index,
        args.shard_count,
        args.candidate_capacity,
    )
    output = args.output or Path(__file__).with_name(
        f"k6_p{args.p}_coefficient_sieve_fast.json"
    )
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
