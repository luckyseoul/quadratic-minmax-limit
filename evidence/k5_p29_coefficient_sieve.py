#!/usr/bin/env python3
"""Exact coefficient, endpoint, and QVAR sieve for the p=29, k=5 stratum.

The square directions form a cyclic set of size 15.  We scan one
representative of every direction-subset orbit and restore exact full counts
with the orbit sizes.  For each nonzero cubic-kernel scalar, all low-energy
type tuples are first filtered by the constant equation, then by the two
degree-one coefficient equations, and finally by Boolean reconstruction.
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "evidence"))

from k5_cubic_energy_barrier import cubic_profile_types  # noqa: E402

P = 29
PARTITIONS = (
    (19, 19, 19, 19, 29),
    (19, 19, 19, 21, 27),
    (19, 19, 19, 23, 25),
    (19, 19, 19, 24, 24),
    (19, 19, 21, 21, 25),
    (19, 19, 21, 23, 23),
    (19, 21, 21, 21, 23),
    (21, 21, 21, 21, 21),
)


@lru_cache(maxsize=None)
def field_ctx(p: int):
    q = p * p

    def irreducible(a, b):
        return all((x * x - a * x - b) % p for x in range(p))

    ia = ib = None
    for a in range(p):
        for b in range(p):
            if irreducible(a, b):
                ia, ib = a, b
                break
        if ia is not None:
            break

    def multiply(u, v):
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        return (
            (c0 * d0 + c1 * d1 * ib) % p
            + (
                (c0 * d1 + c1 * d0 + c1 * d1 * ia) % p
            ) * p
        )

    def power(value, exponent):
        result = 1
        while exponent:
            if exponent & 1:
                result = multiply(result, value)
            value = multiply(value, value)
            exponent >>= 1
        return result

    def character(value):
        return 0 if value == 0 else (1 if power(value, (q - 1) // 2) == 1 else -1)

    def trace(value):
        return (2 * (value % p) + ia * (value // p)) % p

    return q, multiply, character, trace


def square_directions(p: int):
    q, multiply, character, trace = field_ctx(p)
    directions = []
    seen = set()
    for generator in range(1, q):
        if generator in seen:
            continue
        line = [multiply(t, generator) for t in range(1, p)]
        seen.update(line)
        annihilator = next(
            value
            for value in range(1, q)
            if trace(multiply(value, generator)) == 0
        )
        coordinate = np.asarray(
            [trace(multiply(annihilator, x)) for x in range(q)],
            dtype=np.int64,
        )
        form = (trace(multiply(annihilator, 1)), trace(multiply(annihilator, p)))
        if character(generator) == 1:
            directions.append((coordinate, form))
    return directions


def kernel_modp(matrix: np.ndarray, p: int) -> list[np.ndarray]:
    matrix = np.asarray(matrix, dtype=np.int64).copy() % p
    row = 0
    pivots = []
    for column in range(matrix.shape[1]):
        pivot = next(
            (r for r in range(row, matrix.shape[0]) if matrix[r, column]), None
        )
        if pivot is None:
            continue
        matrix[[row, pivot]] = matrix[[pivot, row]]
        matrix[row] *= pow(int(matrix[row, column]), p - 2, p)
        matrix[row] %= p
        for r in range(matrix.shape[0]):
            if r != row and matrix[r, column]:
                matrix[r] -= matrix[r, column] * matrix[row]
                matrix[r] %= p
        pivots.append(column)
        row += 1
    free = [column for column in range(matrix.shape[1]) if column not in pivots]
    basis = []
    for free_column in free:
        vector = np.zeros(matrix.shape[1], dtype=np.int64)
        vector[free_column] = 1
        for r, column in reversed(list(enumerate(pivots))):
            vector[column] = -(matrix[r] @ vector) % p
        basis.append(vector)
    return basis


def homogeneous_matrix(forms, degree: int, p: int) -> np.ndarray:
    return np.asarray(
        [
            [
                pow(form[0], exponent, p)
                * pow(form[1], degree - exponent, p)
                % p
                for form in forms
            ]
            for exponent in range(degree + 1)
        ],
        dtype=np.int64,
    )


def projective_key(value: int, p: int) -> tuple[int, int]:
    x, y = value % p, value // p
    if x:
        inverse = pow(x, p - 2, p)
        return 1, y * inverse % p
    return 0, 1


def cyclic_direction_permutation(forms, p: int) -> tuple[int, ...]:
    q, multiply, character, _trace = field_ctx(p)
    representatives = [
        coefficient_y + ((-coefficient_x) % p) * p
        for coefficient_x, coefficient_y in forms
    ]
    by_line = {
        projective_key(representative, p): index
        for index, representative in enumerate(representatives)
    }
    for multiplier in range(2, q):
        if character(multiplier) != 1:
            continue
        permutation = tuple(
            by_line[projective_key(multiply(multiplier, representative), p)]
            for representative in representatives
        )
        seen = {0}
        value = permutation[0]
        while value not in seen:
            seen.add(value)
            value = permutation[value]
        if len(seen) == len(forms):
            return permutation
    raise RuntimeError("no cyclic square-direction multiplier found")


def subset_orbits(n: int, k: int, permutation) -> list[tuple[tuple[int, ...], int]]:
    unseen = set(itertools.combinations(range(n), k))
    orbits = []
    while unseen:
        representative = min(unseen)
        orbit = set()
        current = representative
        while current not in orbit:
            orbit.add(current)
            current = tuple(sorted(permutation[index] for index in current))
        unseen.difference_update(orbit)
        orbits.append((representative, len(orbit)))
    return orbits


def primitive_element(p: int) -> int:
    q, multiply, _character, _trace = field_ctx(p)
    for candidate in range(2, q):
        value = 1
        seen = set()
        for _ in range(q - 1):
            seen.add(value)
            value = multiply(value, candidate)
        if len(seen) == q - 1:
            return candidate
    raise RuntimeError("no primitive field element")


def quartic_direction_signs(forms, p: int) -> list[int]:
    q, multiply, character, _trace = field_ctx(p)
    quartic_real = np.zeros(q, dtype=np.int8)
    value = 1
    generator = primitive_element(p)
    for exponent in range(q - 1):
        quartic_real[value] = (1, 0, -1, 0)[exponent % 4]
        value = multiply(value, generator)
    signs = []
    for coefficient_x, coefficient_y in forms:
        direction = coefficient_y + ((-coefficient_x) % p) * p
        if character(direction) != 1 or quartic_real[direction] not in (-1, 1):
            raise RuntimeError("invalid square-direction quartic sign")
        signs.append(int(quartic_real[direction]))
    return signs


def scan_p29() -> dict:
    p = P
    midpoint = (p - 1) // 2
    total = (p * p - 1) // 8
    required_rho_constant_sum = (-(5 + 1) * pow(2, p - 2, p)) % p
    required_depressed_constant_sum = (
        required_rho_constant_sum - 5 * midpoint
    ) % p
    types = cubic_profile_types(p)
    minimum = min(record[0] for record in types)
    relevant = [record for record in types if record[0] <= total - 4 * minimum]
    choices = defaultdict(list)
    for energy, leading, linear, constant in relevant:
        choices[(leading, energy)].append((linear, constant))
    assignments = sorted(
        {
            assignment
            for partition in PARTITIONS
            for assignment in itertools.permutations(partition)
        }
    )

    square = square_directions(p)
    coordinates = [coordinate for coordinate, _form in square]
    forms = [form for _coordinate, form in square]
    permutation = cyclic_direction_permutation(forms, p)
    orbits = subset_orbits(len(forms), 5, permutation)
    quartic_signs = quartic_direction_signs(forms, p)

    template_cache = {}
    candidate_histogram = Counter()
    normalized_quartic_abs_histogram = Counter()
    total_type_tuples = total_candidates = endpoint_branches = 0
    boolean_representatives = 0

    def templates(leading):
        key = tuple(map(int, leading))
        if key in template_cache:
            return template_cache[key]
        linear_rows = []
        energy_rows = []
        raw_count = 0
        for assignment in assignments:
            records = np.asarray(
                list(
                    itertools.product(
                        *[
                            choices[(int(value), energy)]
                            for value, energy in zip(leading, assignment)
                        ]
                    )
                ),
                dtype=np.int64,
            )
            raw_count += len(records)
            constant_ok = (
                np.sum(records[:, :, 1], axis=1) % p
                == required_depressed_constant_sum
            )
            if np.any(constant_ok):
                linear_rows.append(records[constant_ok])
                energy_rows.append(
                    np.tile(assignment, (int(np.count_nonzero(constant_ok)), 1))
                )
        result = (
            np.concatenate(linear_rows),
            np.concatenate(energy_rows),
            raw_count,
        )
        template_cache[key] = result
        return result

    for subset, orbit_size in orbits:
        selected = [forms[index] for index in subset]
        selected_coordinates = [coordinates[index] for index in subset]
        top_kernel = kernel_modp(homogeneous_matrix(selected, 3, p), p)
        if len(top_kernel) != 1 or np.any(top_kernel[0] == 0):
            raise RuntimeError("unexpected cubic top kernel")
        top = top_kernel[0]
        linear_matrix = homogeneous_matrix(selected, 1, p)
        subset_candidates = 0

        for scalar in range(1, p):
            leading = scalar * top % p
            records, energies, raw_count = templates(leading)
            total_type_tuples += orbit_size * raw_count
            coefficient_ok = np.all(
                records[:, :, 0] @ linear_matrix.T % p == 0, axis=1
            )
            surviving_records = records[coefficient_ok]
            surviving_energies = energies[coefficient_ok]
            subset_candidates += len(surviving_records)
            s = np.arange(p, dtype=np.int64)
            for record, energy_assignment in zip(
                surviving_records, surviving_energies
            ):
                polynomial = (
                    leading[:, None] * s[None, :] ** 3
                    + record[:, 0, None] * s[None, :]
                    + record[:, 1, None]
                ) % p
                centered = np.where(
                    polynomial <= midpoint, polynomial, polynomial - p
                )
                replacements = np.sum(centered, axis=1) // p
                choices_at_endpoint = [
                    list(
                        itertools.combinations(
                            np.where(polynomial[j] == midpoint)[0],
                            int(replacements[j]),
                        )
                    )
                    for j in range(5)
                ]
                for endpoint_sets in itertools.product(*choices_at_endpoint):
                    endpoint_branches += orbit_size
                    profiles = centered.copy()
                    for j, endpoint_set in enumerate(endpoint_sets):
                        if endpoint_set:
                            profiles[j, list(endpoint_set)] -= p
                    point_sum = sum(
                        profiles[j][selected_coordinates[j]] for j in range(5)
                    )
                    if np.all(
                        (point_sum == midpoint) | (point_sum == -midpoint - 1)
                    ):
                        boolean_representatives += orbit_size
                        normalized_quartic = sum(
                            quartic_signs[index] * energy
                            for index, energy in zip(subset, energy_assignment)
                        )
                        normalized_quartic_abs_histogram[
                            abs(normalized_quartic)
                        ] += orbit_size

        candidate_histogram[subset_candidates] += orbit_size
        total_candidates += orbit_size * subset_candidates

    quartic_threshold = Fraction(3 * total, 8)
    quartic_second_moment = (
        Fraction(
            sum(
                value * value * count
                for value, count in normalized_quartic_abs_histogram.items()
            ),
            boolean_representatives,
        )
        if boolean_representatives
        else None
    )
    return {
        "p": p,
        "n_square_directions": len(forms),
        "n_direction_subsets": sum(size for _subset, size in orbits),
        "n_cyclic_subset_orbits": len(orbits),
        "cyclic_orbit_size_histogram": {
            str(size): count
            for size, count in sorted(Counter(size for _subset, size in orbits).items())
        },
        "normalized_total_T": total,
        "minimum_cubic_b": minimum,
        "energy_partitions": [list(partition) for partition in PARTITIONS],
        "required_depressed_constant_sum_mod_p": required_depressed_constant_sum,
        "total_type_tuples_before_coefficient_sieve": total_type_tuples,
        "coefficient_candidate_histogram": {
            str(count): frequency
            for count, frequency in sorted(candidate_histogram.items())
        },
        "total_coefficient_candidates": total_candidates,
        "endpoint_branches": endpoint_branches,
        "boolean_representatives_mod_translation": boolean_representatives,
        "eps_plus_count_including_translations": boolean_representatives * p * p,
        "normalized_quartic_abs_histogram": {
            str(value): count
            for value, count in sorted(normalized_quartic_abs_histogram.items())
        },
        "E_B2": str(quartic_second_moment) if quartic_second_moment is not None else None,
        "normalized_QVAR_threshold": str(quartic_threshold),
        "clears_QVAR": (
            quartic_second_moment >= quartic_threshold
            if quartic_second_moment is not None
            else True
        ),
        "k5_empty": boolean_representatives == 0,
    }


def main() -> dict:
    report = scan_p29()
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
