#!/usr/bin/env python3
"""Coefficient-compatibility sieve for the p=31, k=5 stratum."""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "evidence"))

from e1_gmin_m4_prop15588 import directions, field_ctx  # noqa: E402
from k4_p3mod4_coefficient_sieve import kernel_modp  # noqa: E402
from k5_cubic_energy_barrier import cubic_profile_types  # noqa: E402
from k5_p37_coefficient_sieve import homogeneous_matrix  # noqa: E402


PARTITIONS = (
    (23, 23, 23, 24, 27),
    (23, 23, 24, 24, 26),
    (24, 24, 24, 24, 24),
)


def quartic_direction_signs(p: int, forms) -> list[int]:
    q, multiply, character, _trace = field_ctx(p)
    generator = next(
        candidate
        for candidate in range(2, q)
        if _is_primitive(candidate, q, multiply)
    )
    quartic_real = np.zeros(q, dtype=np.int8)
    value = 1
    for exponent in range(q - 1):
        quartic_real[value] = (1, 0, -1, 0)[exponent % 4]
        value = multiply(value, generator)
    signs = []
    for coefficient_x, coefficient_y in forms:
        direction = coefficient_y + ((-coefficient_x) % p) * p
        if character(direction) != 1 or quartic_real[direction] not in (-1, 1):
            raise RuntimeError("failed to recover a square-direction quartic sign")
        signs.append(int(quartic_real[direction]))
    return signs


def _is_primitive(candidate: int, q: int, multiply) -> bool:
    value = 1
    seen = set()
    for _ in range(q - 1):
        seen.add(value)
        value = multiply(value, candidate)
    return len(seen) == q - 1


def scan_p31() -> dict:
    p = 31
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
    square, _ = directions(p)
    coordinates = [coordinate for coordinate, _form in square]
    forms = [form for _coordinate, form in square]
    quartic_signs = quartic_direction_signs(p, forms)
    subsets = list(itertools.combinations(range(len(forms)), 5))
    candidate_histogram = Counter()
    scalar_pattern_histogram = Counter()
    total_scalar_patterns = total_type_tuples = total_candidates = 0
    total_endpoint_branches = boolean_representatives = 0
    normalized_quartic_histogram = Counter()

    for subset in subsets:
        selected = [forms[index] for index in subset]
        top_kernel = kernel_modp(homogeneous_matrix(selected, 3, p), p)
        if len(top_kernel) != 1 or np.any(top_kernel[0] == 0):
            raise RuntimeError("unexpected cubic top kernel")
        top = top_kernel[0]
        linear_matrix = homogeneous_matrix(selected, 1, p)
        selected_coordinates = [coordinates[index] for index in subset]
        subset_patterns = subset_candidates = 0

        for scalar in range(1, p):
            leading = scalar * top % p
            compatible_assignments = [
                assignment
                for assignment in assignments
                if all(
                    choices[(int(value), energy)]
                    for value, energy in zip(leading, assignment)
                )
            ]
            if not compatible_assignments:
                continue
            subset_patterns += 1
            for assignment in compatible_assignments:
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
                total_type_tuples += len(records)
                constant_ok = (
                    np.sum(records[:, :, 1], axis=1) % p
                    == required_depressed_constant_sum
                )
                linear = records[constant_ok, :, 0]
                if len(linear):
                    coefficient_ok = np.all(
                        linear @ linear_matrix.T % p == 0, axis=1
                    )
                    surviving = records[constant_ok][coefficient_ok]
                    subset_candidates += len(surviving)
                    s = np.arange(p, dtype=np.int64)
                    for record in surviving:
                        linear_coefficients = record[:, 0]
                        constants = record[:, 1]
                        polynomial = (
                            leading[:, None] * s[None, :] ** 3
                            + linear_coefficients[:, None] * s[None, :]
                            + constants[:, None]
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
                            for j in range(5)
                        ]
                        for endpoint_sets in itertools.product(*endpoint_choices):
                            total_endpoint_branches += 1
                            profiles = centered.copy()
                            for j, endpoint_set in enumerate(endpoint_sets):
                                if endpoint_set:
                                    profiles[j, list(endpoint_set)] -= p
                            if np.any(np.sum(profiles, axis=1)):
                                raise RuntimeError("profile zero-sum audit failed")
                            if any(
                                np.sum(profiles[j] ** 2) != 2 * p * assignment[j]
                                for j in range(5)
                            ):
                                raise RuntimeError("profile energy audit failed")
                            point_sum = sum(
                                profiles[j][selected_coordinates[j]]
                                for j in range(5)
                            )
                            if np.all(
                                (point_sum == midpoint)
                                | (point_sum == -midpoint - 1)
                            ):
                                boolean_representatives += 1
                                normalized_quartic = sum(
                                    quartic_signs[index] * energy
                                    for index, energy in zip(subset, assignment)
                                )
                                normalized_quartic_histogram[normalized_quartic] += 1

        scalar_pattern_histogram[subset_patterns] += 1
        candidate_histogram[subset_candidates] += 1
        total_scalar_patterns += subset_patterns
        total_candidates += subset_candidates

    quartic_second_moment = Fraction(
        sum(value * value * count for value, count in normalized_quartic_histogram.items()),
        boolean_representatives,
    )
    quartic_threshold = Fraction(3 * total, 8)
    return {
        "p": p,
        "n_square_directions": len(forms),
        "n_direction_subsets": len(subsets),
        "normalized_total_T": total,
        "minimum_cubic_b": minimum,
        "energy_partitions": [list(partition) for partition in PARTITIONS],
        "relevant_type_histogram": {
            str(energy): count
            for energy, count in sorted(Counter(row[0] for row in relevant).items())
        },
        "required_depressed_constant_sum_mod_p": required_depressed_constant_sum,
        "scalar_pattern_histogram": {
            str(count): frequency
            for count, frequency in sorted(scalar_pattern_histogram.items())
        },
        "total_scalar_patterns": total_scalar_patterns,
        "total_type_tuples_before_constant_and_linear_sieves": total_type_tuples,
        "coefficient_candidate_histogram": {
            str(count): frequency
            for count, frequency in sorted(candidate_histogram.items())
        },
        "total_coefficient_candidates": total_candidates,
        "endpoint_branches": total_endpoint_branches,
        "boolean_representatives_mod_translation": boolean_representatives,
        "eps_plus_count_including_translations": boolean_representatives * p * p,
        "normalized_quartic_histogram": {
            str(value): count
            for value, count in sorted(normalized_quartic_histogram.items())
        },
        "E_B2": str(quartic_second_moment),
        "normalized_QVAR_threshold": str(quartic_threshold),
        "clears_QVAR": quartic_second_moment >= quartic_threshold,
        "boolean_endpoint_search_needed": False,
        "k5_empty": boolean_representatives == 0,
    }


def main() -> dict:
    report = scan_p31()
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
