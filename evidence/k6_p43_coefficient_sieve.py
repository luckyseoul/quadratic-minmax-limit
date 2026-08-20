#!/usr/bin/env python3
"""Exact coefficient, endpoint, and QVAR sieve for the p=43, k=6 stratum."""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "evidence"))

from k5_p29_coefficient_sieve import (  # noqa: E402
    cyclic_direction_permutation,
    homogeneous_matrix,
    kernel_modp,
    quartic_direction_signs,
    square_directions,
    subset_orbits,
)

P = 43
PARTITION = (38, 38, 38, 38, 38, 41)


def relevant_quartic_types(p: int, cutoff: int) -> list[tuple[int, int, int, int, int]]:
    midpoint = (p - 1) // 2
    s = np.arange(p, dtype=np.int64)
    centered = np.where(s <= midpoint, s, s - p)
    shifted = np.stack(
        [centered[(s + constant) % p] for constant in range(p)]
    )
    shifted_sq = shifted * shifted
    records = []

    for leading_value in range(1, p):
        leading = np.repeat(leading_value, p * p)
        quadratic = np.repeat(np.arange(p, dtype=np.int64), p)
        linear = np.tile(np.arange(p, dtype=np.int64), p)
        values = (
            leading[:, None] * s[None, :] ** 4
            + quadratic[:, None] * s[None, :] ** 2
            + linear[:, None] * s[None, :]
        ) % p
        counts = np.zeros((p * p, p), dtype=np.int16)
        rows = np.repeat(np.arange(p * p), p)
        np.add.at(counts, (rows, values.ravel()), 1)
        standard_sum = counts @ shifted.T
        standard_energy = counts @ shifted_sq.T
        endpoint_counts = np.stack(
            [counts[:, (midpoint - constant) % p] for constant in range(p)],
            axis=1,
        )
        replacements = standard_sum // p
        valid = (
            (standard_sum % p == 0)
            & (replacements >= 0)
            & (replacements <= endpoint_counts)
        )
        energy = (standard_energy + p * replacements) // (2 * p)
        row_indices, constants = np.nonzero(valid & (energy <= cutoff))
        records.extend(
            (
                int(energy[row, constant]),
                leading_value,
                int(quadratic[row]),
                int(linear[row]),
                int(constant),
            )
            for row, constant in zip(row_indices, constants)
        )
    return sorted(records)


def scan_p43() -> dict:
    p = P
    midpoint = (p - 1) // 2
    total = (p * p - 1) // 8
    required_rho_constant_sum = (-(6 + 1) * pow(2, p - 2, p)) % p
    required_depressed_constant_sum = (
        required_rho_constant_sum - 6 * midpoint
    ) % p
    types = relevant_quartic_types(p, cutoff=41)
    type_histogram = Counter(record[0] for record in types)
    choices = defaultdict(list)
    for energy, leading, quadratic, linear, constant in types:
        choices[(leading, energy)].append((quadratic, linear, constant))
    if type_histogram != Counter({38: 42, 41: 42}) or any(
        len(records) != 1 for records in choices.values()
    ):
        raise RuntimeError("unexpected p=43 relevant quartic types")

    assignments = sorted(set(itertools.permutations(PARTITION)))
    square = square_directions(p)
    coordinates = [coordinate for coordinate, _form in square]
    forms = [form for _coordinate, form in square]
    permutation = cyclic_direction_permutation(forms, p)
    orbits = subset_orbits(len(forms), 6, permutation)
    quartic_signs = quartic_direction_signs(forms, p)

    candidate_histogram = Counter()
    normalized_quartic_abs_histogram = Counter()
    total_candidates = endpoint_branches = boolean_representatives = 0

    for orbit_index, (subset, orbit_size) in enumerate(orbits, start=1):
        selected = [forms[index] for index in subset]
        selected_coordinates = [coordinates[index] for index in subset]
        top_kernel = kernel_modp(homogeneous_matrix(selected, 4, p), p)
        if len(top_kernel) != 1 or np.any(top_kernel[0] == 0):
            raise RuntimeError("unexpected quartic top kernel")
        top = top_kernel[0]
        quadratic_matrix = homogeneous_matrix(selected, 2, p)
        linear_matrix = homogeneous_matrix(selected, 1, p)
        subset_candidates = 0

        for scalar in range(1, p):
            leading = scalar * top % p
            for assignment in assignments:
                records = [
                    choices[(int(value), energy)][0]
                    for value, energy in zip(leading, assignment)
                ]
                quadratic = np.asarray([record[0] for record in records])
                linear = np.asarray([record[1] for record in records])
                constants = np.asarray([record[2] for record in records])
                if np.any(quadratic_matrix @ quadratic % p):
                    continue
                if np.any(linear_matrix @ linear % p):
                    continue
                if int(np.sum(constants)) % p != required_depressed_constant_sum:
                    continue
                subset_candidates += 1

                s = np.arange(p, dtype=np.int64)
                polynomial = (
                    leading[:, None] * s[None, :] ** 4
                    + quadratic[:, None] * s[None, :] ** 2
                    + linear[:, None] * s[None, :]
                    + constants[:, None]
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
                    for j in range(6)
                ]
                for endpoint_sets in itertools.product(*choices_at_endpoint):
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
                        normalized_quartic = sum(
                            quartic_signs[index] * energy
                            for index, energy in zip(subset, assignment)
                        )
                        normalized_quartic_abs_histogram[
                            abs(normalized_quartic)
                        ] += orbit_size

        candidate_histogram[subset_candidates] += orbit_size
        total_candidates += orbit_size * subset_candidates
        if orbit_index % 500 == 0:
            print(
                f"orbit {orbit_index}/{len(orbits)} candidates={total_candidates}",
                flush=True,
            )

    quartic_threshold = Fraction(3 * total, 8)
    quartic_moment = (
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
        "energy_partition": list(PARTITION),
        "relevant_type_histogram": {
            str(energy): count for energy, count in sorted(type_histogram.items())
        },
        "required_depressed_constant_sum_mod_p": required_depressed_constant_sum,
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
        "E_B2": str(quartic_moment) if quartic_moment is not None else None,
        "normalized_QVAR_threshold": str(quartic_threshold),
        "clears_QVAR": bool(
            quartic_moment >= quartic_threshold if quartic_moment is not None else True
        ),
        "k6_empty": boolean_representatives == 0,
    }


def main() -> dict:
    report = scan_p43()
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
