#!/usr/bin/env python3
"""Exact coefficient, endpoint, and QVAR sieve for the p=41, k=6 stratum."""
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
    square_directions,
    subset_orbits,
)
from k5_p23_coefficient_sieve import quartic_kernel  # noqa: E402
from k6_p43_coefficient_sieve import relevant_quartic_types  # noqa: E402

P = 41


def scan_prime(p: int, minimum: int) -> dict:
    midpoint = (p - 1) // 2
    total = (p * p - 1) // 8
    cutoff = total - 5 * minimum
    required_rho_constant_sum = (-(6 + 1) * pow(2, p - 2, p)) % p
    required_depressed_constant_sum = (
        required_rho_constant_sum - 6 * midpoint
    ) % p
    types = relevant_quartic_types(p, cutoff)
    type_histogram = Counter(record[0] for record in types)
    choices = defaultdict(list)
    for energy, leading, quadratic, linear, constant in types:
        choices[(leading, energy)].append((quadratic, linear, constant))
    energy_values = sorted(type_histogram)
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

    square = square_directions(p)
    coordinates = [coordinate for coordinate, _form in square]
    forms = [form for _coordinate, form in square]
    permutation = cyclic_direction_permutation(forms, p)
    orbits = subset_orbits(len(forms), 6, permutation)
    kernel_real, kernel_imag = quartic_kernel(p)

    template_cache = {}
    candidate_histogram = Counter()
    quartic_abs_sq_histogram = Counter()
    total_type_tuples = total_candidates = endpoint_branches = 0
    boolean_representatives = 0

    def templates(leading):
        key = tuple(map(int, leading))
        if key in template_cache:
            return template_cache[key]
        record_rows = []
        energy_rows = []
        raw_count = 0
        for assignment in assignments:
            profile_choices = [
                choices[(int(value), energy)]
                for value, energy in zip(leading, assignment)
            ]
            if not all(profile_choices):
                continue
            records = np.asarray(
                list(itertools.product(*profile_choices)), dtype=np.int64
            )
            raw_count += len(records)
            constant_ok = (
                np.sum(records[:, :, 2], axis=1) % p
                == required_depressed_constant_sum
            )
            if np.any(constant_ok):
                record_rows.append(records[constant_ok])
                energy_rows.append(
                    np.tile(assignment, (int(np.count_nonzero(constant_ok)), 1))
                )
        if record_rows:
            result = (
                np.concatenate(record_rows),
                np.concatenate(energy_rows),
                raw_count,
            )
        else:
            result = (
                np.empty((0, 6, 3), dtype=np.int64),
                np.empty((0, 6), dtype=np.int64),
                raw_count,
            )
        template_cache[key] = result
        return result

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
            records, energies, raw_count = templates(leading)
            total_type_tuples += orbit_size * raw_count
            if not len(records):
                continue
            coefficient_ok = np.all(
                records[:, :, 0] @ quadratic_matrix.T % p == 0, axis=1
            ) & np.all(records[:, :, 1] @ linear_matrix.T % p == 0, axis=1)
            surviving_records = records[coefficient_ok]
            surviving_energies = energies[coefficient_ok]
            subset_candidates += len(surviving_records)

            s = np.arange(p, dtype=np.int64)
            for record, energy_assignment in zip(
                surviving_records, surviving_energies
            ):
                polynomial = (
                    leading[:, None] * s[None, :] ** 4
                    + record[:, 0, None] * s[None, :] ** 2
                    + record[:, 1, None] * s[None, :]
                    + record[:, 2, None]
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
                        negative = (point_sum == -midpoint - 1).astype(np.int64)
                        quartic_real = int(negative @ kernel_real @ negative)
                        quartic_imag = int(negative @ kernel_imag @ negative)
                        quartic_abs_sq_histogram[
                            quartic_real * quartic_real
                            + quartic_imag * quartic_imag
                        ] += orbit_size

        candidate_histogram[subset_candidates] += orbit_size
        total_candidates += orbit_size * subset_candidates
        if orbit_index % 500 == 0:
            print(
                f"orbit {orbit_index}/{len(orbits)} candidates={total_candidates}",
                flush=True,
            )

    quartic_threshold = Fraction(3 * p * p * (p * p - 1), 16)
    quartic_moment = (
        Fraction(
            sum(value * count for value, count in quartic_abs_sq_histogram.items()),
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
        "normalized_total_T": total,
        "minimum_quartic_b": minimum,
        "energy_partitions": [list(partition) for partition in partitions],
        "relevant_type_histogram": {
            str(energy): count for energy, count in sorted(type_histogram.items())
        },
        "total_type_tuples_before_coefficient_sieve": total_type_tuples,
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
        "E_abs_Zpsi_sq": str(quartic_moment) if quartic_moment is not None else None,
        "QVAR_threshold": str(quartic_threshold),
        "clears_QVAR": bool(
            quartic_moment >= quartic_threshold if quartic_moment is not None else True
        ),
        "k6_empty": boolean_representatives == 0,
    }


def main() -> dict:
    report = scan_prime(P, minimum=34)
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
