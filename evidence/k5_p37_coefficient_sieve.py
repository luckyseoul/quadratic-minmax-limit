#!/usr/bin/env python3
"""Coefficient-compatibility sieve for the p=37, k=5 stratum.

The only cubic energy partition is 34+34+34+34+35.  After the unique global
translation that depresses all five cubics, the exact low-energy type census
also forces the b=35 profile's constant to be 18 (rather than 19).  We then
exhaust all direction five-subsets, nonzero cubic-kernel scalars, and the
remaining three choices for every linear coefficient.
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "evidence"))

from e1_gmin_m4_prop15588 import directions  # noqa: E402
from k4_p3mod4_coefficient_sieve import kernel_modp  # noqa: E402
from k5_cubic_energy_barrier import cubic_profile_types  # noqa: E402


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


def scan_p37() -> dict:
    p = 37
    midpoint = (p - 1) // 2
    total = (p * p - 1) // 8
    types = cubic_profile_types(p)
    minimum = min(record[0] for record in types)
    relevant = [record for record in types if record[0] <= total - 4 * minimum]
    type_histogram = Counter(record[0] for record in relevant)

    linear_choices = defaultdict(list)
    for energy, leading, linear, constant in relevant:
        # sum_j (d_j+midpoint) = -(k+1)/2 = -3 (mod p).
        # Four b=34 constants are zero, so the b=35 constant must be 18.
        if energy == 35 and constant != 18:
            continue
        linear_choices[(energy, leading)].append(linear)

    leading34 = {leading for energy, leading in linear_choices if energy == 34}
    leading35 = {leading for energy, leading in linear_choices if energy == 35}
    if not (
        type_histogram == Counter({35: 72, 34: 36})
        and len(leading34) == len(leading35) == 12
        and leading34.isdisjoint(leading35)
        and all(len(values) == 3 for values in linear_choices.values())
    ):
        raise RuntimeError("unexpected p=37 low cubic types")

    square, _ = directions(p)
    forms = [form for _coordinate, form in square]
    subsets = list(itertools.combinations(range(len(forms)), 5))
    candidate_histogram = Counter()
    scalar_pattern_histogram = Counter()
    total_scalar_patterns = 0
    total_coefficient_candidates = 0

    for subset in subsets:
        selected = [forms[index] for index in subset]
        top_kernel = kernel_modp(homogeneous_matrix(selected, 3, p), p)
        if len(top_kernel) != 1 or np.any(top_kernel[0] == 0):
            raise RuntimeError("unexpected cubic top kernel")
        top = top_kernel[0]
        linear_matrix = homogeneous_matrix(selected, 1, p)
        subset_patterns = subset_candidates = 0

        for scalar in range(1, p):
            leading = scalar * top % p
            energies = [
                34 if int(value) in leading34
                else 35 if int(value) in leading35
                else None
                for value in leading
            ]
            if energies.count(34) != 4 or energies.count(35) != 1:
                continue
            subset_patterns += 1
            choices = [
                linear_choices[(energy, int(value))]
                for energy, value in zip(energies, leading)
            ]
            candidates = np.asarray(list(itertools.product(*choices)), dtype=np.int64)
            compatible = np.all(candidates @ linear_matrix.T % p == 0, axis=1)
            subset_candidates += int(np.count_nonzero(compatible))

        scalar_pattern_histogram[subset_patterns] += 1
        candidate_histogram[subset_candidates] += 1
        total_scalar_patterns += subset_patterns
        total_coefficient_candidates += subset_candidates

    return {
        "p": p,
        "n_square_directions": len(forms),
        "n_direction_subsets": len(subsets),
        "normalized_total_T": total,
        "minimum_cubic_b": minimum,
        "energy_partition": [34, 34, 34, 34, 35],
        "relevant_type_histogram_before_constant_sieve": {
            str(energy): count for energy, count in sorted(type_histogram.items())
        },
        "constant_sieve": {
            "required_rho_constant_sum_mod_p": 34,
            "surviving_b35_constant": 18,
            "discarded_b35_constant": 19,
        },
        "scalar_pattern_histogram": {
            str(count): frequency
            for count, frequency in sorted(scalar_pattern_histogram.items())
        },
        "total_scalar_patterns": total_scalar_patterns,
        "coefficient_candidate_histogram": {
            str(count): frequency
            for count, frequency in sorted(candidate_histogram.items())
        },
        "total_coefficient_candidates": total_coefficient_candidates,
        "boolean_endpoint_search_needed": total_coefficient_candidates != 0,
        "k5_empty": total_coefficient_candidates == 0,
    }


def main() -> dict:
    report = scan_p37()
    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
