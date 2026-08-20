#!/usr/bin/env python3
"""Exact depressed-quintic profile probe for the k=7 Max+ stratum.

For seven active directions, the degree-five coefficient level has a
one-dimensional kernel.  Translation acts isomorphically on the
two-dimensional degree-four kernel, so every translation orbit has a unique
representative whose profiles are

    a*s**5 + c*s**3 + d*s**2 + e*s + f.

This first-stage probe enumerates every zero-sum odd lift of those reduced
polynomials, records its normalized profile energy, and audits the universal
Vandermonde ranks at all remaining coefficient levels.  It intentionally
stops before the coupled seven-profile Boolean sieve: the output determines
whether that sieve is finite and how much elimination is available.
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "evidence"))

from k5_p29_coefficient_sieve import (  # noqa: E402
    homogeneous_matrix,
    kernel_modp,
    square_directions,
)


def relevant_quintic_types(
    p: int, cutoff: int | None = None
) -> list[tuple[int, int, int, int, int, int]]:
    """Enumerate depressed quintic lifts as (energy,a,c,d,e,f)."""
    midpoint = (p - 1) // 2
    s = np.arange(p, dtype=np.int64)
    centered = np.where(s <= midpoint, s, s - p)
    shifted = np.stack(
        [centered[(s + constant) % p] for constant in range(p)]
    )
    shifted_sq = shifted * shifted
    records: list[tuple[int, int, int, int, int, int]] = []

    cubic = np.repeat(np.arange(p, dtype=np.int64), p * p)
    quadratic = np.tile(np.repeat(np.arange(p, dtype=np.int64), p), p)
    linear = np.tile(np.arange(p, dtype=np.int64), p * p)
    row_count = p**3
    rows = np.repeat(np.arange(row_count), p)
    powers = np.stack([s, s**2, s**3, s**5])

    for leading in range(1, p):
        values = (
            leading * powers[3][None, :]
            + cubic[:, None] * powers[2][None, :]
            + quadratic[:, None] * powers[1][None, :]
            + linear[:, None] * powers[0][None, :]
        ) % p
        counts = np.zeros((row_count, p), dtype=np.int16)
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
        if cutoff is not None:
            valid &= energy <= cutoff
        row_indices, constants = np.nonzero(valid)
        records.extend(
            (
                int(energy[row, constant]),
                leading,
                int(cubic[row]),
                int(quadratic[row]),
                int(linear[row]),
                int(constant),
            )
            for row, constant in zip(row_indices, constants)
        )
    return sorted(records)


def energy_partitions(energies: list[int], total: int) -> list[tuple[int, ...]]:
    return [
        part
        for part in itertools.combinations_with_replacement(energies, 7)
        if sum(part) == total
    ]


def rank_audit(p: int) -> dict:
    """Audit rank d+1 and kernel dimension 7-d-1 for every 7-subset."""
    forms = [form for _coordinate, form in square_directions(p)]
    histograms = {degree: Counter() for degree in range(1, 6)}
    top_full_support = True
    n_subsets = 0
    for subset in itertools.combinations(range(len(forms)), 7):
        selected = [forms[index] for index in subset]
        n_subsets += 1
        for degree in range(1, 6):
            kernel = kernel_modp(homogeneous_matrix(selected, degree, p), p)
            histograms[degree][len(kernel)] += 1
            if degree == 5 and (
                len(kernel) != 1 or np.any(np.asarray(kernel[0]) == 0)
            ):
                top_full_support = False
    return {
        "n_direction_subsets": n_subsets,
        "kernel_dimension_histograms": {
            str(degree): {str(k): v for k, v in sorted(hist.items())}
            for degree, hist in histograms.items()
        },
        "expected_kernel_dimensions": {
            str(degree): 7 - degree - 1 for degree in range(1, 6)
        },
        "top_kernel_full_support": top_full_support,
        "translation_removes_degree_four_kernel": (
            all(hist == Counter({7 - degree - 1: n_subsets})
                for degree, hist in histograms.items())
            and top_full_support
        ),
    }


def scan_prime(p: int) -> dict:
    if p < 13:
        raise ValueError("k=7 arithmetic probe requires p>=13")
    total = (p * p - 1) // 8
    all_types = relevant_quintic_types(p)
    all_histogram = Counter(record[0] for record in all_types)
    minimum = min(all_histogram)
    cutoff = total - 6 * minimum
    relevant = [record for record in all_types if record[0] <= cutoff]
    relevant_histogram = Counter(record[0] for record in relevant)
    partitions = energy_partitions(sorted(relevant_histogram), total)
    return {
        "p": p,
        "algorithm": "exact depressed-quintic lift enumeration",
        "normal_form": "a*s^5+c*s^3+d*s^2+e*s+f",
        "normalized_total_T": total,
        "minimum_profile_energy": minimum,
        "relevant_cutoff": cutoff,
        "all_type_histogram": {
            str(energy): count for energy, count in sorted(all_histogram.items())
        },
        "relevant_type_histogram": {
            str(energy): count
            for energy, count in sorted(relevant_histogram.items())
        },
        "energy_partitions": [list(partition) for partition in partitions],
        "empty_by_energy": not partitions,
        "rank_audit": rank_audit(p),
    }


def main() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("p", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = scan_prime(args.p)
    output = args.output or Path(__file__).with_name(
        f"k7_p{args.p}_quintic_profile_probe.json"
    )
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
