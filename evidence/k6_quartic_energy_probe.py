#!/usr/bin/env python3
"""Exact depressed-quartic profile energies for the k=6 attack."""
from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np


def quartic_profile_energies(p: int, leading_batch: int = 4) -> Counter:
    midpoint = (p - 1) // 2
    s = np.arange(p, dtype=np.int64)
    centered = np.where(s <= midpoint, s, s - p)
    shifted = np.stack(
        [centered[(s + constant) % p] for constant in range(p)]
    )
    shifted_sq = shifted * shifted
    histogram = Counter()

    for leading_start in range(1, p, leading_batch):
        leading_values = np.arange(
            leading_start, min(p, leading_start + leading_batch), dtype=np.int64
        )
        leading = np.repeat(leading_values, p * p)
        quadratic = np.tile(np.repeat(np.arange(p, dtype=np.int64), p), len(leading_values))
        linear = np.tile(np.arange(p, dtype=np.int64), len(leading_values) * p)
        values = (
            leading[:, None] * s[None, :] ** 4
            + quadratic[:, None] * s[None, :] ** 2
            + linear[:, None] * s[None, :]
        ) % p
        counts = np.zeros((len(leading), p), dtype=np.int16)
        rows = np.repeat(np.arange(len(leading)), p)
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
        energy = standard_energy + p * replacements
        if np.any(energy[valid] % (2 * p)):
            raise RuntimeError("quartic profile energy not divisible by 2p")
        normalized = energy[valid] // (2 * p)
        values_b, frequencies = np.unique(normalized, return_counts=True)
        histogram.update(
            {
                int(value): int(frequency)
                for value, frequency in zip(values_b, frequencies)
            }
        )
    return histogram


def scan_prime(p: int, leading_batch: int = 4) -> dict:
    histogram = quartic_profile_energies(p, leading_batch)
    minimum = min(histogram)
    total = (p * p - 1) // 8
    maximum_relevant = total - 5 * minimum
    relevant = sorted(value for value in histogram if value <= maximum_relevant)
    partitions = [
        list(partition)
        for partition in itertools.combinations_with_replacement(relevant, 6)
        if sum(partition) == total
    ]
    return {
        "p": p,
        "normalized_total_T": total,
        "minimum_quartic_b": minimum,
        "six_minima_minus_total": 6 * minimum - total,
        "maximum_relevant_b": maximum_relevant,
        "relevant_type_histogram": {
            str(value): histogram[value] for value in relevant
        },
        "energy_partitions": partitions,
        "nonzero_top_scalar_empty_by_energy": not partitions,
    }


def main() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument("primes", nargs="+", type=int)
    parser.add_argument("--leading-batch", type=int, default=4)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {}
    for p in args.primes:
        report[str(p)] = scan_prime(p, args.leading_batch)
        print(json.dumps(report[str(p)], indent=2), flush=True)
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    main()
