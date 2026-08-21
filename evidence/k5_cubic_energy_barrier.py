#!/usr/bin/env python3
"""Exact finite certificate for emptiness of the k=5 stratum at p>=41.

For a genuinely cubic active profile, translate the input to write its
residue polynomial as

    f(s) = a s^3 + c s + d,  a != 0.

Translation preserves its value multiplicities and hence its lift energy.
The centered lift can only be changed at residue (p-1)/2: replacing that
value by -(p+1)/2 lowers the sum by p and raises the energy by p.  This script
exhausts (a,c,d), performs exactly the number of replacements required for
zero sum, and records b=||h||^2/(2p).
"""
from __future__ import annotations

import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
FINITE_PRIMES = (41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97)
EXPECTED_MINIMA = {
    41: 43,
    43: 45,
    47: 58,
    53: 77,
    59: 97,
    61: 99,
    67: 129,
    71: 144,
    73: 153,
    79: 181,
    83: 210,
    89: 244,
    97: 288,
}


def cubic_profile_types(p: int) -> list[tuple[int, int, int, int]]:
    """Return all admissible depressed-cubic types (b,a,c,d)."""
    midpoint = (p - 1) // 2
    s = np.arange(p, dtype=np.int64)
    leading = np.repeat(np.arange(1, p, dtype=np.int64), p)
    linear = np.tile(np.arange(p, dtype=np.int64), p - 1)
    values = (
        leading[:, None] * s[None, :] ** 3
        + linear[:, None] * s[None, :]
    ) % p

    histogram = np.zeros(((p - 1) * p, p), dtype=np.int16)
    rows = np.repeat(np.arange(len(leading)), p)
    np.add.at(histogram, (rows, values.ravel()), 1)

    centered = np.where(s <= midpoint, s, s - p)
    shifted = np.stack(
        [centered[(s + constant) % p] for constant in range(p)]
    )
    sums = histogram @ shifted.T
    energies = histogram @ (shifted * shifted).T
    endpoint_counts = np.stack(
        [histogram[:, (midpoint - constant) % p] for constant in range(p)],
        axis=1,
    )
    replacements = sums // p
    valid = (
        (sums % p == 0)
        & (replacements >= 0)
        & (replacements <= endpoint_counts)
    )
    normalized = (energies + p * replacements) // (2 * p)

    row_index, constants = np.nonzero(valid)
    return sorted(
        (
            int(normalized[row, constant]),
            int(leading[row]),
            int(linear[row]),
            int(constant),
        )
        for row, constant in zip(row_index, constants)
    )


def energy_partitions(energies: list[int], total: int) -> list[list[int]]:
    return [
        list(partition)
        for partition in itertools.combinations_with_replacement(energies, 5)
        if sum(partition) == total
    ]


def scan_prime(p: int) -> dict:
    types = cubic_profile_types(p)
    minimum = min(record[0] for record in types)
    total = (p * p - 1) // 8
    maximum_relevant = total - 4 * minimum
    relevant = [record for record in types if record[0] <= maximum_relevant]
    histogram = Counter(record[0] for record in relevant)
    partitions = energy_partitions(sorted(histogram), total)
    result = {
        "minimum_cubic_b": minimum,
        "normalized_total_T": total,
        "five_minima_exceed_total": 5 * minimum > total,
        "maximum_relevant_b": maximum_relevant,
        "relevant_type_histogram": {
            str(energy): count for energy, count in sorted(histogram.items())
        },
        "energy_partitions": partitions,
        "k5_empty_by_energy": not partitions,
    }
    if p == 43:
        result["relevant_types"] = [
            {"b": b, "leading": a, "linear": c, "constant": d}
            for b, a, c, d in relevant
        ]
    return result


def main() -> dict:
    finite = {str(p): scan_prime(p) for p in FINITE_PRIMES}
    report = {
        "scope": "k=5 and every prime p>=41",
        "normalization": "b=||h||^2/(2p), T=(p^2-1)/8",
        "finite_exact_range": finite,
        "analytic_range": {
            "first_prime": 101,
            "reason": "general activity barrier p>4k^2 with k=5",
        },
        "proved": all(
            row["minimum_cubic_b"] == EXPECTED_MINIMA[int(p)]
            and row["k5_empty_by_energy"]
            for p, row in finite.items()
        ),
    }
    p43 = finite["43"]
    assert p43["relevant_type_histogram"] == {"45": 28}
    assert p43["energy_partitions"] == []
    assert report["proved"]
    output = REPO / "evidence" / "k5_cubic_energy_barrier.json"
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
