#!/usr/bin/env python3
"""Exact k=4 closure for p=1 (mod 4) below the analytic cutoff.

At p=29,37 the coefficient sieve has no candidates.  At p=13,17 it
regenerates the complete eps=+1 k=4 families and evaluates the Gaussian-
integer quartic pair statistic exactly.
"""
from __future__ import annotations

import itertools
import json
import sys
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "evidence"))

from e1_gmin_m4_prop15588 import directions, field_ctx  # noqa: E402
from e1_gmin_m4_prop15589 import quadratic_profile_min_b  # noqa: E402
from k4_p3mod4_coefficient_sieve import (  # noqa: E402
    homogeneous_matrix,
    kernel_modp,
    quadratic_types,
    scan_prime,
)


def primitive_root(q: int, mul) -> int:
    for generator in range(2, q):
        x = 1
        seen = set()
        for _ in range(q - 1):
            seen.add(x)
            x = mul(x, generator)
        if len(seen) == q - 1:
            return generator
    raise RuntimeError("no primitive root")


def quartic_kernel_parts(p: int) -> tuple[np.ndarray, np.ndarray]:
    q, mul, _chi, _trace = field_ctx(p)
    generator = primitive_root(q, mul)
    real = np.zeros(q, dtype=np.int8)
    imag = np.zeros(q, dtype=np.int8)
    units = ((1, 0), (0, 1), (-1, 0), (0, -1))
    x = 1
    for exponent in range(q - 1):
        real[x], imag[x] = units[exponent % 4]
        x = mul(x, generator)

    kernel_real = np.zeros((q, q), dtype=np.float32)
    kernel_imag = np.zeros((q, q), dtype=np.float32)
    for a in range(q):
        for b in range(q):
            difference = (a % p - b % p) % p + ((a // p - b // p) % p) * p
            kernel_real[a, b] = real[difference]
            kernel_imag[a, b] = imag[difference]
    return kernel_real, kernel_imag


def generate_k4(p: int) -> tuple[np.ndarray, dict[str, int], dict]:
    q = p * p
    midpoint = (p - 1) // 2
    normalized_total = (q - 1) // 8
    types = quadratic_types(p)
    minimum = quadratic_profile_min_b(p)
    maximum_relevant = normalized_total - 3 * minimum
    relevant = {
        sign: [record for record in records if record["b"] <= maximum_relevant]
        for sign, records in types.items()
    }
    square, _ = directions(p)
    forms = [form for _t_of, form in square]
    coordinates = [t_of for t_of, _form in square]
    subsets = list(itertools.combinations(range(len(square)), 4))
    parameters = np.asarray(list(itertools.product(range(p), repeat=2)), dtype=np.int64)
    chi = np.asarray(
        [0]
        + [
            1 if pow(value, (p - 1) // 2, p) == 1 else -1
            for value in range(1, p)
        ],
        dtype=np.int8,
    )
    required_constant_sum = (-5 * pow(2, p - 2, p)) % p
    s = np.arange(p, dtype=np.int64)
    rows = []
    subset_counts = {}
    coefficient_candidates = endpoint_branches = 0

    for subset in subsets:
        selected_forms = [forms[j] for j in subset]
        selected_coordinates = np.stack([coordinates[j] for j in subset])
        top = kernel_modp(homogeneous_matrix(selected_forms, 2, p), p)[0]
        linear_kernel = kernel_modp(homogeneous_matrix(selected_forms, 1, p), p)
        linear = parameters @ np.asarray(linear_kernel) % p
        before = len(rows)

        for scalar in range(1, p):
            leading = scalar * top % p
            inverse_4a = np.asarray(
                [pow(int(4 * value % p), p - 2, p) for value in leading]
            )
            square_completion = linear * linear * inverse_4a % p
            choices = [relevant[int(chi[value])] for value in leading]
            for combination in itertools.product(*choices):
                if sum(record["b"] for record in combination) != normalized_total:
                    continue
                completed = np.asarray(
                    [record["completed_constant"] for record in combination]
                )
                constants = (square_completion + completed) % p
                valid = (
                    np.sum((constants + midpoint) % p, axis=1) % p
                    == required_constant_sum
                )
                for degree_one, constant in zip(linear[valid], constants[valid]):
                    coefficient_candidates += 1
                    polynomial = (
                        leading[:, None] * s[None, :] ** 2
                        + degree_one[:, None] * s[None, :]
                        + constant[:, None]
                    ) % p
                    centered = np.where(
                        polynomial <= midpoint, polynomial, polynomial - p
                    )
                    flip_choices = [
                        list(
                            itertools.combinations(
                                np.where(polynomial[j] == midpoint)[0],
                                record["endpoint_flips"],
                            )
                        )
                        for j, record in enumerate(combination)
                    ]
                    for flips in itertools.product(*flip_choices):
                        endpoint_branches += 1
                        h = centered.copy()
                        for j, selected in enumerate(flips):
                            if selected:
                                h[j, list(selected)] -= p
                        if np.any(h.sum(axis=1)):
                            raise RuntimeError("profile zero-sum audit failed")
                        total = sum(
                            h[j][selected_coordinates[j]] for j in range(4)
                        )
                        if np.all((total == midpoint) | (total == -midpoint - 1)):
                            rows.append(
                                np.where(total == midpoint, 1, -1).astype(np.int8)
                            )

        if len(rows) > before:
            subset_counts[",".join(map(str, subset))] = len(rows) - before

    return (
        np.stack(rows),
        subset_counts,
        {
            "coefficient_candidates": coefficient_candidates,
            "endpoint_branches": endpoint_branches,
        },
    )


def quartic_moment(p: int, rows: np.ndarray, chunk: int = 5_000) -> dict:
    kernel_real, kernel_imag = quartic_kernel_parts(p)
    sum_abs2 = 0
    histogram = Counter()
    for lo in range(0, len(rows), chunk):
        finite = ((1 - rows[lo : lo + chunk]) // 2).astype(np.float32)
        real = np.sum((finite @ kernel_real) * finite, axis=1).astype(np.int64)
        imag = np.sum((finite @ kernel_imag) * finite, axis=1).astype(np.int64)
        abs2 = real * real + imag * imag
        sum_abs2 += int(abs2.sum())
        values, counts = np.unique(abs2, return_counts=True)
        histogram.update({int(value): int(count) for value, count in zip(values, counts)})
    moment = Fraction(sum_abs2, len(rows))
    threshold = Fraction(3 * p * p * (p * p - 1), 16)
    return {
        "E_abs_Zpsi_sq": str(moment),
        "QVAR_threshold": str(threshold),
        "clears_QVAR": moment >= threshold,
        "abs_Zpsi_sq_histogram": {
            str(value): count for value, count in sorted(histogram.items())
        },
    }


def main() -> dict:
    report = {"empty": {}, "nonempty": {}}
    for p in (29, 37):
        report["empty"][str(p)] = scan_prime(p)
    for p in (13, 17):
        rows, subset_counts, audit = generate_k4(p)
        report["nonempty"][str(p)] = {
            "count_eps_plus": len(rows),
            "n_nonzero_direction_subsets": len(subset_counts),
            "vectors_per_nonzero_subset": sorted(set(subset_counts.values())),
            "generation_audit": audit,
            "quartic": quartic_moment(p, rows),
        }

    if not (
        report["empty"]["29"]["coefficient_candidate_histogram"] == {"0": 1365}
        and report["empty"]["37"]["coefficient_candidate_histogram"] == {"0": 3876}
        and report["nonempty"]["13"]["count_eps_plus"] == 28_392
        and report["nonempty"]["13"]["n_nonzero_direction_subsets"] == 7
        and report["nonempty"]["13"]["vectors_per_nonzero_subset"] == [4_056]
        and report["nonempty"]["13"]["quartic"]["E_abs_Zpsi_sq"] == "8788"
        and report["nonempty"]["17"]["count_eps_plus"] == 62_424
        and report["nonempty"]["17"]["n_nonzero_direction_subsets"] == 27
        and report["nonempty"]["17"]["vectors_per_nonzero_subset"] == [2_312]
        and report["nonempty"]["17"]["quartic"]["E_abs_Zpsi_sq"] == "314432/3"
        and all(
            record["quartic"]["clears_QVAR"]
            for record in report["nonempty"].values()
        )
    ):
        raise RuntimeError("p=1 mod 4 k=4 closure audit failed")

    output = Path(__file__).with_suffix(".json")
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return report


if __name__ == "__main__":
    main()
