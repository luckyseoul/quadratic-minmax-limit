#!/usr/bin/env python3
"""Search exact low-degree relations in the order-21 F_8 boundary data.

The signature scanner records the three nonsquare residues modulo the
reciprocal factors 0x57 and 0x75.  Dividing by X and by the first nonzero
coordinate makes each nonzero triple a point of PG(2,8).  This script moves
the 0x75 coordinates into the 0x57 model via beta -> alpha^-1, then tests
projective coverage, Frobenius-twisted pairings, conics, and affine
low-degree interpolation with a chronological holdout.

All arithmetic is exact.  The search is reconnaissance: a fitted relation
is reported separately from holdout validation and is never called a proof.
"""
from __future__ import annotations

import argparse
import collections
import functools
import itertools
import json
from pathlib import Path


COMMON_MODULUS = 0x57
RECIPROCAL_MODULUS = 0x75
ZERO = (0, 0, 0)


def multiply(left: int, right: int, modulus: int = COMMON_MODULUS) -> int:
    """Multiply packed quotient-field elements over F_2."""
    result = 0
    degree = modulus.bit_length() - 1
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left >> degree:
            left ^= modulus
    return result


def power(value: int, exponent: int, modulus: int = COMMON_MODULUS) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply(result, value, modulus)
        value = multiply(value, value, modulus)
        exponent >>= 1
    return result


def reciprocal_to_common(value: int) -> int:
    """Evaluate the packed polynomial value at alpha^-1 modulo 0x57."""
    alpha_inverse = power(0x2, 62)
    result = 0
    term = 1
    for exponent in range(6):
        if value & (1 << exponent):
            result ^= term
        term = multiply(term, alpha_inverse)
    return result


def normalize(vector: tuple[int, int, int]) -> tuple[int, int, int]:
    pivot = next((value for value in vector if value), None)
    if pivot is None:
        return ZERO
    inverse = power(pivot, 6)
    return tuple(multiply(value, inverse) for value in vector)


def projective_points() -> set[tuple[int, int, int]]:
    return {
        normalize(vector)
        for vector in itertools.product(subfield_elements(), repeat=3)
        if vector != ZERO
    }


def subfield_elements() -> tuple[int, ...]:
    return tuple(value for value in range(64) if power(value, 8) == value)


def dot(
    left: tuple[int, int, int],
    right: tuple[int, int, int],
) -> int:
    return (
        multiply(left[0], right[0])
        ^ multiply(left[1], right[1])
        ^ multiply(left[2], right[2])
    )


def monomial_exponents(n_variables: int, max_degree: int):
    """Yield exponent tuples of total degree 1..max_degree."""
    for degree in range(1, max_degree + 1):
        for bars in itertools.combinations_with_replacement(
            range(n_variables), degree
        ):
            exponents = [0] * n_variables
            for index in bars:
                exponents[index] += 1
            yield tuple(exponents)


def monomial_value(values: tuple[int, ...], exponents: tuple[int, ...]) -> int:
    result = 1
    for value, exponent in zip(values, exponents):
        if exponent:
            result = multiply(result, power(value, exponent))
    return result


def solve_linear(matrix: list[list[int]], rhs: list[int]) -> list[int] | None:
    """Return one exact GF(8) solution, or None if inconsistent."""
    if not matrix:
        return []
    n_columns = len(matrix[0])
    augmented = [row[:] + [value] for row, value in zip(matrix, rhs)]
    pivots: list[int] = []
    pivot_row = 0
    for column in range(n_columns):
        found = next(
            (
                row
                for row in range(pivot_row, len(augmented))
                if augmented[row][column]
            ),
            None,
        )
        if found is None:
            continue
        augmented[pivot_row], augmented[found] = (
            augmented[found],
            augmented[pivot_row],
        )
        inverse = power(augmented[pivot_row][column], 6)
        augmented[pivot_row] = [
            multiply(value, inverse) for value in augmented[pivot_row]
        ]
        for row in range(len(augmented)):
            if row == pivot_row or not augmented[row][column]:
                continue
            scale = augmented[row][column]
            augmented[row] = [
                value ^ multiply(scale, pivot_value)
                for value, pivot_value in zip(
                    augmented[row], augmented[pivot_row]
                )
            ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(augmented):
            break
    if any(
        not any(row[:n_columns]) and row[n_columns] for row in augmented
    ):
        return None
    solution = [0] * n_columns
    for row, column in enumerate(pivots):
        solution[column] = augmented[row][n_columns]
    return solution


def evaluate_polynomial(
    values: tuple[int, ...],
    exponents: list[tuple[int, ...]],
    coefficients: list[int],
) -> int:
    result = 0
    for exponent, coefficient in zip(exponents, coefficients):
        if coefficient:
            result ^= multiply(coefficient, monomial_value(values, exponent))
    return result


def quadratic_value(vector: tuple[int, int, int], coefficients: tuple[int, ...]):
    x, y, z = vector
    terms = (
        multiply(x, x),
        multiply(y, y),
        multiply(z, z),
        multiply(x, y),
        multiply(x, z),
        multiply(y, z),
    )
    result = 0
    for coefficient, term in zip(coefficients, terms):
        result ^= multiply(coefficient, term)
    return result


def canonical_coefficients(coefficients: tuple[int, ...]) -> tuple[int, ...]:
    pivot = next(value for value in coefficients if value)
    inverse = power(pivot, 6)
    return tuple(multiply(value, inverse) for value in coefficients)


def projective_coefficients(field: tuple[int, ...], dimension: int):
    """Yield one normalized representative of every projective coefficient."""
    for pivot in range(dimension):
        for tail in itertools.product(field, repeat=dimension - pivot - 1):
            yield (0,) * pivot + (1,) + tail


def load_rows(path: Path):
    data = json.loads(path.read_text())
    ratio_by_prime = {
        row["p"]: row["component_ratios"]
        for row in data["component_ratio_rows"]
    }
    rows = []
    for row in data["nonsquare_projective_rows"]:
        left = tuple(row["vectors"][0])
        right = normalize(
            tuple(reciprocal_to_common(value) for value in row["vectors"][1])
        )
        ratios = ratio_by_prime[row["p"]]
        left_ratio = None if not ratios[0] else ratios[0][0]
        right_ratio = (
            None
            if not ratios[1]
            else reciprocal_to_common(ratios[1][0])
        )
        rows.append(
            {
                "p": row["p"],
                "left": left,
                "right": right,
                "state": (left_ratio, right_ratio),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--train", type=int, default=300)
    parser.add_argument("--max-degree", type=int, default=4)
    parser.add_argument("--top-conics", type=int, default=12)
    args = parser.parse_args()

    rows = load_rows(args.input)
    points = projective_points()
    left_counts = collections.Counter(row["left"] for row in rows)
    right_counts = collections.Counter(row["right"] for row in rows)
    pair_counts = collections.Counter(
        (row["left"], row["right"]) for row in rows
    )
    print(
        "coverage",
        {
            "rows": len(rows),
            "pg_points": len(points),
            "left_projective": len(set(left_counts) - {ZERO}),
            "right_projective": len(set(right_counts) - {ZERO}),
            "left_zero_rows": left_counts[ZERO],
            "right_zero_rows": right_counts[ZERO],
            "distinct_pairs": len(pair_counts),
            "missing_left_points": len(points - set(left_counts)),
            "missing_right_points": len(points - set(right_counts)),
        },
    )
    state_counts = collections.Counter(row["state"] for row in rows)
    print("states", dict(sorted(state_counts.items(), key=str)))

    twisted = []
    for left_twist in range(3):
        for right_twist in range(3):
            for permutation in itertools.permutations(range(3)):
                zero_primes = []
                for row in rows:
                    left = tuple(
                        power(value, 1 << left_twist) for value in row["left"]
                    )
                    right = tuple(
                        power(row["right"][index], 1 << right_twist)
                        for index in permutation
                    )
                    if dot(left, right) == 0:
                        zero_primes.append(row["p"])
                twisted.append(
                    (
                        len(zero_primes),
                        left_twist,
                        right_twist,
                        permutation,
                        zero_primes[:8],
                    )
                )
    print("best_twisted_pairings", sorted(twisted)[:12])

    field = subfield_elements()
    conics = []
    multiplication_table = [
        [multiply(left, right) for right in range(64)] for left in range(64)
    ]
    point_quadratic_terms = {}
    for point in points | {ZERO}:
        x, y, z = point
        point_quadratic_terms[point] = (
            multiply(x, x),
            multiply(y, y),
            multiply(z, z),
            multiply(x, y),
            multiply(x, z),
            multiply(y, z),
        )
    for canonical in projective_coefficients(field, 6):
        zero_vectors = {
            point
            for point, terms in point_quadratic_terms.items()
            if not functools.reduce(
                int.__xor__,
                (
                    multiplication_table[coefficient][term]
                    for coefficient, term in zip(canonical, terms)
                ),
                0,
            )
        }
        zero_points = {
            point for point in zero_vectors if point != ZERO
        }
        if not zero_points:
            continue
        joint_zero_primes = [
            row["p"]
            for row in rows
            if row["left"] in zero_vectors and row["right"] in zero_vectors
        ]
        conics.append(
            (
                len(joint_zero_primes),
                len(zero_points),
                sum(1 for value in canonical if value),
                canonical,
                joint_zero_primes[:8],
            )
        )
    print("best_common_quadrics")
    for result in sorted(conics)[: args.top_conics]:
        print(result)
    conic_size_counts = collections.Counter(result[1] for result in conics)
    print("quadratic_zero_set_size_counts", dict(sorted(conic_size_counts.items())))
    print("best_nondegenerate_conics")
    for result in sorted(item for item in conics if item[1] == 9)[
        : args.top_conics
    ]:
        print(result)

    train_count = min(args.train, len(rows))
    for degree in range(1, args.max_degree + 1):
        exponents = list(monomial_exponents(6, degree))
        train_matrix = [
            [
                monomial_value(row["left"] + row["right"], exponent)
                for exponent in exponents
            ]
            for row in rows[:train_count]
        ]
        solution = solve_linear(train_matrix, [1] * train_count)
        if solution is None:
            print(
                "affine_fit",
                {
                    "degree": degree,
                    "features": len(exponents),
                    "train": train_count,
                    "consistent": False,
                },
            )
            continue
        train_failures = []
        holdout_failures = []
        for index, row in enumerate(rows):
            value = evaluate_polynomial(
                row["left"] + row["right"], exponents, solution
            )
            if value != 1:
                (train_failures if index < train_count else holdout_failures).append(
                    row["p"]
                )
        print(
            "affine_fit",
            {
                "degree": degree,
                "features": len(exponents),
                "train": train_count,
                "consistent": True,
                "nonzero_coefficients": sum(value != 0 for value in solution),
                "train_failures": train_failures,
                "holdout_failure_count": len(holdout_failures),
                "holdout_failure_examples": holdout_failures[:12],
                "terms": (
                    [
                        {
                            "coefficient": hex(coefficient),
                            "exponents": exponent,
                        }
                        for exponent, coefficient in zip(exponents, solution)
                        if coefficient
                    ]
                    if sum(value != 0 for value in solution) <= 12
                    else None
                ),
            },
        )


if __name__ == "__main__":
    main()
