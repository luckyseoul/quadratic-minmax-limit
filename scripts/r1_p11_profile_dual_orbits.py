#!/usr/bin/env python3
"""Build the finite-field data for a p=11 profile/glue theta counter.

The ten-dimensional glue-dual code is

    C = K_1 + K_2 X + K_3 X^2 + K_4 X^3,

where K_d is the kernel of the six d-th powers of the square projective
directions.  Plane translations act triangularly on these four coefficient
blocks.  This script constructs exact bases and audits the orbit counts used
to reduce 11^10 dual words before a GPU theta calculation.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from math import comb
from pathlib import Path
from itertools import product

P = 11
DIRECTIONS = ((0, 1), (1, 0), (1, 8), (1, 2), (1, 9), (1, 3))


def inv(value: int) -> int:
    return pow(value % P, P - 2, P)


def rref(matrix: list[list[int]]) -> tuple[list[list[int]], list[int]]:
    rows = [[value % P for value in row] for row in matrix]
    if not rows:
        return rows, []
    width = len(rows[0])
    pivots: list[int] = []
    pivot_row = 0
    for column in range(width):
        selected = next(
            (index for index in range(pivot_row, len(rows)) if rows[index][column]),
            None,
        )
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        scale = inv(rows[pivot_row][column])
        rows[pivot_row] = [(scale * value) % P for value in rows[pivot_row]]
        for index, row in enumerate(rows):
            if index == pivot_row or not row[column]:
                continue
            multiple = row[column]
            rows[index] = [
                (value - multiple * pivot_value) % P
                for value, pivot_value in zip(row, rows[pivot_row])
            ]
        pivots.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return rows, pivots


def nullspace(matrix: list[list[int]]) -> list[list[int]]:
    reduced, pivots = rref(matrix)
    width = len(matrix[0])
    free = [column for column in range(width) if column not in pivots]
    basis: list[list[int]] = []
    for free_column in free:
        vector = [0] * width
        vector[free_column] = 1
        for row_index, pivot_column in enumerate(pivots):
            vector[pivot_column] = -reduced[row_index][free_column] % P
        basis.append(vector)
    return basis


def rank(matrix: list[list[int]]) -> int:
    return len(rref(matrix)[1])


def direction_power_matrix(degree: int) -> list[list[int]]:
    """Columns are coefficient vectors of (a X + b Y)^degree."""
    return [
        [
            comb(degree, x_degree)
            * pow(a, x_degree, P)
            * pow(b, degree - x_degree, P)
            % P
            for a, b in DIRECTIONS
        ]
        for x_degree in range(degree + 1)
    ]


def kernel_bases() -> dict[int, list[list[int]]]:
    """Return each K_d as a list of six-coordinate basis vectors."""
    return {
        degree: nullspace(direction_power_matrix(degree))
        for degree in range(1, 5)
    }


def linear_combination(basis: list[list[int]], coordinates: tuple[int, ...]) -> list[int]:
    return [
        sum(coefficient * vector[index] for coefficient, vector in zip(coordinates, basis))
        % P
        for index in range(6)
    ]


def coordinates_in_basis(basis: list[list[int]], vector: list[int]) -> tuple[int, ...]:
    """Solve basis-columns * coordinates = vector over F_11."""
    dimension = len(basis)
    augmented = [
        [basis[column][row] for column in range(dimension)] + [vector[row] % P]
        for row in range(6)
    ]
    reduced, pivots = rref(augmented)
    if pivots[:dimension] != list(range(dimension)):
        raise ArithmeticError("basis coordinate solve is singular")
    if any(not any(row[:dimension]) and row[-1] for row in reduced):
        raise ArithmeticError("vector is outside the claimed kernel")
    return tuple(reduced[index][-1] for index in range(dimension))


def multiply_by_translation(
    basis_source: list[list[int]],
    basis_target: list[list[int]],
    translation: tuple[int, int],
) -> list[list[int]]:
    """Matrix for c_j -> t_j(translation)c_j between adjacent kernels."""
    u, v = translation
    values = [(a * u + b * v) % P for a, b in DIRECTIONS]
    columns = []
    for source_vector in basis_source:
        product = [value * coefficient % P for value, coefficient in zip(values, source_vector)]
        columns.append(coordinates_in_basis(basis_target, product))
    return [
        [columns[column][row] for column in range(len(columns))]
        for row in range(len(basis_target))
    ]


def projective_representatives(dimension: int) -> list[tuple[int, ...]]:
    representatives: list[tuple[int, ...]] = []
    for first_nonzero in range(dimension):
        suffix_length = dimension - first_nonzero - 1
        for suffix in product(range(P), repeat=suffix_length):
            representatives.append(
                (0,) * first_nonzero + (1,) + tuple(suffix)
            )
    expected = (P**dimension - 1) // (P - 1)
    if len(representatives) != expected:
        raise ArithmeticError("projective representative census failed")
    return representatives


def translation_image(
    bases: dict[int, list[list[int]]],
    degree: int,
    source_coordinates: tuple[int, ...],
) -> list[list[int]]:
    """Columns are multiplication by the two coordinate translations."""
    source = linear_combination(bases[degree], source_coordinates)
    columns: list[tuple[int, ...]] = []
    for u, v in ((1, 0), (0, 1)):
        values = [(a * u + b * v) % P for a, b in DIRECTIONS]
        product_vector = [
            value * coefficient % P
            for value, coefficient in zip(values, source)
        ]
        columns.append(coordinates_in_basis(bases[degree - 1], product_vector))
    return [
        [columns[column][row] for column in range(2)]
        for row in range(len(bases[degree - 1]))
    ]


def column_complement(matrix: list[list[int]], ambient_dimension: int) -> list[tuple[int, ...]]:
    columns = [
        tuple(matrix[row][column] for row in range(ambient_dimension))
        for column in range(len(matrix[0]))
    ]
    independent: list[tuple[int, ...]] = []
    for column in columns:
        candidate = independent + [column]
        if rank([[value[row] for value in candidate] for row in range(ambient_dimension)]) > len(independent):
            independent.append(column)
    image_rank = len(independent)
    complement: list[tuple[int, ...]] = []
    for index in range(ambient_dimension):
        unit = tuple(1 if row == index else 0 for row in range(ambient_dimension))
        candidate = independent + complement + [unit]
        if rank([[value[row] for value in candidate] for row in range(ambient_dimension)]) > len(independent) + len(complement):
            complement.append(unit)
    if len(complement) != ambient_dimension - image_rank:
        raise ArithmeticError("failed to construct a translation-image complement")
    return complement


def polynomial_histogram(coefficients: tuple[int, int, int, int]) -> tuple[int, ...]:
    counts = [0] * P
    for value in range(P):
        power = value
        image = 0
        for coefficient in coefficients:
            image = (image + coefficient * power) % P
            power = power * value % P
        counts[image] += 1
    return tuple(counts)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    bases = kernel_bases()
    dimensions = {degree: len(basis) for degree, basis in bases.items()}
    if dimensions != {1: 4, 2: 3, 3: 2, 4: 1}:
        raise ArithmeticError(f"wrong glue dimensions: {dimensions}")

    multiplication_ranks: dict[str, int] = {}
    for degree in range(2, 5):
        columns_by_translation = []
        for translation in ((1, 0), (0, 1)):
            matrix = multiply_by_translation(
                bases[degree], bases[degree - 1], translation
            )
            columns_by_translation.extend(
                [[matrix[row][column] for row in range(len(matrix))] for column in range(len(matrix[0]))]
            )
        multiplication_ranks[str(degree)] = rank(
            [[column[row] for column in columns_by_translation] for row in range(len(bases[degree - 1]))]
        )

    histogram_counts = Counter(
        polynomial_histogram((c1, c2, c3, c4))
        for c1 in range(P)
        for c2 in range(P)
        for c3 in range(P)
        for c4 in range(P)
    )
    if sum(histogram_counts.values()) != P**4:
        raise ArithmeticError("quartic histogram census lost coefficient tuples")

    # Translation/scalar representative census by highest nonzero degree.
    h4_image = translation_image(bases, 4, (1,))
    h4_image_rank = rank(h4_image)
    h3_projective = projective_representatives(2)
    h3_image_ranks = Counter(
        rank(translation_image(bases, 3, coordinates))
        for coordinates in h3_projective
    )
    h2_projective = projective_representatives(3)
    h2_image_ranks = Counter(
        rank(translation_image(bases, 2, coordinates))
        for coordinates in h2_projective
    )
    if h4_image_rank != 2 or h3_image_ranks != {2: P + 1}:
        raise ArithmeticError("top translation action is not free")
    h4_representatives = P**7
    h3_representatives = (P + 1) * P**5
    projective_k2 = (P**3 - 1) // (P - 1)
    h2_representatives = sum(
        count * P ** (4 - image_rank)
        for image_rank, count in h2_image_ranks.items()
    )
    h2_weighted_codewords = projective_k2 * (P - 1) * P**4
    h1_representatives = (P**4 - 1) // (P - 1)
    reconstructed = (
        1
        + h4_representatives * P**2 * (P - 1)
        + h3_representatives * P**2 * (P - 1)
        + h2_weighted_codewords
        + h1_representatives * (P - 1)
    )
    if reconstructed != P**10:
        raise ArithmeticError("translation/scalar strata do not reconstruct 11^10")

    report = {
        "experiment": "r1_p11_profile_dual_orbits",
        "status": "exact_finite_field_data_built",
        "p": P,
        "directions": [list(direction) for direction in DIRECTIONS],
        "kernel_dimensions": dimensions,
        "kernel_bases": {str(degree): basis for degree, basis in bases.items()},
        "adjacent_translation_product_span_ranks": multiplication_ranks,
        "quartic_coefficient_tuples": P**4,
        "quartic_value_distribution_types": len(histogram_counts),
        "largest_value_distribution_fibre": max(histogram_counts.values()),
        "translation_scalar_representatives": {
            "highest_degree_4": h4_representatives,
            "highest_degree_3": h3_representatives,
            "highest_degree_2": h2_representatives,
            "highest_degree_1": h1_representatives,
            "degree_2_projective_classes": projective_k2,
            "degree_2_codewords_after_weighting": h2_weighted_codewords,
        },
        "translation_image_ranks": {
            "highest_degree_4": h4_image_rank,
            "highest_degree_3": dict(sorted(h3_image_ranks.items())),
            "highest_degree_2": dict(sorted(h2_image_ranks.items())),
        },
        "total_nonzero_translation_scalar_representatives": (
            h4_representatives
            + h3_representatives
            + h2_representatives
            + h1_representatives
        ),
        "codeword_count_reconstructed": reconstructed,
        "expected_codeword_count": P**10,
    }
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
