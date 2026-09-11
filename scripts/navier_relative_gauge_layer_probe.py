#!/usr/bin/env python3
"""Exact connected-layer scout for relative-gauge conditional densities.

For random two-block signings this computes, in rational arithmetic, the
coefficients of t^4, t^6 and t^8 in log E cosh(t Q_Y) for every projective
relative gauge.  The block-local and rectangular normalizers are gauge
invariant, so collisions in the displayed coefficients are collisions in the
corresponding conditional log-density layers.

It is a finite witness search only.  A witness with equal H4 and unequal H6
proves that the mixed four-cycle layer alone is not a closed state for that
exact block; it does not establish an all-orders composition theorem.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import random
from fractions import Fraction
from math import factorial


def random_symmetric(order: int, generator: random.Random) -> list[list[int]]:
    matrix = [[0] * order for _ in range(order)]
    for row in range(order):
        for column in range(row + 1, order):
            value = generator.choice((-1, 1))
            matrix[row][column] = matrix[column][row] = value
    return matrix


def random_rectangular(rows: int, columns: int, generator: random.Random) -> list[list[int]]:
    return [[generator.choice((-1, 1)) for _ in range(columns)] for _ in range(rows)]


def projective_signs(order: int):
    for tail in itertools.product((-1, 1), repeat=order - 1):
        yield (1,) + tail


def log_even_partition_coefficients(matrix: list[list[int]], maximum_layer: int = 4) -> tuple[Fraction, ...]:
    order = len(matrix)
    values = []
    for spin in projective_signs(order):
        values.append(sum(
            matrix[row][column] * spin[row] * spin[column]
            for row in range(order) for column in range(row + 1, order)
        ))
    count = len(values)
    series = [Fraction(1)]
    for layer in range(1, maximum_layer + 1):
        series.append(Fraction(sum(value ** (2 * layer) for value in values), count * factorial(2 * layer)))
    logarithm = [Fraction(0)] * (maximum_layer + 1)
    for layer in range(1, maximum_layer + 1):
        logarithm[layer] = series[layer] - sum(
            Fraction(previous, layer) * logarithm[previous] * series[layer - previous]
            for previous in range(1, layer)
        )
    return tuple(logarithm[1:])


def switched_block(
    left: list[list[int]], right: list[list[int]], cross: list[list[int]],
    alpha: tuple[int, ...], beta: tuple[int, ...], tau: int,
) -> list[list[int]]:
    left_order, right_order = len(left), len(right)
    output = [[0] * (left_order + right_order) for _ in range(left_order + right_order)]
    for row in range(left_order):
        for column in range(row + 1, left_order):
            output[row][column] = output[column][row] = alpha[row] * left[row][column] * alpha[column]
    for row in range(right_order):
        for column in range(row + 1, right_order):
            value = tau * beta[row] * right[row][column] * beta[column]
            output[left_order + row][left_order + column] = value
            output[left_order + column][left_order + row] = value
    for row in range(left_order):
        for column in range(right_order):
            output[row][left_order + column] = output[left_order + column][row] = cross[row][column]
    return output


def matrix_bits(left: list[list[int]], right: list[list[int]], cross: list[list[int]]) -> str:
    values = []
    for matrix in (left, right):
        values.extend(matrix[row][column] for row in range(len(matrix)) for column in range(row + 1, len(matrix)))
    values.extend(value for row in cross for value in row)
    return ''.join('+' if value > 0 else '-' for value in values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed', type=int, required=True)
    parser.add_argument('--order', type=int, default=5)
    arguments = parser.parse_args()
    if arguments.order < 3 or arguments.order > 6:
        raise ValueError('this exact scout supports 3 <= order <= 6')
    generator = random.Random(arguments.seed)
    left = random_symmetric(arguments.order, generator)
    right = random_symmetric(arguments.order, generator)
    cross = random_rectangular(arguments.order, arguments.order, generator)
    seen_h4: dict[Fraction, tuple[Fraction, Fraction, tuple[tuple[int, ...], tuple[int, ...], int]]] = {}
    seen_h46: dict[tuple[Fraction, Fraction], tuple[Fraction, tuple[tuple[int, ...], tuple[int, ...], int]]] = {}
    h4_h6_collision = None
    h46_h8_collision = None
    layer_count = 0
    for alpha in projective_signs(arguments.order):
        for beta in projective_signs(arguments.order):
            for tau in (-1, 1):
                layers = log_even_partition_coefficients(switched_block(left, right, cross, alpha, beta, tau))
                h4, h6, h8 = layers[1:4]
                gauge = (alpha, beta, tau)
                previous = seen_h4.get(h4)
                if previous is not None and previous[0] != h6 and h4_h6_collision is None:
                    h4_h6_collision = (previous[2], gauge, h4, previous[0], h6)
                else:
                    seen_h4.setdefault(h4, (h6, h8, gauge))
                previous_46 = seen_h46.get((h4, h6))
                if previous_46 is not None and previous_46[0] != h8 and h46_h8_collision is None:
                    h46_h8_collision = (previous_46[1], gauge, h4, h6, previous_46[0], h8)
                else:
                    seen_h46.setdefault((h4, h6), (h8, gauge))
                layer_count += 1
    bits = matrix_bits(left, right, cross)
    record = {
        'kind': 'finite_exact_relative_gauge_connected_layer_scout',
        'seed': arguments.seed,
        'block_order': arguments.order,
        'gauges_checked': layer_count,
        'input_sign_digest': hashlib.sha256(bits.encode()).hexdigest(),
        'distinct_h4': len(seen_h4),
        'distinct_h4_h6': len(seen_h46),
        'h4_h6_collision': h4_h6_collision is not None,
        'h4_h6_witness': str(h4_h6_collision) if h4_h6_collision else None,
        'h46_h8_collision': h46_h8_collision is not None,
        'h46_h8_witness': str(h46_h8_collision) if h46_h8_collision else None,
    }
    print(json.dumps(record, sort_keys=True))


if __name__ == '__main__':
    main()
