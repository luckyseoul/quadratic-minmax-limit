"""Exact arithmetic in ``(Z/4Z)[X]/(Phi_21)`` for W2 probes."""
from __future__ import annotations


MODULUS = 4
# Phi_21 = x^12-x^11+x^9-x^8+x^6-x^4+x^3-x+1.
PHI21 = [1, -1, 0, 1, -1, 0, 1, 0, -1, 1, 0, -1, 1]
DECOMPOSITION_GROUP = (1, 2, 4, 8, 16, 11)


def ring_reduce(coefficients: list[int]) -> list[int]:
    values = [value % MODULUS for value in coefficients]
    if len(values) < 13:
        values.extend([0] * (13 - len(values)))
    for exponent in range(len(values) - 1, 11, -1):
        leading = values[exponent] % MODULUS
        if not leading:
            continue
        shift = exponent - 12
        for index, coefficient in enumerate(PHI21[:-1]):
            values[shift + index] = (
                values[shift + index] - leading * coefficient
            ) % MODULUS
    return (values + [0] * 12)[:12]


def ring_multiply(left: list[int], right: list[int]) -> list[int]:
    product = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            product[i + j] += left_value * right_value
    return ring_reduce(product)


def ring_automorphism(value: list[int], exponent: int) -> list[int]:
    result = [0] * (11 * exponent + 1)
    for index, coefficient in enumerate(value):
        result[index * exponent] = coefficient
    return ring_reduce(result)


def decomposition_norm(value: list[int]) -> list[int]:
    result = [1] + [0] * 11
    for exponent in DECOMPOSITION_GROUP:
        result = ring_multiply(result, ring_automorphism(value, exponent))
    return result
