#!/usr/bin/env python3
"""Exact phase certificate for the d=21 Bose generator.

For p == 29 (mod 84), let psi be the order-14 character of F_p^* induced
by a primitive element of F_(p^2), and put c=4/27.  The coefficients

    b_e = #{x != 0,1 : ind(c (x-1)^3/x) == e (mod 14)}

encode the parity of the affine generator line.  Their odd coefficients are
the seven nonsquare support counts q_a.  If H_7 and H_14 are the associated
order-7 and order-14 Jacobi sums, then

    sigma_4(H_7) - H_14 = 2 zeta^4 O(zeta),

where O contains the odd b_e.  This script verifies the identity, constructs
U=H_14/sigma_4(H_7) modulo four, and certifies that its local cubic norm is
-1.  Thus Tr_F8/F2((U-1)/2)=1, the phase bit missing from the RDS lift.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

ROOT = Path(os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from gf2x_ntl import field_primitive  # noqa: E402
from io_atomic import write_json_atomic  # noqa: E402
from w2_d21_generator_ratio import (  # noqa: E402
    field_power,
    generator_record,
    is_prime,
)


def reduce_phi7(coefficients: list[int], modulus: int) -> list[int]:
    values = [value % modulus for value in coefficients]
    for exponent in range(len(values) - 1, 5, -1):
        leading = values[exponent] % modulus
        if leading:
            for offset in range(6):
                values[exponent - 6 + offset] = (
                    values[exponent - 6 + offset] - leading
                ) % modulus
    return (values + [0] * 6)[:6]


def ring_multiply(left: list[int], right: list[int], modulus: int) -> list[int]:
    product = [0] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            product[i + j] += left_value * right_value
    return reduce_phi7(product, modulus)


def ring_frobenius(value: list[int], exponent: int, modulus: int) -> list[int]:
    result = [0] * (5 * exponent + 1)
    for index, coefficient in enumerate(value):
        result[index * exponent] = coefficient
    return reduce_phi7(result, modulus)


def ring_inverse_mod_two(value: list[int]) -> list[int]:
    for mask in range(64):
        candidate = [(mask >> index) & 1 for index in range(6)]
        if ring_multiply(value, candidate, 2) == [1, 0, 0, 0, 0, 0]:
            return candidate
    raise ArithmeticError("nonunit modulo Phi_7 and two")


def ring_inverse_mod_four(value: list[int]) -> list[int]:
    inverse = ring_inverse_mod_two([coefficient & 1 for coefficient in value])
    product = ring_multiply(value, inverse, 4)
    correction = [(-coefficient) % 4 for coefficient in product]
    correction[0] = (2 - product[0]) % 4
    return ring_multiply(inverse, correction, 4)


def gf8_multiply(left: int, right: int, polynomial: int) -> int:
    result = 0
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left >> 3:
            left ^= polynomial
    return result


def gf8_evaluate(coefficients: list[int], polynomial: int) -> int:
    result = 0
    term = 1
    for coefficient in coefficients:
        if coefficient & 1:
            result ^= term
        term = gf8_multiply(term, 2, polynomial)
    return result


def gf8_trace(value: int, polynomial: int) -> int:
    square = gf8_multiply(value, value, polynomial)
    fourth = gf8_multiply(square, square, polynomial)
    return value ^ square ^ fourth


def crt_class(mod14: int, mod6: int) -> int:
    return next(
        value
        for value in range(42)
        if value % 14 == mod14 and value % 6 == mod6
    )


def record(p: int) -> dict:
    if p % 84 != 29 or not is_prime(p):
        raise ValueError(f"p={p} must be prime and 29 mod 84")
    generator = generator_record(p)
    q, mul, _add, _chi, _frob, _norm, ia, ib = field_ctx(p)
    omega = field_primitive(p, ia, ib)
    norm_omega = mul(omega, field_power(mul, omega, p)) % p
    scalar_log = {}
    value = 1
    for exponent in range(p - 1):
        scalar_log[value] = exponent
        value = value * norm_omega % p
    if value != 1 or len(scalar_log) != p - 1:
        raise AssertionError(f"p={p}: norm(omega) is not primitive")

    constant = 4 * pow(27, p - 2, p) % p
    coefficients = [0] * 14
    for x in range(2, p):
        rational_value = (
            constant * pow(x - 1, 3, p) * pow(x, p - 2, p)
        ) % p
        coefficients[scalar_log[rational_value] % 14] += 1

    line_counts = generator["line_component_class_counts"][1]
    line_word = [0] * 42
    for component in (0, 1):
        for order21_class, count in enumerate(line_counts[component]):
            line_word[2 * order21_class + component] = count & 1
    bridge_checks = []
    for exponent in range(14):
        if exponent & 1:
            expected = line_word[crt_class(exponent, 1)]
        else:
            expected = line_word[crt_class(exponent, 2)] ^ int(exponent == 0)
        bridge_checks.append((coefficients[exponent] & 1) == expected)

    order7_coefficients = [
        coefficients[index] + coefficients[index + 7] for index in range(7)
    ]
    order14_differences = [
        coefficients[index] - coefficients[index + 7] for index in range(7)
    ]
    # xi=-zeta^4.  Convert sum D_j xi^j to the zeta basis.
    order14_zeta = [0] * 7
    for exponent, coefficient in enumerate(order14_differences):
        order14_zeta[4 * exponent % 7] = (-1) ** exponent * coefficient

    h7 = reduce_phi7(order7_coefficients, 4)
    h14 = reduce_phi7(order14_zeta, 4)
    sigma4_h7 = ring_frobenius(h7, 4, 4)
    quotient = ring_multiply(h14, ring_inverse_mod_four(sigma4_h7), 4)
    first_digit_numerators = quotient[:]
    first_digit_numerators[0] = (first_digit_numerators[0] - 1) % 4
    if any(value & 1 for value in first_digit_numerators):
        raise AssertionError(f"p={p}: Jacobi quotient is not one modulo two")
    first_digit = [value // 2 for value in first_digit_numerators]
    local_norm = ring_multiply(
        ring_multiply(quotient, ring_frobenius(quotient, 2, 4), 4),
        ring_frobenius(quotient, 4, 4),
        4,
    )
    local_traces = {
        hex(polynomial): gf8_trace(
            gf8_evaluate(first_digit, polynomial), polynomial
        )
        for polynomial in (0xB, 0xD)
    }

    odd_weight = sum(coefficients[index] & 1 for index in range(1, 14, 2))
    collapsed_parity = [value & 1 for value in order7_coefficients]
    l_size = (p - 1) // 14
    square_norm = sum(value * value for value in order7_coefficients)
    difference_norm = sum(value * value for value in order14_differences)
    signed_cross = sum(
        (-1) ** index
        * order7_coefficients[index]
        * order14_differences[index]
        for index in range(7)
    )
    return {
        "p": p,
        "field_polynomial": [ia, ib],
        "primitive_element": omega,
        "norm_primitive": norm_omega,
        "constant_4_over_27": constant,
        "coefficient_counts_order_14": coefficients,
        "coefficient_parities_order_14": [value & 1 for value in coefficients],
        "order7_coefficient_counts": order7_coefficients,
        "order14_difference_counts": order14_differences,
        "order14_zeta_coefficients": order14_zeta,
        "odd_coefficient_weight": odd_weight,
        "odd_weight_one_mod_four": odd_weight % 4 == 1,
        "line_bridge_all_classes": all(bridge_checks),
        "line_bridge_checks": bridge_checks,
        "order7_parity_is_monomial": sum(collapsed_parity) == 1,
        "order7_parity_support": collapsed_parity.index(1),
        "jacobi_quotient_mod_4": quotient,
        "jacobi_first_digit_mod_2": first_digit,
        "local_norm_mod_4": local_norm,
        "local_trace_by_factor": local_traces,
        "local_trace_one": set(local_traces.values()) == {1},
        "norm_checks": {
            "sum_C_squared": square_norm,
            "sum_C_squared_expected": (p * p + 2 * p + 4) // 7,
            "sum_D_squared": difference_norm,
            "sum_D_squared_expected": (6 * p + 1) // 7,
            "signed_cross": signed_cross,
            "signed_cross_mod_8_expected": (4 * ((p - 29) // 84) + 3) % 8,
            "all_hold": (
                square_norm == (p * p + 2 * p + 4) // 7
                and difference_norm == (6 * p + 1) // 7
                and signed_cross % 8
                == (4 * ((p - 29) // 84) + 3) % 8
            ),
        },
        "subgroup_class_size": l_size,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primes = list(dict.fromkeys(map(int, args.primes.split(","))))
    rows = [record(p) for p in primes]
    for row in rows:
        print(
            f"p={row['p']} weight={row['odd_coefficient_weight']} "
            f"norm4={row['local_norm_mod_4']} "
            f"trace={row['local_trace_by_factor']}",
            flush=True,
        )
    write_json_atomic(
        args.output,
        {
            "description": "order-14 Jacobi phase and local 2-adic norm",
            "proof_target": "d=21 generator state trace equals one",
            "rows": rows,
            "failures": [
                row
                for row in rows
                if not (
                    row["odd_weight_one_mod_four"]
                    and row["line_bridge_all_classes"]
                    and row["order7_parity_is_monomial"]
                    and row["local_norm_mod_4"] == [3, 0, 0, 0, 0, 0]
                    and row["local_trace_one"]
                    and row["norm_checks"]["all_hold"]
                )
            ],
        },
    )


if __name__ == "__main__":
    main()
