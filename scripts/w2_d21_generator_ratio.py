#!/usr/bin/env python3
"""Direct order-21 character table of the two-line Bose generator.

The generator is the parity multiset of ``t*sigma`` and ``1+t*sigma`` for
``t in F_p``.  Only its exponent classes modulo 21 are needed here, so this
avoids the full ``p^2`` inverse-log table and evaluates the order character
directly on the two lines.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

ROOT = Path(os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15613 import _finv  # noqa: E402
from gf2x_ntl import field_primitive  # noqa: E402
from io_atomic import write_json_atomic  # noqa: E402


ORDER = 21
FACTORS = (0x57, 0x75)


def is_prime(value: int) -> bool:
    if value < 2 or value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def remainder_bits(value: int, modulus: int) -> int:
    modulus_degree = modulus.bit_length() - 1
    while value and value.bit_length() - 1 >= modulus_degree:
        value ^= modulus << (value.bit_length() - 1 - modulus_degree)
    return value


def multiply_bits(left: int, right: int, modulus: int) -> int:
    result = 0
    modulus_degree = modulus.bit_length() - 1
    while right:
        if right & 1:
            result ^= left
        right >>= 1
        left <<= 1
        if left >> modulus_degree:
            left ^= modulus
    return result


def power_bits(value: int, exponent: int, modulus: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply_bits(result, value, modulus)
        value = multiply_bits(value, value, modulus)
        exponent >>= 1
    return result


def field_power(mul, value: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = mul(result, value)
        value = mul(value, value)
        exponent >>= 1
    return result


def generator_record(p: int) -> dict:
    if p % 84 != 29 or not is_prime(p):
        raise ValueError(f"p={p} must be prime and 29 mod 84")
    started = time.perf_counter()
    q, mul, add, _chi, _frob, _norm, ia, ib = field_ctx(p)
    sigma = p
    omega = field_primitive(p, ia, ib)
    omega_inverse = _finv(mul, q, omega)
    generator = mul(omega, omega)
    orbit_length = (q - 1) // 2
    root = field_power(mul, generator, orbit_length // ORDER)
    root_log = {}
    value = 1
    for residue in range(ORDER):
        root_log[value] = residue
        value = mul(value, root)
    if value != 1 or len(root_log) != ORDER:
        raise AssertionError(f"p={p}: invalid order-21 root table")

    component_bits = [0, 0]

    def include(point: int) -> None:
        if point == 0:
            return
        a, b = point % p, point // p
        norm = (a * a + ia * a * b - ib * b * b) % p
        square = pow(norm, (p - 1) // 2, p) == 1
        base = point if square else mul(omega_inverse, point)
        character = field_power(mul, base, orbit_length // ORDER)
        component = 0 if square else 1
        component_bits[component] ^= 1 << root_log[character]

    for t in range(p):
        line_point = mul(t, sigma)
        include(line_point)
        include(add(1, line_point))

    factor_rows = []
    for factor in FACTORS:
        square = remainder_bits(component_bits[0], factor)
        nonsquare = remainder_bits(component_bits[1], factor)
        normalized_nonsquare = multiply_bits(
            nonsquare, power_bits(0x2, 62, factor), factor
        )
        ratio = None
        if nonsquare:
            ratio = multiply_bits(
                multiply_bits(0x2, square, factor),
                power_bits(nonsquare, 62, factor),
                factor,
            )
        factor_rows.append(
            {
                "factor_hex": hex(factor),
                "square_remainder_hex": hex(square),
                "nonsquare_remainder_hex": hex(nonsquare),
                "nonsquare_over_x_hex": hex(normalized_nonsquare),
                "square_in_f8": power_bits(square, 8, factor) == square,
                "nonsquare_over_x_in_f8": (
                    power_bits(normalized_nonsquare, 8, factor)
                    == normalized_nonsquare
                ),
                "ratio_hex": None if ratio is None else hex(ratio),
            }
        )

    return {
        "p": p,
        "field_polynomial": [ia, ib],
        "sigma": sigma,
        "primitive_element": omega,
        "component_polynomial_hex": [hex(value) for value in component_bits],
        "component_support": [
            [j for j in range(ORDER) if value >> j & 1]
            for value in component_bits
        ],
        "component_support_mod_3_7": [
            [[j % 3, j % 7] for j in range(ORDER) if value >> j & 1]
            for value in component_bits
        ],
        "factors": factor_rows,
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primes = list(dict.fromkeys(map(int, args.primes.split(","))))
    rows = []
    for p in primes:
        row = generator_record(p)
        rows.append(row)
        print(
            f"p={p} ratios="
            f"{[item['ratio_hex'] for item in row['factors']]} "
            f"seconds={row['elapsed_seconds']:.3f}",
            flush=True,
        )
    write_json_atomic(
        args.output,
        {
            "order": ORDER,
            "factor_hex": [hex(factor) for factor in FACTORS],
            "algorithm": "direct two-line order-character evaluation",
            "rows": rows,
        },
    )


if __name__ == "__main__":
    main()
