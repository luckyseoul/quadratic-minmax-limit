#!/usr/bin/env python3
"""Probe local mod-four norms of the six-line d=21 boundary obstruction.

The binary square/nonsquare component words are reductions of honest
order-21 cyclotomic sums.  This script retains their exponent-class counts
modulo four and takes the decomposition-group norm for <2> in Q(zeta_21).
It is reconnaissance for a local obstruction analogous to the Jacobi norm
that proves the Bose-generator phase; no observed norm pattern is promoted
to a theorem here.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from gf2x_ntl import (  # noqa: E402
    field_primitive,
    selected_line_bins,
    selected_line_counts_mod4,
)
from io_atomic import write_json_atomic  # noqa: E402
from w2_d21_mod4_ring import (  # noqa: E402
    DECOMPOSITION_GROUP,
    decomposition_norm,
    ring_multiply,
    ring_reduce,
)


ORDER = 21
MODULUS = 4


def is_prime(value: int) -> bool:
    if value < 2 or value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def field_power(multiply, value: int, exponent: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply(result, value)
        value = multiply(value, value)
        exponent >>= 1
    return result


def row_levels(p: int, a: int) -> list[int]:
    half = (p - 1) // 2
    enter = (half - a) % p
    leave = (p - 1 - a) % p
    return [enter, leave, (-enter) % p, (-leave) % p]


def component_value(counts: np.ndarray) -> list[int]:
    return ring_reduce([int(value) for value in counts])


def shift(value: list[int], exponent: int) -> list[int]:
    return ring_reduce([0] * exponent + value)


def add(left: list[int], right: list[int], sign: int = 1) -> list[int]:
    return [
        (left[index] + sign * right[index]) % MODULUS
        for index in range(12)
    ]


def record(p: int) -> dict:
    if p % 84 != 29 or not is_prime(p):
        raise ValueError(f"p={p} must be prime and 29 mod 84")
    q, mul, _add, _chi, _frob, _norm, ia, ib = field_ctx(p)
    sigma = p
    sigma_inverse = field_power(mul, sigma, q - 2)
    omega = field_primitive(p, ia, ib)
    generator = mul(omega, omega)
    half = (p - 1) // 2
    pole_t = (sigma_inverse // p) * pow(half, p - 2, p) % p
    rows = [row_levels(p, a) for a in range(1, 4)]
    levels = sorted(set(level for row in rows for level in row))
    level_index = {level: index for index, level in enumerate(levels)}
    parity, _ = selected_line_bins(
        p, ia, ib, generator, omega, sigma, pole_t, levels, [ORDER],
        force_wide=True,
    )
    counts, _ = selected_line_counts_mod4(
        p, ia, ib, generator, omega, sigma, pole_t, levels, [ORDER]
    )
    if not np.array_equal(parity, counts & 1):
        raise AssertionError(f"p={p}: mod-four counts do not reduce to parity")

    differences = []
    for a, row in enumerate(rows, 1):
        components = []
        for component in range(2):
            total = np.zeros(ORDER, dtype=np.uint8)
            for level in row:
                total = (total + counts[component, level_index[level]]) & 3
            components.append(component_value(total))
        square, nonsquare = components
        shifted_nonsquare = shift(nonsquare, 11)
        order21_sum = add(square, shifted_nonsquare)
        order42_sum = add(square, shifted_nonsquare, -1)
        differences.append(
            {
                "a": a,
                "square_mod4": square,
                "nonsquare_mod4": nonsquare,
                "square_norm_mod4": decomposition_norm(square),
                "nonsquare_norm_mod4": decomposition_norm(nonsquare),
                "order21_sum_mod4": order21_sum,
                "order42_sum_mod4": order42_sum,
                "order21_norm_mod4": decomposition_norm(order21_sum),
                "order42_norm_mod4": decomposition_norm(order42_sum),
            }
        )
    return {
        "p": p,
        "field_polynomial": [ia, ib],
        "primitive_element": omega,
        "pole_t": pole_t,
        "differences": differences,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primes = list(dict.fromkeys(map(int, args.primes.split(","))))
    rows = []
    for p in primes:
        row = record(p)
        rows.append(row)
        compact = [
            {
                "a": item["a"],
                "N21": item["order21_norm_mod4"],
                "N42": item["order42_norm_mod4"],
            }
            for item in row["differences"]
        ]
        print(f"p={p} {json.dumps(compact, separators=(',', ':'))}", flush=True)
    write_json_atomic(
        args.output,
        {
            "description": "d=21 boundary component and character norms mod four",
            "decomposition_group": list(DECOMPOSITION_GROUP),
            "warning": "reconnaissance only; no norm pattern is a proof",
            "rows": rows,
        },
    )


if __name__ == "__main__":
    main()
