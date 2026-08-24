#!/usr/bin/env python3
"""Parallel exact census of oriented W2 boundary common divisors.

For the first boundary differences let G be the gcd of Phi_d with every raw
square/nonsquare component.  The Bose generator pair is unimodular on this
layer, so G is the one-sided common divisor of the underlying line words.
The actual Aut obstruction is exactly the reversible core gcd(G, G*).
"""
from __future__ import annotations

import argparse
import collections
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from gf2x_ntl import gcd_bits  # noqa: E402
from io_atomic import write_json_atomic  # noqa: E402
from w2_boundary_joint_norm import cyclotomic_bits  # noqa: E402
from w2_boundary_line_direct import record  # noqa: E402
from w2_translated_antipodal_norm_scan import is_prime  # noqa: E402


def cyclic_reciprocal_bits(value: int, order: int) -> int:
    """Apply X -> X^-1 modulo X^order+1 to packed coefficients."""
    result = 0
    while value:
        bit = value & -value
        exponent = bit.bit_length() - 1
        result |= 1 << (0 if exponent == 0 else order - exponent)
        value ^= bit
    return result


def ordinary_reciprocal_bits(value: int) -> int:
    """Reverse coefficients through the polynomial's actual degree."""
    degree = value.bit_length() - 1
    return sum(
        ((value >> exponent) & 1) << (degree - exponent)
        for exponent in range(degree + 1)
    )


def classify_row(row: dict) -> dict:
    factor_cache = {1: 0b11}
    raw_clearances: collections.Counter = collections.Counter()
    oriented_degree_counts: collections.Counter = collections.Counter()
    oriented_nonunits = []
    reversible_exceptions = []
    records = 0
    max_order = 0
    for clearance in row["clearances"]:
        if clearance["layer"] != "scalar":
            continue
        records += 1
        order = clearance["order"]
        max_order = max(max_order, order)
        factor = cyclotomic_bits(order, factor_cache)
        oriented = factor
        first_reversible_clearance = None
        reversible_core = factor
        for trace in clearance["residual_trace"]:
            for component_hex in trace["raw_component_hex"]:
                oriented = gcd_bits(oriented, int(component_hex, 16))
            reversible_core = gcd_bits(
                oriented, cyclic_reciprocal_bits(oriented, order)
            )
            if reversible_core == 1 and first_reversible_clearance is None:
                first_reversible_clearance = trace["a"]

        expected_clearance = clearance["first_difference_a_clearing_order"]
        expected_residual = int(clearance["residual_hex"], 16)
        if first_reversible_clearance != expected_clearance:
            raise AssertionError(
                f"p={row['p']} d={order}: oriented/raw clearance mismatch"
            )
        if reversible_core != expected_residual:
            raise AssertionError(
                f"p={row['p']} d={order}: oriented/raw residual mismatch"
            )
        raw_clearances[first_reversible_clearance] += 1

        degree = oriented.bit_length() - 1
        oriented_degree_counts[degree] += 1
        if oriented != 1:
            oriented_nonunits.append(
                {
                    "p": row["p"],
                    "order": order,
                    "oriented_gcd_hex": hex(oriented),
                    "oriented_gcd_degree": degree,
                    "reciprocal_hex": hex(ordinary_reciprocal_bits(oriented)),
                    "reversible_core_hex": hex(reversible_core),
                    "reversible_core_degree": (
                        reversible_core.bit_length() - 1
                    ),
                }
            )
        if reversible_core != 1:
            reversible_exceptions.append(oriented_nonunits[-1])

    return {
        "p": row["p"],
        "records": records,
        "max_order": max_order,
        "raw_clearances": dict(raw_clearances),
        "oriented_degree_counts": dict(oriented_degree_counts),
        "oriented_nonunits": oriented_nonunits,
        "reversible_exceptions": reversible_exceptions,
        "line_elapsed_seconds": row["elapsed_seconds"],
    }


def scan_prime(p: int, orders: list[int], prefix: int) -> dict:
    return classify_row(record(p, orders, prefix, False, True))


def encode_counter(counter: collections.Counter) -> dict[str, int]:
    return {
        str(key): counter[key]
        for key in sorted(counter, key=lambda value: (value is None, value or 0))
    }


def oddpart(value: int) -> int:
    while value % 2 == 0:
        value //= 2
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=5)
    parser.add_argument("--stop", type=int, required=True)
    parser.add_argument("--orders")
    parser.add_argument("--max-order", type=int, default=4095)
    parser.add_argument("--prefix", type=int, default=3)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.prefix < 1 or args.workers < 1:
        raise ValueError("prefix and workers must be positive")
    if not args.orders and (args.max_order < 3 or args.max_order % 2 == 0):
        raise ValueError("max-order must be odd and at least three")

    orders = (
        sorted(set(map(int, args.orders.split(","))))
        if args.orders
        else list(range(3, args.max_order + 1, 2))
    )
    if any(order < 3 or order % 2 == 0 for order in orders):
        raise ValueError("orders must be odd and at least three")
    primes = [
        p
        for p in range(max(5, args.start), args.stop + 1)
        if p % 12 == 5 and is_prime(p)
        and (
            not args.orders
            or any(
                ((p + 1) // 2 * oddpart(p - 1)) % order == 0
                and ((p + 1) // 2) % order != 0
                for order in orders
            )
        )
    ]
    started = time.perf_counter()
    rows = []
    if args.workers == 1:
        for index, p in enumerate(primes, 1):
            rows.append(scan_prime(p, orders, args.prefix))
            if index % 25 == 0 or index == len(primes):
                print(f"[{index}/{len(primes)}] p={p}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(scan_prime, p, orders, args.prefix): p
                for p in primes
            }
            for index, future in enumerate(as_completed(futures), 1):
                rows.append(future.result())
                if index % 25 == 0 or index == len(primes):
                    print(f"[{index}/{len(primes)}]", flush=True)

    rows.sort(key=lambda row: row["p"])
    raw_clearances: collections.Counter = collections.Counter()
    oriented_degrees: collections.Counter = collections.Counter()
    for row in rows:
        raw_clearances.update(row["raw_clearances"])
        oriented_degrees.update(row["oriented_degree_counts"])
    result = {
        "range": [args.start, args.stop],
        "congruence_class": "p == 5 (mod 12)",
        "prefix": args.prefix,
        "requested_orders": orders,
        "max_requested_order": max(orders),
        "algorithm": "parallel direct lines; exact oriented gcd and reversible core",
        "n_primes": len(rows),
        "records": sum(row["records"] for row in rows),
        "max_order_tested": max((row["max_order"] for row in rows), default=0),
        "raw_clearance_counts": encode_counter(raw_clearances),
        "oriented_gcd_degree_counts": encode_counter(oriented_degrees),
        "oriented_nonunits": [
            item for row in rows for item in row["oriented_nonunits"]
        ],
        "reversible_exceptions": [
            item for row in rows for item in row["reversible_exceptions"]
        ],
        "line_cpu_seconds": sum(
            row["line_elapsed_seconds"] for row in rows
        ),
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(args.output, result)
    print(
        f"done primes={result['n_primes']} records={result['records']} "
        f"oriented_nonunits={len(result['oriented_nonunits'])} "
        f"reversible_exceptions={len(result['reversible_exceptions'])} "
        f"seconds={result['elapsed_seconds']:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
