#!/usr/bin/env python3
"""Exact O(p log p) selected-line attack on the W2 scalar obstruction.

The orbit scanners visit all ``p^2-1`` finite nonzero points.  For a fixed
prefix only four affine lines per difference matter.  This program enumerates
those lines directly and recovers exponent classes with small character
tables in the native NTL bridge.

It also checks a useful repeated-root reformulation.  If

    D(Z) = D_sq(Z^2) + Z D_ns(Z^2),

then in characteristic two ``D'(Z)=D_ns(Z^2)``.  Hence an odd irreducible
factor divides both raw components exactly when it is a repeated factor of
the unsplit exponent polynomial D.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import _finv  # noqa: E402
from gf2x_ntl import available as ntl_available  # noqa: E402
from gf2x_ntl import selected_line_bins  # noqa: E402
from io_atomic import write_json_atomic  # noqa: E402
from w2_low_order_atomic_gpu import (  # noqa: E402
    cyclotomic_bits,
    oddpart,
    polynomial_bits,
    reciprocal_bits,
)
from w2_translated_antipodal_norm_scan import f2_gcd_bits, is_prime  # noqa: E402


def row_levels(p: int, a: int) -> list[int]:
    half = (p - 1) // 2
    enter = (half - a) % p
    leave = (p - 1 - a) % p
    return [enter, leave, (-enter) % p, (-leave) % p]


def spread_even_bits(value: int) -> int:
    """Substitute X^2 into a packed binary polynomial."""
    result = 0
    while value:
        bit = value & -value
        result |= 1 << (2 * (bit.bit_length() - 1))
        value ^= bit
    return result


def record(
    p: int,
    requested_orders: list[int],
    prefix: int,
    include_projective: bool,
    trace_full_prefix: bool = False,
) -> dict:
    started = time.perf_counter()
    q, mul, _add, chi, _frob, _norm, ia, ib = field_ctx(p)
    sigma = next(value for value in range(1, q) if chi(value) == -1)
    sigma_inverse = _finv(mul, q, sigma)
    omega = _primitive(mul, q)
    generator = mul(omega, omega)
    projective_order = (p + 1) // 2
    scalar_order = oddpart(p - 1)
    ambient_order = projective_order * scalar_order
    orders = [
        order
        for order in requested_orders
        if ambient_order % order == 0
        and (include_projective or projective_order % order != 0)
    ]
    if not orders:
        return {
            "p": p,
            "ambient_oddpart": ambient_order,
            "projective_order": projective_order,
            "scalar_order": scalar_order,
            "orders_tested": [],
            "clearances": [],
            "all_orders_cleared_by_differences": True,
            "elapsed_seconds": time.perf_counter() - started,
        }

    rows = [row_levels(p, a) for a in range(1, prefix + 1)]
    levels = sorted(set(level for row in rows for level in row))
    level_index = {level: index for index, level in enumerate(levels)}
    half = (p - 1) // 2
    pole_t = (sigma_inverse // p) * pow(half, p - 2, p) % p
    line_bins, offsets = selected_line_bins(
        p,
        ia,
        ib,
        generator,
        omega,
        sigma,
        pole_t,
        levels,
        orders,
    )

    factors = {order: cyclotomic_bits(order) for order in orders}
    common = dict(factors)
    clearance = {order: None for order in orders}
    traces = {order: [] for order in orders}
    for a, row in enumerate(rows, 1):
        for order, offset in zip(orders, offsets):
            if common[order] == 1 and not trace_full_prefix:
                continue
            was_uncleared = common[order] != 1
            component_values = []
            component_polynomials = []
            for component in range(2):
                delta = np.zeros(order, dtype=np.uint8)
                for level in row:
                    delta ^= line_bins[
                        component, level_index[level], offset : offset + order
                    ]
                value = polynomial_bits(delta)
                component_polynomials.append(value)
                component_values.append(hex(value))

            square, nonsquare = component_polynomials
            component_common = f2_gcd_bits(factors[order], square)
            component_common = f2_gcd_bits(component_common, nonsquare)
            full = spread_even_bits(square) | (spread_even_bits(nonsquare) << 1)
            derivative = spread_even_bits(nonsquare)
            repeated_common = f2_gcd_bits(factors[order], full)
            repeated_common = f2_gcd_bits(repeated_common, derivative)
            if component_common != repeated_common:
                raise AssertionError(
                    f"p={p} a={a} d={order}: component/repeated-root mismatch"
                )

            for value in component_polynomials:
                common[order] = f2_gcd_bits(common[order], value)
                common[order] = f2_gcd_bits(
                    common[order],
                    reciprocal_bits(
                        np.array(
                            [(value >> j) & 1 for j in range(order)],
                            dtype=np.uint8,
                        )
                    ),
                )
            traces[order].append(
                {
                    "a": a,
                    "raw_component_hex": component_values,
                    "repeated_root_hex": hex(repeated_common),
                    "running_aut_hex": hex(common[order]),
                }
            )
            if was_uncleared and common[order] == 1:
                clearance[order] = a
        if (
            not trace_full_prefix
            and all(value == 1 for value in common.values())
        ):
            break

    clearances = [
        {
            "order": order,
            "layer": (
                "projective" if projective_order % order == 0 else "scalar"
            ),
            "first_difference_a_clearing_order": clearance[order],
            "residual_degree": common[order].bit_length() - 1,
            "residual_hex": hex(common[order]),
            "residual_trace": traces[order],
        }
        for order in orders
    ]
    return {
        "p": p,
        "ambient_oddpart": ambient_order,
        "projective_order": projective_order,
        "scalar_order": scalar_order,
        "orders_tested": orders,
        "prefix": prefix,
        "trace_full_prefix": trace_full_prefix,
        "n_selected_levels": len(levels),
        "pole_t": pole_t,
        "clearances": clearances,
        "all_orders_cleared_by_differences": all(
            value == 1 for value in common.values()
        ),
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=5)
    parser.add_argument("--stop", type=int, required=True)
    parser.add_argument("--orders")
    parser.add_argument("--max-order", type=int, default=255)
    parser.add_argument("--prefix", type=int, default=3)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--include-projective", action="store_true")
    parser.add_argument("--trace-full-prefix", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.prefix < 1 or args.workers < 1:
        raise ValueError("prefix and workers must be positive")
    requested_orders = (
        sorted(set(map(int, args.orders.split(","))))
        if args.orders
        else list(range(3, args.max_order + 1, 2))
    )
    if any(order < 3 or order % 2 == 0 for order in requested_orders):
        raise ValueError("orders must be odd and at least three")
    if not ntl_available():
        raise RuntimeError("the direct selected-line scan requires the NTL bridge")

    primes = [
        p
        for p in range(max(5, args.start), args.stop + 1)
        if p % 12 == 5 and is_prime(p)
    ]
    started = time.perf_counter()
    rows = []
    if args.workers == 1:
        for index, p in enumerate(primes, 1):
            row = record(
                p,
                requested_orders,
                args.prefix,
                args.include_projective,
                args.trace_full_prefix,
            )
            rows.append(row)
            print(
                f"[{index}/{len(primes)}] p={p} "
                f"cleared={row['all_orders_cleared_by_differences']} "
                f"seconds={row['elapsed_seconds']:.3f}",
                flush=True,
            )
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {
                pool.submit(
                    record,
                    p,
                    requested_orders,
                    args.prefix,
                    args.include_projective,
                    args.trace_full_prefix,
                ): p
                for p in primes
            }
            for index, future in enumerate(as_completed(futures), 1):
                row = future.result()
                rows.append(row)
                print(
                    f"[{index}/{len(primes)}] p={row['p']} "
                    f"cleared={row['all_orders_cleared_by_differences']} "
                    f"seconds={row['elapsed_seconds']:.3f}",
                    flush=True,
                )
    rows.sort(key=lambda row: row["p"])
    result = {
        "range": [args.start, args.stop],
        "congruence_class": "p == 5 (mod 12)",
        "requested_orders": requested_orders,
        "prefix": args.prefix,
        "include_projective": args.include_projective,
        "trace_full_prefix": args.trace_full_prefix,
        "algorithm": "direct selected affine lines; O(p log p) per order",
        "n_primes": len(rows),
        "uncleared": [
            {"p": row["p"], **item}
            for row in rows
            for item in row["clearances"]
            if item["residual_hex"] != "0x1"
        ],
        "rows": rows,
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(args.output, result)
    print(
        f"done primes={len(rows)} uncleared={len(result['uncleared'])}",
        flush=True,
    )


if __name__ == "__main__":
    main()
