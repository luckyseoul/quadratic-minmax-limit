#!/usr/bin/env python3
"""Parallel direct-line scan of the W2 boundary cross-norm Gram criterion."""
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

from gf2x_ntl import cyclic_star_product_bits, gcd_bits  # noqa: E402
from io_atomic import write_json_atomic  # noqa: E402
from w2_boundary_joint_norm import cyclotomic_bits  # noqa: E402
from w2_boundary_line_direct import record  # noqa: E402
from w2_translated_antipodal_norm_scan import is_prime  # noqa: E402


def classify_row(row: dict) -> dict:
    factor_cache = {1: 0b11}
    gram_clearances: collections.Counter = collections.Counter()
    raw_clearances: collections.Counter = collections.Counter()
    exceptions = []
    raw_exceptions = []
    records = 0
    max_order = 0
    for clearance in row["clearances"]:
        if clearance["layer"] != "scalar":
            continue
        records += 1
        order = clearance["order"]
        max_order = max(max_order, order)
        raw_clearance = clearance["first_difference_a_clearing_order"]
        raw_clearances[raw_clearance] += 1
        if raw_clearance is None:
            raw_exceptions.append(
                {
                    "p": row["p"],
                    "order": order,
                    "raw_residual_hex": clearance["residual_hex"],
                    "raw_residual_degree": clearance["residual_degree"],
                }
            )

        factor = cyclotomic_bits(order, factor_cache)
        gram_gcd = factor
        gram_clearance = None
        vectors = []
        for trace in clearance["residual_trace"]:
            current = [int(value, 16) for value in trace["raw_component_hex"]]
            for previous in vectors:
                forward = 0
                backward = 0
                for component, previous_component in zip(current, previous):
                    forward ^= cyclic_star_product_bits(
                        component, previous_component, order
                    )
                    backward ^= cyclic_star_product_bits(
                        previous_component, component, order
                    )
                gram_gcd = gcd_bits(gram_gcd, forward)
                gram_gcd = gcd_bits(gram_gcd, backward)
            diagonal = 0
            for component in current:
                diagonal ^= cyclic_star_product_bits(
                    component, component, order
                )
            gram_gcd = gcd_bits(gram_gcd, diagonal)
            vectors.append(current)
            if gram_gcd == 1 and gram_clearance is None:
                gram_clearance = trace["a"]

        gram_clearances[gram_clearance] += 1
        if gram_clearance is None:
            exceptions.append(
                {
                    "p": row["p"],
                    "order": order,
                    "raw_clearance": raw_clearance,
                    "gram_gcd_hex": hex(gram_gcd),
                    "gram_gcd_degree": gram_gcd.bit_length() - 1,
                }
            )
    return {
        "p": row["p"],
        "records": records,
        "max_order": max_order,
        "raw_clearances": dict(raw_clearances),
        "gram_clearances": dict(gram_clearances),
        "raw_exceptions": raw_exceptions,
        "gram_exceptions": exceptions,
        "line_elapsed_seconds": row["elapsed_seconds"],
    }


def scan_prime(p: int, orders: list[int], prefix: int) -> dict:
    return classify_row(record(p, orders, prefix, False, True))


def encode_counter(counter: collections.Counter) -> dict[str, int]:
    return {
        str(key): counter[key]
        for key in sorted(counter, key=lambda value: (value is None, value or 0))
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=5)
    parser.add_argument("--stop", type=int, required=True)
    parser.add_argument("--max-order", type=int, default=4095)
    parser.add_argument("--prefix", type=int, default=3)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.prefix < 1 or args.workers < 1:
        raise ValueError("prefix and workers must be positive")
    if args.max_order < 3 or args.max_order % 2 == 0:
        raise ValueError("max-order must be odd and at least three")

    orders = list(range(3, args.max_order + 1, 2))
    primes = [
        p
        for p in range(max(5, args.start), args.stop + 1)
        if p % 12 == 5 and is_prime(p)
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
    gram_clearances: collections.Counter = collections.Counter()
    for row in rows:
        raw_clearances.update(row["raw_clearances"])
        gram_clearances.update(row["gram_clearances"])
    result = {
        "range": [args.start, args.stop],
        "congruence_class": "p == 5 (mod 12)",
        "prefix": args.prefix,
        "max_requested_order": args.max_order,
        "algorithm": "parallel direct selected lines plus F2 cross-norm Gram",
        "n_primes": len(rows),
        "records": sum(row["records"] for row in rows),
        "max_order_tested": max((row["max_order"] for row in rows), default=0),
        "raw_clearance_counts": encode_counter(raw_clearances),
        "gram_clearance_counts": encode_counter(gram_clearances),
        "raw_exceptions": [
            item for row in rows for item in row["raw_exceptions"]
        ],
        "gram_exceptions": [
            item for row in rows for item in row["gram_exceptions"]
        ],
        "line_cpu_seconds": sum(
            row["line_elapsed_seconds"] for row in rows
        ),
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(args.output, result)
    print(
        f"done primes={result['n_primes']} records={result['records']} "
        f"raw_exceptions={len(result['raw_exceptions'])} "
        f"gram_exceptions={len(result['gram_exceptions'])} "
        f"seconds={result['elapsed_seconds']:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
