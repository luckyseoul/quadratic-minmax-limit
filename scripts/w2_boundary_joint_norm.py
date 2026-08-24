#!/usr/bin/env python3
"""Test a joint Hermitian-norm certificate on exact W2 line traces.

For packed line-difference polynomials P_{a,c} in F2[X]/(X^d-1), form

    N_A(X) = sum_{a <= A, c} P_{a,c}(X) P_{a,c}(X^-1).

If gcd(N_A, Phi_d) is one, the component polynomials cannot vanish together
at any primitive d-th root.  Thus a unit joint norm is a sufficient scalar
clearance certificate, potentially replacing six separate nonvanishings by
one intersection-parity identity.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT / "src"))

from gf2x_ntl import (  # noqa: E402
    available as ntl_available,
    cyclic_star_product_bits,
    gcd_bits,
)
from io_atomic import write_json_atomic  # noqa: E402


def exact_divide_bits(numerator: int, denominator: int) -> int:
    """Exact polynomial division in F2[X], using packed coefficients."""
    quotient = 0
    denominator_degree = denominator.bit_length() - 1
    while numerator and numerator.bit_length() - 1 >= denominator_degree:
        shift = numerator.bit_length() - 1 - denominator_degree
        quotient |= 1 << shift
        numerator ^= denominator << shift
    if numerator:
        raise ArithmeticError("nonzero remainder in cyclotomic factorization")
    return quotient


def divisors(value: int) -> list[int]:
    lower = []
    upper = []
    divisor = 1
    while divisor * divisor <= value:
        if value % divisor == 0:
            lower.append(divisor)
            if divisor * divisor != value:
                upper.append(value // divisor)
        divisor += 1
    return lower + list(reversed(upper))


def cyclotomic_bits(order: int, cache: dict[int, int]) -> int:
    """Return Phi_order modulo two from X^order+1 factorization."""
    if order not in cache:
        result = (1 << order) | 1
        for divisor in divisors(order)[:-1]:
            result = exact_divide_bits(
                result, cyclotomic_bits(divisor, cache)
            )
        cache[order] = result
    return cache[order]


def sorted_counter(counter: collections.Counter) -> dict[str, int]:
    return {
        str(key): counter[key]
        for key in sorted(counter, key=lambda value: (value is None, value or 0))
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not ntl_available():
        raise RuntimeError("the joint-norm test requires the NTL bridge")

    started = time.perf_counter()
    input_bytes = args.input.read_bytes()
    data = json.loads(input_bytes)
    factors: dict[int, int] = {1: 0b11}
    raw_clearances: collections.Counter = collections.Counter()
    norm_clearances: collections.Counter = collections.Counter()
    common_norm_clearances: collections.Counter = collections.Counter()
    gram_clearances: collections.Counter = collections.Counter()
    transitions: collections.Counter = collections.Counter()
    exceptions = []
    common_norm_exceptions = []
    gram_exceptions = []
    records = 0
    tested_orders = set()

    for row in data["rows"]:
        for clearance in row["clearances"]:
            if clearance["layer"] != "scalar":
                continue
            records += 1
            order = clearance["order"]
            tested_orders.add(order)
            if order not in factors:
                factors[order] = cyclotomic_bits(order, factors)
            factor = factors[order]
            raw_clearance = clearance["first_difference_a_clearing_order"]
            raw_clearances[raw_clearance] += 1
            joint_norm = 0
            norm_clearance = None
            joint_gcd = factor
            common_norm_gcd = factor
            common_norm_clearance = None
            gram_gcd = factor
            gram_clearance = None
            component_vectors = []
            for trace in clearance["residual_trace"]:
                trace_norm = 0
                components = [
                    int(component_hex, 16)
                    for component_hex in trace["raw_component_hex"]
                ]
                for component in components:
                    trace_norm ^= cyclic_star_product_bits(
                        component, component, order
                    )
                joint_norm ^= trace_norm
                joint_gcd = gcd_bits(factor, joint_norm)
                if joint_gcd == 1 and norm_clearance is None:
                    norm_clearance = trace["a"]
                common_norm_gcd = gcd_bits(common_norm_gcd, trace_norm)
                if common_norm_gcd == 1 and common_norm_clearance is None:
                    common_norm_clearance = trace["a"]
                for previous in component_vectors:
                    forward = 0
                    backward = 0
                    for component, previous_component in zip(
                        components, previous
                    ):
                        forward ^= cyclic_star_product_bits(
                            component, previous_component, order
                        )
                        backward ^= cyclic_star_product_bits(
                            previous_component, component, order
                        )
                    gram_gcd = gcd_bits(gram_gcd, forward)
                    gram_gcd = gcd_bits(gram_gcd, backward)
                gram_gcd = gcd_bits(gram_gcd, trace_norm)
                component_vectors.append(components)
                if gram_gcd == 1 and gram_clearance is None:
                    gram_clearance = trace["a"]

            norm_clearances[norm_clearance] += 1
            common_norm_clearances[common_norm_clearance] += 1
            gram_clearances[gram_clearance] += 1
            transitions[(raw_clearance, norm_clearance)] += 1
            if norm_clearance is None:
                exceptions.append(
                    {
                        "p": row["p"],
                        "order": order,
                        "raw_clearance": raw_clearance,
                        "trace_length": len(clearance["residual_trace"]),
                        "joint_norm_hex": hex(joint_norm),
                        "joint_gcd_hex": hex(joint_gcd),
                        "joint_gcd_degree": joint_gcd.bit_length() - 1,
                    }
                )
            if common_norm_clearance is None:
                common_norm_exceptions.append(
                    {
                        "p": row["p"],
                        "order": order,
                        "raw_clearance": raw_clearance,
                        "trace_length": len(clearance["residual_trace"]),
                        "common_norm_gcd_hex": hex(common_norm_gcd),
                        "common_norm_gcd_degree": (
                            common_norm_gcd.bit_length() - 1
                        ),
                    }
                )
            if gram_clearance is None:
                gram_exceptions.append(
                    {
                        "p": row["p"],
                        "order": order,
                        "raw_clearance": raw_clearance,
                        "trace_length": len(clearance["residual_trace"]),
                        "gram_gcd_hex": hex(gram_gcd),
                        "gram_gcd_degree": gram_gcd.bit_length() - 1,
                    }
                )

    result = {
        "input": str(args.input),
        "input_sha256": hashlib.sha256(input_bytes).hexdigest(),
        "criterion": "gcd(Phi_d, sum_{a,c} P_{a,c}(X)P_{a,c}(X^-1)) == 1",
        "records": records,
        "unit_by_available_trace": records - len(exceptions),
        "exceptions": exceptions,
        "common_norm_unit_by_available_trace": (
            records - len(common_norm_exceptions)
        ),
        "common_norm_exceptions": common_norm_exceptions,
        "gram_unit_by_available_trace": records - len(gram_exceptions),
        "gram_exceptions": gram_exceptions,
        "raw_clearance_counts": sorted_counter(raw_clearances),
        "joint_norm_clearance_counts": sorted_counter(norm_clearances),
        "common_norm_clearance_counts": sorted_counter(
            common_norm_clearances
        ),
        "gram_clearance_counts": sorted_counter(gram_clearances),
        "raw_to_joint_counts": {
            f"{raw}->{norm}": count
            for (raw, norm), count in sorted(
                transitions.items(), key=lambda item: str(item[0])
            )
        },
        "distinct_orders": len(tested_orders),
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(args.output, result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
