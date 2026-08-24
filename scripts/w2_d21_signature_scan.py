#!/usr/bin/env python3
"""Exact mask census for the reciprocal sextics of Phi_21 over F_2.

For each eligible prime, bit 2*(a-1)+c records whether the raw component
Delta W_(a,c) is divisible by the indicated irreducible factor.  The full
mask 0x3f is the one-sided flat-family event; simultaneous full masks for
0x57 and its reciprocal 0x75 would be the actual W2 obstruction.
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

from io_atomic import write_json_atomic  # noqa: E402
from w2_boundary_line_direct import record  # noqa: E402
from w2_translated_antipodal_norm_scan import is_prime  # noqa: E402


FACTORS = (0x57, 0x75)
ORDER = 21
FULL_MASK = (1 << 6) - 1


def remainder_bits(value: int, modulus: int) -> int:
    """Polynomial remainder over F_2 for packed little-endian bits."""
    modulus_degree = modulus.bit_length() - 1
    while value and value.bit_length() - 1 >= modulus_degree:
        value ^= modulus << (value.bit_length() - 1 - modulus_degree)
    return value


def multiply_bits(left: int, right: int, modulus: int) -> int:
    """Multiply two elements of the packed degree-six quotient field."""
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


def prime_signature(p: int) -> dict:
    row = record(p, [ORDER], 3, False, True)
    if len(row["clearances"]) != 1:
        raise AssertionError(f"p={p}: missing d=21 scalar trace")
    trace = row["clearances"][0]["residual_trace"]
    if len(trace) != 3:
        raise AssertionError(f"p={p}: incomplete three-difference trace")
    masks = [0, 0]
    remainders = [[], []]
    component_hex = []
    for a_index, item in enumerate(trace):
        component_hex.append(item["raw_component_hex"])
        for component, value_hex in enumerate(item["raw_component_hex"]):
            value = int(value_hex, 16)
            bit = 1 << (2 * a_index + component)
            for factor_index, factor in enumerate(FACTORS):
                remainder = remainder_bits(value, factor)
                if component == 0:
                    remainders[factor_index].append([remainder])
                else:
                    remainders[factor_index][-1].append(remainder)
                if remainder == 0:
                    masks[factor_index] |= bit
    ratios = []
    nonsquare_projective = []
    nonsquare_zero_violations = []
    subfield_violations = []
    for factor_index, factor in enumerate(FACTORS):
        factor_ratios = set()
        inverse_x = power_bits(0x2, 62, factor)
        for difference_index, (square, nonsquare) in enumerate(
            remainders[factor_index], 1
        ):
            normalized_nonsquare = multiply_bits(
                nonsquare, inverse_x, factor
            )
            if power_bits(square, 8, factor) != square:
                subfield_violations.append(
                    [factor_index, difference_index, "square"]
                )
            if power_bits(normalized_nonsquare, 8, factor) != normalized_nonsquare:
                subfield_violations.append(
                    [factor_index, difference_index, "nonsquare_over_x"]
                )
            if nonsquare == 0:
                if square != 0:
                    nonsquare_zero_violations.append(
                        [factor_index, difference_index]
                    )
                continue
            factor_ratios.add(
                multiply_bits(
                    multiply_bits(0x2, square, factor),
                    power_bits(nonsquare, 62, factor),
                    factor,
                )
            )
        ratios.append(sorted(factor_ratios))
        normalized_vector = [
            multiply_bits(nonsquare, inverse_x, factor)
            for _square, nonsquare in remainders[factor_index]
        ]
        pivot = next((value for value in normalized_vector if value), None)
        if pivot is not None:
            pivot_inverse = power_bits(pivot, 62, factor)
            normalized_vector = [
                multiply_bits(value, pivot_inverse, factor)
                for value in normalized_vector
            ]
        nonsquare_projective.append(normalized_vector)
    return {
        "p": p,
        "masks": masks,
        "popcounts": [mask.bit_count() for mask in masks],
        "remainders": remainders,
        "component_ratios": ratios,
        "nonsquare_projective": nonsquare_projective,
        "nonsquare_zero_violations": nonsquare_zero_violations,
        "subfield_violations": subfield_violations,
        "component_hex": component_hex,
        "elapsed_seconds": row["elapsed_seconds"],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=29)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--primes")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--retain-popcount", type=int, default=4)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.workers < 1:
        raise ValueError("workers must be positive")
    if not 0 <= args.retain_popcount <= 6:
        raise ValueError("retain-popcount must lie in [0,6]")

    if args.primes:
        primes = sorted(set(map(int, args.primes.split(","))))
        if any(p % 84 != 29 or not is_prime(p) for p in primes):
            raise ValueError("every explicit prime must be prime and 29 mod 84")
    else:
        if args.stop is None:
            raise ValueError("--stop is required unless --primes is supplied")
        primes = [
            p
            for p in range(max(29, args.start), args.stop + 1)
            if p % 84 == 29 and is_prime(p)
        ]
    started = time.perf_counter()
    rows = []
    if args.workers == 1:
        for index, p in enumerate(primes, 1):
            rows.append(prime_signature(p))
            if index % 25 == 0 or index == len(primes):
                print(f"[{index}/{len(primes)}] p={p}", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(prime_signature, p): p for p in primes}
            for index, future in enumerate(as_completed(futures), 1):
                rows.append(future.result())
                if index % 25 == 0 or index == len(primes):
                    print(f"[{index}/{len(primes)}]", flush=True)
    rows.sort(key=lambda row: row["p"])

    mask_pairs = collections.Counter(tuple(row["masks"]) for row in rows)
    remainder_pairs = [collections.Counter(), collections.Counter()]
    remainder_pairs_by_difference = [
        [collections.Counter() for _ in range(3)] for _ in range(2)
    ]
    for row in rows:
        for factor_index in range(2):
            remainder_pairs[factor_index].update(
                tuple(pair) for pair in row["remainders"][factor_index]
            )
            for difference_index, pair in enumerate(
                row["remainders"][factor_index]
            ):
                remainder_pairs_by_difference[factor_index][
                    difference_index
                ][tuple(pair)] += 1
    retained = [
        row for row in rows if max(row["popcounts"]) >= args.retain_popcount
    ]
    result = {
        "range": None if args.primes else [args.start, args.stop],
        "explicit_primes": primes if args.primes else None,
        "congruence_class": "p == 29 (mod 84)",
        "order": ORDER,
        "factor_hex": [hex(factor) for factor in FACTORS],
        "full_mask_hex": hex(FULL_MASK),
        "n_primes": len(rows),
        "mask_pair_counts": {
            f"{left:#04x},{right:#04x}": count
            for (left, right), count in sorted(mask_pairs.items())
        },
        "component_ratio_counts": [
            {
                ("none" if not ratio else ",".join(hex(x) for x in ratio)): count
                for ratio, count in sorted(
                    collections.Counter(
                        tuple(row["component_ratios"][factor_index])
                        for row in rows
                    ).items()
                )
            }
            for factor_index in range(2)
        ],
        "component_ratio_pair_counts": {
            "|".join(
                "none" if not ratio else ",".join(hex(x) for x in ratio)
                for ratio in pair
            ): count
            for pair, count in sorted(
                collections.Counter(
                    tuple(tuple(ratio) for ratio in row["component_ratios"])
                    for row in rows
                ).items()
            )
        },
        "component_ratio_rows": [
            {
                "p": row["p"],
                "component_ratios": row["component_ratios"],
            }
            for row in rows
        ],
        "nonsquare_projective_pair_counts": {
            "|".join(
                ",".join(hex(value) for value in vector)
                for vector in pair
            ): count
            for pair, count in sorted(
                collections.Counter(
                    tuple(tuple(vector) for vector in row["nonsquare_projective"])
                    for row in rows
                ).items()
            )
        },
        "nonsquare_projective_rows": [
            {
                "p": row["p"],
                "vectors": row["nonsquare_projective"],
            }
            for row in rows
        ],
        "nonsquare_zero_without_square_zero": [
            {
                "p": row["p"],
                "violations": row["nonsquare_zero_violations"],
            }
            for row in rows
            if row["nonsquare_zero_violations"]
        ],
        "multiple_component_ratios": [
            {
                "p": row["p"],
                "component_ratios": row["component_ratios"],
            }
            for row in rows
            if any(len(ratio) > 1 for ratio in row["component_ratios"])
        ],
        "subfield_violations": [
            {"p": row["p"], "violations": row["subfield_violations"]}
            for row in rows
            if row["subfield_violations"]
        ],
        "component_remainder_pair_counts": [
            {
                f"{square:#04x},{nonsquare:#04x}": count
                for (square, nonsquare), count in sorted(counter.items())
            }
            for counter in remainder_pairs
        ],
        "component_remainder_pair_counts_by_difference": [
            [
                {
                    f"{square:#04x},{nonsquare:#04x}": count
                    for (square, nonsquare), count in sorted(counter.items())
                }
                for counter in factor_counters
            ]
            for factor_counters in remainder_pairs_by_difference
        ],
        "one_sided_full": [
            row for row in rows if FULL_MASK in row["masks"]
        ],
        "reciprocal_pair_failures": [
            row
            for row in rows
            if row["masks"] == [FULL_MASK, FULL_MASK]
        ],
        "retain_popcount": args.retain_popcount,
        "retained": retained,
        "line_cpu_seconds": sum(row["elapsed_seconds"] for row in rows),
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(args.output, result)
    print(
        f"done primes={len(rows)} mask_pairs={len(mask_pairs)} "
        f"one_sided_full={len(result['one_sided_full'])} "
        f"reciprocal_failures={len(result['reciprocal_pair_failures'])} "
        f"seconds={result['elapsed_seconds']:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
