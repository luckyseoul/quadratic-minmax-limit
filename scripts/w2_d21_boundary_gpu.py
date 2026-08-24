#!/usr/bin/env python3
"""GPU search for a failure of the first-three d=21 boundary criterion.

For p == 29 (mod 84), the proved Jacobi-norm lemma makes the nonsquare Bose
generator residue nonzero at both reciprocal sextics of Phi_21.  Therefore a
boundary scalar vanishes at a sextic exactly when its nonsquare four-line word
does.  This kernel evaluates an order-42 character directly on the twelve
selected affine lines; no p^2 orbit, primitive element, or inverse table is
constructed.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT / "src"))

from io_atomic import write_json_atomic  # noqa: E402


ORDER = 21
FACTORS = (0x57, 0x75)


KERNEL_SOURCE = r"""
__device__ __forceinline__ unsigned int mul_mod(
    unsigned int left, unsigned int right, unsigned int p) {
    return (unsigned int)(((unsigned long long)left * right) % p);
}

__device__ __forceinline__ unsigned int pow_mod(
    unsigned int value, unsigned int exponent, unsigned int p) {
    unsigned int result = 1U;
    while (exponent) {
        if (exponent & 1U) result = mul_mod(result, value, p);
        value = mul_mod(value, value, p);
        exponent >>= 1U;
    }
    return result;
}

__device__ __forceinline__ unsigned long long field_mul(
    unsigned long long left,
    unsigned long long right,
    unsigned int p,
    unsigned int ib) {
    const unsigned int a0 = (unsigned int)(left % p);
    const unsigned int a1 = (unsigned int)(left / p);
    const unsigned int b0 = (unsigned int)(right % p);
    const unsigned int b1 = (unsigned int)(right / p);
    const unsigned int c0 = (
        mul_mod(a0, b0, p) +
        mul_mod(mul_mod(a1, b1, p), ib, p)) % p;
    const unsigned int c1 = (
        mul_mod(a0, b1, p) + mul_mod(a1, b0, p)) % p;
    return (unsigned long long)c0 + (unsigned long long)p * c1;
}

__device__ __forceinline__ unsigned long long field_pow(
    unsigned long long value,
    unsigned long long exponent,
    unsigned int p,
    unsigned int ib) {
    unsigned long long result = 1ULL;
    while (exponent) {
        if (exponent & 1ULL) result = field_mul(result, value, p, ib);
        value = field_mul(value, value, p, ib);
        exponent >>= 1ULL;
    }
    return result;
}

__device__ __forceinline__ unsigned long long field_inverse(
    unsigned long long value, unsigned int p, unsigned int ib) {
    const unsigned int a = (unsigned int)(value % p);
    const unsigned int b = (unsigned int)(value / p);
    const unsigned int norm = (
        mul_mod(a, a, p) + p -
        mul_mod(mul_mod(b, b, p), ib, p)) % p;
    if (!norm) return 0ULL;
    const unsigned int inverse_norm = pow_mod(norm, p - 2U, p);
    const unsigned int inverse_a = mul_mod(a, inverse_norm, p);
    const unsigned int inverse_b = mul_mod(b ? p - b : 0U, inverse_norm, p);
    return (unsigned long long)inverse_a + (unsigned long long)p * inverse_b;
}

extern "C" __global__ void d21_selected_lines(
    unsigned int p,
    unsigned int ib,
    unsigned int pole_t,
    unsigned long long character_exponent,
    unsigned long long sigma_character_inverse,
    unsigned int n_slots,
    const unsigned int* levels,
    const unsigned long long* roots,
    unsigned int* output) {
    const unsigned long long count = (unsigned long long)n_slots * p;
    const unsigned long long stride =
        (unsigned long long)gridDim.x * blockDim.x;
    for (unsigned long long index =
             (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
         index < count;
         index += stride) {
        const unsigned int slot = (unsigned int)(index / p);
        const unsigned int z = (unsigned int)(index - (unsigned long long)slot * p);
        const unsigned int level = levels[slot];
        const unsigned int image_a = mul_mod(ib, level, p);
        const unsigned int image_b = z;
        unsigned int denominator_a = mul_mod(pole_t, image_a, p);
        denominator_a = denominator_a ? denominator_a - 1U : p - 1U;
        const unsigned int denominator_b = mul_mod(pole_t, image_b, p);
        const unsigned long long denominator =
            (unsigned long long)denominator_a +
            (unsigned long long)p * denominator_b;
        const unsigned long long denominator_inverse =
            field_inverse(denominator, p, ib);
        if (!denominator_inverse) continue;
        const unsigned long long image =
            (unsigned long long)image_a + (unsigned long long)p * image_b;
        const unsigned long long point =
            field_mul(image, denominator_inverse, p, ib);
        const unsigned long long half_character =
            field_pow(point, character_exponent, p, ib);
        // The 21st power is the quadratic character.  Retain only nonsquares.
        if (field_pow(half_character, 21ULL, p, ib) != p - 1U) continue;
        unsigned long long character =
            field_mul(half_character, half_character, p, ib);
        character = field_mul(character, sigma_character_inverse, p, ib);
        unsigned int residue = 21U;
        for (unsigned int candidate = 0; candidate < 21U; ++candidate) {
            if (character == roots[candidate]) {
                residue = candidate;
                break;
            }
        }
        if (residue < 21U) atomicXor(&output[slot * 21U + residue], 1U);
    }
}
"""


def field_multiply(left: int, right: int, p: int, ib: int) -> int:
    a0, a1 = left % p, left // p
    b0, b1 = right % p, right // p
    return (
        (a0 * b0 + (a1 * b1 % p) * ib) % p
        + p * ((a0 * b1 + a1 * b0) % p)
    )


def field_power(value: int, exponent: int, p: int, ib: int) -> int:
    result = 1
    while exponent:
        if exponent & 1:
            result = field_multiply(result, value, p, ib)
        value = field_multiply(value, value, p, ib)
        exponent >>= 1
    return result


def remainder_bits(value: int, modulus: int) -> int:
    degree = modulus.bit_length() - 1
    while value and value.bit_length() - 1 >= degree:
        value ^= modulus << (value.bit_length() - 1 - degree)
    return value


def polynomial_bits(coefficients: np.ndarray) -> int:
    return sum(int(value & 1) << index for index, value in enumerate(coefficients))


def field_setup(p: int) -> tuple[int, int, np.ndarray]:
    ib = next(
        value
        for value in range(2, p)
        if pow(value, (p - 1) // 2, p) == p - 1
    )
    sigma = p
    q = p * p
    character_exponent = (q - 1) // 42
    sigma_character = field_power(sigma, (q - 1) // 21, p, ib)
    sigma_character_inverse = field_power(sigma_character, 20, p, ib)
    root = None
    for b in range(1, 16):
        for a in range(1, 16):
            candidate = field_power(a + p * b, (q - 1) // 21, p, ib)
            if (
                field_power(candidate, 21, p, ib) == 1
                and field_power(candidate, 7, p, ib) != 1
                and field_power(candidate, 3, p, ib) != 1
            ):
                root = candidate
                break
        if root is not None:
            break
    if root is None:
        raise ArithmeticError(f"p={p}: failed to find a primitive 21st root")
    roots = []
    value = 1
    for _ in range(21):
        roots.append(value)
        value = field_multiply(value, root, p, ib)
    return ib, sigma_character_inverse, np.asarray(roots, dtype=np.uint64)


def levels(p: int) -> np.ndarray:
    half = (p - 1) // 2
    result = []
    for a in range(1, 4):
        enter = (half - a) % p
        leave = (p - 1 - a) % p
        result.extend((enter, leave, (-enter) % p, (-leave) % p))
    if len(set(result)) != 12:
        raise AssertionError(f"p={p}: selected lines are not distinct")
    return np.asarray(result, dtype=np.uint32)


def line_bins_for_levels(
    p: int,
    ib: int,
    pole_t: int,
    sigma_character_inverse: int,
    roots: np.ndarray,
    kernel,
    line_levels: np.ndarray,
) -> np.ndarray:
    import cupy as cp

    d_levels = cp.asarray(line_levels)
    d_roots = cp.asarray(roots)
    output = cp.zeros((len(line_levels), 21), dtype=cp.uint32)
    threads = 256
    blocks = (len(line_levels) * p + threads - 1) // threads
    kernel(
        (blocks,),
        (threads,),
        (
            np.uint32(p),
            np.uint32(ib),
            np.uint32(pole_t),
            np.uint64((p * p - 1) // 42),
            np.uint64(sigma_character_inverse),
            np.uint32(len(line_levels)),
            d_levels,
            d_roots,
            output,
        ),
    )
    return (cp.asnumpy(output) & 1).astype(np.uint8)


def record(p: int, kernel) -> dict:
    started = time.perf_counter()
    ib, sigma_character_inverse, roots = field_setup(p)
    pole_t = (-2 * pow(ib, p - 2, p)) % p
    line_bins = line_bins_for_levels(
        p,
        ib,
        pole_t,
        sigma_character_inverse,
        roots,
        kernel,
        levels(p),
    )
    remainders = [[], []]
    polynomials = []
    for a in range(3):
        delta = np.bitwise_xor.reduce(line_bins[4 * a : 4 * a + 4], axis=0)
        polynomial = polynomial_bits(delta)
        polynomials.append(hex(polynomial))
        for factor_index, factor in enumerate(FACTORS):
            remainders[factor_index].append(remainder_bits(polynomial, factor))
    zero_triples = [all(value == 0 for value in side) for side in remainders]
    return {
        "p": p,
        "ib": ib,
        "pole_t": pole_t,
        "nonsquare_component_hex": polynomials,
        "remainders": remainders,
        "zero_triples": zero_triples,
        "simultaneous_zero_triples": all(zero_triples),
        "elapsed_seconds": time.perf_counter() - started,
    }


def staged_record(p: int, kernel) -> dict:
    """Stop as soon as one difference is nonzero at either sextic."""
    started = time.perf_counter()
    ib, sigma_character_inverse, roots = field_setup(p)
    pole_t = (-2 * pow(ib, p - 2, p)) % p
    all_levels = levels(p)
    remainders = [[], []]
    polynomials = []
    global_zero_differences = []
    for difference in range(3):
        line_bins = line_bins_for_levels(
            p,
            ib,
            pole_t,
            sigma_character_inverse,
            roots,
            kernel,
            all_levels[4 * difference : 4 * difference + 4],
        )
        delta = np.bitwise_xor.reduce(line_bins, axis=0)
        polynomial = polynomial_bits(delta)
        polynomials.append(hex(polynomial))
        factor_remainders = [
            remainder_bits(polynomial, factor) for factor in FACTORS
        ]
        for factor_index, value in enumerate(factor_remainders):
            remainders[factor_index].append(value)
        globally_zero = all(value == 0 for value in factor_remainders)
        global_zero_differences.append(globally_zero)
        if not globally_zero:
            break
    simultaneous = len(global_zero_differences) == 3 and all(
        global_zero_differences
    )
    return {
        "p": p,
        "ib": ib,
        "pole_t": pole_t,
        "nonsquare_component_hex": polynomials,
        "remainders": remainders,
        "global_zero_differences": global_zero_differences,
        "differences_evaluated": len(global_zero_differences),
        "simultaneous_zero_triples": simultaneous,
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int)
    parser.add_argument("--stop", type=int)
    parser.add_argument("--primes")
    parser.add_argument("--retain", action="store_true")
    parser.add_argument("--staged-counterexample", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.primes:
        primes = sorted(set(map(int, args.primes.split(","))))
    else:
        if args.start is None or args.stop is None:
            raise ValueError("start and stop are required without explicit primes")
        primes = [
            int(p)
            for p in sp.primerange(max(29, args.start), args.stop + 1)
            if p % 84 == 29
        ]
    import cupy as cp

    kernel = cp.RawKernel(KERNEL_SOURCE, "d21_selected_lines")
    rows = []
    retained = []
    started = time.perf_counter()
    for index, p in enumerate(primes, 1):
        row = staged_record(p, kernel) if args.staged_counterexample else record(p, kernel)
        one_sided = (not args.staged_counterexample) and any(row["zero_triples"])
        if args.retain or one_sided or row["simultaneous_zero_triples"]:
            retained.append(row)
        rows.append(
            {
                "p": p,
                "simultaneous_zero_triples": row["simultaneous_zero_triples"],
                "elapsed_seconds": row["elapsed_seconds"],
                **(
                    {
                        "global_zero_differences": row["global_zero_differences"],
                        "differences_evaluated": row["differences_evaluated"],
                    }
                    if args.staged_counterexample
                    else {"zero_triples": row["zero_triples"]}
                ),
            }
        )
        noteworthy = row["simultaneous_zero_triples"] or (
            not args.staged_counterexample and any(row["zero_triples"])
        )
        if noteworthy or index % 25 == 0 or index == len(primes):
            print(
                f"[{index}/{len(primes)}] p={p} "
                + (
                    f"global_prefix={row['global_zero_differences']} "
                    if args.staged_counterexample
                    else f"zero={row['zero_triples']} "
                )
                + f"seconds={row['elapsed_seconds']:.4f}",
                flush=True,
            )
        if row["simultaneous_zero_triples"]:
            print(f"COUNTEREXAMPLE p={p}", flush=True)
            break
    result = {
        "range": None if args.primes else [args.start, args.stop],
        "explicit_primes": primes if args.primes else None,
        "congruence_class": "p == 29 (mod 84)",
        "criterion": "three nonsquare residues vanish at both reciprocal sextics",
        "search_mode": "staged-counterexample" if args.staged_counterexample else "full",
        "n_primes_requested": len(primes),
        "n_primes_completed": len(rows),
        "one_sided_zero_triples": [
            row
            for row in retained
            if not args.staged_counterexample and any(row["zero_triples"])
        ],
        "counterexamples": [
            row for row in retained if row["simultaneous_zero_triples"]
        ],
        "rows": rows,
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(args.output, result)
    print(
        f"done completed={len(rows)} one_sided={len(result['one_sided_zero_triples'])} "
        f"counterexamples={len(result['counterexamples'])} "
        f"seconds={result['elapsed_seconds']:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
