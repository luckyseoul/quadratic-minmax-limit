#!/usr/bin/env python3
"""GPU mod-four lifts of rare d=21 boundary-prefix records.

The counterexample kernel only needs parity and therefore uses atomic XOR.
For local 2-adic reconnaissance this companion kernel replaces that final
operation by atomic addition and retains the selected-line class counts
modulo four.  Reducing the resulting order-21 vector modulo ``Phi_21`` gives
a twelve-coordinate element of ``(Z/4Z)[zeta_21]``.  When its parity is zero,
the coefficientwise half-lift is also recorded.

Only the nonsquare component is needed: the proved d=21 Bose-generator lemma
makes its generator factor nonzero at both reciprocal sextics.
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
sys.path.insert(0, str(ROOT / "scripts"))

from io_atomic import write_json_atomic  # noqa: E402
from w2_d21_mod4_ring import decomposition_norm, ring_reduce  # noqa: E402
from w2_d21_boundary_gpu import (  # noqa: E402
    FACTORS,
    KERNEL_SOURCE,
    field_setup,
    levels,
    polynomial_bits,
    remainder_bits,
)


ORDER = 21
ATOMIC_XOR = "atomicXor(&output[slot * 21U + residue], 1U);"
ATOMIC_ADD = "atomicAdd(&output[slot * 21U + residue], 1U);"
if KERNEL_SOURCE.count(ATOMIC_XOR) != 1:
    raise RuntimeError("unexpected d21 kernel atomic-XOR source")
MOD4_KERNEL_SOURCE = KERNEL_SOURCE.replace(ATOMIC_XOR, ATOMIC_ADD)


def line_counts_mod4(
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
    output = cp.zeros((len(line_levels), ORDER), dtype=cp.uint32)
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
    return (cp.asnumpy(output) & 3).astype(np.uint8)


def record(p: int, kernel) -> dict:
    started = time.perf_counter()
    ib, sigma_character_inverse, roots = field_setup(p)
    pole_t = (-2 * pow(ib, p - 2, p)) % p
    all_levels = levels(p)
    differences = []
    for difference in range(3):
        counts = line_counts_mod4(
            p,
            ib,
            pole_t,
            sigma_character_inverse,
            roots,
            kernel,
            all_levels[4 * difference : 4 * difference + 4],
        )
        class_counts = np.asarray(
            np.sum(counts.astype(np.uint16), axis=0) & 3,
            dtype=np.uint8,
        )
        quotient = ring_reduce([int(value) for value in class_counts])
        parity_polynomial = polynomial_bits(class_counts & 1)
        factor_remainders = [
            remainder_bits(parity_polynomial, factor) for factor in FACTORS
        ]
        parity_zero = all(value == 0 for value in factor_remainders)
        quotient_even = all(value % 2 == 0 for value in quotient)
        if parity_zero != quotient_even:
            raise AssertionError(
                f"p={p} difference={difference + 1}: CRT/parity mismatch"
            )
        differences.append(
            {
                "a": difference + 1,
                "nonsquare_class_counts_mod4": [int(value) for value in class_counts],
                "quotient_mod4": quotient,
                "decomposition_norm_mod4": decomposition_norm(quotient),
                "factor_remainders_mod2": factor_remainders,
                "parity_zero": parity_zero,
                "half_lift_mod2": (
                    [value // 2 for value in quotient] if quotient_even else None
                ),
            }
        )
    prefix = [item["parity_zero"] for item in differences]
    return {
        "p": p,
        "ib": ib,
        "pole_t": pole_t,
        "global_zero_differences": prefix,
        "deep_prefix": prefix[:2] == [True, True],
        "simultaneous_zero_triples": all(prefix),
        "differences": differences,
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    primes = list(dict.fromkeys(map(int, args.primes.split(","))))
    invalid = [p for p in primes if p % 84 != 29 or not sp.isprime(p)]
    if invalid:
        raise ValueError(f"primes must be prime and 29 mod 84: {invalid}")

    import cupy as cp

    kernel = cp.RawKernel(MOD4_KERNEL_SOURCE, "d21_selected_lines")
    rows = []
    started = time.perf_counter()
    for index, p in enumerate(primes, 1):
        row = record(p, kernel)
        rows.append(row)
        print(
            f"[{index}/{len(primes)}] p={p} "
            f"prefix={row['global_zero_differences']} "
            f"seconds={row['elapsed_seconds']:.4f}",
            flush=True,
        )
    result = {
        "description": "atomic GPU mod-four nonsquare d=21 boundary lifts",
        "primes": primes,
        "n_primes": len(primes),
        "all_deep_prefixes": all(row["deep_prefix"] for row in rows),
        "counterexamples": [row for row in rows if row["simultaneous_zero_triples"]],
        "rows": rows,
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(args.output, result)
    print(
        f"done primes={len(rows)} deep={result['all_deep_prefixes']} "
        f"counterexamples={len(result['counterexamples'])} "
        f"seconds={result['elapsed_seconds']:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
