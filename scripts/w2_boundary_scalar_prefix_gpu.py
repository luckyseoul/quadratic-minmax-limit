#!/usr/bin/env python3
"""Exact all-orders scalar test for the first boundary differences.

For a boundary difference ``delta c`` and either multiplicative point orbit,

    Delta W_i(alpha) = Delta c(alpha) Gamma_i(alpha).

The Bose relative-difference-set norm says that the two ``Gamma_i`` cannot
both vanish at a nonprincipal odd root.  Thus ``Delta c(alpha)=0`` exactly
when both raw four-line residues ``Delta W_i(alpha)`` vanish.  Testing both
components and their reciprocals therefore avoids reconstructing the content.

This program folds the selected affine-line pullbacks directly modulo the
complete ambient odd order H.  A running NTL gcd with

    R_s(X^M),  M=(p+1)/2,  s=oddpart(p-1),  H=M*s,

tests every scalar-layer root order at once; there is no cyclotomic-order
cutoff and no FFT or floating-point arithmetic.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import _finv  # noqa: E402
from gf2x_ntl import available as ntl_available  # noqa: E402
from gf2x_ntl import field_two_orbits, gcd_bits  # noqa: E402
from io_atomic import write_json_atomic  # noqa: E402
from w2_endpoint_norm_gpu import scalar_factor_bits  # noqa: E402
from w2_low_order_atomic_gpu import (  # noqa: E402
    inverse_and_legendre_tables,
    oddpart,
    polynomial_bits,
    reciprocal_bits,
)
from w2_translated_antipodal_norm_scan import is_prime  # noqa: E402


KERNEL_SOURCE = r"""
extern "C" __global__
void selected_line_residue_bins(
    const unsigned int* points,
    unsigned long long orbit_length,
    unsigned int p,
    unsigned int ia,
    unsigned int ib,
    unsigned int sinv_a,
    unsigned int sinv_b,
    unsigned int pole_t,
    const unsigned int* inverse,
    const unsigned int* level_slots,
    unsigned int n_levels,
    unsigned int order,
    unsigned int* output) {
    const unsigned long long stride =
        (unsigned long long)gridDim.x * blockDim.x;
    for (unsigned long long j =
             (unsigned long long)blockIdx.x * blockDim.x + threadIdx.x;
         j < orbit_length;
         j += stride) {
        const unsigned int point = points[j];
        const unsigned long long a = point % p;
        const unsigned long long b = point / p;
        const unsigned long long t = pole_t;
        const unsigned long long e = (t * a + p - 1U) % p;
        const unsigned long long f = (t * b) % p;
        const unsigned int positive = (unsigned int)(
            (e * e + (unsigned long long)ia * e * f) % p);
        const unsigned int negative = (unsigned int)(
            ((unsigned long long)ib * f % p) * f % p);
        const unsigned int norm = (positive + p - negative) % p;
        if (norm == 0U) continue;
        const unsigned long long inv_norm = inverse[norm];
        const unsigned long long inv_a =
            ((e + (unsigned long long)ia * f) * inv_norm) % p;
        const unsigned long long inv_b = ((p - f) * inv_norm) % p;
        const unsigned int image_a = (unsigned int)(
            (a * inv_a + (unsigned long long)ib * b * inv_b) % p);
        const unsigned int image_b = (unsigned int)(
            (a * inv_b + b * inv_a +
             (unsigned long long)ia * b * inv_b) % p);
        const unsigned int level = (unsigned int)(
            ((unsigned long long)sinv_a * image_b +
             (unsigned long long)sinv_b * image_a +
             (unsigned long long)sinv_b * image_b * ia) % p);
        const unsigned int slot = level_slots[level];
        if (slot >= n_levels) continue;
        atomicXor(
            &output[(unsigned long long)slot * order +
                    (unsigned int)(j % order)],
            1U);
    }
}
"""


def difference_levels(p: int, a: int) -> tuple[list[int], list[int]]:
    half = (p - 1) // 2
    r_enter = (half - a) % p
    r_leave = (p - 1 - a) % p
    row = [r_enter, r_leave, (-r_enter) % p, (-r_leave) % p]
    levels = sorted(set(row))
    slots = {level: index for index, level in enumerate(levels)}
    return levels, [slots[level] for level in row]


def record(p: int, prefix: int, kernel) -> dict:
    import cupy as cp

    started = time.perf_counter()
    q, mul, add, chi, _frob, _norm, ia, ib = field_ctx(p)
    sig = next(value for value in range(1, q) if chi(value) == -1)
    sinv = _finv(mul, q, sig)
    sinv_a, sinv_b = sinv % p, sinv // p
    omega = _primitive(mul, q)
    square, nonsquare = field_two_orbits(
        p, ia, ib, mul(omega, omega), omega
    )
    scalar_order = oddpart(p - 1)
    projective_order = (p + 1) // 2
    ambient_order = projective_order * scalar_order
    if ambient_order != oddpart((p * p - 1) // 2):
        raise AssertionError("ambient odd-order decomposition failed")

    scalar_factor = scalar_factor_bits(projective_order, scalar_order)
    if scalar_order == 1:
        return {
            "p": p,
            "ambient_oddpart": ambient_order,
            "projective_order": projective_order,
            "scalar_order": scalar_order,
            "prefix": prefix,
            "scalar_layer_vacuous": True,
            "first_difference_a_clearing_scalar_layer": 0,
            "residual_trace": [],
            "elapsed_seconds": time.perf_counter() - started,
        }

    inverse, _legendre = inverse_and_legendre_tables(p)
    half = (p - 1) // 2
    pole_t = sinv_b * pow(half, p - 2, p) % p

    d_inverse = cp.asarray(inverse)
    threads = 256
    blocks = min(2048, (len(square) + threads - 1) // threads)

    common = scalar_factor
    residual_trace = []
    clearance = None
    n_selected_levels = 0
    for a in range(1, prefix + 1):
        levels, slots = difference_levels(p, a)
        n_selected_levels += len(levels)
        level_slots = np.full(p, np.uint32(0xFFFFFFFF), dtype=np.uint32)
        level_slots[levels] = np.arange(len(levels), dtype=np.uint32)
        d_level_slots = cp.asarray(level_slots)
        component_weights = []
        for points in (square, nonsquare):
            d_points = cp.asarray(points)
            output = cp.zeros((len(levels), ambient_order), dtype=cp.uint32)
            kernel(
                (blocks,),
                (threads,),
                (
                    d_points,
                    np.uint64(len(points)),
                    np.uint32(p),
                    np.uint32(ia),
                    np.uint32(ib),
                    np.uint32(sinv_a),
                    np.uint32(sinv_b),
                    np.uint32(pole_t),
                    d_inverse,
                    d_level_slots,
                    np.uint32(len(levels)),
                    np.uint32(ambient_order),
                    output,
                ),
            )
            line_bins = (cp.asnumpy(output) & 1).astype(np.uint8)
            delta = np.bitwise_xor.reduce(
                line_bins[slots], axis=0
            )
            component_weights.append(int(np.count_nonzero(delta)))
            common = gcd_bits(common, polynomial_bits(delta))
            if common != 1:
                common = gcd_bits(common, reciprocal_bits(delta))
            del d_points, output, line_bins, delta
            if common == 1:
                break
        del d_level_slots
        degree = common.bit_length() - 1
        residual_trace.append(
            {
                "a": a,
                "component_weights": component_weights,
                "residual_scalar_degree": degree,
                "residual_scalar_hex": hex(common) if degree <= 4096 else None,
            }
        )
        if common == 1:
            clearance = a
            break

    return {
        "p": p,
        "ambient_oddpart": ambient_order,
        "projective_order": projective_order,
        "scalar_order": scalar_order,
        "prefix": prefix,
        "n_selected_levels": n_selected_levels,
        "pole_t": pole_t,
        "scalar_layer_vacuous": False,
        "first_difference_a_clearing_scalar_layer": clearance,
        "residual_scalar_degree": common.bit_length() - 1,
        "residual_scalar_hex": (
            hex(common) if common.bit_length() <= 4097 else None
        ),
        "residual_trace": residual_trace,
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=5)
    parser.add_argument("--stop", type=int, required=True)
    parser.add_argument("--prefix", type=int, default=3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.prefix < 1:
        raise ValueError("prefix must be positive")
    if not ntl_available():
        raise RuntimeError("the exact all-orders scan requires the NTL bridge")

    import cupy as cp

    kernel = cp.RawKernel(KERNEL_SOURCE, "selected_line_residue_bins")
    primes = [
        p
        for p in range(max(5, args.start), args.stop + 1)
        if p % 12 == 5 and is_prime(p)
    ]
    started = time.perf_counter()
    rows = []
    for index, p in enumerate(primes, 1):
        row = record(p, args.prefix, kernel)
        rows.append(row)
        print(
            f"[{index}/{len(primes)}] p={p} "
            f"clear={row['first_difference_a_clearing_scalar_layer']} "
            f"residual_degree={row.get('residual_scalar_degree', 0)} "
            f"seconds={row['elapsed_seconds']:.3f}",
            flush=True,
        )
    result = {
        "range": [args.start, args.stop],
        "congruence_class": "p == 5 (mod 12)",
        "prefix": args.prefix,
        "order_scope": "all roots of R_s(X^M); no order cutoff",
        "gcd_backend": "ntl",
        "n_primes": len(rows),
        "failures": [
            {
                "p": row["p"],
                "residual_scalar_degree": row.get("residual_scalar_degree", 0),
                "residual_scalar_hex": row.get("residual_scalar_hex"),
            }
            for row in rows
            if row["first_difference_a_clearing_scalar_layer"] is None
        ],
        "rows": rows,
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(args.output, result)
    print(
        f"done primes={len(rows)} failures={len(result['failures'])}",
        flush=True,
    )


if __name__ == "__main__":
    main()
