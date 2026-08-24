#!/usr/bin/env python3
"""Atomic GPU census of translated-boundary four-line differences.

For boundary contents c_a, the difference c_{a+1}+c_a is recovered from
four affine-line pullbacks.  If their running Aut-gcd is one at order d,
then the genuine boundary family collectively clears every primitive-order-d
Aut orbit, without evaluating any interval halfspace.
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
from gf2x_ntl import field_orbits  # noqa: E402
from io_atomic import write_json_atomic  # noqa: E402
from w2_low_order_atomic_gpu import (  # noqa: E402
    cyclic_crosscorrelation,
    cyclotomic_bits,
    gamma_bins,
    inverse_and_legendre_tables,
    oddpart,
    polynomial_bits,
    reciprocal_bits,
)
from w2_translated_antipodal_norm_scan import f2_gcd_bits, is_prime  # noqa: E402


KERNEL_SOURCE = r"""
extern "C" __global__
void line_residue_bins(
    const unsigned int* points,
    unsigned long long orbit_length,
    unsigned int p,
    unsigned int ia,
    unsigned int ib,
    unsigned int sinv_a,
    unsigned int sinv_b,
    unsigned int pole_t,
    const unsigned int* inverse,
    const unsigned int* orders,
    const unsigned int* order_offsets,
    unsigned int n_orders,
    unsigned int total_order,
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
        const unsigned int base = level * total_order;
        for (unsigned int index = 0; index < n_orders; ++index) {
            const unsigned int d = orders[index];
            atomicXor(
                &output[base + order_offsets[index] + (unsigned int)(j % d)],
                1U);
        }
    }
}
"""


def record(p: int, requested_orders: list[int], kernel) -> dict:
    import cupy as cp

    started = time.perf_counter()
    q, mul, add, chi, _frob, _norm, ia, ib = field_ctx(p)
    sig = next(value for value in range(1, q) if chi(value) == -1)
    sinv = _finv(mul, q, sig)
    sinv_a, sinv_b = sinv % p, sinv // p
    omega = _primitive(mul, q)
    square, nonsquare, logarithm = field_orbits(
        p, ia, ib, mul(omega, omega), omega
    )
    scalar_order = oddpart(p - 1)
    projective_order = (p + 1) // 2
    ambient_order = scalar_order * projective_order
    orders = [d for d in requested_orders if ambient_order % d == 0]
    offsets: list[int] = []
    total_order = 0
    for d in orders:
        offsets.append(total_order)
        total_order += d
    if not orders:
        return {"p": p, "orders_tested": [], "clearances": []}

    gamma = gamma_bins(
        p, sig, mul, add, logarithm, orders, offsets, total_order
    )
    inverse, _legendre = inverse_and_legendre_tables(p)
    half = (p - 1) // 2
    pole_t = sinv_b * pow(half, p - 2, p) % p
    d_inverse = cp.asarray(inverse)
    d_orders = cp.asarray(orders, dtype=cp.uint32)
    d_offsets = cp.asarray(offsets, dtype=cp.uint32)
    line_bins = np.empty((2, p, total_order), dtype=np.uint8)
    threads = 256
    blocks = min(2048, (len(square) + threads - 1) // threads)
    for component, points in enumerate((square, nonsquare)):
        d_points = cp.asarray(points)
        output = cp.zeros((p, total_order), dtype=cp.uint32)
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
                d_orders,
                d_offsets,
                np.uint32(len(orders)),
                np.uint32(total_order),
                output,
            ),
        )
        line_bins[component] = (cp.asnumpy(output) & 1).astype(np.uint8)

    common = {d: cyclotomic_bits(d) for d in orders}
    clearance = {d: None for d in orders}
    residual_trace = {d: [] for d in orders}
    for a in range(1, max(1, half - 1)):
        r_enter = (half - a) % p
        r_leave = (p - 1 - a) % p
        levels = (r_enter, r_leave, (-r_enter) % p, (-r_leave) % p)
        for d, offset in zip(orders, offsets):
            if common[d] == 1:
                continue
            delta = np.zeros(d, dtype=np.uint8)
            for component in range(2):
                line_delta = np.zeros(d, dtype=np.uint8)
                for level in levels:
                    line_delta ^= line_bins[
                        component, level, offset : offset + d
                    ]
                delta ^= cyclic_crosscorrelation(
                    line_delta, gamma[component, offset : offset + d]
                )
            value = polynomial_bits(delta)
            common[d] = f2_gcd_bits(common[d], value)
            common[d] = f2_gcd_bits(common[d], reciprocal_bits(delta))
            residual_trace[d].append(
                {
                    "a": a,
                    "degree": common[d].bit_length() - 1,
                    "hex": hex(common[d]),
                }
            )
            if common[d] == 1:
                clearance[d] = a
        if all(value == 1 for value in common.values()):
            break

    clearances = [
        {
            "order": d,
            "layer": "projective" if projective_order % d == 0 else "scalar",
            "first_difference_a_clearing_order": clearance[d],
            "residual_degree": common[d].bit_length() - 1,
            "residual_hex": hex(common[d]),
            "residual_trace": residual_trace[d],
        }
        for d in orders
    ]
    return {
        "p": p,
        "orders_tested": orders,
        "boundary_u": half,
        "pole_t": pole_t,
        "clearances": clearances,
        "all_orders_cleared_by_differences": all(value == 1 for value in common.values()),
        "elapsed_seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=5)
    parser.add_argument("--stop", type=int, required=True)
    parser.add_argument("--orders")
    parser.add_argument("--max-order", type=int, default=63)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    orders = (
        sorted(set(map(int, args.orders.split(","))))
        if args.orders
        else list(range(3, args.max_order + 1, 2))
    )
    if any(d < 3 or d % 2 == 0 for d in orders):
        raise ValueError("orders must be odd and at least three")
    import cupy as cp

    kernel = cp.RawKernel(KERNEL_SOURCE, "line_residue_bins")
    primes = [
        p for p in range(max(5, args.start), args.stop + 1)
        if p % 12 == 5 and is_prime(p)
    ]
    started = time.perf_counter()
    rows = []
    for index, p in enumerate(primes, 1):
        row = record(p, orders, kernel)
        rows.append(row)
        print(
            f"[{index}/{len(primes)}] p={p} "
            f"cleared={row.get('all_orders_cleared_by_differences')} "
            f"seconds={row.get('elapsed_seconds', 0):.3f}",
            flush=True,
        )
    result = {
        "range": [args.start, args.stop],
        "requested_orders": orders,
        "n_primes": len(rows),
        "uncleared": [
            {"p": row["p"], **item}
            for row in rows for item in row.get("clearances", [])
            if item["residual_hex"] != "0x1"
        ],
        "rows": rows,
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(args.output, result)
    print(f"done primes={len(rows)} uncleared={len(result['uncleared'])}", flush=True)


if __name__ == "__main__":
    main()
