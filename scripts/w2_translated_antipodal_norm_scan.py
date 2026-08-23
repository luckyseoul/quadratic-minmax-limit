#!/usr/bin/env python3
"""Factorization-free scan of translated antipodal pole norms.

Write N=(p^2-1)/2=2^a m with m odd and g=(X^m+1)/(X+1).
For a finite W-vector w, its square/nonsquare cyclic autocorrelations,
folded modulo m and added over F2, evaluate at a root zeta of g as

    W_square(zeta) W_square(zeta^-1)
      + W_nonsquare(zeta) W_nonsquare(zeta^-1).

If w=c(D) gamma, the corresponding expression for gamma is exactly 1.
Indeed gamma restricted to F_{p^2}^* is R+H, where R is the Bose affine
(p+1,p-1,p,1) relative difference set and H is its trace-zero F_p^*-coset;
reducing the group-ring identity modulo 2 and projecting to square shifts
gives norm 1 at every nonprincipal odd character.  Thus this folded norm is
c(X)c(X^-1) modulo g.  It is coprime to g exactly when c is.  This scanner
applies the test to antipodal pole pairs based at the translated halfspaces
s=-a.  Both pole endpoints lie in U precisely on

    1 <= a < u <= (p-1)/2,   a+u > (p-1)/2.

The parallel edge u=(p-1)/2-r is scanned one representative at a time
and also through the xor of the norms with a<(p-1)/4.
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15611 import _v2  # noqa: E402
from e1_gmin_m4_prop15613 import _finv  # noqa: E402
from w2_pole_fourier_fast import (  # noqa: E402
    named_z_without_conference,
)


def poly_bits(coefficients: np.ndarray) -> int:
    """Pack little-endian F2 coefficients into a Python integer."""
    packed = np.packbits(coefficients, bitorder="little")
    return int.from_bytes(packed.tobytes(), "little")


def f2_gcd_bits(a: int, b: int) -> int:
    """Polynomial gcd over F2 using Python's word-parallel integers."""
    while b:
        if not a:
            return b
        shift = a.bit_length() - b.bit_length()
        if shift < 0:
            a, b = b, a
        else:
            a ^= b << shift
    return a


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    return all(p % d for d in range(3, math.isqrt(p) + 1, 2))


def folded_autocorrelation(v: np.ndarray, points: np.ndarray, m: int) -> np.ndarray:
    """Return cyclic autocorrelation over F2, folded modulo X^m+1."""
    seq = np.asarray(v[points], dtype=np.float64)
    spectrum = np.fft.fft(seq)
    corr = np.rint(np.fft.ifft(spectrum * np.conj(spectrum)).real)
    corr = corr.astype(np.int64) & 1
    return np.bitwise_xor.reduce(corr.reshape(-1, m), axis=0).astype(np.uint8)


def pole_action_data(
    p: int,
    t: int,
    ia: int,
    ib: int,
    inv_fp: np.ndarray,
    legendre: np.ndarray,
):
    """Vectorized permutation/signs for x -> x/(t*x-1).

    The field basis satisfies omega^2=ia*omega+ib.  For d=e+f*omega,
    d^-1=((e+ia*f)-f*omega)/Norm(d), where
    Norm(d)=e^2+ia*e*f-ib*f^2.  Since the norm lies in F_p, one inverse
    lookup replaces an exponentiation in F_{p^2} at every point.
    """
    q = p * p
    points = np.arange(q, dtype=np.int64)
    a = points % p
    b = points // p
    e = (t * a - 1) % p
    f = (t * b) % p
    norm = (e * e + ia * e * f - ib * f * f) % p
    zero = norm == 0
    inv_norm = inv_fp[norm]
    inv_a = ((e + ia * f) * inv_norm) % p
    inv_b = ((-f) * inv_norm) % p
    image_a = (a * inv_a + b * inv_b * ib) % p
    image_b = (a * inv_b + b * inv_a + b * inv_b * ia) % p
    image = image_a + p * image_b

    perm = np.empty(q + 1, dtype=np.int64)
    perm[0] = 1 + int(inv_fp[t])
    perm[1:] = 1 + image
    perm[1:][zero] = 0

    switch = np.empty(q + 1, dtype=np.int8)
    switch[0] = 1
    switch[1:] = legendre[norm]
    switch[1:][zero] = 1
    return perm, switch


def apply_pole(z: np.ndarray, bits: np.ndarray, perm, switch) -> tuple[bool, np.ndarray]:
    """Apply precomputed pole data and return U-membership and the W-vector."""
    y = switch * z[perm]
    ybits = ((1 - y) // 2).astype(np.uint8)
    in_u = int(ybits[0]) == 1 and int(ybits[1]) == 0
    diff = (bits ^ ybits) & 1
    wfn = diff[1:].copy()
    if diff[0]:
        wfn ^= 1
    return in_u, wfn


def record(job: tuple[int, int, int, bool]) -> dict:
    p, max_edges, max_a, stop_after_first_unit = job
    t0 = time.time()
    _z0, _b0, q, mul, add, chi, lam = named_z_without_conference(p)
    q2, _mul2, _add2, _chi2, _frob, _norm, ia, ib = field_ctx(p)
    if q2 != q:
        raise AssertionError("field contexts disagree")
    inv_fp = np.zeros(p, dtype=np.int64)
    legendre = np.zeros(p, dtype=np.int8)
    for value in range(1, p):
        inv_fp[value] = pow(value, p - 2, p)
        legendre[value] = 1 if pow(value, (p - 1) // 2, p) == 1 else -1
    sig = next(value for value in range(1, q) if chi(value) == -1)
    sinv = _finv(mul, q, sig)
    ell = np.fromiter((mul(sinv, x) // p for x in range(q)), dtype=np.int32)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    ncoord = (q - 1) // 2
    oddpart = ncoord >> _v2(ncoord)
    square_points = []
    point = 1
    for _ in range(ncoord):
        square_points.append(point)
        point = mul(gen, point)
    square_points = np.asarray(square_points, dtype=np.int64)
    nonsquare_points = np.fromiter(
        (mul(omega, point) for point in square_points),
        dtype=np.int64,
        count=ncoord,
    )
    # For odd m, (X^m+1)/(X+1)=1+X+...+X^(m-1).
    g = (1 << oddpart) - 1

    half = (p - 1) // 2
    common_gcd = g
    edge_rows = []
    first_unit_pair = None
    for offset in range(min(max_edges, max(0, (half - 1) // 2))):
        u = half - offset
        actions = []
        for opposite_u in (u, p - u):
            pole_t = lam * pow(opposite_u, p - 2, p) % p
            actions.append(
                pole_action_data(p, pole_t, ia, ib, inv_fp, legendre)
            )
        aggregate = np.zeros(oddpart, dtype=np.uint8)
        n_rows = 0
        unit_pairs = []
        a_stop = (half + 1) // 2
        if max_a > 0:
            a_stop = min(a_stop, max_a + 1)
        for a in range(offset + 1, a_stop):
            s = p - a
            shifted = (ell - s) % p
            z = np.empty(q + 1, dtype=np.int8)
            z[0] = -1
            z[1:] = np.where(shifted <= half, 1, -1).astype(np.int8)
            bits = ((1 - z) // 2).astype(np.uint8)
            pair = None
            for opposite_u, (perm, switch) in zip((u, p - u), actions):
                in_u, wfn = apply_pole(z, bits, perm, switch)
                if not in_u:
                    raise AssertionError(
                        f"translated antipodal endpoint missed U: "
                        f"p={p}, a={a}, u={opposite_u}"
                    )
                pair = wfn.copy() if pair is None else pair ^ wfn
            norm = folded_autocorrelation(pair, square_points, oddpart)
            norm ^= folded_autocorrelation(pair, nonsquare_points, oddpart)
            aggregate ^= norm
            n_rows += 1
            norm_gcd = f2_gcd_bits(poly_bits(norm), g)
            if norm_gcd == 1:
                unit_pairs.append([a, u])
                if first_unit_pair is None:
                    first_unit_pair = [a, u]
                if stop_after_first_unit:
                    break
        aggregate_bits = poly_bits(aggregate)
        aggregate_gcd = f2_gcd_bits(aggregate_bits, g)
        common_gcd = f2_gcd_bits(common_gcd, aggregate_bits)
        edge_rows.append(
            {
                "offset": offset,
                "u": u,
                "n_half_edge_rows": n_rows,
                "aggregate_norm_weight": int(aggregate.sum()),
                "aggregate_gcd_degree": aggregate_gcd.bit_length() - 1,
                "common_gcd_degree": common_gcd.bit_length() - 1,
                "unit_pairs": unit_pairs,
            }
        )
        if common_gcd == 1 and first_unit_pair is not None:
            break
    return {
        "p": p,
        "max_edges": max_edges,
        "max_a": max_a,
        "stopped_after_first_unit": stop_after_first_unit,
        "first_unit_pair": first_unit_pair,
        "aggregate_unit_ideal": common_gcd == 1,
        "final_common_gcd_degree": common_gcd.bit_length() - 1,
        "edges": edge_rows,
        "elapsed_seconds": time.time() - t0,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=5)
    ap.add_argument("--stop", type=int, required=True)
    ap.add_argument("--max-edges", type=int, default=8)
    ap.add_argument(
        "--max-a",
        type=int,
        default=0,
        help="scan only a <= this value (0 means the full half-edge)",
    )
    ap.add_argument("--stop-after-first-unit", action="store_true")
    ap.add_argument("--workers", type=int, default=min(16, os.cpu_count() or 1))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    primes = [p for p in range(max(5, args.start), args.stop + 1) if is_prime(p)]
    t0 = time.time()
    jobs = [
        (p, args.max_edges, args.max_a, args.stop_after_first_unit)
        for p in primes
    ]
    if args.workers == 1:
        rows = list(map(record, jobs))
    else:
        with ProcessPoolExecutor(
            max_workers=args.workers,
            mp_context=mp.get_context("fork"),
        ) as ex:
            rows = list(ex.map(record, jobs, chunksize=1))
    out = {
        "range": [args.start, args.stop],
        "max_edges": args.max_edges,
        "max_a": args.max_a,
        "stopped_after_first_unit": args.stop_after_first_unit,
        "n_primes": len(rows),
        "no_unit_pair": [row["p"] for row in rows if row["first_unit_pair"] is None],
        "aggregate_failures": [
            row["p"] for row in rows if not row["aggregate_unit_ideal"]
        ],
        "rows": rows,
        "elapsed_seconds": time.time() - t0,
    }
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"p={args.start}..{args.stop} primes={len(rows)} "
        f"no-unit={out['no_unit_pair']} "
        f"aggregate-failures={out['aggregate_failures']} "
        f"seconds={out['elapsed_seconds']:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
