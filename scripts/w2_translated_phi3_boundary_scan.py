#!/usr/bin/env python3
"""Factor-free Phi_3 scan for the translated antipodal boundary edge.

For the valid translated family, put m=(p-1)/2 and use the antipodal pole
parameters u=m,m+1, based at s=-a.  Both endpoints are in U for
1<=a<(m+1)/2.  This scanner asks for the first a whose content is nonzero
at a root of Phi_3=X^2+X+1.

No polynomial factorization or FFT is needed.  Exponents on each of the
square/nonsquare point orbits are reduced modulo 3, giving two residues in
F4.  The Bose relative-difference-set normalization makes their Hermitian
norm exactly c(zeta)c(zeta^-1), which is 0 or 1.
"""
from __future__ import annotations

import argparse
import json
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
from e1_gmin_m4_prop15613 import _finv  # noqa: E402
from w2_pole_fourier_fast import named_z_without_conference  # noqa: E402
from w2_translated_antipodal_norm_scan import (  # noqa: E402
    apply_pole,
    is_prime,
    pole_action_data,
)


def gf4_mul(x: int, y: int) -> int:
    """Multiply bit-basis residues modulo Z^2+Z+1."""
    a, b = x & 1, (x >> 1) & 1
    c, d = y & 1, (y >> 1) & 1
    return (a * c ^ b * d) | ((a * d ^ b * c ^ b * d) << 1)


def gf4_star(x: int) -> int:
    """Apply Z -> Z^-1=Z^2=Z+1."""
    a, b = x & 1, (x >> 1) & 1
    return (a ^ b) | (b << 1)


def phi3_residue(v: np.ndarray, points: np.ndarray) -> int:
    """Evaluate sum_j v(points[j]) Z^j in F2[Z]/Phi_3."""
    seq = v[points]
    parity = [int(seq[r::3].sum()) & 1 for r in range(3)]
    # Z^2=Z+1.
    return (parity[0] ^ parity[2]) | ((parity[1] ^ parity[2]) << 1)


def phi3_norm(pair: np.ndarray, square: np.ndarray, nonsquare: np.ndarray) -> int:
    square_value = phi3_residue(pair, square)
    nonsquare_value = phi3_residue(pair, nonsquare)
    value = gf4_mul(square_value, gf4_star(square_value))
    value ^= gf4_mul(nonsquare_value, gf4_star(nonsquare_value))
    if value not in (0, 1):
        raise AssertionError(f"Hermitian norm escaped F2: {value}")
    return value


def record(job: tuple[int, int]) -> dict:
    p, max_a = job
    t0 = time.time()
    _z, _bits, q, mul, _add, chi, lam = named_z_without_conference(p)
    q2, _mul2, _add2, _chi2, _frob, _norm, ia, ib = field_ctx(p)
    if q2 != q:
        raise AssertionError("field contexts disagree")

    sig = next(value for value in range(1, q) if chi(value) == -1)
    sinv = _finv(mul, q, sig)
    field_points = np.arange(q, dtype=np.int64)
    point_a = field_points % p
    point_b = field_points // p
    sinv_a, sinv_b = sinv % p, sinv // p
    ell = (
        sinv_a * point_b
        + sinv_b * point_a
        + sinv_b * point_b * ia
    ) % p

    inv_fp = np.zeros(p, dtype=np.int64)
    legendre = np.zeros(p, dtype=np.int8)
    for value in range(1, p):
        inv_fp[value] = pow(value, p - 2, p)
        legendre[value] = 1 if pow(value, (p - 1) // 2, p) == 1 else -1

    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    ncoord = (q - 1) // 2
    square = np.empty(ncoord, dtype=np.int64)
    point = 1
    for j in range(ncoord):
        square[j] = point
        point = mul(gen, point)
    nonsquare = np.fromiter(
        (mul(omega, point) for point in square),
        dtype=np.int64,
        count=ncoord,
    )

    half = (p - 1) // 2
    actions = []
    for u in (half, p - half):
        pole_t = lam * pow(u, p - 2, p) % p
        actions.append(pole_action_data(p, pole_t, ia, ib, inv_fp, legendre))

    valid_stop = half // 2 + 1
    scan_stop = valid_stop if max_a == 0 else min(valid_stop, max_a + 1)
    bad_prefix = []
    first_good_a = None
    for a in range(1, scan_stop):
        shifted = (ell + a) % p
        z = np.empty(q + 1, dtype=np.int8)
        z[0] = -1
        z[1:] = np.where(shifted <= half, 1, -1).astype(np.int8)
        bits = ((1 - z) // 2).astype(np.uint8)
        pair = None
        for perm, switch in actions:
            in_u, wfn = apply_pole(z, bits, perm, switch)
            if not in_u:
                raise AssertionError(f"boundary endpoint outside U: p={p}, a={a}")
            pair = wfn.copy() if pair is None else pair ^ wfn
        good = phi3_norm(pair, square, nonsquare) == 1
        bad_prefix.append(not good)
        if good:
            first_good_a = a
            break

    return {
        "p": p,
        "boundary_u": half,
        "max_a": max_a,
        "first_good_a": first_good_a,
        "tested_bad_prefix": bad_prefix,
        "elapsed_seconds": time.time() - t0,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=5)
    ap.add_argument("--stop", type=int, required=True)
    ap.add_argument("--max-a", type=int, default=64)
    ap.add_argument("--workers", type=int, default=min(16, os.cpu_count() or 1))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    primes = [p for p in range(max(5, args.start), args.stop + 1) if is_prime(p)]
    jobs = [(p, args.max_a) for p in primes]
    t0 = time.time()
    if args.workers == 1:
        rows = list(map(record, jobs))
    else:
        with ProcessPoolExecutor(
            max_workers=args.workers,
            mp_context=mp.get_context("fork"),
        ) as ex:
            rows = list(ex.map(record, jobs, chunksize=1))
    failures = [row["p"] for row in rows if row["first_good_a"] is None]
    finite = [row["first_good_a"] for row in rows if row["first_good_a"] is not None]
    out = {
        "range": [args.start, args.stop],
        "max_a": args.max_a,
        "n_primes": len(rows),
        "failures": failures,
        "max_first_good_a": max(finite) if finite else None,
        "rows": rows,
        "elapsed_seconds": time.time() - t0,
    }
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"p={args.start}..{args.stop} primes={len(rows)} "
        f"max-first-good-a={out['max_first_good_a']} "
        f"failures={failures} seconds={out['elapsed_seconds']:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
