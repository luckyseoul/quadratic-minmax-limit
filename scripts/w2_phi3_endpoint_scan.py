#!/usr/bin/env python3
"""Search for a fixed-endpoint failure at the W2 Phi_3 gate.

Phi_3=X^2+X+1 is a singleton factor-orbit for every odd prime p != 3.
This scanner skips the rest of g and asks how many consecutive normalized
pole parameters u=(p+1)/2+j are Phi_3-bad before the first good endpoint.
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import named_gamma  # noqa: E402
from w2_pole_fourier_fast import (  # noqa: E402
    fourier_factor_mask,
    fourier_factor_residues,
    named_z_without_conference,
    switched_wfn,
)


PHI3 = [[1, 1, 1]]


def is_prime(p: int) -> bool:
    if p < 2:
        return False
    if p % 2 == 0:
        return p == 2
    for d in range(3, math.isqrt(p) + 1, 2):
        if p % d == 0:
            return False
    return True


def endpoint_record(job: tuple[int, int]) -> dict:
    p, max_offsets = job
    z, bits, q, mul, add, chi, lam = named_z_without_conference(p)
    omega = _primitive(mul, q)
    gen = mul(omega, omega)
    ncoord = (q - 1) // 2
    square_points = []
    point = 1
    for _ in range(ncoord):
        square_points.append(point)
        point = mul(gen, point)
    nonsquare_points = [mul(omega, point) for point in square_points]

    gamma, _q2, _mul2, _b = named_gamma(p)
    gamma_square = fourier_factor_residues(gamma, square_points, PHI3)[0]
    gamma_nonsquare = fourier_factor_residues(gamma, nonsquare_points, PHI3)[0]
    if gamma_square == 0 and gamma_nonsquare == 0:
        raise AssertionError(f"Phi_3 gamma projection vanished at p={p}")
    use_square = [gamma_square != 0]

    lo = (p + 1) // 2
    bad = []
    first_good = None
    for offset in range(min(max_offsets, (p - 1) // 2)):
        u = lo + offset
        t = lam * pow(u, p - 2, p) % p
        in_u, wfn = switched_wfn(p, t, z, bits, q, mul, add, chi)
        if not in_u:
            raise AssertionError(f"endpoint outside U at p={p}, offset={offset}")
        mask = fourier_factor_mask(
            wfn, square_points, nonsquare_points, PHI3, use_square
        )
        is_bad = mask == 1
        bad.append(is_bad)
        if not is_bad:
            first_good = offset
            break
    return {
        "p": p,
        "pole_lambda": lam,
        "first_good_offset": first_good,
        "tested_bad_prefix": bad,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=5)
    ap.add_argument("--stop", type=int, required=True)
    ap.add_argument("--max-offsets", type=int, default=8)
    ap.add_argument("--workers", type=int, default=min(16, os.cpu_count() or 1))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    primes = [
        p for p in range(max(5, args.start), args.stop + 1)
        if p != 3 and is_prime(p)
    ]
    t0 = time.time()
    jobs = [(p, args.max_offsets) for p in primes]
    if args.workers == 1:
        rows = list(map(endpoint_record, jobs))
    else:
        with ProcessPoolExecutor(
            max_workers=args.workers,
            mp_context=mp.get_context("fork"),
        ) as ex:
            rows = list(ex.map(endpoint_record, jobs, chunksize=1))
    failures = [row["p"] for row in rows if row["first_good_offset"] is None]
    finite = [row["first_good_offset"] for row in rows if row["first_good_offset"] is not None]
    out = {
        "range": [args.start, args.stop],
        "max_offsets": args.max_offsets,
        "n_primes": len(rows),
        "max_first_good_offset": max(finite) if finite else None,
        "failures": failures,
        "rows": rows,
        "elapsed_seconds": time.time() - t0,
    }
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"p={args.start}..{args.stop} primes={len(rows)} "
        f"max-first-good={out['max_first_good_offset']} "
        f"failures={failures} seconds={out['elapsed_seconds']:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
