#!/usr/bin/env python3
"""Frobenius reduction for the normalized Phi_3 endpoint at p=5 mod 12.

For the endpoint u=(p+1)/2, the pole parameter lies in F_p.  The named
halfspace and its pole image are therefore Frobenius invariant.  On the
nonsquare multiplicative orbit, Frobenius sends exponent class j modulo 3
to 2-j.  Thus the endpoint Phi_3 residue is controlled by the parity of one
class count C_0=C_2; class C_1 is automatically even.

This scanner verifies the exact symmetries, the surviving parity, and the
normalized content value without factoring the ambient polynomial g.
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

ROOT = Path(
    os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1])
)
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import named_gamma  # noqa: E402
from w2_pole_fourier_fast import named_z_without_conference  # noqa: E402
from w2_translated_antipodal_norm_scan import (  # noqa: E402
    apply_pole,
    is_prime,
    pole_action_data,
)
from w2_translated_phi3_boundary_scan import (  # noqa: E402
    gf4_mul,
    gf4_star,
    phi3_residue,
)


def record(p: int) -> dict:
    if p % 12 != 5:
        raise ValueError("this reduction requires p == 5 (mod 12)")
    t0 = time.time()
    z, bits, q, mul, _add, chi, lam = named_z_without_conference(p)
    q2, _mul2, _add2, _chi2, _frob, _norm, ia, ib = field_ctx(p)
    if q2 != q:
        raise AssertionError("field contexts disagree")

    inv_fp = np.zeros(p, dtype=np.int64)
    legendre = np.zeros(p, dtype=np.int8)
    for value in range(1, p):
        inv_fp[value] = pow(value, p - 2, p)
        legendre[value] = (
            1 if pow(value, (p - 1) // 2, p) == 1 else -1
        )

    # u=1/2 in F_p, so t=lambda/u=2*lambda also lies in F_p.  Hence its
    # projective action and quadratic switch commute with p-Frobenius.
    u = (p + 1) // 2
    pole_t = lam * pow(u, p - 2, p) % p
    perm, switch = pole_action_data(
        p, pole_t, ia, ib, inv_fp, legendre
    )
    in_u, wfn = apply_pole(z, bits, perm, switch)
    if not in_u:
        raise AssertionError(f"normalized endpoint missed U at p={p}")

    points = np.arange(q, dtype=np.int64)
    point_a = points % p
    point_b = points // p
    # If omega^2=ia*omega+ib, then omega^p=ia-omega.
    conjugates = (
        (point_a + ia * point_b) % p + ((-point_b) % p) * p
    )
    named_z_frobenius_mismatches = int(
        np.count_nonzero(z[1:] != z[1 + conjugates])
    )
    endpoint_frobenius_mismatches = int(
        np.count_nonzero(wfn != wfn[conjugates])
    )

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

    gamma, qg, _mulg, _b = named_gamma(p)
    if qg != q:
        raise AssertionError("gamma field disagrees")
    gamma_residues = [
        phi3_residue(gamma, square),
        phi3_residue(gamma, nonsquare),
    ]
    gamma_class_counts = [
        int(gamma[nonsquare[j::3]].sum()) for j in range(3)
    ]
    expected_gamma_class_counts = [
        (p + 1) // 6,
        (7 * p - 11) // 6,
        (p + 1) // 6,
    ]
    endpoint_residues = [
        phi3_residue(wfn, square),
        phi3_residue(wfn, nonsquare),
    ]
    active = [i for i, value in enumerate(gamma_residues) if value]
    if active != [1]:
        raise AssertionError(
            f"expected nonsquare Phi_3 component at p={p}: "
            f"{gamma_residues}"
        )
    content_value = gf4_mul(
        endpoint_residues[1], gf4_star(gamma_residues[1])
    )

    class_counts = [
        int(wfn[nonsquare[j::3]].sum()) for j in range(3)
    ]
    c0, c1, c2 = class_counts
    return {
        "p": p,
        "endpoint_u": u,
        "pole_t": pole_t,
        "named_z_frobenius_mismatches": named_z_frobenius_mismatches,
        "endpoint_frobenius_mismatches": endpoint_frobenius_mismatches,
        "gamma_phi3_residues": gamma_residues,
        "gamma_nonsquare_exponent_class_counts": gamma_class_counts,
        "expected_gamma_class_counts": expected_gamma_class_counts,
        "gamma_class_count_formula": (
            gamma_class_counts == expected_gamma_class_counts
        ),
        "endpoint_phi3_residues": endpoint_residues,
        "normalized_content_value": content_value,
        "nonsquare_exponent_class_counts": class_counts,
        "frobenius_c0_equals_c2": c0 == c2,
        "frobenius_c1_even": c1 % 2 == 0,
        "surviving_c0_odd": c0 % 2 == 1,
        "observed_c1_mod4": c1 % 4,
        "observed_nonsquare_weight_mod4": (c0 + c1 + c2) % 4,
        "elapsed_seconds": time.time() - t0,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=5)
    ap.add_argument("--stop", type=int, required=True)
    ap.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    primes = [
        p for p in range(max(5, args.start), args.stop + 1)
        if p % 12 == 5 and is_prime(p)
    ]
    t0 = time.time()
    if args.workers == 1:
        rows = list(map(record, primes))
    else:
        with ProcessPoolExecutor(
            max_workers=args.workers,
            mp_context=mp.get_context("fork"),
        ) as ex:
            rows = list(ex.map(record, primes, chunksize=1))
    failures = {
        "named_z_frobenius": [
            row["p"] for row in rows
            if row["named_z_frobenius_mismatches"]
        ],
        "endpoint_frobenius": [
            row["p"] for row in rows
            if row["endpoint_frobenius_mismatches"]
        ],
        "c0_equals_c2": [
            row["p"] for row in rows
            if not row["frobenius_c0_equals_c2"]
        ],
        "c1_even": [
            row["p"] for row in rows if not row["frobenius_c1_even"]
        ],
        "c0_odd": [
            row["p"] for row in rows if not row["surviving_c0_odd"]
        ],
        "content_not_one": [
            row["p"] for row in rows
            if row["normalized_content_value"] != 1
        ],
        "gamma_class_count_formula": [
            row["p"] for row in rows
            if not row["gamma_class_count_formula"]
        ],
    }
    out = {
        "range": [args.start, args.stop],
        "congruence_class": "p == 5 (mod 12)",
        "n_primes": len(rows),
        "failures": failures,
        "all_observed_c1_mod4_two": all(
            row["observed_c1_mod4"] == 2 for row in rows
        ),
        "all_observed_nonsquare_weight_mod4_zero": all(
            row["observed_nonsquare_weight_mod4"] == 0 for row in rows
        ),
        "rows": rows,
        "elapsed_seconds": time.time() - t0,
    }
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(
        f"p={args.start}..{args.stop}, p=5 mod12 primes={len(rows)} "
        f"failures={failures} seconds={out['elapsed_seconds']:.3f}",
        flush=True,
    )


if __name__ == "__main__":
    main()
