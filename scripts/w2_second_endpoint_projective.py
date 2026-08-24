#!/usr/bin/env python3
"""Projective-line parity theorem for the second normalized W2 endpoint.

For p=5 mod 12 and u=(p+3)/2=3/2 in F_p, every finite nonsquare
F_p-direction has odd endpoint-support parity except the two slopes

    d k^2 = 9/8.

They exist exactly for p=5 mod 24.  The pure-omega direction is always even.
Thus the projective quotient is either R_M+X^r, or

    R_M + X^r + X^(r+a) + X^(r-a),  M=(p+1)/2.

At nontrivial M-th roots the latter is a monomial times Phi_3(X^a).
Endpoint zero already clears all of R_M, so these possible second-endpoint
zeros cannot obstruct the two-endpoint W2 certificate.
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
from w2_pole_fourier_fast import (  # noqa: E402
    named_z_without_conference,
)
from w2_translated_antipodal_norm_scan import (  # noqa: E402
    apply_pole,
    is_prime,
    pole_action_data,
)


def f2_gcd_bits(left: int, right: int) -> int:
    while right:
        while left.bit_length() >= right.bit_length():
            left ^= right << (left.bit_length() - right.bit_length())
        left, right = right, left
    return left


def record(p: int) -> dict:
    if p % 12 != 5:
        raise ValueError("this endpoint reduction requires p == 5 (mod 12)")
    t0 = time.time()
    z, bits, q, mul, add, chi, lam = named_z_without_conference(p)
    q2, _mul2, _add2, _chi2, _frob, _norm, ia, d = field_ctx(p)
    if q2 != q or ia != 0:
        raise AssertionError("expected the omega^2=d coordinate field")

    projective_order = (p + 1) // 2
    u = projective_order + 1
    pole_t = lam * pow(u, p - 2, p) % p
    inv_fp = np.zeros(p, dtype=np.int64)
    legendre = np.zeros(p, dtype=np.int8)
    for value in range(1, p):
        inv_fp[value] = pow(value, p - 2, p)
        legendre[value] = (
            1 if pow(value, (p - 1) // 2, p) == 1 else -1
        )
    perm, switch = pole_action_data(
        p, pole_t, ia, d, inv_fp, legendre
    )
    in_u, wfn = apply_pole(z, bits, perm, switch)
    if not in_u:
        raise AssertionError(f"second endpoint missed U at p={p}")

    omega = _primitive(mul, q)
    generator = mul(omega, omega)
    orbit_length = (q - 1) // 2
    nonsquare = np.empty(orbit_length, dtype=np.int64)
    point = omega
    for j in range(orbit_length):
        nonsquare[j] = point
        point = mul(point, generator)
    if point != omega:
        raise AssertionError("nonsquare orbit did not close")

    point_a = nonsquare % p
    point_b = nonsquare // p
    directions = np.where(
        point_a == 0,
        p,
        point_b * inv_fp[point_a] % p,
    )
    projective_parities = np.bitwise_xor.reduce(
        wfn[nonsquare].reshape(p - 1, projective_order), axis=0
    )
    representative_directions = directions[:projective_order]
    if len(np.unique(representative_directions)) != projective_order:
        raise AssertionError("projective representatives are not distinct")

    pure_indices = np.flatnonzero(representative_directions == p)
    if len(pure_indices) != 1:
        raise AssertionError("pure-omega direction is not unique")
    pure_index = int(pure_indices[0])
    pure_parity = int(projective_parities[pure_index])

    finite_nonsquare_slopes = []
    for k in range(p):
        n = (1 - d * k * k) % p
        if pow(n, (p - 1) // 2, p) == p - 1:
            finite_nonsquare_slopes.append(k)
    observed_even_slopes = sorted(
        int(k)
        for k in finite_nonsquare_slopes
        if not int(
            projective_parities[
                int(np.flatnonzero(representative_directions == k)[0])
            ]
        )
    )

    rhs = 9 * pow(8, p - 2, p) % p
    predicted_even_slopes = sorted(
        k for k in finite_nonsquare_slopes if d * k * k % p == rhs
    )
    expected_n_exceptions = 2 if p % 24 == 5 else 0
    if len(predicted_even_slopes) != expected_n_exceptions:
        raise AssertionError(f"exceptional-slope count failed at p={p}")

    expected_parities = np.ones(projective_order, dtype=np.uint8)
    expected_parities[pure_index] = 0
    exceptional_indices = []
    for k in predicted_even_slopes:
        index = int(np.flatnonzero(representative_directions == k)[0])
        exceptional_indices.append(index)
        expected_parities[index] = 0
    quotient_mismatches = int(
        np.count_nonzero(projective_parities ^ expected_parities)
    )

    all_ones = (1 << projective_order) - 1
    quotient_bits = sum(
        int(value) << j for j, value in enumerate(projective_parities)
    )
    observed_bad = f2_gcd_bits(all_ones, quotient_bits)
    relative_exception = None
    exception_indices_symmetric = True
    trinomial_bad = 1
    if exceptional_indices:
        relative_exception = (exceptional_indices[0] - pure_index) % projective_order
        exception_indices_symmetric = (
            (exceptional_indices[1] - pure_index) % projective_order
            == (-relative_exception) % projective_order
        )
        trinomial = (
            1
            ^ (1 << relative_exception)
            ^ (1 << ((-relative_exception) % projective_order))
        )
        trinomial_bad = f2_gcd_bits(all_ones, trinomial)

    return {
        "p": p,
        "p_mod_24": p % 24,
        "endpoint_u": u,
        "pole_t": pole_t,
        "projective_order": projective_order,
        "n_finite_nonsquare_directions": len(finite_nonsquare_slopes),
        "pure_omega_index": pure_index,
        "pure_omega_parity": pure_parity,
        "observed_even_finite_slopes": observed_even_slopes,
        "predicted_even_finite_slopes": predicted_even_slopes,
        "exceptional_indices": exceptional_indices,
        "relative_exception_index": relative_exception,
        "exception_indices_symmetric": exception_indices_symmetric,
        "quotient_mismatches": quotient_mismatches,
        "observed_projective_bad_degree": observed_bad.bit_length() - 1,
        "predicted_projective_bad_degree": trinomial_bad.bit_length() - 1,
        "projective_bad_polynomial_match": observed_bad == trinomial_bad,
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
        ) as executor:
            rows = list(executor.map(record, primes, chunksize=1))
    failures = {
        "pure_omega_parity": [row["p"] for row in rows if row["pure_omega_parity"]],
        "finite_slope_formula": [
            row["p"] for row in rows
            if row["observed_even_finite_slopes"]
            != row["predicted_even_finite_slopes"]
        ],
        "exception_symmetry": [
            row["p"] for row in rows if not row["exception_indices_symmetric"]
        ],
        "projective_quotient": [
            row["p"] for row in rows if row["quotient_mismatches"]
        ],
        "projective_bad_polynomial": [
            row["p"] for row in rows
            if not row["projective_bad_polynomial_match"]
        ],
    }
    result = {
        "range": [args.start, args.stop],
        "congruence_class": "p == 5 (mod 12)",
        "n_primes": len(rows),
        "failures": failures,
        "rows": rows,
        "elapsed_seconds": time.time() - t0,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"primes={len(rows)} failures={failures} "
        f"elapsed={result['elapsed_seconds']:.3f}s"
    )


if __name__ == "__main__":
    main()
