#!/usr/bin/env python3
"""Factorization-free oriented-content test for two normalized W2 endpoints.

For w=c(D)gamma, split the multiplicative Fourier transform over the square
and nonsquare point-orbits.  The exact Bose generator norm gives

  Gamma_sq(X) Gamma_sq(X^-1) + Gamma_ns(X) Gamma_ns(X^-1) = 1 mod g.

Therefore the two cross-correlations reconstruct c itself modulo g, not only
its reciprocal norm.  Frobenius invariance makes the bad factors p-invariant;
gcd(g,c,c*) is consequently the exact product of complete reciprocal Aut
orbits missed by an endpoint.  Two endpoints certify W2 iff their two such
gcds are coprime.
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

from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import named_gamma  # noqa: E402
from w2_pole_fourier_fast import (  # noqa: E402
    named_z_without_conference,
    switched_wfn,
)
from w2_translated_antipodal_norm_scan import is_prime  # noqa: E402


def oddpart(value: int) -> int:
    while value % 2 == 0:
        value //= 2
    return value


def f2_gcd_bits(left: int, right: int) -> int:
    while right:
        while left.bit_length() >= right.bit_length():
            left ^= right << (left.bit_length() - right.bit_length())
        left, right = right, left
    return left


def poly_bits(coeffs: np.ndarray) -> int:
    return sum(int(value) << i for i, value in enumerate(coeffs))


def reciprocal_bits(coeffs: np.ndarray) -> int:
    """C(X^-1) modulo X^H+1, where H=len(coeffs)."""
    reciprocal = np.concatenate((coeffs[:1], coeffs[:0:-1]))
    return poly_bits(reciprocal)


def folded_crosscorrelation(
    values: np.ndarray,
    points: np.ndarray,
    gamma_spectrum: np.ndarray,
    orbit_length: int,
    odd_length: int,
) -> np.ndarray:
    sequence = np.asarray(values[points], dtype=np.float64)
    spectrum = np.fft.rfft(sequence)
    correlation = np.rint(
        np.fft.irfft(spectrum * np.conj(gamma_spectrum), n=orbit_length)
    ).astype(np.int64) & 1
    return np.bitwise_xor.reduce(
        correlation.reshape(-1, odd_length), axis=0
    ).astype(np.uint8)


def record(p: int) -> dict:
    if p % 12 != 5:
        raise ValueError("this endpoint reduction requires p == 5 (mod 12)")
    t0 = time.time()
    z, bits, q, mul, add, chi, lam = named_z_without_conference(p)
    gamma, gamma_q, _gamma_mul, _b = named_gamma(p)
    if gamma_q != q:
        raise AssertionError("gamma field disagrees")

    omega = _primitive(mul, q)
    generator = mul(omega, omega)
    orbit_length = (q - 1) // 2
    ambient_order = oddpart(orbit_length)
    projective_order = (p + 1) // 2
    p_minus_oddpart = oddpart(p - 1)
    if ambient_order != projective_order * p_minus_oddpart:
        raise AssertionError("odd-part decomposition failed")

    square_points = np.empty(orbit_length, dtype=np.int64)
    point = 1
    for j in range(orbit_length):
        square_points[j] = point
        point = mul(point, generator)
    if point != 1:
        raise AssertionError("square orbit did not close")
    nonsquare_points = np.fromiter(
        (mul(omega, int(point)) for point in square_points),
        dtype=np.int64,
        count=orbit_length,
    )

    gamma_square = np.asarray(gamma[square_points], dtype=np.float64)
    gamma_nonsquare = np.asarray(gamma[nonsquare_points], dtype=np.float64)
    gamma_square_spectrum = np.fft.rfft(gamma_square)
    gamma_nonsquare_spectrum = np.fft.rfft(gamma_nonsquare)

    # Executable check of the dual-vector identity used to reconstruct c.
    gamma_norm = np.rint(
        np.fft.irfft(
            gamma_square_spectrum * np.conj(gamma_square_spectrum),
            n=orbit_length,
        )
    ).astype(np.int64) & 1
    gamma_norm ^= np.rint(
        np.fft.irfft(
            gamma_nonsquare_spectrum * np.conj(gamma_nonsquare_spectrum),
            n=orbit_length,
        )
    ).astype(np.int64) & 1
    gamma_norm = np.bitwise_xor.reduce(
        gamma_norm.reshape(-1, ambient_order), axis=0
    ).astype(np.uint8)
    expected_norm = np.zeros(ambient_order, dtype=np.uint8)
    expected_norm[0] = 1
    gamma_norm_mismatches = int(np.count_nonzero(gamma_norm ^ expected_norm))
    if gamma_norm_mismatches:
        raise AssertionError(f"generator norm mismatch at p={p}")

    ambient_g = (1 << ambient_order) - 1
    endpoints = []
    aut_bad_polynomials = []
    for offset in (0, 1):
        u = projective_order + offset
        pole_t = lam * pow(u, p - 2, p) % p
        in_u, wfn = switched_wfn(p, pole_t, z, bits, q, mul, add, chi)
        if not in_u:
            raise AssertionError(f"endpoint missed U at p={p}, offset={offset}")

        content = folded_crosscorrelation(
            wfn,
            square_points,
            gamma_square_spectrum,
            orbit_length,
            ambient_order,
        )
        content ^= folded_crosscorrelation(
            wfn,
            nonsquare_points,
            gamma_nonsquare_spectrum,
            orbit_length,
            ambient_order,
        )
        content_bits = poly_bits(content)
        content_gcd = f2_gcd_bits(ambient_g, content_bits)
        aut_bad = f2_gcd_bits(content_gcd, reciprocal_bits(content))
        aut_bad_polynomials.append(aut_bad)
        aut_bad_degree = aut_bad.bit_length() - 1
        endpoints.append(
            {
                "offset": offset,
                "u": u,
                "pole_t": pole_t,
                "content_weight_mod_g": int(np.count_nonzero(content)),
                "content_gcd_degree": content_gcd.bit_length() - 1,
                "aut_bad_degree": aut_bad_degree,
                "aut_bad_hex": hex(aut_bad) if aut_bad_degree <= 4096 else None,
                "endpoint_w2_witness": aut_bad == 1,
            }
        )

    common_aut_bad = f2_gcd_bits(
        aut_bad_polynomials[0], aut_bad_polynomials[1]
    )
    common_degree = common_aut_bad.bit_length() - 1
    return {
        "p": p,
        "ambient_oddpart": ambient_order,
        "projective_order": projective_order,
        "p_minus_oddpart": p_minus_oddpart,
        "gamma_norm_mismatches": gamma_norm_mismatches,
        "endpoints": endpoints,
        "common_aut_bad_degree": common_degree,
        "common_aut_bad_hex": (
            hex(common_aut_bad) if common_degree <= 4096 else None
        ),
        "two_endpoint_w2_witness": common_aut_bad == 1,
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
    result = {
        "range": [args.start, args.stop],
        "congruence_class": "p == 5 (mod 12)",
        "n_primes": len(rows),
        "n_first_endpoint_w2": sum(
            row["endpoints"][0]["endpoint_w2_witness"] for row in rows
        ),
        "n_second_endpoint_w2": sum(
            row["endpoints"][1]["endpoint_w2_witness"] for row in rows
        ),
        "n_two_endpoint_w2": sum(
            row["two_endpoint_w2_witness"] for row in rows
        ),
        "first_endpoint_failures": [
            row["p"] for row in rows
            if not row["endpoints"][0]["endpoint_w2_witness"]
        ],
        "second_endpoint_failures": [
            row["p"] for row in rows
            if not row["endpoints"][1]["endpoint_w2_witness"]
        ],
        "two_endpoint_failures": [
            row["p"] for row in rows if not row["two_endpoint_w2_witness"]
        ],
        "rows": rows,
        "elapsed_seconds": time.time() - t0,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"primes={len(rows)} first={result['n_first_endpoint_w2']} "
        f"second={result['n_second_endpoint_w2']} "
        f"pair={result['n_two_endpoint_w2']} "
        f"pair_failures={result['two_endpoint_failures']} "
        f"elapsed={result['elapsed_seconds']:.3f}s"
    )


if __name__ == "__main__":
    main()
