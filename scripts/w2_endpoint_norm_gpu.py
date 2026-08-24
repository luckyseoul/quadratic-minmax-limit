#!/usr/bin/env python3
"""GPU census of several normalized W2 endpoints without factorization.

For an endpoint vector ``w=c(D) gamma``, the sum of the square- and
nonsquare-orbit autocorrelations equals ``c(X)c(X^-1)`` at every
nonprincipal odd character.  Its gcd with

    g = 1 + X + ... + X^(H-1),
    H = oddpart((p^2-1)/2),

is the reciprocal *union* of factors missed in either orientation.  The
oriented cross-correlation with the exact Bose generator separately recovers
``c(X)``.  Hence

    A(c) = gcd(g, c(X), c(X^-1))

is the exact reciprocal Aut-bad polynomial, while the norm gcd is a stronger
unit-content diagnostic.  The gcd of several A(c)'s tests the endpoints
collectively for W2; the gcd of their norm polynomials tests the stronger
collective unit-ideal condition.

The orbit sequences are folded modulo H *before* autocorrelation.  This is
the same calculation in F2[X]/(X^H+1), but cuts the transform length by the
full 2-part of (p^2-1)/2.  All requested endpoints and both point-orbits are
then transformed as one GPU batch.  No irreducible factorization is used.

The output separately intersects each obstruction with the projective
factor R_M and scalar factor R_s(X^M), where

    M=(p+1)/2, s=oddpart(p-1), H=M*s.

This makes the scanner useful for testing endpoint-pair, longer-prefix, and
projective-versus-scalar hypotheses in one pass.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
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
from gf2x_ntl import available as ntl_available  # noqa: E402
from gf2x_ntl import gcd_bits as ntl_gcd_bits  # noqa: E402
from io_atomic import write_json_atomic  # noqa: E402
from w2_pole_fourier_fast import (  # noqa: E402
    named_z_without_conference,
)
from w2_translated_antipodal_norm_scan import (  # noqa: E402
    apply_pole,
    f2_gcd_bits,
    is_prime,
    pole_action_data,
    poly_bits,
)


POLYNOMIAL_GCD = f2_gcd_bits


def oddpart(value: int) -> int:
    while value % 2 == 0:
        value //= 2
    return value


def scalar_factor_bits(projective_order: int, scalar_order: int) -> int:
    """R_s(X^M)=sum_(j=0)^(s-1) X^(M*j) as packed F2 coefficients."""
    return sum(1 << (projective_order * j) for j in range(scalar_order))


def f2_multiply_bits(left: int, right: int) -> int:
    """Carryless multiplication of packed F2 polynomials."""
    result = 0
    while right:
        low_bit = right & -right
        result ^= left << (low_bit.bit_length() - 1)
        right ^= low_bit
    return result


def fold_sequence(values: np.ndarray, points: np.ndarray, order: int) -> np.ndarray:
    sequence = np.asarray(values[points], dtype=np.uint8)
    if sequence.size % order:
        raise AssertionError("orbit length is not a multiple of folded order")
    return np.bitwise_xor.reduce(sequence.reshape(-1, order), axis=0)


def spectral_certificates(
    batch: np.ndarray,
    gamma_batch: np.ndarray,
    backend: str,
    precision: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return oriented contents, endpoint norms, and the generator norm."""
    order = batch.shape[1]
    if backend == "cpu":
        spectrum = np.fft.rfft(batch.astype(np.float64), axis=1)
        gamma_spectrum = np.fft.rfft(
            gamma_batch.astype(np.float64), axis=1
        )
        oriented_spectrum = (
            spectrum[0::2] * np.conj(gamma_spectrum[0])
            + spectrum[1::2] * np.conj(gamma_spectrum[1])
        )
        norm_spectrum = (
            spectrum[0::2] * np.conj(spectrum[0::2])
            + spectrum[1::2] * np.conj(spectrum[1::2])
        )
        gamma_norm_spectrum = (
            gamma_spectrum[0] * np.conj(gamma_spectrum[0])
            + gamma_spectrum[1] * np.conj(gamma_spectrum[1])
        )
        transforms = [
            np.fft.irfft(oriented_spectrum, n=order, axis=1),
            np.fft.irfft(norm_spectrum, n=order, axis=1),
            np.fft.irfft(gamma_norm_spectrum, n=order),
        ]
        return tuple(
            (np.rint(values).astype(np.int64) & 1).astype(np.uint8)
            for values in transforms
        )

    import cupy as cp

    memory_pool = cp.get_default_memory_pool()
    if memory_pool.get_limit() == 0:
        memory_pool.set_limit(fraction=0.75)
    device_dtype = cp.float32 if precision == "single" else cp.float64
    device_batch = cp.asarray(batch, dtype=device_dtype)
    device_gamma = cp.asarray(gamma_batch, dtype=device_dtype)
    spectrum = cp.fft.rfft(device_batch, axis=1)
    gamma_spectrum = cp.fft.rfft(device_gamma, axis=1)
    oriented_spectrum = (
        spectrum[0::2] * cp.conj(gamma_spectrum[0])
        + spectrum[1::2] * cp.conj(gamma_spectrum[1])
    )
    norm_spectrum = (
        spectrum[0::2] * cp.conj(spectrum[0::2])
        + spectrum[1::2] * cp.conj(spectrum[1::2])
    )
    gamma_norm_spectrum = (
        gamma_spectrum[0] * cp.conj(gamma_spectrum[0])
        + gamma_spectrum[1] * cp.conj(gamma_spectrum[1])
    )
    device_transforms = [
        cp.fft.irfft(oriented_spectrum, n=order, axis=1),
        cp.fft.irfft(norm_spectrum, n=order, axis=1),
        cp.fft.irfft(gamma_norm_spectrum, n=order),
    ]
    result = tuple(
        cp.asnumpy(
            (cp.rint(values).astype(cp.int64) & 1).astype(cp.uint8)
        )
        for values in device_transforms
    )
    del (
        device_batch,
        device_gamma,
        spectrum,
        gamma_spectrum,
        oriented_spectrum,
        norm_spectrum,
        gamma_norm_spectrum,
        device_transforms,
    )
    return result


def reciprocal_bits(coefficients: np.ndarray) -> int:
    reciprocal = np.concatenate((coefficients[:1], coefficients[:0:-1]))
    return poly_bits(reciprocal)


def record(p: int, offsets: list[int], backend: str, precision: str) -> dict:
    if p % 12 != 5:
        raise ValueError("this endpoint reduction requires p == 5 (mod 12)")
    started = time.perf_counter()
    setup_started = started
    z, bits, q, mul, add, chi, lam = named_z_without_conference(p)
    q2, _mul2, _add2, _chi2, _frob, _norm, ia, ib = field_ctx(p)
    if q2 != q:
        raise AssertionError("field contexts disagree")

    omega = _primitive(mul, q)
    generator = mul(omega, omega)
    orbit_length = (q - 1) // 2
    ambient_order = oddpart(orbit_length)
    projective_order = (p + 1) // 2
    scalar_order = oddpart(p - 1)
    if ambient_order != projective_order * scalar_order:
        raise AssertionError("projective/scalar order decomposition failed")

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
    inv_fp = np.zeros(p, dtype=np.int64)
    legendre = np.zeros(p, dtype=np.int8)
    for value in range(1, p):
        inv_fp[value] = pow(value, p - 2, p)
        legendre[value] = (
            1 if pow(value, (p - 1) // 2, p) == 1 else -1
        )
    gamma, gamma_q, _gamma_mul, _b = named_gamma(p)
    if gamma_q != q:
        raise AssertionError("gamma field disagrees")
    folded_gamma = np.stack(
        (
            fold_sequence(gamma, square_points, ambient_order),
            fold_sequence(gamma, nonsquare_points, ambient_order),
        )
    )
    setup_seconds = time.perf_counter() - setup_started

    endpoint_started = time.perf_counter()
    folded_rows = []
    endpoint_metadata = []
    for offset in offsets:
        u = projective_order + offset
        if not projective_order <= u < p:
            raise ValueError(f"offset {offset} selects invalid endpoint u={u}")
        pole_t = lam * pow(u, p - 2, p) % p
        perm, switch = pole_action_data(
            p, pole_t, ia, ib, inv_fp, legendre
        )
        in_u, wfn = apply_pole(z, bits, perm, switch)
        if not in_u:
            raise AssertionError(f"endpoint missed U at p={p}, offset={offset}")
        folded_rows.append(fold_sequence(wfn, square_points, ambient_order))
        folded_rows.append(fold_sequence(wfn, nonsquare_points, ambient_order))
        endpoint_metadata.append({"offset": offset, "u": u, "pole_t": pole_t})
        del perm, switch, wfn
    endpoint_seconds = time.perf_counter() - endpoint_started

    transform_started = time.perf_counter()
    contents, norms, gamma_norm = spectral_certificates(
        np.stack(folded_rows), folded_gamma, backend, precision
    )
    transform_seconds = time.perf_counter() - transform_started

    gcd_started = time.perf_counter()
    ambient_g = (1 << ambient_order) - 1
    projective_factor = (1 << projective_order) - 1
    scalar_factor = scalar_factor_bits(projective_order, scalar_order)
    if POLYNOMIAL_GCD(projective_factor, scalar_factor) != 1:
        raise AssertionError("projective and scalar factors are not coprime")

    expected_gamma_norm = np.zeros(ambient_order, dtype=np.uint8)
    expected_gamma_norm[0] = 1
    gamma_norm_mismatches = int(
        np.count_nonzero(gamma_norm ^ expected_gamma_norm)
    )
    if gamma_norm_mismatches:
        raise AssertionError(f"generator norm mismatch at p={p}")

    aut_bad_polynomials = []
    norm_bad_polynomials = []
    for endpoint_index, metadata in enumerate(endpoint_metadata):
        content = contents[endpoint_index]
        content_bad = POLYNOMIAL_GCD(ambient_g, poly_bits(content))
        if content_bad == 1:
            reciprocal_bad = aut_bad = norm_bad = 1
        else:
            reciprocal_bad = POLYNOMIAL_GCD(
                ambient_g, reciprocal_bits(content)
            )
            aut_bad = POLYNOMIAL_GCD(content_bad, reciprocal_bad)
            norm_bad = POLYNOMIAL_GCD(
                ambient_g, poly_bits(norms[endpoint_index])
            )
        aut_bad_polynomials.append(aut_bad)
        norm_bad_polynomials.append(norm_bad)
        if aut_bad == 1:
            projective_bad = scalar_bad = 1
        else:
            projective_bad = POLYNOMIAL_GCD(aut_bad, projective_factor)
            scalar_bad = POLYNOMIAL_GCD(aut_bad, scalar_factor)
        if aut_bad != f2_multiply_bits(projective_bad, scalar_bad):
            raise AssertionError("bad polynomial did not split by layer")
        metadata.update(
            {
                "aut_bad_degree": aut_bad.bit_length() - 1,
                "norm_bad_degree": norm_bad.bit_length() - 1,
                "content_gcd_degree": content_bad.bit_length() - 1,
                "projective_bad_degree": projective_bad.bit_length() - 1,
                "scalar_bad_degree": scalar_bad.bit_length() - 1,
                "endpoint_w2_witness": aut_bad == 1,
                "unit_content_witness": norm_bad == 1,
                "aut_bad_hex": hex(aut_bad) if aut_bad.bit_length() <= 4097 else None,
            }
        )

    prefixes = []
    common = ambient_g
    common_norm = ambient_g
    for endpoint_count, (aut_bad, norm_bad) in enumerate(
        zip(aut_bad_polynomials, norm_bad_polynomials), 1
    ):
        if common != 1:
            common = POLYNOMIAL_GCD(common, aut_bad)
        if common_norm != 1:
            common_norm = POLYNOMIAL_GCD(common_norm, norm_bad)
        if common == 1:
            common_projective = common_scalar = 1
        else:
            common_projective = POLYNOMIAL_GCD(common, projective_factor)
            common_scalar = POLYNOMIAL_GCD(common, scalar_factor)
        prefixes.append(
            {
                "n_endpoints": endpoint_count,
                "offsets": offsets[:endpoint_count],
                "common_aut_bad_degree": common.bit_length() - 1,
                "common_projective_bad_degree": common_projective.bit_length() - 1,
                "common_scalar_bad_degree": common_scalar.bit_length() - 1,
                "collective_w2_witness": common == 1,
                "common_norm_bad_degree": common_norm.bit_length() - 1,
                "collective_unit_ideal_witness": common_norm == 1,
                "common_aut_bad_hex": (
                    hex(common) if common.bit_length() <= 4097 else None
                ),
            }
        )
    gcd_seconds = time.perf_counter() - gcd_started
    return {
        "p": p,
        "orbit_length": orbit_length,
        "ambient_oddpart": ambient_order,
        "fold_factor": orbit_length // ambient_order,
        "projective_order": projective_order,
        "scalar_order": scalar_order,
        "gamma_norm_mismatches": gamma_norm_mismatches,
        "endpoints": endpoint_metadata,
        "prefixes": prefixes,
        "timing_seconds": {
            "setup": setup_seconds,
            "endpoints_and_fold": endpoint_seconds,
            "batched_transform": transform_seconds,
            "polynomial_gcd": gcd_seconds,
            "total": time.perf_counter() - started,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=5)
    parser.add_argument("--stop", type=int, required=True)
    parser.add_argument("--offsets", default="0,1,2,3")
    parser.add_argument("--backend", choices=("gpu", "cpu"), default="gpu")
    parser.add_argument(
        "--precision", choices=("double", "single"), default="double"
    )
    parser.add_argument(
        "--gcd-backend", choices=("auto", "ntl", "python"), default="auto"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    global POLYNOMIAL_GCD
    use_ntl = args.gcd_backend != "python" and ntl_available()
    if args.gcd_backend == "ntl" and not use_ntl:
        raise RuntimeError("--gcd-backend=ntl requested but NTL bridge is unavailable")
    POLYNOMIAL_GCD = ntl_gcd_bits if use_ntl else f2_gcd_bits
    offsets = [int(value) for value in args.offsets.split(",")]
    if not offsets or offsets[0] != 0 or offsets != sorted(set(offsets)):
        raise ValueError("offsets must be distinct, increasing, and start at zero")
    primes = [
        p
        for p in range(max(5, args.start), args.stop + 1)
        if p % 12 == 5 and is_prime(p)
    ]
    started = time.perf_counter()
    rows = []
    for index, p in enumerate(primes, 1):
        row = record(p, offsets, args.backend, args.precision)
        rows.append(row)
        prefix_summary = ",".join(
            f"{item['n_endpoints']}:{item['common_aut_bad_degree']}"
            for item in row["prefixes"]
        )
        print(
            f"[{index}/{len(primes)}] p={p} H={row['ambient_oddpart']} "
            f"prefix_degrees={prefix_summary} "
            f"seconds={row['timing_seconds']['total']:.3f}",
            flush=True,
        )

    max_prefix = len(offsets)
    result = {
        "range": [args.start, args.stop],
        "congruence_class": "p == 5 (mod 12)",
        "offsets": offsets,
        "backend": args.backend,
        "precision": args.precision,
        "gcd_backend": "ntl" if use_ntl else "python",
        "n_primes": len(rows),
        "prefix_failures": {
            str(n): [
                row["p"]
                for row in rows
                if not row["prefixes"][n - 1]["collective_w2_witness"]
            ]
            for n in range(1, max_prefix + 1)
        },
        "unit_ideal_prefix_failures": {
            str(n): [
                row["p"]
                for row in rows
                if not row["prefixes"][n - 1]["collective_unit_ideal_witness"]
            ]
            for n in range(1, max_prefix + 1)
        },
        "scalar_prefix_failures": {
            str(n): [
                row["p"]
                for row in rows
                if row["prefixes"][n - 1]["common_scalar_bad_degree"]
            ]
            for n in range(1, max_prefix + 1)
        },
        "projective_prefix_failures": {
            str(n): [
                row["p"]
                for row in rows
                if row["prefixes"][n - 1]["common_projective_bad_degree"]
            ]
            for n in range(1, max_prefix + 1)
        },
        "rows": rows,
        "elapsed_seconds": time.perf_counter() - started,
    }
    write_json_atomic(args.output, result)
    print(
        f"done primes={len(rows)} prefix_failures={result['prefix_failures']} "
        f"elapsed={result['elapsed_seconds']:.3f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
