#!/usr/bin/env python3
"""Test the normalized p=5 mod 12 endpoint on reciprocal Aut-orbits.

The projective-torus theorem clears every factor whose roots have order
dividing M=(p+1)/2.  For the remaining endpoint factors, Frobenius
invariance makes the bad-factor set invariant under X -> X^p.  Hence W2
can fail only if both a p-orbit and its reciprocal p-orbit are bad.

This scanner avoids factoring the enormous ambient repetition polynomial.
The exact two-component Fourier norm first isolates the reciprocal closure
of the bad factors.  Only that usually tiny gcd is factored, after which
direct Fourier evaluation determines which member(s) are actually bad.
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
import sympy as sp

ROOT = Path(
    os.environ.get("QML_REPO_ROOT", Path(__file__).resolve().parents[1])
)
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15604 import _primitive  # noqa: E402
from e1_gmin_m4_prop15613 import named_gamma  # noqa: E402
from e1_gmin_m4_prop15617 import _poly_gcd  # noqa: E402
from w2_pole_fourier_fast import (  # noqa: E402
    fourier_factor_mask,
    fourier_factor_residues,
    named_z_without_conference,
    switched_wfn,
)
from w2_translated_antipodal_norm_scan import (  # noqa: E402
    f2_gcd_bits,
    folded_autocorrelation,
    is_prime,
    poly_bits,
)


def oddpart(value: int) -> int:
    while value % 2 == 0:
        value //= 2
    return value


def bits_to_coeffs(value: int) -> list[int]:
    return [(value >> i) & 1 for i in range(value.bit_length())]


def factor_squarefree_bits(value: int) -> list[list[int]]:
    x = sp.symbols("x")
    poly = sp.Poly(
        sum(((value >> i) & 1) * x**i for i in range(value.bit_length())),
        x,
        modulus=2,
    )
    _unit, raw = sp.factor_list(poly, modulus=2)
    factors = []
    for fac, exponent in raw:
        if exponent != 1:
            raise AssertionError(f"candidate norm gcd is not square-free: {fac}")
        factors.append([int(v) & 1 for v in reversed(fac.all_coeffs())])
    factors.sort(key=lambda fac: (len(fac), fac))
    return factors


def substitute_power_mod(
    poly: list[int], power: int, modulus: list[int]
) -> list[int]:
    substituted = [0] * ((len(poly) - 1) * power + 1)
    for i, coeff in enumerate(poly):
        if coeff:
            substituted[i * power] = 1
    # Binary long division; modulus is monic.
    while len(substituted) >= len(modulus):
        if substituted[-1]:
            shift = len(substituted) - len(modulus)
            for j, coeff in enumerate(modulus):
                substituted[shift + j] ^= coeff
        while len(substituted) > 1 and not substituted[-1]:
            substituted.pop()
    return substituted


def candidate_aut_orbits(
    p: int,
    ambient_order: int,
    candidate: list[int],
    factors: list[list[int]],
) -> tuple[list[list[int]], list[int], list[int]]:
    """Aut orbits inside an invariant factor of R_ambient_order."""
    neighbours = {i: set() for i in range(len(factors))}
    p_map = []
    inverse_map = []
    for power, mapping in ((p, p_map), (ambient_order - 1, inverse_map)):
        for i, fac in enumerate(factors):
            transformed = substitute_power_mod(fac, power, candidate)
            hits = [
                j for j, other in enumerate(factors)
                if _poly_gcd(transformed, other) != [1]
            ]
            if len(hits) != 1:
                raise AssertionError(
                    f"candidate transform p={p}, i={i}, power={power}: {hits}"
                )
            mapping.append(hits[0])
            neighbours[i].add(hits[0])
            neighbours[hits[0]].add(i)

    orbits = []
    unseen = set(range(len(factors)))
    while unseen:
        todo = [min(unseen)]
        orbit = set()
        while todo:
            i = todo.pop()
            if i in orbit:
                continue
            orbit.add(i)
            todo.extend(neighbours[i] - orbit)
        unseen -= orbit
        orbits.append(sorted(orbit))
    return orbits, p_map, inverse_map


def record(p: int, repair_offsets: int = 8) -> dict:
    if p % 12 != 5:
        raise ValueError("this endpoint reduction requires p == 5 (mod 12)")
    t0 = time.time()
    z, bits, q, mul, add, chi, lam = named_z_without_conference(p)
    u = (p + 1) // 2
    pole_t = lam * pow(u, p - 2, p) % p
    in_u, wfn = switched_wfn(p, pole_t, z, bits, q, mul, add, chi)
    if not in_u:
        raise AssertionError(f"normalized endpoint missed U at p={p}")

    omega = _primitive(mul, q)
    generator = mul(omega, omega)
    ncoord = (q - 1) // 2
    ambient_order = oddpart(ncoord)
    projective_order = (p + 1) // 2
    p_minus_oddpart = oddpart(p - 1)
    if ambient_order != projective_order * p_minus_oddpart:
        raise AssertionError("odd-part decomposition failed")

    square_points = np.empty(ncoord, dtype=np.int64)
    point = 1
    for j in range(ncoord):
        square_points[j] = point
        point = mul(point, generator)
    nonsquare_points = np.fromiter(
        (mul(omega, int(point)) for point in square_points),
        dtype=np.int64,
        count=ncoord,
    )

    # gamma*gamma^-1=1+G+F_p^* makes this combined autocorrelation exactly
    # c(X)c(X^-1) at every nonprincipal odd-order character.
    norm = folded_autocorrelation(wfn, square_points, ambient_order)
    norm ^= folded_autocorrelation(wfn, nonsquare_points, ambient_order)
    ambient_g = (1 << ambient_order) - 1
    norm_gcd = f2_gcd_bits(poly_bits(norm), ambient_g)

    row = {
        "p": p,
        "u": u,
        "pole_t": pole_t,
        "ambient_oddpart": ambient_order,
        "projective_order": projective_order,
        "p_minus_oddpart": p_minus_oddpart,
        "norm_gcd_degree": norm_gcd.bit_length() - 1,
        "unit_content": norm_gcd == 1,
    }
    if norm_gcd == 1:
        row.update(
            {
                "candidate_factor_degrees": [],
                "candidate_aut_orbits": [],
                "bad_factor_indices": [],
                "bad_aut_orbits": [],
                "endpoint_w2_witness": True,
                "repair_rows": [],
                "first_repair_offset": 0,
                "two_endpoint_w2_witness": True,
                "elapsed_seconds": time.time() - t0,
            }
        )
        return row

    candidate = bits_to_coeffs(norm_gcd)
    factors = factor_squarefree_bits(norm_gcd)
    orbits, p_map, inverse_map = candidate_aut_orbits(
        p, ambient_order, candidate, factors
    )

    gamma, gamma_q, _gamma_mul, _b = named_gamma(p)
    if gamma_q != q:
        raise AssertionError("gamma field disagrees")
    gamma_square = fourier_factor_residues(
        gamma, square_points.tolist(), factors
    )
    gamma_nonsquare = fourier_factor_residues(
        gamma, nonsquare_points.tolist(), factors
    )
    if any(a == 0 and b == 0 for a, b in zip(gamma_square, gamma_nonsquare)):
        raise AssertionError("gamma vanished on both Fourier components")
    use_square = [value != 0 for value in gamma_square]
    bad_mask = fourier_factor_mask(
        wfn,
        square_points.tolist(),
        nonsquare_points.tolist(),
        factors,
        use_square,
    )
    bad = [i for i in range(len(factors)) if (bad_mask >> i) & 1]
    bad_orbits = [
        orbit for orbit in orbits
        if all((bad_mask >> i) & 1 for i in orbit)
    ]

    repair_rows = []
    first_repair_offset = None
    for offset in range(1, min(repair_offsets, (p - 1) // 2)):
        repair_u = u + offset
        repair_t = lam * pow(repair_u, p - 2, p) % p
        repair_in_u, repair_wfn = switched_wfn(
            p, repair_t, z, bits, q, mul, add, chi
        )
        if not repair_in_u:
            raise AssertionError(
                f"repair endpoint missed U at p={p}, offset={offset}"
            )
        repair_mask = fourier_factor_mask(
            repair_wfn,
            square_points.tolist(),
            nonsquare_points.tolist(),
            factors,
            use_square,
        )
        clears = all(
            any(not ((repair_mask >> i) & 1) for i in orbit)
            for orbit in bad_orbits
        )
        repair_rows.append(
            {
                "offset": offset,
                "u": repair_u,
                "pole_t": repair_t,
                "bad_factor_indices": [
                    i for i in range(len(factors)) if (repair_mask >> i) & 1
                ],
                "clears_endpoint_bad_orbits": clears,
            }
        )
        if clears and first_repair_offset is None:
            first_repair_offset = offset
            break

    # Frobenius invariance is the reason each p-orbit is all bad or all
    # good.  Keep it as an executable assertion, not an assumed symmetry.
    if any(((bad_mask >> i) & 1) != ((bad_mask >> p_map[i]) & 1)
           for i in range(len(factors))):
        raise AssertionError(f"bad-factor mask is not p-invariant at p={p}")
    if any(inverse_map[inverse_map[i]] != i for i in range(len(factors))):
        raise AssertionError(f"candidate reciprocal map is not involutive at p={p}")

    row.update(
        {
            "candidate_factor_degrees": [len(fac) - 1 for fac in factors],
            "candidate_aut_orbits": orbits,
            "candidate_p_map": p_map,
            "candidate_inverse_map": inverse_map,
            "bad_factor_indices": bad,
            "bad_aut_orbits": bad_orbits,
            "endpoint_w2_witness": not bad_orbits,
            "repair_rows": repair_rows,
            "first_repair_offset": first_repair_offset,
            "two_endpoint_w2_witness": not bad_orbits or first_repair_offset is not None,
            "elapsed_seconds": time.time() - t0,
        }
    )
    return row


def record_job(job: tuple[int, int]) -> dict:
    return record(*job)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", type=int, default=5)
    ap.add_argument("--stop", type=int, required=True)
    ap.add_argument("--workers", type=int, default=min(8, os.cpu_count() or 1))
    ap.add_argument("--repair-offsets", type=int, default=8)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    primes = [
        p for p in range(max(5, args.start), args.stop + 1)
        if p % 12 == 5 and is_prime(p)
    ]
    t0 = time.time()
    jobs = [(p, args.repair_offsets) for p in primes]
    if args.workers == 1:
        rows = list(map(record_job, jobs))
    else:
        with ProcessPoolExecutor(
            max_workers=args.workers,
            mp_context=mp.get_context("fork"),
        ) as ex:
            rows = list(ex.map(record_job, jobs, chunksize=1))
    result = {
        "range": [args.start, args.stop],
        "congruence_class": "p == 5 (mod 12)",
        "n_primes": len(rows),
        "n_unit_endpoints": sum(row["unit_content"] for row in rows),
        "n_endpoint_w2_witnesses": sum(
            row["endpoint_w2_witness"] for row in rows
        ),
        "n_two_endpoint_w2_witnesses": sum(
            row["two_endpoint_w2_witness"] for row in rows
        ),
        "failures": [row["p"] for row in rows if not row["endpoint_w2_witness"]],
        "repair_failures": [
            row["p"] for row in rows if not row["two_endpoint_w2_witness"]
        ],
        "rows": rows,
        "elapsed_seconds": time.time() - t0,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        f"primes={len(rows)} units={result['n_unit_endpoints']} "
        f"W2={result['n_endpoint_w2_witnesses']} "
        f"two={result['n_two_endpoint_w2_witnesses']} "
        f"failures={result['failures']} "
        f"repair_failures={result['repair_failures']} "
        f"elapsed={result['elapsed_seconds']:.3f}s"
    )


if __name__ == "__main__":
    main()
