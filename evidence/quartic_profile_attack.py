#!/usr/bin/env python3
"""Diagnostics for the exceptional quartic statistic in profile coordinates.

This is a lab script, not a numbered proposition and not a general floor proof.

For p == 3 (mod 4), a quartic character psi of F_{p^2} is trivial on
F_p^*.  If sigma_L(s) is the line-sum profile of y on the affine lines
parallel to L and eps=y_inf, then direct pair counting gives

    Z_psi(y) = sum_L psi(g_L) a_L(y),
    a_L(y)   = (1/4) sum_s (sigma_L(s)-eps)^2.

The a_L are nonnegative integers.  Orthogonality of the mean-zero ridge
functions also gives the pointwise conservation law

    sum_L a_L(y) = p(p^2-1)/4.

The script verifies these identities against direct difference-pair sums and
reports the quartic variance for persisted Max+ arrays.  For p == 1 (mod 4),
psi|F_p^* is quadratic, so the line-energy identity is deliberately skipped;
the direct quartic variance is still reported.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path

import numpy as np


def field_context(p: int):
    q = p * p

    def is_irreducible(a: int, b: int) -> bool:
        return all((x * x - a * x - b) % p != 0 for x in range(p))

    irr_a = irr_b = None
    for a in range(p):
        for b in range(p):
            if is_irreducible(a, b):
                irr_a, irr_b = a, b
                break
        if irr_a is not None:
            break
    if irr_a is None or irr_b is None:
        raise RuntimeError("no irreducible quadratic found")

    def mul(u: int, v: int) -> int:
        c0, c1 = u % p, u // p
        d0, d1 = v % p, v // p
        return (c0 * d0 + c1 * d1 * irr_b) % p + (
            (c0 * d1 + c1 * d0 + c1 * d1 * irr_a) % p
        ) * p

    def trace(x: int) -> int:
        return (2 * (x % p) + irr_a * (x // p)) % p

    def add(u: int, v: int) -> int:
        return (u % p + v % p) % p + ((u // p + v // p) % p) * p

    return q, mul, trace, add


def primitive_root(q: int, mul) -> int:
    for g in range(2, q):
        x = 1
        seen = set()
        for _ in range(q - 1):
            seen.add(x)
            x = mul(x, g)
        if len(seen) == q - 1 and x == 1:
            return g
    raise RuntimeError("no primitive root")


def quartic_character(p: int, mul) -> np.ndarray:
    q = p * p
    psi = np.zeros(q, dtype=np.complex128)
    generator = primitive_root(q, mul)
    x = 1
    for exponent in range(q - 1):
        psi[x] = (1j) ** (exponent % 4)
        x = mul(x, generator)
    return psi


def quartic_kernel(p: int, add, psi: np.ndarray) -> np.ndarray:
    q = p * p
    neg = np.array([((-x) % p) + ((-(x // p)) % p) * p for x in range(q)])
    kernel = np.empty((q, q), dtype=np.complex128)
    for a in range(q):
        for b in range(q):
            kernel[a, b] = psi[add(a, int(neg[b]))]
    np.fill_diagonal(kernel, 0)
    return kernel


def projective_directions(p: int, mul, trace):
    q = p * p
    seen: set[int] = set()
    out = []
    for g in range(1, q):
        if g in seen:
            continue
        line = [mul(t, g) for t in range(1, p)]
        seen.update(line)
        annihilator = next(c for c in range(1, q) if trace(mul(c, g)) == 0)
        t_of = np.array([trace(mul(annihilator, x)) for x in range(q)], dtype=np.int16)
        out.append((g, t_of))
    if len(out) != p + 1:
        raise RuntimeError(f"expected {p + 1} directions, got {len(out)}")
    return out


def finite_and_eps(rows: np.ndarray, q: int, finite_only: bool):
    if finite_only:
        if rows.shape[1] != q:
            raise ValueError(f"finite-only rows must have {q} columns")
        return rows, np.ones(len(rows), dtype=np.int8)
    if rows.shape[1] != q + 1:
        raise ValueError(f"full rows must have {q + 1} columns")
    return rows[:, 1:], rows[:, 0].astype(np.int8)


def direct_values(y_finite: np.ndarray, kernel: np.ndarray, chunk: int) -> np.ndarray:
    out = np.empty(len(y_finite), dtype=np.complex128)
    for lo in range(0, len(y_finite), chunk):
        hi = min(len(y_finite), lo + chunk)
        d = ((1 - np.asarray(y_finite[lo:hi], dtype=np.int8)) // 2).astype(np.float64)
        out[lo:hi] = np.einsum("bi,ij,bj->b", d, kernel, d, optimize=True)
    return out


def line_energy_values(
    p: int,
    y_finite: np.ndarray,
    eps: np.ndarray,
    directions,
    psi: np.ndarray,
    chunk: int,
):
    total = np.empty(len(y_finite), dtype=np.int64)
    signed = np.empty(len(y_finite), dtype=np.complex128)
    for lo in range(0, len(y_finite), chunk):
        hi = min(len(y_finite), lo + chunk)
        y = np.asarray(y_finite[lo:hi], dtype=np.int16)
        e = eps[lo:hi].astype(np.int16)
        a = np.empty((len(y), len(directions)), dtype=np.int64)
        weights = np.empty(len(directions), dtype=np.complex128)
        for j, (g, t_of) in enumerate(directions):
            sigma = np.empty((len(y), p), dtype=np.int16)
            for s in range(p):
                sigma[:, s] = y[:, t_of == s].sum(axis=1)
            delta = sigma - e[:, None]
            sq = np.sum(delta.astype(np.int64) ** 2, axis=1)
            if np.any(sq % 4):
                raise RuntimeError("profile energy is not integral")
            a[:, j] = sq // 4
            weights[j] = psi[g]
        total[lo:hi] = a.sum(axis=1)
        signed[lo:hi] = a @ weights
    return signed, total


def moment_record(values: np.ndarray, p: int) -> dict:
    q = p * p
    abs2 = np.abs(values) ** 2
    threshold = 3 * q * (q - 1) / 16
    mean = complex(values.mean())
    variance = float(abs2.mean())
    return {
        "count": int(len(values)),
        "mean": [float(mean.real), float(mean.imag)],
        "E_abs_Zpsi_sq": variance,
        "threshold": float(threshold),
        "ratio_to_threshold": variance / threshold,
        "lambda_exceptional": 32 * variance / (q * (q - 1)),
        "min_abs": float(np.min(np.abs(values))),
        "max_abs": float(np.max(np.abs(values))),
    }


def spherical_benchmark_record(p: int) -> dict:
    """Exact radius-sqrt(q+1) spherical comparison in dim (q+1)/2."""
    q = p * p
    sphere = Fraction(q * (q - 1) * (q + 1), 4 * (q + 5))
    threshold = Fraction(3 * q * (q - 1), 16)
    gap = sphere - threshold
    return {
        "E_abs_Zpsi_sq": str(sphere),
        "QVAR_threshold": str(threshold),
        "sphere_minus_threshold": str(gap),
        "sphere_clears_QVAR": gap > 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("p", type=int)
    parser.add_argument("array", type=Path)
    parser.add_argument("--finite-only", action="store_true")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--chunk", type=int, default=20_000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    rows = np.load(args.array, mmap_mode="r")
    if args.limit is not None:
        rows = rows[: args.limit]
    q, mul, trace, add = field_context(args.p)
    y_finite, eps = finite_and_eps(rows, q, args.finite_only)
    psi = quartic_character(args.p, mul)
    kernel = quartic_kernel(args.p, add, psi)
    direct = direct_values(y_finite, kernel, args.chunk)

    report = {
        "p": args.p,
        "array": str(args.array),
        "finite_only": args.finite_only,
        "direct": moment_record(direct, args.p),
        "spherical_benchmark": spherical_benchmark_record(args.p),
    }
    report["direct"]["excess_over_spherical_benchmark"] = (
        report["direct"]["E_abs_Zpsi_sq"]
        - float(Fraction(report["spherical_benchmark"]["E_abs_Zpsi_sq"]))
    )
    if args.p % 4 == 3:
        directions = projective_directions(args.p, mul, trace)
        via_profiles, total_energy = line_energy_values(
            args.p, y_finite, eps, directions, psi, args.chunk
        )
        target_total = args.p * (q - 1) // 4
        report["profile_identity"] = {
            "max_abs_error": float(np.max(np.abs(via_profiles - direct))),
            "target_total_energy": target_total,
            "total_energy_min": int(total_energy.min()),
            "total_energy_max": int(total_energy.max()),
            "conservation_holds": bool(np.all(total_energy == target_total)),
            "moment": moment_record(via_profiles, args.p),
        }

    rendered = json.dumps(report, indent=2)
    print(rendered)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
