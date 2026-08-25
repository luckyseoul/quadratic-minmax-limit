#!/usr/bin/env python3
"""Finite audit of the two-root half-conic obstruction in Prop. 15.638."""
from __future__ import annotations

import argparse
import json
from itertools import product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def legendre_table(p: int) -> tuple[int, ...]:
    table = [0] * p
    for x in range(1, p):
        table[x] = 1 if pow(x, (p - 1) // 2, p) == 1 else -1
    return tuple(table)


def first_nonsquare(p: int, eta: tuple[int, ...]) -> int:
    return next(x for x in range(2, p) if eta[x] == -1)


def selected_half(p: int, eta: tuple[int, ...]) -> tuple[int | None, ...]:
    nu = first_nonsquare(p, eta)
    affine = tuple(x for x in range(p) if eta[(x * x - nu) % p] == 1)
    # N(1,0)=1, so infinity belongs to the eta(N)=+1 half.
    return affine + (None,)


def quadratic_value(
    p: int, coeffs: tuple[int, int, int], point: int | None
) -> int:
    A, B, C = coeffs
    if point is None:
        return A % p
    return (A * point * point + B * point + C) % p


def audit_prime(p: int) -> dict:
    eta = legendre_table(p)
    selected = selected_half(p, eta)
    forbidden = []
    checked = 0
    for coeffs in product(range(p), repeat=3):
        if coeffs == (0, 0, 0):
            continue
        checked += 1
        signs = tuple(
            eta[quadratic_value(p, coeffs, point)] for point in selected
        )
        if signs.count(0) == 2 and all(sign >= 0 for sign in signs):
            forbidden.append(coeffs)
    return {
        "p": p,
        "selected_directions": len(selected),
        "quadratic_forms_checked": checked,
        "forbidden_two_root_nonnegative_forms": len(forbidden),
        "examples": forbidden[:5],
        "checks": len(selected) == (p + 1) // 2 and not forbidden,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--primes",
        nargs="+",
        type=int,
        default=[11, 13, 17, 19, 23, 29, 31, 37, 41, 43],
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "evidence" / "r1_next_shell_half_conic_11_43.json",
    )
    args = parser.parse_args()
    rows = [audit_prime(p) for p in args.primes]
    out = {
        "role": "finite audit only; the Hasse character-sum argument is proof",
        "all_checks": all(row["checks"] for row in rows),
        "rows": rows,
    }
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(f"half-conic audit: {out['all_checks']}")
    print(f"  wrote {args.output}")


if __name__ == "__main__":
    main()
