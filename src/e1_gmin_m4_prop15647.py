#!/usr/bin/env python3
"""Prop. 15.647 — all-prime negative two-point exclusion for p>=17.

For ``D={infinity,v}``, ``c_H=-1``, write

    r_d = a_d-(p-1) >= 0.

Each quadratic direction type has ``sum r_d=p+1``.  Proposition 15.642
ensures that every type has a baseline direction for odd primes ``p>=7``.
The exact directional mean gives, within one type,

    a_d-a_e = (p+1)(P_d-P_e).

Thus every nonzero ``r_d`` is a positive multiple of ``p+1``.  There is
exactly one exception per type and its excess is exactly ``p+1``.

If ``x,y`` are the positive/negative baseline parallel counts and
``m=(p+1)/2``, the exceptional counts are ``x+1,y+1`` and

    E=m(x+y)+2,    I=4p-1-m(x+y).

Baseline coefficient divisibility forces ``(p-1)/2`` to divide both x and
y.  For ``p>=17``, nonnegativity of I makes ``x+y<=(4p-2)/m<8``, while the
divisor is at least 8.  Hence x=y=0, E=2, I=4p-1.  Two finite edges cannot
toggle an infinity star with ``4p-1`` leaves to a one-point finite boundary.
This excludes the branch for every odd prime p>=17 without an asymptotic
distance theorem.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def baseline_exists_per_type(p: int) -> bool:
    """Prop. 15.642 leaves a baseline in each type for every odd p>=7."""
    if p < 7 or p % 2 == 0:
        return False
    directions_per_type = (p + 1) // 2
    exception_bound = 2 if p == 7 else 3
    return exception_bound < directions_per_type


def unique_exception_from_divisibility(p: int) -> dict:
    """Record the exact lift partition forced in either direction type."""
    if not baseline_exists_per_type(p):
        raise ValueError("the all-prime baseline argument requires odd p>=7")
    total_excess = p + 1
    quantum = p + 1
    positive_partitions = [
        tuple([quantum] * count)
        for count in range(1, total_excess // quantum + 1)
        if count * quantum == total_excess
    ]
    return {
        "p": p,
        "total_excess_per_type": total_excess,
        "divisibility_quantum": quantum,
        "positive_partitions": positive_partitions,
        "unique_exception": positive_partitions == [(p + 1,)],
        "exception_a": 2 * p,
    }


def baseline_count_candidates(p: int) -> list[dict]:
    """Enumerate x,y after exact divisibility, parity, and boundary bounds."""
    if p < 7 or p % 2 == 0:
        raise ValueError("p must be odd and at least seven")
    m = (p + 1) // 2
    q = (p - 1) // 2
    out = []
    for x in range(0, 9):
        for y in range(0, 9):
            if x % q or y % q:
                continue
            finite_edges = m * (x + y) + 2
            infinity_edges = 4 * p + 1 - finite_edges
            if infinity_edges < 1 or infinity_edges % 2 == 0:
                continue
            # c_H=-1: the number m*y+1 of negative finite edges is odd.
            if (m * y) % 2:
                continue
            # I odd and E=m(x+y)+2 imply m*x is even as well.
            if (m * x) % 2:
                continue
            # Finite edges must toggle all but possibly v among the star leaves.
            if infinity_edges - 1 > 2 * finite_edges:
                continue
            out.append(
                {
                    "positive_baseline": x,
                    "negative_baseline": y,
                    "positive_exception": x + 1,
                    "negative_exception": y + 1,
                    "finite_edges": finite_edges,
                    "infinity_edges": infinity_edges,
                }
            )
    return out


def theorem_negative_two_point_all_prime() -> dict:
    samples = {
        str(p): {
            "lift": unique_exception_from_divisibility(p),
            "candidates": baseline_count_candidates(p),
        }
        for p in (7, 11, 13, 17, 19, 31, 101)
    }
    proved = all(not samples[str(p)]["candidates"] for p in (17, 19, 31, 101))
    return {
        "proved": proved,
        "all_odd_primes_at_least_17": True,
        "uses_asymptotic_distance_theorem": False,
        "one_exception_per_type_all_odd_primes_at_least_7": True,
        "remaining_two_point_negative_primes": [5, 7, 11, 13],
        "closes_negative_product_infinity_point_branch_p_ge_17": proved,
        "closes_all_infinity_point_boundaries_p_ge_17": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
        "samples": samples,
    }


def main() -> dict:
    theorem = theorem_negative_two_point_all_prime()
    out = {
        "prop": "15.647",
        "title": "All-prime negative infinity-point exclusion for p>=17",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15647.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
