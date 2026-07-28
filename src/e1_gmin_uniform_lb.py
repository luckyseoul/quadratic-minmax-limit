#!/usr/bin/env python3
"""
Uniform lower bound candidates for g_min vs bi-tight threshold.

Target: prove g_min > -(p-2)/(p(2p-1)) for all primes p>=5.
Candidate LB: g_min >= -(p-2)/(2 p^2)  (strictly above threshold for p>2).
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def thresh(p: int) -> Fraction:
    return Fraction(-(p - 2), p * (2 * p - 1))


def candidate_lb(p: int) -> Fraction:
    return Fraction(-(p - 2), 2 * p * p)


def main() -> None:
    known = {
        3: Fraction(-1, 3),
        5: Fraction(-3, 65),
        7: Fraction(-109, 2863),
    }
    rows = []
    for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]:
        t = thresh(p)
        lb = candidate_lb(p)
        row = {
            "p": p,
            "threshold": float(t),
            "candidate_lb": float(lb),
            "lb_gt_threshold": lb > t,
            "known_gmin": float(known[p]) if p in known else None,
            "known_ge_lb": (known[p] >= lb) if p in known else None,
            "known_gt_threshold": (known[p] > t) if p in known else None,
        }
        rows.append(row)
        print(
            f"p={p}: lb={float(lb):.6f} thresh={float(t):.6f} "
            f"lb>thresh={lb > t} known={row['known_gmin']}"
        )

    # Algebra: lb > thresh iff (p-2)/(2p^2) < (p-2)/(p(2p-1)) for p>2
    # iff 2p^2 > p(2p-1) iff 2p > 2p-1 iff 0 > -1. Always true for p>2.
    algebra_ok = all(
        candidate_lb(p) > thresh(p) for p in range(3, 200) if p % 2 == 1
    )
    out = {
        "candidate": "g_min >= -(p-2)/(2 p^2)",
        "algebra_lb_beats_threshold_all_odd_p_ge_3": algebra_ok,
        "rows": rows,
        "status": (
            "Candidate LB beats bi-tight threshold for all odd p>2 by algebra. "
            "Holds for certified g_min at p=5,7; FAILS to be a valid LB at p=3 "
            "(as required: bi-tight exists at p=3). "
            "OPEN: prove g_min >= -(p-2)/(2p^2) for all primes p>=5."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_uniform_lb.json"
    path.write_text(json.dumps(out, indent=2))
    print("wrote", path)
    print(out["status"])


if __name__ == "__main__":
    main()
