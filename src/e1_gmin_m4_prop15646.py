#!/usr/bin/env python3
"""Prop. 15.646 — the negative two-point normal form is impossible.

Inside Proposition 15.644's sufficiently-large-prime normal form, the two
exceptional parallel counts are ``(U,V)=(3,1)`` or ``(1,3)`` according to
positive/negative quadratic direction type.  Every baseline direction has
two parallel edges.

The additive inter-fibre identity in a baseline direction is

    K_st = -eps_d (a_s+a_t),     sum_s a_s=0.

Hence its signed transverse-edge sum is ``sum K_st=0``.  On the other hand,
the total signed sum of all finite edges is ``U-V`` because the two edges in
each positive and negative baseline direction cancel globally.  Removing
the two parallel edges from a baseline direction of type ``eps`` leaves
signed transverse sum ``U-V-2*eps``.  This equals 4 in a negative baseline
when ``(U,V)=(3,1)``, and -4 in a positive baseline when ``(U,V)=(1,3)``.
Both contradict zero.  Thus the branch is empty for all sufficiently large
odd primes.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def additive_matrix_total(deviations: tuple[int, ...], eps: int) -> int:
    """Return ``sum_{s<t} -eps*(a_s+a_t)`` exactly."""
    if eps not in (-1, 1):
        raise ValueError("eps must be +1 or -1")
    if sum(deviations) != 0:
        raise ValueError("deviations must sum to zero")
    return sum(
        -eps * (deviations[s] + deviations[t])
        for s in range(len(deviations))
        for t in range(s + 1, len(deviations))
    )


def normal_form_signed_audit(
    p: int, positive_exception: int, negative_exception: int
) -> dict:
    """Audit one exceptional-count split from Proposition 15.644."""
    if p < 3 or p % 2 == 0:
        raise ValueError("p must be odd and at least three")
    if (positive_exception, negative_exception) not in ((3, 1), (1, 3)):
        raise ValueError("normal-form exceptional counts must be (3,1) or (1,3)")
    baseline_per_type = (p - 1) // 2
    total_finite_signed_sum = (
        2 * baseline_per_type
        - 2 * baseline_per_type
        + positive_exception
        - negative_exception
    )
    transverse_by_type = {
        "+1": total_finite_signed_sum - 2,
        "-1": total_finite_signed_sum + 2,
    }
    contradictory_type = (
        -1 if (positive_exception, negative_exception) == (3, 1) else 1
    )
    return {
        "p": p,
        "baseline_directions_per_type": baseline_per_type,
        "positive_negative_exception_counts": [positive_exception, negative_exception],
        "total_finite_signed_sum": total_finite_signed_sum,
        "baseline_transverse_signed_sum": transverse_by_type,
        "K_required_transverse_signed_sum": 0,
        "contradictory_baseline_type": contradictory_type,
        "contradiction": baseline_per_type >= 1
        and transverse_by_type[f"{contradictory_type:+d}"] != 0,
    }


def theorem_negative_branch_exclusion() -> dict:
    samples = {
        str(p): [
            normal_form_signed_audit(p, 3, 1),
            normal_form_signed_audit(p, 1, 3),
        ]
        for p in (5, 7, 31, 101, 1009)
    }
    proved = all(row["contradiction"] for rows in samples.values() for row in rows)
    return {
        "proved": proved,
        "all_sufficiently_large_odd_primes": True,
        "conditional_input": "Proposition 15.644 normal form",
        "identity": "sum_{s<t} K_st=-(p-1) eps_d sum_s a_s=0",
        "exception_splits_excluded": [[3, 1], [1, 3]],
        "closes_negative_product_infinity_point_branch": proved,
        "closes_positive_product_infinity_point_branch": False,
        "closes_all_infinity_point_boundaries": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
        "samples": samples,
    }


def main() -> dict:
    theorem = theorem_negative_branch_exclusion()
    out = {
        "prop": "15.646",
        "title": "Signed transverse sums exclude the negative normal form",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15646.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
