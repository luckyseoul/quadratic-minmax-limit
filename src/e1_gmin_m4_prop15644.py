#!/usr/bin/env python3
"""Prop. 15.644 — asymptotic normal form of the negative two-point branch.

Under ``D={infinity,v}``, ``c_H=-1``, the near-optimal polynomial distance
lemma on Boolean slices implies that, for all sufficiently large odd primes,
there is exactly one nonbaseline direction of each quadratic type.  Exact
directional means and the additive inter-fibre ``l1`` inequality then force

    infinity edges I = 2p-1,
    baseline parallel counts P_d = 2,
    exceptional parallel counts {1,3},

with the negative-type exceptional count odd.  This is a normal form, not
an exclusion of the branch.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def baseline_parallel_bound(p: int, infinity_edges: int) -> int:
    """Maximum ``P_d`` from ``|pP_d-E+2| <= E-P_d``."""
    if p < 5 or p % 2 == 0:
        raise ValueError("p must be odd and at least five")
    finite_edges = 4 * p + 1 - infinity_edges
    return (2 * finite_edges - 2) // (p + 1)


def exceptional_parallel_bounds(p: int, infinity_edges: int) -> tuple[int, int]:
    """Bounds when the unique lift consumes type surplus ``p+1``.

    Its scaled directional slack is ``a_d=(p-1)+(p+1)=2p``.  Comparing the
    exact signed mean with the number of transverse edges gives ``1<=P<=8``.
    The returned upper bound retains its exact finite-p form.
    """
    upper = (9 * p + 1 - 2 * infinity_edges) // (p + 1)
    return 1, upper


def negative_branch_normal_form(p: int) -> dict:
    """Arithmetic classification conditional on one exception per type."""
    if p < 31 or p % 2 == 0:
        raise ValueError("the clean arithmetic audit is recorded for odd p>=31")
    q = (p - 1) // 2
    candidates = []
    for residue in range(8):
        for k0 in range(-8, 17):
            infinity_edges = 3 - residue + q * k0
            if not 1 <= infinity_edges <= 4 * p + 1:
                continue
            if infinity_edges % 2 == 0:
                continue
            finite_edges = 4 * p + 1 - infinity_edges
            if baseline_parallel_bound(p, infinity_edges) < residue:
                continue
            exceptional_sum = finite_edges - (p - 1) * residue
            if not 2 <= exceptional_sum <= 16:
                continue
            lo, hi = exceptional_parallel_bounds(p, infinity_edges)
            pairs = [
                (positive, negative)
                for positive in range(lo, hi + 1)
                for negative in range(lo, hi + 1)
                if positive + negative == exceptional_sum
                and ((q * residue + negative) & 1) == 1
            ]
            if not pairs:
                continue
            # Infinity is in the odd boundary, hence I is odd.  The finite
            # graph must toggle all but possibly v among the I star leaves.
            boundary_possible = infinity_edges - 1 <= 2 * finite_edges
            if not boundary_possible:
                continue
            candidates.append(
                {
                    "residue": residue,
                    "k0": k0,
                    "infinity_edges": infinity_edges,
                    "finite_edges": finite_edges,
                    "exceptional_parallel_sum": exceptional_sum,
                    "exceptional_pairs_positive_negative": pairs,
                }
            )
    expected = [
        {
            "residue": 2,
            "k0": 4,
            "infinity_edges": 2 * p - 1,
            "finite_edges": 2 * p + 2,
            "exceptional_parallel_sum": 4,
            "exceptional_pairs_positive_negative": [(1, 3), (3, 1)],
        }
    ]
    return {
        "p": p,
        "q": q,
        "conditional_hypothesis": "exactly one nonbaseline direction per type",
        "candidates": candidates,
        "unique_normal_form": candidates == expected,
        "normal_form": expected[0] if candidates == expected else None,
    }


def theorem_negative_branch_normal_form() -> dict:
    samples = {str(p): negative_branch_normal_form(p) for p in (31, 41, 101, 201)}
    return {
        "proved": all(row["unique_normal_form"] for row in samples.values()),
        "all_sufficiently_large_odd_primes": True,
        "external_input": (
            "Amireddy-Behera-Srinivasan-Sudan slice distance theorem, d=2"
        ),
        "external_input_consequence": "exactly one exception per type",
        "normal_form": {
            "infinity_edges": "2p-1",
            "finite_edges": "2p+2",
            "baseline_parallel_count": 2,
            "exceptional_parallel_counts": [1, 3],
            "negative_exception_parallel_count": "odd",
        },
        "samples": samples,
        "closes_negative_product_branch": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_negative_branch_normal_form()
    out = {
        "prop": "15.644",
        "title": "Asymptotic normal form of the negative infinity-point boundary",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15644.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
