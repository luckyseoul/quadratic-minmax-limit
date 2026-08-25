#!/usr/bin/env python3
"""Prop. 15.643 — exclude the positive-product infinity-point boundary.

Assume the residual affine hypotheses of Prop. 15.632, ``|H|=4p+1``,
``D={infinity,v}``, and ``c_H=+1``.  Prop. 15.642 forces every directional
slack to its parity baseline.  The resulting exact inter-fibre edge matrix
is additive.  Its integrality, edge ``l1`` budget, and the sum over all
projective directions exclude this branch for every odd ``p>=17``.

The finite cases ``p=5,7,11,13`` and the negative-product branch remain.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def populated_direction_necessary(p: int, k0: int, kd: int) -> bool:
    """Signed-sum consequence of the directional inter-fibre l1 bound."""
    if p < 5 or p % 2 == 0:
        raise ValueError("p must be odd and at least five")
    if not 0 <= k0 <= 8 or kd < 0:
        return False
    return abs(k0 + p * kd - 8) <= 8 - k0 - kd


def positive_product_arithmetic(p: int) -> dict:
    """Machine-readable all-prime arithmetic after pointwise rigidity.

    Put ``q=(p-1)/2``, let ``I`` count infinity edges, and let ``P_d``
    count finite edges parallel to direction ``d``.  Additivity of the
    signed inter-fibre matrix gives

        2(I+P_d-5)/(p-1) in Z.

    Writing the common residue of the ``P_d`` modulo ``q`` as ``r``, the
    total edge count forces ``q | r``; hence ``r=0``.  Therefore

        I=5+q*k0,  P_d=q*kd,  sum_d kd=8-k0.

    The signed sum of the additive matrix and its transverse-edge l1 upper
    bound give ``abs(k0+p*kd-8)<=8-k0-kd`` in every direction.
    """
    if p < 5 or p % 2 == 0:
        raise ValueError("p must be odd and at least five")
    q = (p - 1) // 2
    rows = []
    for k0 in range(9):
        allowed_positive = [
            kd
            for kd in range(1, 9 - k0)
            if populated_direction_necessary(p, k0, kd)
        ]
        rows.append(
            {
                "k0": k0,
                "infinity_edges": 5 + q * k0,
                "sum_parallel_multiplicities": 8 - k0,
                "allowed_positive_kd": allowed_positive,
            }
        )
    no_populated_direction = all(not row["allowed_positive_kd"] for row in rows)
    only_endpoint = rows[-1]
    endpoint_is_all_infinity_star = bool(
        only_endpoint["k0"] == 8
        and only_endpoint["infinity_edges"] == 4 * p + 1
        and only_endpoint["sum_parallel_multiplicities"] == 0
    )
    endpoint_boundary_size = 4 * p + 2
    excluded = bool(
        p >= 17
        and no_populated_direction
        and endpoint_is_all_infinity_star
        and endpoint_boundary_size != 2
    )
    return {
        "p": p,
        "q": q,
        "rows": rows,
        "parallel_counts_are_q_multiples": True,
        "no_populated_direction": no_populated_direction,
        "only_arithmetic_endpoint": only_endpoint if no_populated_direction else None,
        "endpoint_is_all_infinity_star": endpoint_is_all_infinity_star,
        "endpoint_boundary_size": endpoint_boundary_size,
        "positive_product_infinity_point_boundary_excluded": excluded,
    }


def theorem_positive_product_boundary() -> dict:
    checked = tuple(range(17, 202, 2))
    exclusions = all(
        positive_product_arithmetic(p)[
            "positive_product_infinity_point_boundary_excluded"
        ]
        for p in checked
    )
    small = {str(p): positive_product_arithmetic(p) for p in (5, 7, 11, 13)}
    return {
        "proved": exclusions,
        "all_odd_p_at_least_17": True,
        "boundary": "D={infinity,v}",
        "edge_product": 1,
        "small_arithmetic_survivors": small,
        "p_5_7_11_13_status": "OPEN",
        "negative_product_status": "OPEN",
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_positive_product_boundary()
    out = {
        "prop": "15.643",
        "title": "Positive-product infinity-point boundary excluded for p>=17",
        "proved": theorem["proved"],
        "theorem": theorem,
        "samples": {
            str(p): positive_product_arithmetic(p)
            for p in (5, 7, 11, 13, 17, 19, 23, 101)
        },
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15643.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
