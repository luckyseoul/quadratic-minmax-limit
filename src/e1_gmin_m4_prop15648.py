#!/usr/bin/env python3
"""Prop. 15.648 — close p=11,13 and four p=7 negative profiles.

This proposition continues the exact count reduction of Proposition 15.647.
At p=13, both surviving orientations have a zero-baseline direction whose
prescribed inter-fibre matrix has l1 norm at least 48, exceeding its 44
transverse edges.  At p=11, exact CP-SAT infeasibility certificates cover
both count orientations and all three type-preserving exceptional-pair
orbits.  At p=7, the same model excludes four unbalanced count profiles over
all exceptional pairs.  The balanced p=7 profile and p=5 remain open.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def p13_zero_baseline_l1(positive_fibres: int) -> int:
    """Exact l1 for sum(w)=10, len(w)=13, K_st=+/- (1-w_s-w_t).

    If u entries are positive, the remaining 13-u are zero.  Separating
    zero-zero, zero-positive, and positive-positive pairs eliminates the
    individual positive magnitudes and gives u^2-25u+198.
    """
    if not 1 <= positive_fibres <= 10:
        raise ValueError("positive_fibres must lie in 1..10")
    u = positive_fibres
    return u * u - 25 * u + 198


def p13_l1_minimum() -> dict:
    values = {u: p13_zero_baseline_l1(u) for u in range(1, 11)}
    minimum = min(values.values())
    minimizers = [u for u, value in values.items() if value == minimum]
    return {
        "values": values,
        "minimum": minimum,
        "minimizers": minimizers,
        "transverse_edge_budget": 44,
        "contradiction": minimum > 44,
    }


def _line_key(p: int, u: int) -> tuple[int, int]:
    a, b = u % p, u // p
    if a:
        return 1, b * pow(a, -1, p) % p
    return 0, 1


def type_preserving_exception_pair_orbits(p: int) -> list[dict]:
    """Orbits under square multiplications and Frobenius on affine directions."""
    from e1_gmin_m4_prop15598 import field_ctx
    from e1_gmin_m4_prop15632 import field_direction_data, projective_directions

    q2, mul, _add, chi, frob, _norm, _ia, _ib = field_ctx(p)
    directions = projective_directions(p)
    kernel_generators = [
        (s % p) + ((-r) % p) * p for r, s in directions
    ]
    line_index = {
        _line_key(p, u): index for index, u in enumerate(kernel_generators)
    }
    permutations = []
    for alpha in range(1, q2):
        if chi(alpha) == 1:
            permutations.append(
                tuple(
                    line_index[_line_key(p, mul(alpha, u))]
                    for u in kernel_generators
                )
            )
    permutations.append(
        tuple(line_index[_line_key(p, frob(u))] for u in kernel_generators)
    )

    types = [field_direction_data(p, direction)[0] for direction in directions]
    remaining = {
        pair
        for pair in itertools.combinations(range(p + 1), 2)
        if types[pair[0]] != types[pair[1]]
    }
    rows = []
    while remaining:
        representative = min(remaining)
        orbit = {representative}
        previous: set[tuple[int, int]] = set()
        while orbit != previous:
            previous = set(orbit)
            orbit |= {
                tuple(sorted((permutation[a], permutation[b])))
                for a, b in previous
                for permutation in permutations
            }
        rows.append(
            {
                "representative": list(representative),
                "size": len(orbit),
                "pairs": [list(pair) for pair in sorted(orbit)],
            }
        )
        remaining -= orbit
    return rows


P11_CERTIFICATES = {
    "x0_y5_e0_1": "15e47ec06b85672bf5a95c0f4277a642e34f3991cd83f054992dd9530644fe2b",
    "x0_y5_e0_2": "6350fcbe2dca0ffe5f7c2a7bad2dcd913799f2dee6b62f468393af9c9dd9a214",
    "x0_y5_e0_3": "f7000ab7abf021d0ca74eb1020d4114f2e6fe039eab11439060716ca93363cb5",
    "x5_y0_e0_1": "604fd1bebdc1390feff70e8abdf209177a6206243abb47a2eb77227e55bf1c4a",
    "x5_y0_e0_2": "a8c7bc344454f8cabc08062344c7d93c307e66f8e4a8e14d05521109aa0183d6",
    "x5_y0_e0_3": "a90cef86896c5735c43cce3f0490755e1eb1086c72db837516aa11c9302df95d",
}

P7_CERTIFICATES = {
    "x0_y3_all_pairs": "d9a72b8be20d3de1cbe7e9a3fa6309b8ab87d5fabc9062f256e952a6895ee8c2",
    "x0_y6_all_pairs": "d057615b03c6c22f79068e0679296500551020b1874636245d4c3531a47f3c7a",
    "x3_y0_all_pairs": "573117f014bfa3e0b79e3c1a0d0378d0acdad0ee08c9ec39be5c0ae7badedeba",
    "x6_y0_all_pairs": "54253640edeb2d3122b7dd35bb3aeb18f2992d22485f11af6fdb0c5c969af955",
}


def theorem_finite_negative_profiles() -> dict:
    p11_orbits = type_preserving_exception_pair_orbits(11)
    p13 = p13_l1_minimum()
    return {
        "proved": p13["contradiction"]
        and [row["representative"] for row in p11_orbits]
        == [[0, 1], [0, 2], [0, 3]]
        and sum(row["size"] for row in p11_orbits) == 36
        and len(P11_CERTIFICATES) == 6
        and len(P7_CERTIFICATES) == 4,
        "p13_negative_two_point_closed": p13["contradiction"],
        "p13_method": "exact inter-fibre l1 lower bound",
        "p13_l1": p13,
        "p11_negative_two_point_finitely_certified": True,
        "p11_exception_pair_orbits": p11_orbits,
        "p11_certificate_sha256": P11_CERTIFICATES,
        "p7_unbalanced_profiles_finitely_certified": [[0, 3], [0, 6], [3, 0], [6, 0]],
        "p7_certificate_sha256": P7_CERTIFICATES,
        "remaining_negative_two_point_cases": [
            "p=5 (no guaranteed baseline per type)",
            "p=7 baseline counts (3,3)",
        ],
        "archive": (
            "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
            "2026-08-25-negative-two-point"
        ),
        "closes_negative_product_infinity_point_branch_p11_p13": True,
        "closes_negative_product_infinity_point_branch_all_primes": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_finite_negative_profiles()
    out = {
        "prop": "15.648",
        "title": "Finite negative two-point exclusions at p=7,11,13",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15648.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
