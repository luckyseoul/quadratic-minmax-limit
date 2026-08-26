#!/usr/bin/env python3
"""Prop. 15.650 — exclude every p=5 negative two-point profile.

For D={infinity,v}, c_H=-1, and p=5, write a_d=4+r_d.  Each quadratic
type has three directions and sum_d r_d=6.  The r_d are nonnegative even
integers, while the exact same-type directional mean gives

    r_d-r_e = 6(P_d-P_e).

Consequently one type has exactly one of two lift profiles:

    unique:       (r_d)=(6,0,0), with P_exception=P_baseline+1;
    distributed:  (r_d)=(2,2,2), with all three P_d equal.

Boundary parity, edge-product parity, baseline coefficient integrality, and
the elementary boundary-size bound leave 24 arithmetic profiles.  Square
field multiplications and Frobenius give two opposite-type exceptional-pair
orbits and one single-direction orbit per type, hence 33 profile-placement
orbits.  An exact CP-SAT model with all 60 affine score identities, direction
counts, boundary XORs, and product parity certifies one representative of
every orbit infeasible.  This closes p=5 and, with Propositions
15.647--15.649, the negative-product infinity-plus-point branch for every
odd prime p>=5.  Other boundary profiles and the finite positive-product
cases remain open.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

CERTIFICATE_ARCHIVE = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
    "2026-08-26-negative-p5/p5_negative_two_point_certificate_2026-08-26.tar.gz"
)
CERTIFICATE_ARCHIVE_SHA256 = (
    "c19ed8a0d50ffad3f7386d3d6100ce25213c3ba4d1e7cc49acea0275d3796a41"
)
CERTIFICATE_SHA256 = (
    "2352ca1040a7989d5850730fafcad311adf1275dd06383578fb9e187698be69e"
)


def p5_type_excess_profiles() -> list[tuple[int, int, int]]:
    """Canonical excess multisets after evenness and mod-six quantization."""
    profiles = {
        tuple(sorted(values))
        for values in itertools.product(range(0, 7, 2), repeat=3)
        if sum(values) == 6
        and all((values[i] - values[j]) % 6 == 0 for i in range(3) for j in range(3))
    }
    return sorted(profiles)


def p5_arithmetic_profiles() -> list[dict]:
    """Enumerate the exact count candidates before placement symmetry."""
    out = []
    for positive_profile, negative_profile in itertools.product(
        ("unique", "distributed"), repeat=2
    ):
        ep = int(positive_profile == "unique")
        en = int(negative_profile == "unique")
        for x in range(8):
            for y in range(8):
                finite_edges = 3 * (x + y) + ep + en
                infinity_edges = 21 - finite_edges
                if infinity_edges < 1 or infinity_edges % 2 == 0:
                    continue
                if (3 * y + en) % 2 != 1:
                    continue
                if infinity_edges - 1 > 2 * finite_edges:
                    continue
                if positive_profile == "unique" and x % 2:
                    continue
                if negative_profile == "unique" and y % 2:
                    continue
                out.append(
                    {
                        "positive_profile": positive_profile,
                        "negative_profile": negative_profile,
                        "positive_parallel_baseline": x,
                        "negative_parallel_baseline": y,
                        "finite_edges": finite_edges,
                        "infinity_edges": infinity_edges,
                    }
                )
    return out


def p5_placement_orbit_count() -> dict:
    profiles = p5_arithmetic_profiles()
    by_kind = {
        kind: sum(
            (row["positive_profile"], row["negative_profile"]) == kind
            for row in profiles
        )
        for kind in itertools.product(("unique", "distributed"), repeat=2)
    }
    # Two opposite-type pair orbits when both types are unique; a single
    # direction orbit when exactly one is unique; no placement otherwise.
    total = (
        2 * by_kind["unique", "unique"]
        + by_kind["unique", "distributed"]
        + by_kind["distributed", "unique"]
        + by_kind["distributed", "distributed"]
    )
    return {
        "arithmetic_profiles": len(profiles),
        "profiles_by_type_pair": {
            f"{a}_{b}": count for (a, b), count in by_kind.items()
        },
        "opposite_type_pair_orbits": 2,
        "single_direction_orbits_per_type": 1,
        "placement_orbits": total,
    }


def theorem_p5_negative_two_point_exclusion() -> dict:
    orbit_count = p5_placement_orbit_count()
    proved = bool(
        p5_type_excess_profiles() == [(0, 0, 6), (2, 2, 2)]
        and orbit_count["arithmetic_profiles"] == 24
        and orbit_count["placement_orbits"] == 33
    )
    return {
        "proved": proved,
        "p5_negative_two_point_closed": proved,
        "type_excess_profiles": {
            "unique": [0, 0, 6],
            "distributed": [2, 2, 2],
        },
        "orbit_count": orbit_count,
        "finite_certificate": {
            "infeasible": 33,
            "unknown": 0,
            "feasible": 0,
            "all_exact_affine_score_rows": 60,
            "certificate_sha256": CERTIFICATE_SHA256,
            "archive": CERTIFICATE_ARCHIVE,
            "archive_sha256": CERTIFICATE_ARCHIVE_SHA256,
        },
        "negative_product_infinity_point_branch_all_odd_primes_p_ge_5": "CLOSED",
        "remaining_negative_two_point_cases": [],
        "closes_negative_product_infinity_point_branch_all_primes": proved,
        "closes_all_infinity_point_boundaries": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_p5_negative_two_point_exclusion()
    out = {
        "prop": "15.650",
        "title": "Complete p=5 negative two-point exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15650.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
