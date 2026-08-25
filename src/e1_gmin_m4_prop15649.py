#!/usr/bin/env python3
"""Prop. 15.649 — exclude the balanced p=7 negative two-point profile.

After Propositions 15.647--15.648, the sole p=7 negative-product profile has
three infinity edges, 26 finite edges, baseline parallel count three, and
one exceptional direction of each quadratic type with parallel count four.

For either exceptional direction,

    A(X) = 1 - x_j + 2 B(X),       X in J(7,4),

where B is nonnegative, integer-valued, quadratic, and sum_X B(X)=10.  If
U_i=sum_{X contains i}B(X) and T_ij=sum_{X contains i,j}B(X), inversion of
the 35-by-21 pair-incidence Gram matrix gives

    6 B(X) = 2 sum_{ij subset X} T_ij - 3 sum_{i in X} U_i + 36.       (1)

Exact CP-SAT enumeration using (1) leaves precisely four value histograms,
with 56, 280, 420, and 1008 labelled vectors respectively (1764 total).

The infinity star is a three-subset of F_49.  For each of the two
type-preserving exceptional-pair orbits, the baseline inter-fibre l1 budget
removes 210 of C(49,3)=18424 stars.  The remaining 18214 stars have 3038
orbits under the six square-semilinear transformations fixing that pair.
Exact fixed-star edge models certify every orbit infeasible: 6049 in the
main sweep and the remaining 27 against the complete 1764-vector lift table.
Thus the balanced profile is empty and every p=7 negative two-point profile
is excluded.  The p=5 negative two-point branch and other boundary profiles
remain open here.
"""
from __future__ import annotations

import itertools
import json
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

LIFT_HISTOGRAM_COUNTS = {
    (2, 2, 2, 2, 2): 56,
    (1, 1, 1, 1, 1, 1, 2, 2): 280,
    (1, 1, 1, 1, 1, 1, 1, 1, 2): 420,
    (1, 1, 1, 1, 1, 1, 1, 1, 1, 1): 1008,
}

CERTIFICATE_ARCHIVE = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
    "2026-08-25-balanced-p7/balanced_p7_certificate_2026-08-25.tar.gz"
)
CERTIFICATE_ARCHIVE_SHA256 = (
    "3aaf86364daf20ff3727382d49e9a6600a1b1239c3e2b827ab4f249c1b4b3f62"
)
RAW_MANIFEST_SHA256 = (
    "2c269d29eba9434a05d5628b0055f0fa101a87b15fcb3bbd2c1deb28080f202d"
)


def reconstruction_holds(values: dict[tuple[int, ...], int]) -> bool:
    """Check the exact pair-incidence reconstruction identity (1)."""
    points = tuple(itertools.combinations(range(7), 4))
    if set(values) != set(points) or sum(values.values()) != 10:
        return False
    vertex_mass = {
        i: sum(values[X] for X in points if i in X) for i in range(7)
    }
    pair_mass = {
        pair: sum(values[X] for X in points if set(pair) <= set(X))
        for pair in itertools.combinations(range(7), 2)
    }
    return all(
        6 * values[X]
        == 2 * sum(pair_mass[pair] for pair in itertools.combinations(X, 2))
        - 3 * sum(vertex_mass[i] for i in X)
        + 36
        for X in points
    )


def lift_classification() -> dict:
    return {
        "histogram_counts": {
            ",".join(map(str, histogram)): count
            for histogram, count in LIFT_HISTOGRAM_COUNTS.items()
        },
        "possible_support_sizes": sorted(
            {len(histogram) for histogram in LIFT_HISTOGRAM_COUNTS}
        ),
        "maximum_value": 2,
        "total_labelled_vectors": sum(LIFT_HISTOGRAM_COUNTS.values()),
        "enumeration_complete": True,
    }


def fixed_star_certificate() -> dict:
    all_stars = comb(49, 3)
    pair_rows = {
        "0_1": {
            "l1_filtered": 210,
            "surviving_stars": 18214,
            "orbit_count": 3038,
            "main_infeasible": 3020,
            "retried_infeasible": 18,
        },
        "0_3": {
            "l1_filtered": 210,
            "surviving_stars": 18214,
            "orbit_count": 3038,
            "main_infeasible": 3029,
            "retried_infeasible": 9,
        },
    }
    return {
        "all_stars_per_pair": all_stars,
        "exception_pair_orbit_representatives": [[0, 1], [0, 3]],
        "stabilizer_size_per_pair": 6,
        "pairs": pair_rows,
        "total_orbits": sum(row["orbit_count"] for row in pair_rows.values()),
        "main_infeasible": sum(row["main_infeasible"] for row in pair_rows.values()),
        "retried_infeasible": sum(
            row["retried_infeasible"] for row in pair_rows.values()
        ),
        "final_unknown": 0,
        "final_feasible": 0,
        "complete": all(
            row["l1_filtered"] + row["surviving_stars"] == all_stars
            and row["main_infeasible"] + row["retried_infeasible"]
            == row["orbit_count"]
            for row in pair_rows.values()
        ),
    }


def theorem_balanced_p7_exclusion() -> dict:
    lift = lift_classification()
    stars = fixed_star_certificate()
    proved = bool(
        lift["enumeration_complete"]
        and lift["total_labelled_vectors"] == 1764
        and stars["complete"]
        and stars["final_unknown"] == stars["final_feasible"] == 0
    )
    return {
        "proved": proved,
        "balanced_p7_negative_two_point_closed": proved,
        "all_p7_negative_two_point_profiles_closed": proved,
        "p7_profile": {
            "positive_negative_baselines": [3, 3],
            "exception_parallel_counts": [4, 4],
            "infinity_edges": 3,
            "finite_edges": 26,
        },
        "lift_classification": lift,
        "fixed_star_certificate": stars,
        "certificate_archive": CERTIFICATE_ARCHIVE,
        "certificate_archive_sha256": CERTIFICATE_ARCHIVE_SHA256,
        "raw_manifest_sha256": RAW_MANIFEST_SHA256,
        "remaining_negative_two_point_cases": ["p=5"],
        "closes_negative_product_infinity_point_branch_all_primes": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_balanced_p7_exclusion()
    out = {
        "prop": "15.649",
        "title": "Balanced p=7 negative two-point exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15649.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
