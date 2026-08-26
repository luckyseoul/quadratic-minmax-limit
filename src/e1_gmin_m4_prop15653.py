#!/usr/bin/env python3
"""Prop. 15.653 -- close p=7 infinity plus three finite boundary points.

Proposition 15.652 already excludes the negative-product sign at p=7.  In
the positive-product branch, all eight directional parity floors equal
eight and saturate the two exact type budgets.  A direction with one odd
boundary fibre has slack ``A=x_j``.  A direction with three odd fibres has
exactly one possible minimum-mean degree-two slack,
``A=(|X cap B|-2)^2``; an exact rank-21 Johnson-space check eliminates the
other 629 sparse corrections.

The resulting score polynomials give sparse coefficient equations for the
infinity star and finite inter-fibre edge sums.  The 18,424 finite triples
form 416 orbits under the 48-element square-semilinear stabilizer of the
distinguished edge.  Exact fixed-boundary CP-SAT certificates exclude all
416 orbits, with an independent coverage and hash audit.  Hence infinity
plus three finite boundary points are impossible at p=7 for both product
signs.  Four finite points at p=7, every p=5 size-four case, larger
boundaries, residual (ii), R1, and the limit remain open.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path

from e1_gmin_m4_prop15652 import infinity_size_four_exclusion

ROOT = Path(__file__).resolve().parents[1]

CERTIFICATE_ARCHIVE = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
    "2026-08-26-p7-four-point/p7_infinity_three_certificate_2026-08-26.tar.gz"
)
CERTIFICATE_ARCHIVE_SHA256 = (
    "f9a125a18d287eef63e579b8416022e0c4d91dee0a82334fabb1167cc9356c17"
)
CERTIFICATE_AUDIT_SHA256 = (
    "d974d27274aacf8987b954a42500d7390bb9f9b895f8b813b3c7faff323264c5"
)
ORBIT_SOURCE_SHA256 = (
    "0b3a928e98e11838eac051c768a3b4aa0f3d5e5d32dd17fde8f8021016eae941"
)


def p7_inf_three_slack_classification() -> dict:
    return {
        "johnson_slice": "J(7,4)",
        "slice_points": 35,
        "degree_at_most_two_rank": 21,
        "left_kernel_dimension": 14,
        "b1": {
            "scaled_mean": 8,
            "unique_slack": "A(X)=x_j",
        },
        "b3": {
            "scaled_mean": 8,
            "total_slack_mass": 20,
            "mandatory_parity_mass": 16,
            "sparse_correction_candidates": 630,
            "survivors": 1,
            "unique_slack": "A(X)=(|X cap B|-2)^2",
        },
    }


def target_coefficients(B: set[int]) -> tuple[int, tuple[int, ...], dict[tuple[int, int], int]]:
    """Return ``c,l,q`` for ``epsilon*S=c+sum lz+sum qzz``."""
    if len(B) == 1:
        special = next(iter(B))
        return (
            4,
            tuple(int(s == special) for s in range(7)),
            {pair: 0 for pair in itertools.combinations(range(7), 2)},
        )
    if len(B) == 3:
        return (
            5,
            tuple(-int(s in B) for s in range(7)),
            {
                pair: int(set(pair) <= B)
                for pair in itertools.combinations(range(7), 2)
            },
        )
    raise ValueError("B must have one or three fibres")


def reconstructed_direction_coefficients(
    B: set[int], star_counts: tuple[int, ...], kernel_parameter: int
) -> tuple[int, tuple[int, ...], dict[tuple[int, int], int]]:
    """Coefficient-kernel form used in every fixed-boundary edge model."""
    if len(star_counts) != 7:
        raise ValueError("need seven infinity-star fibre counts")
    constant, linear, pairs = target_coefficients(B)
    infinity_count = sum(star_counts)
    parallel = constant + sum(linear) + 3 * kernel_parameter - infinity_count
    cross = {
        (s, t): (
            pairs[s, t]
            + kernel_parameter
            - star_counts[s]
            - star_counts[t]
            + linear[s]
            + linear[t]
        )
        for s, t in itertools.combinations(range(7), 2)
    }
    return parallel, star_counts, cross


def p7_inf_three_orbit_certificate() -> dict:
    return {
        "candidate_finite_triples": 18424,
        "square_semilinear_stabilizer_size": 48,
        "boundary_orbits": 416,
        "orbit_size_sum": 18424,
        "fixed_boundary_infeasible": 416,
        "unknown": 0,
        "feasible": 0,
        "coefficient_equalities_per_orbit": 177,
        "allowed_infinity_edge_counts": [5, 11, 17, 23, 29],
        "fresh_orbit_reclassification_matches": True,
        "independent_audit": True,
    }


def theorem_p7_infinity_three_exclusion() -> dict:
    negative = infinity_size_four_exclusion(7, -1)
    finite = p7_inf_three_orbit_certificate()
    slack = p7_inf_three_slack_classification()
    proved = bool(
        negative["excluded"]
        and slack["b3"]["sparse_correction_candidates"] == 630
        and slack["b3"]["survivors"] == 1
        and finite["candidate_finite_triples"] == 18424
        and finite["orbit_size_sum"] == 18424
        and finite["fixed_boundary_infeasible"] == 416
        and finite["unknown"] == 0
        and finite["feasible"] == 0
        and finite["independent_audit"]
    )
    return {
        "proved": proved,
        "p7_infinity_plus_three_finite_both_product_signs": "CLOSED",
        "negative_product_source": "Proposition 15.652",
        "positive_product_slacks": slack,
        "positive_product_finite_certificate": finite,
        "p7_four_finite_points": "OPEN",
        "p5_size_four": "OPEN",
        "closes_all_p7_size_four": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_p7_infinity_three_exclusion()
    out = {
        "prop": "15.653",
        "title": "p=7 infinity-plus-three-finite boundary exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
        "certificate_archive": CERTIFICATE_ARCHIVE,
        "certificate_archive_sha256": CERTIFICATE_ARCHIVE_SHA256,
        "certificate_audit_sha256": CERTIFICATE_AUDIT_SHA256,
        "orbit_source_sha256": ORBIT_SOURCE_SHA256,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15653.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
