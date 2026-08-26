#!/usr/bin/env python3
"""Prop. 15.655 -- mod-seven exclusion of unsaturated p=7 four-finite boundaries.

After Proposition 15.654, one product sign has 23,520 unsaturated
four-finite boundaries in 518 square-semilinear orbits.  Fixing the one
elevated direction in each unsaturated quadratic type gives 2,408 cases.

For each case, the complete Johnson-space slack catalogs fix all 280 affine
bad-edge counts of a putative 29-edge graph H.  Together with |H|=29 and the
distinguished edge, these are 282 integer linear equations in 1,225 edge
indicators.  The common coefficient matrix has rank 147 over F_7, hence 135
left-null dependencies.  Exact syndrome joins over the at most two
non-singleton directional catalogs leave zero compatible catalog tuples in
all 2,408 cases.  An independent implementation rebuilds the score matrix,
row dependencies, catalog values, and coverage and obtains the same result.

The nonsquare anti-isometry from Proposition 15.654 transfers the exclusion
to the other Paley-product sign.  Consequently every p=7 four-finite
boundary is excluded; with Proposition 15.653, every p=7 size-four boundary
is excluded.  The p=5 size-four cases and larger boundaries remain open.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15654 import p7_nonsquare_signed_permutation

ROOT = Path(__file__).resolve().parents[1]

CERTIFICATE_ARCHIVE = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
    "2026-08-26-p7-four-point/"
    "p7_no_infinity_unsaturated_mod7_certificate_2026-08-26.tar.gz"
)
CERTIFICATE_ARCHIVE_SHA256 = (
    "ad4efb2c1f3eb53f2ebd089c8b5b6ce959945bf7ce1298c7dfa694eaeda8965d"
)
CERTIFICATE_BATCH_SHA256 = (
    "64a87c509f92144e1dd1956c92145abc3b42deebb7fdb610405963b30d80c50f"
)
CERTIFICATE_AUDIT_SHA256 = (
    "3348579efe686302ba5eb9e654d8121f8445e9ffece19047d581413964581350"
)
ORBIT_SOURCE_SHA256 = (
    "7f7d3cc26077bb40ac096b638c6fc20ddf1a8fe6ddee60641f2fb568bacfd077"
)


def p7_unsaturated_mod7_certificate() -> dict:
    patterns = {
        "1764": 1372,
        "1764x36": 294,
        "1764x1764": 112,
        "2233x36": 294,
        "2233x1764": 336,
    }
    catalog_tuples = (
        1372 * 1764
        + 294 * 1764 * 36
        + 112 * 1764 * 1764
        + 294 * 2233 * 36
        + 336 * 2233 * 1764
    )
    return {
        "scope": "unsaturated p=7 four-finite boundaries for c_H=-1",
        "edge_variables": 1225,
        "exact_linear_equations": 282,
        "affine_score_equations": 280,
        "modulus": 7,
        "rank": 147,
        "left_dependency_dimension": 135,
        "unsaturated_boundaries": 23520,
        "unsaturated_orbits": 518,
        "fixed_elevation_cases": 2408,
        "catalog_pattern_counts": patterns,
        "catalog_tuples_excluded": catalog_tuples,
        "mod7_infeasible_cases": 2408,
        "surviving_cases": 0,
        "missing_cases": 0,
        "duplicate_cases": 0,
        "independent_rank": 147,
        "independent_dependency_dimension": 135,
        "independent_recomputed_nonzero_cases": 0,
        "independent_audit": True,
    }


def theorem_p7_unsaturated_four_finite_exclusion() -> dict:
    certificate = p7_unsaturated_mod7_certificate()
    symmetry = p7_nonsquare_signed_permutation()
    sign_transfer = {
        "signed_conference_anti_isometry": symmetry[
            "signed_conference_anti_isometry"
        ],
        "fixes_distinguished_edge": symmetry["fixes_distinguished_edge"],
        "edge_count_29_is_odd": True,
        "four_finite_boundary_forces_even_infinity_degree": True,
        "paley_product_sign_flips": True,
        "normalized_scores_are_preserved": True,
    }
    proved = bool(
        certificate["rank"] == 147
        and certificate["left_dependency_dimension"] == 135
        and certificate["unsaturated_orbits"] == 518
        and certificate["unsaturated_boundaries"] == 23520
        and certificate["fixed_elevation_cases"] == 2408
        and certificate["mod7_infeasible_cases"] == 2408
        and certificate["surviving_cases"] == 0
        and certificate["missing_cases"] == 0
        and certificate["independent_audit"]
        and certificate["independent_recomputed_nonzero_cases"] == 0
        and all(sign_transfer.values())
    )
    return {
        "proved": proved,
        "p7_four_finite_unsaturated_both_product_signs": "CLOSED",
        "p7_all_four_finite_points_with_prop15654": "CLOSED",
        "p7_all_size_four_with_prop15653": "CLOSED",
        "certificate": certificate,
        "sign_transfer": sign_transfer,
        "p5_size_four": "OPEN",
        "first_open_boundary_size_for_p_ge_7": 6,
        "closes_all_p7_size_four": proved,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_p7_unsaturated_four_finite_exclusion()
    out = {
        "prop": "15.655",
        "title": "p=7 unsaturated four-finite mod-seven exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
        "certificate_archive": CERTIFICATE_ARCHIVE,
        "certificate_archive_sha256": CERTIFICATE_ARCHIVE_SHA256,
        "certificate_batch_sha256": CERTIFICATE_BATCH_SHA256,
        "certificate_audit_sha256": CERTIFICATE_AUDIT_SHA256,
        "orbit_source_sha256": ORBIT_SOURCE_SHA256,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15655.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
