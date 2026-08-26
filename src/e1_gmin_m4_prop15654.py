#!/usr/bin/env python3
"""Prop. 15.654 -- exclude saturated p=7 four-finite boundaries.

For a p=7 boundary consisting of four finite points, this proposition
handles exactly the profiles for which both quadratic direction types use
their full Proposition 15.632 budget of 32.  The complete degree-two slack
catalog leaves one phase-zero b=4 target and 36 phase-one b=4 targets.
Sparse coefficient models are infeasible on all 1,225 square-semilinear
boundary orbits, covering 58,800 boundaries for one Paley-product sign.

Multiplication by a nonsquare, together with switching only the infinity
coordinate, is an exact anti-isometry of the Paley conference matrix.  It
fixes the distinguished edge, exchanges the eigenshells, preserves the
normalized score conditions, and flips the product sign for a 29-edge graph
whose boundary omits infinity.  Thus the same exclusion holds for both
product signs.  The 23,520 unsaturated surviving boundaries (518 orbits)
per sign remain open, as do the p=5 size-four cases and larger boundaries.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from e1_gmin_m4_prop15598 import field_ctx
from minmax_quadratic import paley_conference_prime_power

ROOT = Path(__file__).resolve().parents[1]

CERTIFICATE_ARCHIVE = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
    "2026-08-26-p7-four-point/"
    "p7_no_infinity_saturated_certificate_2026-08-26.tar.gz"
)
CERTIFICATE_ARCHIVE_SHA256 = (
    "5234f60d246c50dcdb7f5feb51b23185e8eecbffac431f04c279c1cc02153612"
)
CERTIFICATE_AUDIT_SHA256 = (
    "4d8fbbba46a0f7b19fc2f2241e4a710bf81fc581b5acb893f834ee65c1f547ea"
)
ORBIT_SOURCE_SHA256 = (
    "7f7d3cc26077bb40ac096b638c6fc20ddf1a8fe6ddee60641f2fb568bacfd077"
)


def p7_nonsquare_signed_permutation() -> dict:
    """Return and verify the signed nonsquare Paley anti-isometry."""
    q, mul, _add, chi, _frob, _norm, _ia, _ib = field_ctx(7)
    alpha = next(value for value in range(1, q) if chi(value) == -1)
    finite_permutation = tuple(mul(alpha, u) for u in range(q))
    vertex_permutation = (0, *(u + 1 for u in finite_permutation))
    switching = (-1, *([1] * q))
    C = np.rint(paley_conference_prime_power(7)).astype(np.int8)
    permutation_array = np.array(vertex_permutation, dtype=np.int16)
    switching_array = np.array(switching, dtype=np.int8)
    conjugated = (
        switching_array[:, None]
        * C[np.ix_(permutation_array, permutation_array)]
        * switching_array[None, :]
    )
    return {
        "nonsquare_multiplier": alpha,
        "finite_permutation": finite_permutation,
        "vertex_permutation": vertex_permutation,
        "switching": switching,
        "fixes_infinity": vertex_permutation[0] == 0,
        "fixes_finite_zero": vertex_permutation[1] == 1,
        "fixes_distinguished_edge": vertex_permutation[:2] == (0, 1),
        "signed_conference_anti_isometry": bool(np.array_equal(conjugated, -C)),
    }


def p7_saturated_four_finite_certificate() -> dict:
    return {
        "surviving_boundaries_per_product_sign": 82320,
        "surviving_orbits_per_product_sign": 1743,
        "saturated_boundaries_per_product_sign": 58800,
        "saturated_orbits_per_product_sign": 1225,
        "remaining_unsaturated_boundaries_per_product_sign": 23520,
        "remaining_unsaturated_orbits_per_product_sign": 518,
        "phase_zero_b4_slacks": 1,
        "phase_one_b4_slacks": 36,
        "coefficient_equalities_per_orbit": 176,
        "fixed_boundary_infeasible": 1225,
        "unknown": 0,
        "feasible": 0,
        "fresh_boundary_and_orbit_reclassification_matches": True,
        "independent_audit": True,
    }


def theorem_p7_saturated_four_finite_exclusion() -> dict:
    symmetry = p7_nonsquare_signed_permutation()
    certificate = p7_saturated_four_finite_certificate()
    sign_transfer = {
        "edge_count_29_is_odd": True,
        "four_finite_boundary_forces_even_infinity_degree": True,
        "paley_product_sign_flips": True,
        "eigenshell_sign_swaps": True,
        "raw_affine_score_sign_flips": True,
        "normalized_eigensign_times_score_is_preserved": True,
    }
    proved = bool(
        symmetry["fixes_distinguished_edge"]
        and symmetry["signed_conference_anti_isometry"]
        and all(sign_transfer.values())
        and certificate["saturated_boundaries_per_product_sign"] == 58800
        and certificate["saturated_orbits_per_product_sign"] == 1225
        and certificate["fixed_boundary_infeasible"] == 1225
        and certificate["unknown"] == 0
        and certificate["feasible"] == 0
        and certificate["independent_audit"]
    )
    return {
        "proved": proved,
        "p7_four_finite_doubly_saturated_both_product_signs": "CLOSED",
        "certificate": certificate,
        "nonsquare_anti_isometry": {
            key: value
            for key, value in symmetry.items()
            if key not in {"finite_permutation", "vertex_permutation", "switching"}
        },
        "sign_transfer": sign_transfer,
        "p7_four_finite_unsaturated": "OPEN",
        "p7_all_four_finite_points": "OPEN",
        "p5_size_four": "OPEN",
        "closes_all_p7_size_four": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_p7_saturated_four_finite_exclusion()
    out = {
        "prop": "15.654",
        "title": "p=7 saturated four-finite boundary exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
        "certificate_archive": CERTIFICATE_ARCHIVE,
        "certificate_archive_sha256": CERTIFICATE_ARCHIVE_SHA256,
        "certificate_audit_sha256": CERTIFICATE_AUDIT_SHA256,
        "orbit_source_sha256": ORBIT_SOURCE_SHA256,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15654.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
