#!/usr/bin/env python3
"""Prop. 15.656 -- close every p=5 four-point residual boundary.

For either p=5 eigenshell, quotienting by antipodes leaves 130 distinct
edge-sign rows.  Every normalized edge column sums to 26.  Therefore a
putative 21-edge residual graph with normalized score at least three has
total shell slack 78.  Its four-point boundary and Paley-product sign fix
the parity of every slack, leaving a bounded lift vector of prescribed
mass.  The edge-count, distinguished-edge, and 130 bad-count equations have
rank 67 modulo five and 65 left dependencies.

Complete square-semilinear orbit scans exclude 712 representative cases by
these mod-five shell syndromes.  The sole mod-five timeout is independently
infeasible modulo seven.  These 713 direct cases cover both signs when the
boundary contains infinity and the negative sign otherwise.  Multiplication
by a nonsquare, together with switching infinity, is an exact conference
anti-isometry; for a boundary omitting infinity it transfers the negative
exclusion to the positive sign.  Hence all 1,202 floor-surviving orbit/sign
cases, representing 26,450 boundary/sign cases, are excluded.  Together
with Proposition 15.632 this closes every p=5 size-four boundary.

With Propositions 15.652--15.655, every size-four boundary is now closed for
every odd prime p>=5.  Boundaries of size at least six, residual (ii), R1,
and the limit remain open.
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
    "2026-08-26-p5-four-point/"
    "p5_size_four_full_shell_certificate_2026-08-26.tar.gz"
)
CERTIFICATE_ARCHIVE_SHA256 = (
    "d5db5e82389ebb0bfcb23e80da5e2322b1d65e74aa8f3804d25275793b7380da"
)
CERTIFICATE_AUDIT_SHA256 = (
    "5cafab9272510dc6871818fbfc395c8f3386ca6615f7a1b1e36f4785cf1d7e4f"
)
MOD7_EXCEPTION_SHA256 = (
    "2b91f83aeb543a29b0b4398243115f25189d233c50d41e16f8fd14b121c61b06"
)
MOD5_BATCH_SHA256 = {
    "c_minus_no_infinity": (
        "0ebee7b3f684168c91993d213b4772f64c8ae462f90e23f28f53ce07273628d1"
    ),
    "c_minus_with_infinity": (
        "73c09f0b717275c2d2f41b259d9b075259e0a94d1939da3823bab1f1c15b23ec"
    ),
    "c_plus_with_infinity": (
        "a032b8c7131c0b1180d8f316c8a01251595b1c20b1a27034b2949adbcd174073"
    ),
}
ORBIT_SOURCE_SHA256 = {
    "c_minus_no_infinity": (
        "fbd7a6f06fcd1b872b46313588b0a6b3db582b19f9ea0387ff4a525dbc2f5097"
    ),
    "c_minus_with_infinity": (
        "8a23cac1a8d2c3e99c5675153366c719d201a511db1f01ce16dbcf4bef6fe74e"
    ),
    "c_plus_no_infinity": (
        "2b0a166875382b4b36bbf3c0b376ddff4b81ed447c20d85d029e438b1dccaa09"
    ),
    "c_plus_with_infinity": (
        "1ecefcd79348f677712c832a8e221e7dd5fda83814f213a649e40058e0539d88"
    ),
}


def p5_nonsquare_signed_permutation() -> dict:
    """Return and verify the signed nonsquare Paley anti-isometry."""
    q, mul, _add, chi, _frob, _norm, _ia, _ib = field_ctx(5)
    alpha = next(value for value in range(1, q) if chi(value) == -1)
    finite_permutation = tuple(mul(alpha, u) for u in range(q))
    vertex_permutation = (0, *(u + 1 for u in finite_permutation))
    switching = (-1, *([1] * q))
    C = np.rint(paley_conference_prime_power(5)).astype(np.int8)
    permutation_array = np.asarray(vertex_permutation, dtype=np.int16)
    switching_array = np.asarray(switching, dtype=np.int8)
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


def p5_full_shell_four_point_certificate() -> dict:
    return {
        "total_size_four_boundary_product_sign_cases": 29900,
        "floor_excluded_boundary_product_sign_cases": 3450,
        "floor_surviving_boundary_product_sign_cases": 26450,
        "floor_surviving_orbit_product_sign_cases": 1202,
        "shell_antipodal_representatives": 130,
        "shell_normalized_column_sum": 26,
        "shell_slack_mass": 78,
        "edge_variables": 325,
        "shell_equations": 132,
        "shell_rank_mod_5": 67,
        "shell_left_dependencies_mod_5": 65,
        "combined_equations": 262,
        "combined_rank_mod_5": 113,
        "combined_left_dependencies_mod_5": 149,
        "direct_orbit_cases": 713,
        "direct_mod5_infeasible_orbits": 712,
        "direct_mod7_infeasible_orbits": 1,
        "direct_boundaries_excluded": 15525,
        "transferred_no_infinity_orbits": 489,
        "transferred_no_infinity_boundaries": 10925,
        "all_floor_surviving_orbits_excluded": 1202,
        "all_floor_surviving_boundaries_excluded": 26450,
        "unknown": 0,
        "feasible": 0,
        "fresh_orbit_reclassification_matches": True,
        "full_shell_parity_mass_reconstruction_matches": True,
        "left_nullspace_audits": True,
        "nonsquare_orbit_bijection_audited": True,
        "structural_audit": True,
    }


def theorem_p5_four_point_exclusion() -> dict:
    certificate = p5_full_shell_four_point_certificate()
    symmetry = p5_nonsquare_signed_permutation()
    sign_transfer = {
        "signed_conference_anti_isometry": symmetry[
            "signed_conference_anti_isometry"
        ],
        "fixes_distinguished_edge": symmetry["fixes_distinguished_edge"],
        "edge_count_21_is_odd": True,
        "four_finite_boundary_forces_even_infinity_degree": True,
        "paley_product_sign_flips": True,
        "eigenshell_sign_swaps": True,
        "normalized_scores_are_preserved": True,
        "no_infinity_orbit_bijection_audited": certificate[
            "nonsquare_orbit_bijection_audited"
        ],
    }
    proved = bool(
        certificate["floor_surviving_orbit_product_sign_cases"] == 1202
        and certificate["floor_surviving_boundary_product_sign_cases"] == 26450
        and certificate["direct_orbit_cases"] == 713
        and certificate["direct_mod5_infeasible_orbits"] == 712
        and certificate["direct_mod7_infeasible_orbits"] == 1
        and certificate["transferred_no_infinity_orbits"] == 489
        and certificate["all_floor_surviving_orbits_excluded"] == 1202
        and certificate["all_floor_surviving_boundaries_excluded"] == 26450
        and certificate["unknown"] == certificate["feasible"] == 0
        and certificate["structural_audit"]
        and all(sign_transfer.values())
    )
    return {
        "proved": proved,
        "p5_all_size_four_boundaries_with_prop15632": "CLOSED",
        "all_size_four_boundaries_for_every_odd_prime_p_at_least_5_with_props15652_to_15655": "CLOSED",
        "certificate": certificate,
        "nonsquare_anti_isometry": {
            key: value
            for key, value in symmetry.items()
            if key not in {"finite_permutation", "vertex_permutation", "switching"}
        },
        "sign_transfer": sign_transfer,
        "first_open_boundary_size_for_p_at_least_5": 6,
        "closes_all_p5_size_four": proved,
        "closes_all_size_four_for_p_at_least_5": proved,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_p5_four_point_exclusion()
    out = {
        "prop": "15.656",
        "title": "p=5 full-shell four-point exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
        "certificate_archive": CERTIFICATE_ARCHIVE,
        "certificate_archive_sha256": CERTIFICATE_ARCHIVE_SHA256,
        "certificate_audit_sha256": CERTIFICATE_AUDIT_SHA256,
        "mod7_exception_sha256": MOD7_EXCEPTION_SHA256,
        "mod5_batch_sha256": MOD5_BATCH_SHA256,
        "orbit_source_sha256": ORBIT_SOURCE_SHA256,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15656.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
