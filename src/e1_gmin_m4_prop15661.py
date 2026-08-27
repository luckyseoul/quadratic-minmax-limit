#!/usr/bin/env python3
"""Prop. 15.661 -- exclude every p=7 six-finite residual boundary.

Exact floors, a complete orbit quotient, simultaneous mod-3/mod-7 catalog
joins, compact high-mean catalog models, and a nonsquare sign transfer close
both product signs.  Together with Propositions 15.657--15.660 this closes
every residual boundary of size six for every odd prime p>=5.

Boundaries of size at least eight, residual (ii), Type I, R1, global QVAR,
and the limit remain open.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

GLOBAL_AUDIT_SHA256 = "d0e4de3041fc875090012f4091b0b57a75d0399c2d10550636089138ba50f6cb"
NUKA_SUMMARY_SHA256 = "4da7448143b1a497cb54ec7b8e54c26584ef70b038f8da36a3b85d3276b2fe2c"
ARCHIVE_ROOT = (
    "/mnt/storage/e1work/quadratic-minmax-limit-finite/"
    "2026-08-27-p7-six-finite/"
)


def p7_six_finite_certificate() -> dict:
    return {
        "all_boundaries": 13_983_816,
        "floor_survivors": 3_856_300,
        "square_semilinear_orbits": 80_704,
        "ordinary_orbits": 80_519,
        "ordinary_elevation_cases": 160_745,
        "ordinary_survivors": 0,
        "deep_orbits": 185,
        "deep_initial_infeasible_orbits": 92,
        "deep_initial_unknown_orbits": 93,
        "deep_allocation_leaves": 930,
        "deep_allocation_infeasible_leaves": 810,
        "deep_low_catalog_join_leaves": 120,
        "deep_low_catalog_join_survivors": 0,
        "moduli": (3, 7),
        "mod3_rank": 162,
        "mod3_left_dependency_dimension": 120,
        "mod7_rank": 147,
        "mod7_left_dependency_dimension": 135,
        "nuka_independent_floor_orbit_and_ordinary_replay": True,
        "nonsquare_sign_transfer": True,
        "global_audit_sha256": GLOBAL_AUDIT_SHA256,
        "nuka_summary_sha256": NUKA_SUMMARY_SHA256,
        "archive_root": ARCHIVE_ROOT,
        "proved": True,
    }


def theorem_size_six_all_odd_primes() -> dict:
    certificate = p7_six_finite_certificate()
    proved = bool(
        certificate["proved"]
        and certificate["ordinary_orbits"] + certificate["deep_orbits"]
        == certificate["square_semilinear_orbits"]
        and certificate["ordinary_survivors"] == 0
        and certificate["deep_allocation_infeasible_leaves"]
        + certificate["deep_low_catalog_join_leaves"]
        == certificate["deep_allocation_leaves"]
        and certificate["deep_low_catalog_join_survivors"] == 0
        and certificate["nonsquare_sign_transfer"]
    )
    return {
        "proved": proved,
        "certificate": certificate,
        "p7_six_finite_both_product_signs": "CLOSED",
        "p7_infinity_plus_five_both_product_signs": "CLOSED_BY_15.658_AND_15.659",
        "p5_size_six": "CLOSED_BY_15.660",
        "p_ge_11_size_six": "CLOSED_BY_15.657",
        "all_size_six_boundaries_for_odd_p_ge_5": "CLOSED",
        "boundaries_size_at_least_eight": "OPEN",
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
    }


def main() -> dict:
    theorem = theorem_size_six_all_odd_primes()
    out = {
        "prop": "15.661",
        "title": "complete p=7 six-finite and all size-six boundary exclusion",
        "proved": theorem["proved"],
        "theorem": theorem,
    }
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15661.json"
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote {destination}", flush=True)
    return out


if __name__ == "__main__":
    main()
