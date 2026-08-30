#!/usr/bin/env python3
"""Prop. 15.712 -- close the p=17 second all-finite endpoint.

Every profile left by Proposition 15.711 has nine phase-one directions with
``b_d=16``. For a sixteen-point affine boundary this means that each of
those directions has sixteen singleton fibres and is not determined by a
pair of boundary points. Thus the boundary determines at most nine
directions.

Szőnyi's direction theorem says that a noncollinear k-point subset of
``AG(2,p)``, with ``k<=p``, determines at least ``(k+3)/2`` directions.
At ``(k,p)=(16,17)`` the lower bound is ten. Hence the boundary would have
to be collinear, but a sixteen-point affine line has phase-labelled profile
``{0:1,16:8}``, ``{16:9}``, absent from the exact ledger.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15710 import p17_phase_one_b16_global_sign_reduction
from e1_gmin_m4_prop15711 import p17_residue_zero_uniform_mean_exclusion


ROOT = Path(__file__).resolve().parents[1]
P = 17
BOUNDARY_SIZE = 16


def p17_redei_szonyi_direction_endpoint_exclusion() -> dict[str, object]:
    """Exclude all fourteen rows using Szőnyi's sharp direction bound."""
    previous = p17_residue_zero_uniform_mean_exclusion()
    source = p17_phase_one_b16_global_sign_reduction()
    by_index = {
        int(row["census_index"]): row for row in source["surviving_profiles"]
    }
    profiles = [
        by_index[int(index)] for index in previous["remaining_profile_indices"]
    ]
    if len(profiles) != 14 or int(previous["profile_count_after"]) != 14:
        raise ArithmeticError("pre-15.712 p17 ledger changed")
    if any(row["phase_profiles_b"]["1"] != {16: 9} for row in profiles):
        raise ArithmeticError("post-15.711 phase-one b=16 core changed")

    # If b_d=16, the sixteen boundary points occupy sixteen odd fibres.
    # Their total occupancy is also sixteen, so each occupied fibre is a
    # singleton. No pair of boundary points therefore determines d.
    phase_one_nondirections = 9
    total_projective_directions = P + 1
    maximum_determined_directions = total_projective_directions - phase_one_nondirections
    minimum_noncollinear_directions = (BOUNDARY_SIZE + 4) // 2
    if maximum_determined_directions != 9 or minimum_noncollinear_directions != 10:
        raise ArithmeticError("Szőnyi direction-count arithmetic changed")
    if not maximum_determined_directions < minimum_noncollinear_directions:
        raise ArithmeticError("Szőnyi direction contradiction disappeared")

    # A 16-subset of one affine line has b=0 in the line direction and b=16
    # in every other direction. Since all nine phase-one directions have
    # b=16, the line direction would be phase zero.
    collinear_phase_zero_profile = {0: 1, 16: 8}
    collinear_phase_one_profile = {16: 9}
    collinear_profile_present = any(
        row["phase_profiles_b"]["0"] == collinear_phase_zero_profile
        and row["phase_profiles_b"]["1"] == collinear_phase_one_profile
        for row in profiles
    )
    if collinear_profile_present:
        raise ArithmeticError("collinear p17 profile unexpectedly survived")

    return {
        "proposition": "15.712",
        "p": P,
        "boundary_size": BOUNDARY_SIZE,
        "profile_count_before": len(profiles),
        "profiles_excluded_here": len(profiles),
        "profile_count_after": 0,
        "excluded_profile_indices_here": sorted(
            int(index) for index in previous["remaining_profile_indices"]
        ),
        "remaining_profile_indices": [],
        "common_phase_one_profile": {"16": 9},
        "phase_one_b16_directions_are_nondirections": True,
        "phase_one_nondirection_count": phase_one_nondirections,
        "maximum_boundary_direction_count": maximum_determined_directions,
        "szonyi_theorem_hypotheses": {
            "prime_affine_plane": True,
            "point_count_at_most_prime": BOUNDARY_SIZE <= P,
            "noncollinear_direction_lower_bound": "(k+3)/2",
        },
        "minimum_noncollinear_direction_count": minimum_noncollinear_directions,
        "therefore_boundary_collinear": True,
        "collinear_phase_zero_profile": {"0": 1, "16": 8},
        "collinear_phase_one_profile": {"16": 9},
        "collinear_profile_present_in_remainder": collinear_profile_present,
        "direction_theorem_reference": {
            "original": (
                "T. Szőnyi, On the number of directions determined by a set "
                "of points in an affine Galois plane, JCTA 74 (1996), 141-146"
            ),
            "modern_statement": (
                "G. Somlai, A new proof of Rédei's theorem on the number of "
                "directions, Arch. Math. 122 (2024), 575-580"
            ),
            "doi": "10.1007/s00013-024-01979-x",
        },
        "remaining_pair_slack_histogram": {},
        "remaining_residue_pair_histogram": {},
        "p17_second_all_finite_endpoint_closed": True,
        "top_level_gates_changed": False,
        "uses_solver": False,
        "uses_new_arc_classification": False,
        "proved_analytically": True,
    }


def main() -> None:
    theorem = p17_redei_szonyi_direction_endpoint_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15712.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.712: p=17 profiles "
        f"{theorem['profile_count_before']} -> {theorem['profile_count_after']}"
    )


if __name__ == "__main__":
    main()
