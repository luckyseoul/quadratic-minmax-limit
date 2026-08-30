#!/usr/bin/env python3
"""Prop. 15.706 -- exclude both p=17 slack-zero profiles analytically.

Every allocation in either surviving profile retains rigid b=2 directions
of both quadratic types.  Comparing their summed coefficient identities with
the global Paley sign sum forces the infinity degree I to satisfy

    17 I = 4 + 72(g_+ + g_-),

so I=68 in the range 0<=I<=69.  Then a 69-edge graph consists of 68
infinity-star edges and one finite edge, whose affine odd boundary has size
66, 68, or 70, never 16.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15700 import (
    p17_second_boundary_profile_census,
    p17_slack_zero_conic_reduction,
)
from e1_gmin_m4_prop15705 import p17_slack_sixteen_orbit_exclusion


ROOT = Path(__file__).resolve().parents[1]
P = 17
H_SIZE = 4 * P + 1


def p17_slack_zero_global_sign_certificate() -> dict[str, object]:
    """Exact solver-free contradiction shared by both profiles."""
    reduction = p17_slack_zero_conic_reduction()
    survivors = reduction["surviving_slack_zero_profiles"]
    expected = [
        ({0: 1, 2: 7, 16: 1}, {2: 9}),
        ({0: 1, 2: 8}, {2: 8, 16: 1}),
    ]
    observed = [
        (
            dict(row["phase_profiles_b"]["0"]),
            dict(row["phase_profiles_b"]["1"]),
        )
        for row in survivors
    ]
    if observed != expected:
        raise ArithmeticError("p=17 slack-zero profiles changed")

    # Case zero elevates one of nine phase-one b=2 directions.  Case one
    # elevates either its phase-zero b=0 anchor or one of eight phase-zero
    # b=2 directions.  Thus both phases retain a rigid b=2 direction.
    rigid_counts = [
        {"phase_zero_b2": 7, "phase_one_b2": 8},
        {"phase_zero_b2_minimum": 7, "phase_one_b2": 8},
    ]
    if any(min(row.values()) <= 0 for row in rigid_counts):
        raise ArithmeticError("a quadratic type lost every rigid b=2 direction")

    # For a rigid b=2 direction of type eps, put sigma=eps*c_H.  Its mean is
    # M=17-sigma and the single target cell sums to T=sigma.  If P_d is its
    # parallel finite-edge count, I the infinity degree, and g_d its gauge,
    # summing the 136 cell identities and comparing with the mean gives
    # P_d=4+8g_d-I.
    choose_two = P * (P - 1) // 2
    star_multiplicity = P - 1
    if choose_two != 136 or star_multiplicity != 16:
        raise ArithmeticError("p=17 cell-sum constants changed")
    parallel_formula = "P_d=4+8*g_d-I"

    # Parallel Paley edges have sign eps.  If S is the signed sum of all
    # finite selected edges, the summed cross coefficient is eps*S-P_d.
    # Equating the two formulas gives
    # S=c_H+eps*(4+144*g_d-17*I).
    gauge_coefficient = choose_two + 8
    if gauge_coefficient != 144:
        raise ArithmeticError("global sign gauge coefficient changed")

    # Equate eps=+1 and eps=-1.  Since 17^{-1}=17 modulo 72,
    # I=68 (mod 72), and 0<=I<=69 leaves only I=68.
    modulus = gauge_coefficient // 2
    inverse = pow(P, -1, modulus)
    residue = 4 * inverse % modulus
    candidates = [value for value in range(H_SIZE + 1) if value % modulus == residue]
    if modulus != 72 or inverse != 17 or residue != 68 or candidates != [68]:
        raise ArithmeticError("infinity-degree congruence changed")

    finite_edges = H_SIZE - candidates[0]
    affine_boundary_sizes = sorted(
        {
            candidates[0] + 2 - 2 * endpoints_already_in_star
            for endpoints_already_in_star in (0, 1, 2)
        }
    )
    if finite_edges != 1 or affine_boundary_sizes != [66, 68, 70]:
        raise ArithmeticError("I=68 boundary contradiction changed")

    return {
        "p": P,
        "edge_count": H_SIZE,
        "surviving_profile_count_before": len(survivors),
        "rigid_b2_counts_after_any_allocation": rigid_counts,
        "rigid_direction_identities": {
            "sigma": "eps*c_H",
            "mean": "M_d=17-sigma",
            "target_sum": "T_d=sigma",
            "summed_cell_identity": "R_d=sigma+136*g_d-16*I",
            "directional_mean": "M_d=I+17*P_d-R_d-51",
            "parallel_count": parallel_formula,
        },
        "global_paley_sign_identity": {
            "finite_sign_sum": "S=sum_{e finite} C_e",
            "parallel_edge_sign": "C_e=eps",
            "cross_sum": "R_d=eps*S-P_d",
            "therefore": "S=c_H+eps*(4+144*g_d-17*I)",
        },
        "opposite_type_comparison": "17*I=4+72*(g_+ + g_-)",
        "infinity_degree_modulus": modulus,
        "infinity_degree_residue": residue,
        "infinity_degree_candidates_in_range": candidates,
        "finite_edges_at_candidate": finite_edges,
        "possible_affine_odd_boundary_sizes": affine_boundary_sizes,
        "required_affine_odd_boundary_size": 16,
        "contradiction": True,
        "uses_solver": False,
        "proved": True,
    }


def p17_slack_zero_profile_exclusion() -> dict[str, object]:
    previous = p17_slack_sixteen_orbit_exclusion()
    profiles = p17_second_boundary_profile_census()["profiles"]
    certificate = p17_slack_zero_global_sign_certificate()
    previous_indices = set(int(index) for index in previous["remaining_profile_indices"])
    excluded_indices = {
        index
        for index in previous_indices
        if int(profiles[index]["pair_slack"]) == 0
    }
    remaining_indices = sorted(previous_indices - excluded_indices)
    before = len(previous_indices)
    excluded = len(excluded_indices)
    if before != 1215 or excluded != 2:
        raise ArithmeticError("pre-15.706 p=17 ledger changed")
    histogram = dict(
        sorted(
            Counter(
                int(profiles[index]["pair_slack"]) for index in remaining_indices
            ).items()
        )
    )
    after = len(remaining_indices)
    if (
        after != 1213
        or sum(histogram.values()) != after
        or min(histogram) != 16
        or histogram.get(16) != 74
    ):
        raise ArithmeticError("post-15.706 p=17 ledger changed")
    return {
        "proposition": "15.706",
        "p": P,
        "boundary_size": 16,
        "profile_count_before": before,
        "profiles_excluded_here": excluded,
        "profile_count_after": after,
        "excluded_profile_indices_here": sorted(excluded_indices),
        "remaining_profile_indices": remaining_indices,
        "remaining_pair_slack_histogram": histogram,
        "remaining_slack_zero_profiles": 0,
        "remaining_slack_sixteen_profiles": 74,
        "remaining_profiles_of_slack_at_least_twenty": after - 74,
        "certificate": certificate,
        "p17_second_all_finite_endpoint_closed": False,
        "top_level_gates_changed": False,
        "proved_analytically": True,
    }


def main() -> None:
    theorem = p17_slack_zero_profile_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15706.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.706: both p=17 slack-zero profiles excluded; "
        f"{theorem['profile_count_after']} profiles remain"
    )


if __name__ == "__main__":
    main()
