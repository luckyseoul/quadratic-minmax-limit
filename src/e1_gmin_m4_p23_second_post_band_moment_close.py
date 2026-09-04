#!/usr/bin/env python3
r"""Close the exceptional ``p=23, k=112`` second post-band endpoint.

Proposition 15.770 advances the ``p=3 (mod 4)`` band to ``t=q-1`` for
``p>=31``.  Its local ``p+13`` exclusion does not include ``p=23``, but the
equality/global-moment certificate at the preceding endpoint does.

At ``p=23,t=10`` the residue ledger has three branches.  The two old exact
endpoint types again force Proposition 15.752's forbidden mass ``p+9``.  The
new mass-``p-1`` lift is Boolean and its density is absent from the fixed
four-bit catalog.  In the carried sharp branch there are eleven low hard rows
and one high row.  The common row sum and coefficient offset leave only the
``P=4,Q=5,F5`` family.  Eleven low triangle-minus-star rows are still more
than the degrees four and eight of the two global moment forms, so the same
all-``binom(23,5)`` K5 certificate gives a contradiction.

This closes one finite endpoint for every boundary size.  It is not a graph
census, a later-layer theorem, or a global residual-(ii) closure.
"""
from __future__ import annotations

import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from e1_gmin_m4_p23_post_band_moment_close import (
    p23_hard_moment_root_certificate,
    p23_k5_moment_sieve,
    p23_sharp_hard_family_catalog,
    p23_slice_half_mean_classification,
)
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15734 import BRANCH_B2, BRANCH_P3_LAST, residual_even_floor_table
from e1_gmin_m4_prop15751 import exact_four_cube_catalog, profile_density
from e1_gmin_m4_prop15752 import p_plus_nine_local_exclusion
from io_atomic import write_json_atomic


ROOT = Path(__file__).resolve().parents[1]
P = 23
Q = 11
M = 12
LAYER_INDEX = 10
ORIGINAL_K = 112
H_EDGE_COUNT = 113
P3_CARRIED_BRANCH = "carried_sharp_p_minus_3"
P3_NEW_LOCAL_BRANCH = "all_low_p_minus_one_lift"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


@lru_cache(maxsize=1)
def p23_second_post_band_residue_ledger() -> dict[str, object]:
    """Classify every phase-one residue at ``p=23,t=10``."""
    floors = {
        int(boundary): int(value)
        for boundary, value in residual_even_floor_table(P)[
            "phase_one_floors"
        ].items()
    }
    lift_floor = int(
        sharp_integral_quadratic_lift_floor(P)["sharp_scaled_floor"]
    )
    expected_live = {
        M - 3: [(2, "sharp_p_minus_3"), (P - 1, "sharp_p_minus_3")],
        M - 2: [(2, "p_minus_one"), (P - 1, "p_minus_one")],
        M - 1: [(2, "exact"), (P - 1, "exact")],
    }

    rows = []
    for residue in range(M):
        quotient_sum = M + LAYER_INDEX - residue
        if quotient_sum >= M:
            low_quotient = 1
            forced_low_count = 2 * M - quotient_sum
            low_mean = P + 1 + 2 * residue
        else:
            low_quotient = 0
            forced_low_count = M - quotient_sum
            low_mean = 2 * residue

        candidates = []
        live = []
        for boundary, floor in floors.items():
            if floor > low_mean:
                continue
            excess = low_mean - floor
            if excess == 0:
                classification = "exact"
            elif excess < lift_floor:
                classification = "excluded_sub_sharp_lift"
            elif excess == lift_floor:
                classification = "sharp_p_minus_3"
            elif residue == M - 2 and excess == P - 1:
                classification = "p_minus_one"
            else:  # pragma: no cover - guarded by the exact floor table
                raise ArithmeticError(
                    f"unclassified p23 excess at u={residue},b={boundary}"
                )
            candidates.append(
                {
                    "b": boundary,
                    "floor": floor,
                    "excess": excess,
                    "classification": classification,
                }
            )
            if classification != "excluded_sub_sharp_lift":
                live.append((boundary, classification))

        _require(
            live == expected_live.get(residue, []) and forced_low_count > 0,
            f"the p23 second-post-band residue u={residue} changed",
        )
        rows.append(
            {
                "u": residue,
                "quotient_sum": quotient_sum,
                "forced_low_quotient": low_quotient,
                "forced_low_direction_count_at_least": forced_low_count,
                "forced_low_mean": low_mean,
                "candidate_rows": candidates,
                "live_rows": [
                    {"b": boundary, "classification": classification}
                    for boundary, classification in live
                ],
                "excluded_arithmetically": not live,
            }
        )

    isolated = P * P + 1 - 2 * H_EDGE_COUNT
    proved = bool(
        H_EDGE_COUNT == 4 * P + 2 * LAYER_INDEX + 1
        and isolated == 304 > 0
        and lift_floor == P - 3 == 20
        and [row["u"] for row in rows if row["live_rows"]]
        == [M - 3, M - 2, M - 1]
    )
    _require(proved, "the p23 second-post-band residue ledger failed")
    return {
        "p": P,
        "q": Q,
        "m": M,
        "layer_index_t": LAYER_INDEX,
        "original_k": ORIGINAL_K,
        "H_edge_count": H_EDGE_COUNT,
        "guaranteed_isolated_vertices": isolated,
        "phase_one_mean_form": "a_L=2u+24*k_L",
        "phase_one_quotient_sum": "sum k_L=22-u",
        "sharp_lift_floor": lift_floor,
        "arithmetic_surviving_residues": [M - 3, M - 2, M - 1],
        "rows": rows,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p23_p_minus_one_local_exclusion() -> dict[str, object]:
    """Exclude the new scaled-mass-22 lift by the fixed cube catalog."""
    lift = sharp_integral_quadratic_lift_floor(P)
    mass = P - 1
    density = Fraction(mass, 4 * P)
    influence_floor = Fraction(
        (P + 1) * (P - 3), 16 * P * (P - 2)
    )
    total_influence_upper = (P - 1) * density * (1 - density)
    junta_bound = Fraction(
        2 * (P - 1) ** 2 * (P - 2) * (3 * P + 1),
        P * P * (P + 1) * (P - 3),
    )
    catalog = exact_four_cube_catalog()
    densities = sorted(
        {
            profile_density(tuple(row["layer_counts"]), P)
            for row in catalog["profiles"]
        }
    )
    proved = bool(
        lift["proved"]
        and int(lift["H_at_least_two_scaled_floor"]) == P + 1 > mass
        and density == Fraction(11, 46)
        and junta_bound == Fraction(5929, 1058) < 6
        and 5 < Q
        and catalog["proved"]
        and Fraction(P - 3, 4 * P) < density < Fraction(P + 1, 4 * P)
        and density not in densities
    )
    _require(proved, "the p23 mass-p-1 local exclusion failed")
    return {
        "p": P,
        "scaled_mass": mass,
        "H_at_least_two_scaled_floor": P + 1,
        "therefore_height_one_boolean": True,
        "density": str(density),
        "relevant_pair_influence_floor": str(influence_floor),
        "total_influence_upper_bound": str(total_influence_upper),
        "largest_zero_influence_class_complement_bound": str(junta_bound),
        "junta_coordinates_at_most": 5,
        "five_less_than_both_slice_sides": 5 < Q,
        "cube_active_coordinates_at_most": 4,
        "four_bit_density_values": [str(value) for value in densities],
        "target_density_absent": True,
        "finite_slice_or_graph_census_used": False,
        "proved": proved,
    }


def _old_endpoint_branch(branch: str) -> dict[str, object]:
    """Replay either old exact branch at edge count 113."""
    if branch == BRANCH_B2:
        hard_parallel = 4
        minimum_Q = 3
    elif branch == BRANCH_P3_LAST:
        hard_parallel = 3
        minimum_Q = 4
    else:  # pragma: no cover - private API
        raise ValueError("unknown old p23 branch")

    hard_mean = P - 1
    h_times_T = (P + 1) * hard_parallel - 3 * P - hard_mean
    hard_edges = (H_EDGE_COUNT + h_times_T) // 2
    opposite_edges = H_EDGE_COUNT - hard_edges
    next_Q = minimum_Q + 1

    def opposite_mean(parallel: int) -> int:
        return (P + 1) * parallel + h_times_T - 3 * P

    surplus = opposite_edges - M * next_Q
    dependency = p_plus_nine_local_exclusion(P)
    proved = bool(
        2 * hard_edges == H_EDGE_COUNT + h_times_T
        and opposite_mean(minimum_Q) == 8
        and opposite_mean(next_Q) == P + 9
        and surplus == 6
        and dependency["proved"]
    )
    _require(proved, f"the old p23 branch {branch} failed")
    return {
        "branch": branch,
        "hard_parallel_count": hard_parallel,
        "hard_sign_times_global_T": h_times_T,
        "hard_edge_count": hard_edges,
        "opposite_edge_count": opposite_edges,
        "forbidden_minimum_Q": minimum_Q,
        "forbidden_minimum_mass": 8,
        "forced_next_Q": next_Q,
        "forced_next_scaled_mass": P + 9,
        "directions_at_next_Q_at_least": M - surplus,
        "dependency": "Proposition 15.752 p+9 local exclusion",
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p23_carried_sharp_moment_exclusion() -> dict[str, object]:
    """Use eleven hard roots to exclude every carried sharp family."""
    catalog = p23_sharp_hard_family_catalog()
    slice_forms = p23_slice_half_mean_classification()
    root_identities = p23_hard_moment_root_certificate()
    sieve = p23_k5_moment_sieve()
    phase_zero = {
        int(boundary): int(value)
        for boundary, value in residual_even_floor_table(P)[
            "phase_zero_floors"
        ].items()
    }
    lift_floor = int(
        sharp_integral_quadratic_lift_floor(P)["sharp_scaled_floor"]
    )

    low_hard_count = M - 1
    family_rows = []
    for family in catalog["hard_families"]:
        offset = int(family["coefficient_offset"])
        parallel_candidates = [
            value
            for value in range(H_EDGE_COUNT // M + 1)
            if (value - offset) % Q == 0
        ]
        _require(parallel_candidates == [offset], "a carried p23 offset changed")
        hard_parallel = offset
        high_parallel = hard_parallel + 1
        hard_edges = low_hard_count * hard_parallel + high_parallel
        opposite_edges = H_EDGE_COUNT - hard_edges
        h_times_T = (P + 1) * hard_parallel - 5 * P + 4

        def opposite_mean(parallel: int) -> int:
            return (P + 1) * parallel + h_times_T - 3 * P

        forbidden_Q = 8 - hard_parallel
        forced_Q = 9 - hard_parallel
        surplus = opposite_edges - M * forced_Q
        nonzero_rows = [
            [boundary, floor, opposite_mean(forced_Q) - floor]
            for boundary, floor in phase_zero.items()
            if boundary and floor <= opposite_mean(forced_Q)
        ]
        compatible_forms = [
            form["name"]
            for form in slice_forms["global_slice_forms"]
            if (forced_Q - int(form["coefficient_offset"])) % Q == 0
        ]
        family_rows.append(
            {
                **family,
                "hard_parallel_candidates": parallel_candidates,
                "low_hard_direction_count": low_hard_count,
                "unique_high_direction_count": 1,
                "forced_high_parallel_count": high_parallel,
                "hard_edge_count": hard_edges,
                "opposite_edge_count": opposite_edges,
                "hard_sign_times_global_T": h_times_T,
                "forbidden_Q": forbidden_Q,
                "forbidden_scaled_mass": opposite_mean(forbidden_Q),
                "forced_low_Q": forced_Q,
                "forced_low_scaled_mass": opposite_mean(forced_Q),
                "surplus_after_every_opposite_Q_at_least_forced_low": surplus,
                "directions_at_forced_low_Q_at_least": M - surplus,
                "nonzero_boundary_floor_and_lift_rows": nonzero_rows,
                "forced_low_Q_is_boundary_zero": True,
                "compatible_slice_forms": compatible_forms,
                "excluded_by_local_offset": not compatible_forms,
            }
        )

    survivors = [row for row in family_rows if row["compatible_slice_forms"]]
    maximum_form_degree = max(
        int(value) for value in root_identities["form_degrees"].values()
    )
    roots_force_zero = low_hard_count > maximum_form_degree
    moment_survivor_excluded = bool(
        len(survivors) == 1
        and int(survivors[0]["coefficient_offset"]) == 4
        and int(survivors[0]["forced_low_Q"]) == 5
        and survivors[0]["compatible_slice_forms"] == ["F5"]
        and roots_force_zero
        and int(sieve["simultaneous_zero_count"]) == 0
    )
    proved = bool(
        catalog["proved"]
        and slice_forms["proved"]
        and root_identities["proved"]
        and sieve["proved"]
        and [int(row["coefficient_offset"]) for row in family_rows]
        == [2, 4, 3, 5]
        and all(int(row["hard_edge_count"]) == M * int(row["coefficient_offset"]) + 1 for row in family_rows)
        and all(int(row["forbidden_scaled_mass"]) == 12 for row in family_rows)
        and all(int(row["forced_low_scaled_mass"]) == 36 for row in family_rows)
        and all(int(row["surplus_after_every_opposite_Q_at_least_forced_low"]) == 4 for row in family_rows)
        and all(
            row["nonzero_boundary_floor_and_lift_rows"]
            == [[2, 24, 12], [22, 24, 12]]
            for row in family_rows
        )
        and lift_floor == 20
        and moment_survivor_excluded
    )
    _require(proved, "the carried p23 moment branch survived")
    return {
        "branch": P3_CARRIED_BRANCH,
        "changed_premise": (
            "the p23,t=9 all-low sharp family becomes eleven low P rows "
            "and one high P+1 row"
        ),
        "family_ledgers": family_rows,
        "unique_survivor_before_moments": {
            "hard_P": 4,
            "opposite_Q": 5,
            "opposite_form": "F5",
        },
        "low_triangle_minus_star_projective_roots": low_hard_count,
        "maximum_common_form_degree": maximum_form_degree,
        "low_roots_force_G4_and_G8_identically_zero": roots_force_zero,
        "opposite_K5_simultaneous_zero_count": sieve[
            "simultaneous_zero_count"
        ],
        "fixed_five_set_certificate_reused": True,
        "excluded": moment_survivor_excluded,
        "proved": proved,
    }


@lru_cache(maxsize=1)
def p23_second_post_band_moment_close() -> dict[str, object]:
    """Package the complete ``p=23,t=10,k=112`` endpoint proof."""
    residues = p23_second_post_band_residue_ledger()
    old_branches = {
        branch: _old_endpoint_branch(branch)
        for branch in (BRANCH_B2, BRANCH_P3_LAST)
    }
    carried = p23_carried_sharp_moment_exclusion()
    p_minus_one = p23_p_minus_one_local_exclusion()
    branches = {
        **old_branches,
        P3_CARRIED_BRANCH: carried,
        P3_NEW_LOCAL_BRANCH: {
            "residue": M - 2,
            "baseline_boundary_values": [2, P - 1],
            "difference_scaled_mass": P - 1,
            "local_exclusion": p_minus_one,
            "proved": p_minus_one["proved"],
        },
    }
    proved = bool(
        residues["proved"]
        and set(branches)
        == {BRANCH_B2, BRANCH_P3_LAST, P3_CARRIED_BRANCH, P3_NEW_LOCAL_BRANCH}
        and all(row["proved"] for row in branches.values())
    )
    _require(proved, "the p23 second post-band endpoint did not close")
    return {
        "title": "p23 second post-band one-row carry and eleven-root close",
        "result_status": (
            "proved endpoint theorem reusing one fixed finite coefficient certificate"
        ),
        "statement": (
            "the residual-(ii) isolated-chart branch at p=23,t=10,k=112 "
            "is empty for every boundary size"
        ),
        "residue_ledger": residues,
        "branch_exclusions": branches,
        "p23_k112_closed": proved,
        "all_boundary_sizes_excluded": proved,
        "new_graph_or_residual_configuration_census_used": False,
        "fixed_p23_five_set_coefficient_certificate_reused": True,
        "later_layers_closed": False,
        "residual_ii_closed_globally": False,
        "E1_closed": False,
        "quadratic_minmax_limit_closed": False,
        "proved": proved,
    }


def write_evidence(path: Path | None = None) -> Path:
    if path is None:
        path = (
            ROOT
            / "evidence"
            / "e1_gmin_m4_p23_second_post_band_moment_close.json"
        )
    write_json_atomic(path, p23_second_post_band_moment_close())
    return path


def main() -> None:
    path = write_evidence()
    print(json.dumps({"proved": True, "wrote": str(path)}, sort_keys=True))


if __name__ == "__main__":
    main()
