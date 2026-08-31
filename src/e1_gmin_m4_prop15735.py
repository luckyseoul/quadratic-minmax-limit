#!/usr/bin/env python3
r"""Prop. 15.735 -- exclude the first three residual-(ii) shells.

Proposition 15.734 excludes the critical shell ``k=4p`` for every prime
``p>=13`` and every boundary size.  The same isolated-chart argument also
closes the next two even shells

    k=4p+2t,                 t in {1,2}.

Indeed, after adjoining the distinguished repair edge the odd flip graph
``H`` has

    |H|=4p+2t+1.

It has at most ``2|H|`` nonisolated vertices.  Even in the worst case
``(p,t)=(13,2)``, the projective line therefore contains at least

    p^2+1-2|H| = p^2-8p-4t-1 = 56

isolated vertices.  Signed PSL transport sends one to infinity, preserves
both separator inequalities, and gives ``I=0`` with the whole odd-degree
boundary finite.  Handshake then makes every directional odd-fibre count
``b_d`` even.

Put ``q=(p-1)/2`` and ``m=q+1``.  The exact budget in either quadratic type
is now

    2m(m+t).

The parity phase alternates with ``t`` (the phase-one type is
``eps_d=(-1)^t c_H``), but one type is still phase one and the other phase
zero.  In the phase-one type write

    a_d=2u+2m k_d,          sum_d k_d=m+t-u.          (1)

For ``t+1<=u<=m-2``, every direction needs ``k_d>=1`` while the right side
of (1) is below ``m``.  For ``0<=u<=t``, at least one direction has
``k_d=1`` and mean ``p+1+2u``.  The only even-``b`` cells below that mean
are ``b=2`` and ``b=p-1``.  Relative to their explicit parity baselines,
every such cell is a nonzero two-, four-, or six-unit integral quadratic
lift, except ``u=0, p=1 mod 4, b=p-1``.  Proposition 15.688 forbids all
those lifts because ``6<p-3`` for ``p>=13``.  Finally, at ``u=m-1`` there
are ``t+1`` quotient units, so there is a low direction of mean ``p-1``.
The two low cells in the ``p=3 mod 4`` case cannot mix: their equal means
give the same parallel count, while their coefficient offsets are 4 and 3.

Consequently the same three baselines as Proposition 15.734 are exhaustive:

* ``A``: phase-one ``b=2``, target ``4+z_i z_j``, offset 4;
* ``B``: when ``p=1 mod 4``, ``u=0`` and ``b=p-1``, target ``4+z_j``,
  offset 5;
* ``C``: when ``p=3 mod 4``, low ``b=p-1``, target ``4-z_j``, offset 3.

If ``P`` is their common low parallel count, coefficient comparison gives
``q | I+P-C`` for offset ``C``.  With ``I=0``, nonnegative quotient and
opposite-edge count force respectively

    (P,rho,s)=(4,0,4), (5,0,5), (3,0,3),

where ``rho=(P-C)/q`` and ``s=P+rho``.  The hard finite-edge totals are
``mP+t+1`` in A/C and ``mP+t`` in B.  Thus the opposite parallel-count
totals are

    A: 4q+t,          B: 3q+t,          C: 5q+t+1.   (2)

Nonnegativity of an opposite-direction mean forces respectively
``Q>=3,2,4``.  The surplus in (2) above these common minima is
``q+t-3, q+t-2, q+t-3``.  For ``t<=2`` each is nonnegative and strictly
less than ``m=q+1``.  Hence some opposite direction attains the minimum and
has scaled mean ``8,6,8``.  Its phase is zero.  A nonempty even odd-fibre
set costs at least ``p-1>=12``; if it is empty, ``A_d=2B_d`` for a nonzero
nonnegative integral quadratic, and Proposition 15.688 instead costs at
least ``p-3>=10``.  Both cases contradict the displayed mean.

Together with Proposition 15.734, the shells
``k in {4p,4p+2,4p+4}`` are therefore empty for every prime ``p>=13`` and
every boundary size.  The strict surplus argument for branch B stops at
``t=3``; no claim is made here for ``k>=4p+6`` or for ``p<=11``.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15721 import is_prime, signed_relative_flip_transport
from e1_gmin_m4_prop15734 import (
    BRANCH_B2,
    BRANCH_P1_LAST,
    BRANCH_P3_LAST,
    baseline_coefficient_rules,
    critical_residual_exclusion,
    residual_even_floor_table,
)


ROOT = Path(__file__).resolve().parents[1]
NEW_LAYER_INDICES = (1, 2)


def _check_parameters(p: int, t: int) -> None:
    """Validate the uniform theorem range and one of the two new layers."""
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 13
        or not is_prime(p)
    ):
        raise ValueError("need an odd prime p>=13")
    if not isinstance(t, int) or isinstance(t, bool) or t not in NEW_LAYER_INDICES:
        raise ValueError("need layer index t in {1,2}")


def residual_layer_arithmetic(p: int, t: int) -> dict[str, object]:
    """Return the exact shell, type budget, and alternating phase data."""
    _check_parameters(p, t)
    q = (p - 1) // 2
    m = q + 1
    original_k = 4 * p + 2 * t
    edge_count = original_k + 1
    type_budget = 2 * m * (m + t)
    phase_exponent = 2 * p + t - 1
    hard_type = "eps_d=c_H" if t % 2 == 0 else "eps_d=-c_H"
    return {
        "p": p,
        "layer_index_t": t,
        "original_k": original_k,
        "H_edge_count": edge_count,
        "q": q,
        "m": m,
        "ambient_vertex_count": p * p + 1,
        "type_budget": type_budget,
        "type_budget_formula": "2m(m+t)",
        "phase_exponent": phase_exponent,
        "phase_exponent_parity": phase_exponent % 2,
        "hard_type": hard_type,
        "hard_type_formula": "eps_d=(-1)^t*c_H",
        "hard_phase": 1,
        "opposite_phase": 0,
        "proved": bool(
            original_k == 4 * p + 2 * t
            and edge_count == 4 * p + 2 * t + 1
            and type_budget
            == (p + 1) * (edge_count - 3 * p) // 2
            and phase_exponent == (edge_count - 3) // 2
        ),
    }


def isolated_layer_chart(p: int, t: int) -> dict[str, object]:
    """Transport an isolated vertex to infinity at a new residual layer."""
    arithmetic = residual_layer_arithmetic(p, t)
    transport = signed_relative_flip_transport()
    maximum_nonisolated = 2 * int(arithmetic["H_edge_count"])
    guaranteed_isolated = int(arithmetic["ambient_vertex_count"]) - maximum_nonisolated
    checks = {
        "support_bound": maximum_nonisolated == 8 * p + 4 * t + 2,
        "isolated_count": guaranteed_isolated == p * p - 8 * p - 4 * t - 1,
        "isolated_gap_positive": guaranteed_isolated > 0,
        "worst_case_gap_at_least_56": guaranteed_isolated >= 56,
        "signed_transport": bool(transport["proved"]),
        "flip_size_preserved": bool(transport["flip_set_size_preserved"]),
        "boundary_permuted": bool(transport["odd_degree_boundary_is_permuted"]),
        "separation_preserved": bool(
            transport["both_separation_inequalities_preserved"]
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError("the new-layer isolated chart dependency failed")
    return {
        "p": p,
        "layer_index_t": t,
        "ambient_vertex_count": p * p + 1,
        "maximum_nonisolated_vertices": maximum_nonisolated,
        "guaranteed_isolated_vertices": guaranteed_isolated,
        "chosen_vertex_is_isolated_and_outside_D": True,
        "transported_boundary_is_all_finite": True,
        "transported_boundary_size_is_even_by_handshake": True,
        "every_transported_directional_b_is_even": True,
        "transported_infinity_degree_I": 0,
        "transported_H_edge_count": 4 * p + 2 * t + 1,
        "transported_separator_inequalities_preserved": True,
        "signed_transport_dependency": transport,
        "checks": checks,
        "proved": True,
    }


def _low_cell_rows(p: int, u: int) -> list[dict[str, object]]:
    """Audit every even phase-one floor at the low mean ``p+1+2u``."""
    floors = residual_even_floor_table(p)["phase_one_floors"]
    lift_floor = int(sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"])
    low_mean = p + 1 + 2 * u
    rows: list[dict[str, object]] = []
    for raw_b, raw_floor in floors.items():
        b = int(raw_b)
        floor = int(raw_floor)
        if floor > low_mean:
            continue
        excess = low_mean - floor
        explicit_baseline = b in (2, p - 1)
        exact = excess == 0
        forbidden_lift = bool(explicit_baseline and 0 < excess < lift_floor)
        rows.append(
            {
                "b": b,
                "floor": floor,
                "low_mean": low_mean,
                "excess_above_explicit_parity_baseline": excess,
                "explicit_baseline_available": explicit_baseline,
                "exact_baseline": exact,
                "forbidden_nonzero_integral_lift": forbidden_lift,
                "survives": exact or not forbidden_lift,
            }
        )
    return rows


def layer_hard_residue_ledger(p: int, t: int) -> dict[str, object]:
    """Classify all common residues in the new shell's phase-one type."""
    _check_parameters(p, t)
    q = (p - 1) // 2
    m = q + 1
    lift_floor = int(sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"])
    rows: list[dict[str, object]] = []
    for u in range(m):
        quotient_sum = m + t - u
        if u == m - 1:
            low_count_lower_bound = m - (t + 1)
            endpoint_candidates = list(
                residual_even_floor_table(p)["phase_one_cells_at_mean_p_minus_one"]
            )
            excluded = False
            branch = "u=m-1 low-baseline dichotomy"
            reason = (
                "sum k=t+1<m, so at least m-(t+1) directions have k=0 "
                "and mean p-1"
            )
            low_rows: list[dict[str, object]] = []
        elif u > t:
            low_count_lower_bound = 0
            endpoint_candidates = []
            excluded = True
            branch = None
            reason = "every direction needs k>=1 but sum k=m+t-u<m"
            low_rows = []
        else:
            extra_quotient_units = t - u
            low_count_lower_bound = m - extra_quotient_units
            endpoint_candidates = []
            low_rows = _low_cell_rows(p, u)
            survivors = [row for row in low_rows if row["survives"]]
            if p % 4 == 1 and u == 0:
                expected = [(p - 1, p + 1, 0)]
                got = [
                    (
                        int(row["b"]),
                        int(row["floor"]),
                        int(row["excess_above_explicit_parity_baseline"]),
                    )
                    for row in survivors
                ]
                if got != expected:
                    raise ArithmeticError("the exact p=1 mod 4 low branch changed")
                excluded = False
                branch = BRANCH_P1_LAST
                reason = "only the exact b=p-1 linear baseline survives"
            else:
                if survivors:
                    raise ArithmeticError("a forbidden low-residue lift survived")
                excluded = True
                branch = None
                reason = "every low cell is a nonzero lift below p-3"
        rows.append(
            {
                "u": u,
                "common_residue": 2 * u,
                "quotient_sum": quotient_sum,
                "low_direction_count_lower_bound": low_count_lower_bound,
                "low_cell_rows": low_rows,
                "endpoint_low_b_candidates": endpoint_candidates,
                "excluded": excluded,
                "surviving_branch": branch,
                "reason": reason,
            }
        )

    feasible_u = [int(row["u"]) for row in rows if not row["excluded"]]
    expected_u = [0, m - 1] if p % 4 == 1 else [m - 1]
    endpoint_candidates = list(
        residual_even_floor_table(p)["phase_one_cells_at_mean_p_minus_one"]
    )
    if p % 4 == 1:
        possible_branches = [BRANCH_B2, BRANCH_P1_LAST]
        endpoint_homogeneous = endpoint_candidates == [2]
    else:
        possible_branches = [BRANCH_B2, BRANCH_P3_LAST]
        offsets = baseline_coefficient_rules(p)
        endpoint_homogeneous = bool(
            endpoint_candidates == [2, p - 1]
            and int(offsets[BRANCH_B2]["offset"])
            - int(offsets[BRANCH_P3_LAST]["offset"])
            == 1
        )
    proved = bool(
        feasible_u == expected_u
        and endpoint_homogeneous
        and m - (t + 1) > 0
        and 2 * t + 2 <= 6 < lift_floor
    )
    if not proved:
        raise ArithmeticError("the new-layer hard residue classification failed")
    return {
        "p": p,
        "layer_index_t": t,
        "same_type_mean_form": f"a_d=2u+{p + 1}*k_d",
        "quotient_identity": "sum_d k_d=m+t-u",
        "residue_rows": rows,
        "feasible_u": feasible_u,
        "low_lift_excess_upper_bound": 2 * t + 2,
        "nonzero_integral_lift_floor": lift_floor,
        "endpoint_low_direction_count_lower_bound": m - (t + 1),
        "endpoint_low_b_candidates": endpoint_candidates,
        "equal_mean_endpoint_cells_cannot_mix": endpoint_homogeneous,
        "possible_branches": possible_branches,
        "proved": proved,
    }


def layer_branch_exclusion(p: int, t: int, branch: str) -> dict[str, object]:
    """Apply the coefficient residue and opposite small-mean contradiction."""
    _check_parameters(p, t)
    allowed = layer_hard_residue_ledger(p, t)["possible_branches"]
    if not isinstance(branch, str) or branch not in allowed:
        raise ValueError(f"branch must be one of {allowed}")
    chart = isolated_layer_chart(p, t)
    q = (p - 1) // 2
    m = q + 1
    rules = baseline_coefficient_rules(p)
    if branch == BRANCH_B2:
        offset = 4
        hard_edge_increment = t + 1
        opposite_edge_delta = t
        mean_constant = 9
        expected_P = 4
        minimum_Q = 3
        expected_mean = 8
        hard_edge_formula = "mP+t+1"
    elif branch == BRANCH_P1_LAST:
        offset = 5
        hard_edge_increment = t
        opposite_edge_delta = t
        mean_constant = 9
        expected_P = 5
        minimum_Q = 2
        expected_mean = 6
        hard_edge_formula = "mP+t"
    else:
        offset = 3
        hard_edge_increment = t + 1
        opposite_edge_delta = t + 1
        mean_constant = 7
        expected_P = 3
        minimum_Q = 4
        expected_mean = 8
        hard_edge_formula = "mP+t+1"

    parameter_rows: list[dict[str, object]] = []
    for parallel_count in range(9):
        numerator = parallel_count - offset
        congruence = numerator % q == 0
        rho = numerator // q if congruence else None
        rho_nonnegative = rho is not None and rho >= 0
        s = int(rho) + parallel_count if rho_nonnegative else None
        opposite_edges = (
            q * (8 - int(s)) + opposite_edge_delta if s is not None else None
        )
        opposite_nonnegative = opposite_edges is not None and opposite_edges >= 0
        feasible = congruence and rho_nonnegative and opposite_nonnegative
        parameter_rows.append(
            {
                "P": parallel_count,
                "coefficient_numerator_I_plus_P_minus_offset": numerator,
                "coefficient_congruence_holds": congruence,
                "rho": rho,
                "rho_nonnegative": rho_nonnegative,
                "s": s,
                "opposite_finite_edge_count": opposite_edges,
                "opposite_edge_nonnegative": opposite_nonnegative,
                "feasible": feasible,
            }
        )
    feasible_rows = [row for row in parameter_rows if row["feasible"]]
    if [(row["P"], row["rho"], row["s"]) for row in feasible_rows] != [
        (expected_P, 0, expected_P)
    ]:
        raise ArithmeticError("I=0 no longer fixes the new-layer branch parameter")

    s = expected_P
    hard_edges = m * expected_P + hard_edge_increment
    total_edges = 4 * p + 2 * t + 1
    opposite_edges = q * (8 - s) + opposite_edge_delta
    if hard_edges + opposite_edges != total_edges:
        raise ArithmeticError("hard and opposite finite-edge counts do not close")
    previous_mean = (
        (p - 1) * s
        + (p + 1) * (minimum_Q - 1)
        + mean_constant
        - 7 * p
    )
    minimum_mean = (
        (p - 1) * s
        + (p + 1) * minimum_Q
        + mean_constant
        - 7 * p
    )
    surplus = opposite_edges - m * minimum_Q
    forced_minimum_direction = 0 <= surplus < m
    floors = residual_even_floor_table(p)
    least_nonzero_floor = int(floors["least_nonzero_phase_zero_floor"])
    lift_floor = int(sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"])
    nonzero_b_excluded = minimum_mean < least_nonzero_floor
    b0_lift_excluded = 0 < minimum_mean < lift_floor
    excluded = bool(
        int(chart["transported_infinity_degree_I"]) == 0
        and previous_mean < 0
        and forced_minimum_direction
        and minimum_mean == expected_mean
        and nonzero_b_excluded
        and b0_lift_excluded
    )
    if not excluded:
        raise ArithmeticError("a new-layer isolated branch survived")
    return {
        "p": p,
        "layer_index_t": t,
        "branch": branch,
        "branch_baseline": rules[branch],
        "transported_infinity_degree_I": 0,
        "coefficient_offset": offset,
        "parameter_rows_P_0_through_8": parameter_rows,
        "forced_P": expected_P,
        "forced_rho": 0,
        "forced_s": s,
        "hard_finite_edge_count_formula": hard_edge_formula,
        "hard_finite_edge_count": hard_edges,
        "opposite_finite_edge_count": opposite_edges,
        "opposite_parallel_count_sum": opposite_edges,
        "opposite_mean_formula": f"a=(p-1)s+(p+1)Q+{mean_constant}-7p",
        "minimum_parallel_count": minimum_Q,
        "mean_one_parallel_step_below": previous_mean,
        "minimum_direction_mean": minimum_mean,
        "parallel_surplus_above_minimum": surplus,
        "opposite_direction_count": m,
        "surplus_less_than_opposite_direction_count": surplus < m,
        "a_minimum_direction_is_forced": forced_minimum_direction,
        "least_nonzero_phase_zero_b_floor": least_nonzero_floor,
        "nonzero_integral_lift_floor": lift_floor,
        "nonzero_b_excluded": nonzero_b_excluded,
        "b0_positive_lift_excluded": b0_lift_excluded,
        "branch_excluded": True,
        "proved": True,
    }


def residual_layer_exclusion(p: int, t: int) -> dict[str, object]:
    """Exclude one of ``k=4p+2`` and ``k=4p+4`` for all boundaries."""
    arithmetic = residual_layer_arithmetic(p, t)
    chart = isolated_layer_chart(p, t)
    residues = layer_hard_residue_ledger(p, t)
    branches = {
        branch: layer_branch_exclusion(p, t, branch)
        for branch in residues["possible_branches"]
    }
    proved = bool(
        arithmetic["proved"]
        and chart["proved"]
        and residues["proved"]
        and branches
        and all(row["proved"] for row in branches.values())
    )
    return {
        "p": p,
        "layer_index_t": t,
        "original_k": 4 * p + 2 * t,
        "H_edge_count": 4 * p + 2 * t + 1,
        "arithmetic": arithmetic,
        "isolated_outside_chart": chart,
        "hard_residue_branches": residues,
        "branch_exclusions": branches,
        "boundary_size_hypothesis_used": False,
        "all_boundary_sizes_excluded": proved,
        "finite_configuration_search_used": False,
        "residual_ii_layer_excluded": proved,
        "result_status": "proved theorem",
        "proved": proved,
    }


def three_layer_uniform_schema() -> dict[str, object]:
    """Record the threshold inequalities for ``t=0,1,2``."""
    p = 13
    maximum_t = 2
    q = (p - 1) // 2
    m = q + 1
    return {
        "prime_range": "odd primes p>=13",
        "closed_layer_indices_t": [0, 1, 2],
        "closed_even_k": ["4p", "4p+2", "4p+4"],
        "maximum_H_size": "4p+5",
        "worst_isolated_gap": p * p - 8 * p - 4 * maximum_t - 1,
        "minimum_q": q,
        "minimum_m": m,
        "maximum_low_lift_excess": 2 * maximum_t + 2,
        "minimum_integral_lift_floor": p - 3,
        "minimum_endpoint_low_direction_count": m - (maximum_t + 1),
        "branch_A_C_maximum_surplus": q + maximum_t - 3,
        "branch_B_maximum_surplus": q + maximum_t - 2,
        "opposite_direction_count_at_threshold": m,
        "next_layer_not_claimed": (
            "at t=3 the branch-B surplus equals the number of opposite "
            "directions, so the minimum-Q direction is no longer forced"
        ),
        "proved": bool(
            p * p - 8 * p - 4 * maximum_t - 1 > 0
            and 2 * maximum_t + 2 < p - 3
            and m - (maximum_t + 1) > 0
            and q + maximum_t - 3 < m
            and q + maximum_t - 2 < m
        ),
    }


def proposition_15735() -> dict[str, object]:
    """Package Proposition 15.734 with the two new shell exclusions."""
    schema = three_layer_uniform_schema()
    sample_primes = (13, 17, 19, 23, 29, 31, 37, 43)
    critical_rows = {
        str(p): critical_residual_exclusion(p) for p in sample_primes
    }
    new_rows = {
        str(p): {
            str(t): residual_layer_exclusion(p, t) for t in NEW_LAYER_INDICES
        }
        for p in sample_primes
    }
    proved = bool(
        schema["proved"]
        and all(row["proved"] for row in critical_rows.values())
        and all(
            row["proved"]
            for by_layer in new_rows.values()
            for row in by_layer.values()
        )
    )
    return {
        "prop": "15.735",
        "title": "First-three-shell isolated-chart residual exclusion",
        "result_status": "proved theorem",
        "statement": (
            "residual (ii) at k in {4p,4p+2,4p+4} is empty for every "
            "prime p>=13, independently of boundary size"
        ),
        "universal_proof_schema": schema,
        "prop_15734_critical_layer_dependencies": critical_rows,
        "representative_new_layer_ledgers": new_rows,
        "first_three_even_residual_layers_empty_p_ge_13": proved,
        "all_boundary_sizes_in_first_three_layers_closed_p_ge_13": proved,
        "k_eq_4p_plus_2_closed_p_ge_13": proved,
        "k_eq_4p_plus_4_closed_p_ge_13": proved,
        "p_at_most_11_closed": False,
        "k_at_least_4p_plus_6_closed": False,
        "residual_ii_k_ge_4p_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "finite_configuration_search_used": False,
        "remaining_scope": (
            "these shells at p=5,7,11; k>=4p+6; multi-level Type I; "
            "and the quadratic-minmax limit"
        ),
        "proved": proved,
    }


def write_evidence() -> Path:
    """Write the deterministic symbolic certificate when explicitly run."""
    output = ROOT / "evidence" / "e1_gmin_m4_prop15735.json"
    payload = json.dumps(proposition_15735(), indent=2, sort_keys=True) + "\n"
    output.write_text(payload)
    return output


def main() -> None:
    result = proposition_15735()
    if not result["proved"]:
        raise ArithmeticError("Proposition 15.735 shell audit failed")
    path = write_evidence()
    print("Prop 15.735 k in {4p,4p+2,4p+4}, p>=13: excluded")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
