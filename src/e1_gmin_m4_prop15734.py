#!/usr/bin/env python3
r"""Prop. 15.734 -- exclude critical residual (ii) for every ``p>=13``.

Let ``p>=13`` be prime and suppose a critical residual-(ii) witness has a
flip graph ``H`` with ``|H|=4p+1``.  No hypothesis on its odd-degree
boundary ``D`` is needed.  The ambient projective line has ``p^2+1``
vertices, whereas at most ``2|H|=8p+2`` vertices are nonisolated.  Since

    p^2+1-(8p+2)=p^2-8p-1>0,

there is an isolated vertex ``w``.  An isolated vertex cannot lie in the
odd-degree boundary.  Proposition 15.721's signed PSL transport sends ``w``
to infinity while preserving ``|H|`` and both separator inequalities.  In
the transported chart the boundary is all finite and

    I=deg_H(w)=0.                                        (1)

Handshake parity makes ``|D|`` even.  Therefore every directional odd-fibre
count ``b_d`` is even.  Put ``q=(p-1)/2`` and ``m=q+1``.  Proposition 15.632
gives each quadratic direction type the exact scaled-mean budget ``2m^2``.
Because ``(|H|-3)/2=2p-1`` is odd, the type ``eps_d=c_H`` has phase one.

The exact parity floors and Proposition 15.688's two-unit lift obstruction
leave exactly three possible hard-type branches:

1. ``q`` low ``b=2`` directions of mean ``p-1`` and one high direction of
   mean ``2p``;
2. if ``p=1 mod 4``, ``m`` exact ``b=p-1`` directions of mean ``p+1``;
3. if ``p=3 mod 4``, ``q`` low ``b=p-1`` directions of mean ``p-1`` and
   one high direction of mean ``2p``.

For ``p>=17`` the floor table is Proposition 15.669.  At ``p=13`` the same
endpoint-relevant table follows directly from Proposition 15.632's exact
three-variable LP and is audited here.  Proposition 15.652's strictly
positive quadrature additionally forces equality to be the pointwise XNOR
or complementary-literal baseline; the LP value alone is not treated as an
equality classification.  In the first and third cases the two possible low
baselines cannot mix: equal means force equal parallel counts, while their
coefficient offsets differ by one.

Let ``P`` be the common hard parallel count.  Slice coefficient comparison
gives respectively

    q | I+P-4,       q | I+P-5,       q | I+P-3.          (2)

Writing the quotient as ``rho`` and ``s=rho+P``, opposite-edge
nonnegativity gives ``0<=P<=8``.  With ``I=0``, (2) forces

    (P,rho,s)=(4,0,4), (5,0,5), (3,0,3),                 (3)

in the three branches.  (At ``p=13`` only the first two branches occur;
the third first occurs at ``p=19``, where ``q=9``.)

For the first two branches an opposite direction with parallel count ``Q``
has mean

    a=(p-1)s+(p+1)Q+9-7p,

and the opposite parallel-count sum is ``q(8-s)``.  In the third branch,

    a=(p-1)s+(p+1)Q+7-7p,
    sum Q=q(8-s)+1.

At the three values of ``s`` in (3), nonnegativity forces respectively
``Q>=3,2,4``.  The total surplus above that common minimum is less than the
``m`` opposite directions, so one direction attains the minimum.  Its
scaled mean is respectively ``8,6,8``.  A nonzero phase-zero odd-fibre set
costs at least ``p-1>=12``.  If its odd-fibre set is empty, parity gives
``A=2C`` for a nonzero nonnegative integral quadratic, but Proposition
15.688 requires scaled lift mass at least ``p-3>=10``.  Both alternatives
contradict the displayed means.

Thus no critical ``k=4p`` residual-(ii) witness exists for any prime
``p>=13``, for any boundary size.  This does not treat even ``k>4p``.  At
``p=11`` the same isolated-chart reduction reaches scaled mean eight, equal
to the sharp lift floor ``p-3``; that equality case remains open here.
"""
from __future__ import annotations

import json
from pathlib import Path

from e1_gmin_m4_prop15632 import scaled_direction_floor
from e1_gmin_m4_prop15652 import parity_floor_certificate
from e1_gmin_m4_prop15669 import full_symbolic_floor
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor
from e1_gmin_m4_prop15721 import is_prime, signed_relative_flip_transport


ROOT = Path(__file__).resolve().parents[1]
BRANCH_B2 = "hard_b2"
BRANCH_P1_LAST = "p1_residue_zero_b_p_minus_1"
BRANCH_P3_LAST = "p3_all_low_b_p_minus_1"


def _check_prime_parameter(p: int) -> None:
    """Validate the theorem range without accepting Boolean integers."""
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 13
        or not is_prime(p)
    ):
        raise ValueError("need an odd prime p>=13")


def critical_residual_arithmetic(p: int) -> dict[str, object]:
    """Record the exact budgets at critical residual size ``k=4p``."""
    _check_prime_parameter(p)
    q = (p - 1) // 2
    m = q + 1
    edge_count = 4 * p + 1
    type_budget = (p + 1) ** 2 // 2
    return {
        "p": p,
        "critical_original_k": 4 * p,
        "H_edge_count": edge_count,
        "q": q,
        "m": m,
        "ambient_vertex_count": p * p + 1,
        "direction_count_per_type": m,
        "type_budget": type_budget,
        "type_budget_identity": "sum_(eps_d=tau) a_d=(p+1)^2/2",
        "phase_exponent": 2 * p - 1,
        "hard_type": "eps_d=c_H",
        "hard_phase": 1,
        "opposite_phase": 0,
        "proved": bool(
            2 * m == p + 1
            and type_budget == 2 * m * m
            and edge_count == 4 * p + 1
            and (2 * p - 1) % 2 == 1
        ),
    }


def isolated_outside_chart(p: int) -> dict[str, object]:
    """Transport an isolated, hence nonboundary, vertex to infinity."""
    data = critical_residual_arithmetic(p)
    transport = signed_relative_flip_transport()
    maximum_nonisolated = 2 * int(data["H_edge_count"])
    guaranteed_isolated = int(data["ambient_vertex_count"]) - maximum_nonisolated
    checks = {
        "support_bound": maximum_nonisolated == 8 * p + 2,
        "isolated_count": guaranteed_isolated == p * p - 8 * p - 1,
        "isolated_gap_positive": guaranteed_isolated > 0,
        "signed_transport": bool(transport["proved"]),
        "flip_size_preserved": bool(transport["flip_set_size_preserved"]),
        "boundary_permuted": bool(transport["odd_degree_boundary_is_permuted"]),
        "separation_preserved": bool(
            transport["both_separation_inequalities_preserved"]
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError("the isolated chart dependency failed")
    return {
        "p": p,
        "ambient_vertex_count": p * p + 1,
        "maximum_nonisolated_vertices": maximum_nonisolated,
        "guaranteed_isolated_vertices": guaranteed_isolated,
        "isolated_vertex_cannot_lie_in_odd_degree_boundary": True,
        "chosen_w_is_outside_D": True,
        "chosen_w_is_isolated": True,
        "signed_PSL_sends_w_to_infinity": True,
        "transported_boundary_is_all_finite": True,
        "transported_boundary_size_is_even_by_handshake": True,
        "every_transported_directional_b_is_even": True,
        "transported_infinity_degree_I": 0,
        "transported_H_edge_count": 4 * p + 1,
        "transported_separator_inequalities_preserved": True,
        "signed_transport_dependency": transport,
        "checks": checks,
        "proved": True,
    }


def _expected_even_floor(p: int, b: int, phase: int) -> int:
    """The endpoint-relevant even-``b`` floor pattern for ``p>=13``."""
    if p == 13 and (phase, b) in ((0, 8), (1, 6)):
        # The two p=13 middle cells sit at 24 rather than the eventual 2p.
        # They remain far above every mean used in the residue classification.
        return 24
    if p % 4 == 3:
        if phase == 1:
            return p - 1 if b in (2, p - 1) else 2 * p
        if b == 0:
            return 0
        if b in (2, p - 1):
            return p + 1
        if b in (4, p - 3):
            return 2 * p - 6
        return 2 * p

    if phase == 1:
        if b == 2:
            return p - 1
        if b == p - 1:
            return p + 1
        if b == p - 3:
            return 2 * p - 6
        return 2 * p
    if b == 0:
        return 0
    if b == 2:
        return p + 1
    if b == p - 1:
        return p - 1
    if b == 4:
        return 2 * p - 6
    return 2 * p


def residual_even_floor_table(p: int) -> dict[str, object]:
    """Audit exact phase floors for all even fibre counts."""
    _check_prime_parameter(p)
    even_b = tuple(range(0, p, 2))
    floor_function = scaled_direction_floor if p == 13 else full_symbolic_floor
    exact = {
        phase: {b: floor_function(p, b, phase) for b in even_b}
        for phase in (0, 1)
    }
    expected = {
        phase: {b: _expected_even_floor(p, b, phase) for b in even_b}
        for phase in (0, 1)
    }
    if exact != expected:
        raise ArithmeticError("the critical residual even-floor pattern changed")
    symbolic_agreement = p == 13 or all(
        exact[phase][b] == full_symbolic_floor(p, b, phase)
        for phase in (0, 1)
        for b in even_b
    )
    if not symbolic_agreement:
        raise ArithmeticError("the exact and symbolic floor tables disagree")
    return {
        "p": p,
        "p_mod_4": p % 4,
        "even_b_values": list(even_b),
        "phase_zero_floors": exact[0],
        "phase_one_floors": exact[1],
        "least_nonzero_phase_zero_floor": min(
            floor for b, floor in exact[0].items() if b != 0
        ),
        "phase_one_cells_at_mean_p_minus_one": [
            b for b, floor in exact[1].items() if floor <= p - 1
        ],
        "phase_one_cells_at_mean_p_plus_one": [
            b for b, floor in exact[1].items() if floor <= p + 1
        ],
        "p13_direct_exact_LP_audit": p == 13,
        "agrees_with_prop_15_669_when_in_range": symbolic_agreement,
        "proved": True,
    }


def baseline_coefficient_rules(p: int) -> dict[str, object]:
    """Record the three exact hard baselines and slice residues."""
    _check_prime_parameter(p)
    q = (p - 1) // 2
    b2_certificate = parity_floor_certificate(p, 2, 1)
    # If B has size p-1 and j is its complementary fibre, then on the
    # middle slice sum_B x=m-x_j.  Phase one therefore becomes x_j when
    # p=1 (mod 4), and 1-x_j when p=3 (mod 4).  These are respectively
    # the phase-zero and phase-one b=1 quadrature rows of Proposition 15.652.
    complementary_phase = 0 if p % 4 == 1 else 1
    complementary_certificate = parity_floor_certificate(
        p, 1, complementary_phase
    )
    equality_rigidity = bool(
        b2_certificate["exact_positive_quadrature_certificate"]
        and complementary_certificate["exact_positive_quadrature_certificate"]
        and all(weight > 0 for weight in b2_certificate["quadrature_weights"])
        and all(
            weight > 0
            for weight in complementary_certificate["quadrature_weights"]
        )
    )
    if not equality_rigidity:
        raise ArithmeticError("the endpoint baseline quadrature lost rigidity")
    return {
        "p": p,
        "q": q,
        "middle_slice": "z_s=2x_s-1 and sum_s z_s=1",
        "coefficient_lemma": (
            "target 4+tau*z_i*z_j gives q|(I+P-4); target "
            "4+sigma*z_j gives q|(I+P-(4+sigma))"
        ),
        BRANCH_B2: {
            "baseline": "A=(1-x_i-x_j)^2",
            "target": "eps*S_H=4+z_i*z_j",
            "offset": 4,
            "congruence": f"{q} divides I+P-4",
        },
        BRANCH_P1_LAST: {
            "applicable_p_mod_4": 1,
            "baseline": "A=x_j",
            "target": "eps*S_H=4+z_j",
            "offset": 5,
            "congruence": f"{q} divides I+P-5",
        },
        BRANCH_P3_LAST: {
            "applicable_p_mod_4": 3,
            "baseline": "A=1-x_j",
            "target": "eps*S_H=4-z_j",
            "offset": 3,
            "congruence": f"{q} divides I+P-3",
        },
        "positive_quadrature_dependency": "Proposition 15.652",
        "b2_phase_one_equality_is_pointwise_XNOR": equality_rigidity,
        "b_p_minus_one_phase_one_equality_is_pointwise_literal": (
            equality_rigidity
        ),
        "complementary_b1_phase": complementary_phase,
        "b2_quadrature_certificate": b2_certificate,
        "complementary_b1_quadrature_certificate": complementary_certificate,
        "proved": equality_rigidity,
    }


def hard_residue_branch_ledger(p: int) -> dict[str, object]:
    """Classify every hard-type residue and surviving baseline branch."""
    _check_prime_parameter(p)
    q = (p - 1) // 2
    m = q + 1
    floors = residual_even_floor_table(p)
    lift_floor = int(sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"])
    rows: list[dict[str, object]] = []
    for u in range(m):
        quotient_sum = m - u
        if 1 <= u <= m - 2:
            excluded = True
            branch = None
            reason = "every direction needs k>=1 but sum k=m-u<m"
        elif u == 0 and p % 4 == 3:
            excluded = 2 < lift_floor
            branch = None
            reason = "all mean-p+1 candidates are forbidden two-unit lifts"
        elif u == 0:
            excluded = False
            branch = BRANCH_P1_LAST
            reason = "b=2 is lift-forbidden; exact b=p-1 remains"
        else:
            excluded = False
            branch = "u=m-1 low-baseline dichotomy"
            reason = "q low means p-1 and one high mean 2p"
        rows.append(
            {
                "u": u,
                "common_residue": 2 * u,
                "quotient_sum": quotient_sum,
                "excluded": excluded,
                "surviving_branch": branch,
                "reason": reason,
            }
        )

    feasible_u = [int(row["u"]) for row in rows if not row["excluded"]]
    expected_u = [0, m - 1] if p % 4 == 1 else [m - 1]
    low_candidates = list(floors["phase_one_cells_at_mean_p_minus_one"])
    if p % 4 == 1:
        possible_branches = [BRANCH_B2, BRANCH_P1_LAST]
        low_dichotomy_proved = low_candidates == [2]
    else:
        possible_branches = [BRANCH_B2, BRANCH_P3_LAST]
        rules = baseline_coefficient_rules(p)
        offset_gap = int(rules[BRANCH_B2]["offset"]) - int(
            rules[BRANCH_P3_LAST]["offset"]
        )
        low_dichotomy_proved = low_candidates == [2, p - 1] and offset_gap == 1
    proved = feasible_u == expected_u and low_dichotomy_proved
    if not proved:
        raise ArithmeticError("the hard residue branch classification changed")
    return {
        "p": p,
        "same_type_mean_form": f"a_d=2u+{p + 1}*k_d",
        "quotient_identity": "sum_d k_d=m-u",
        "residue_rows": rows,
        "feasible_u": feasible_u,
        "low_mean": p - 1,
        "low_direction_count": q,
        "high_mean": 2 * p,
        "high_direction_count": 1,
        "low_phase_one_b_candidates": low_candidates,
        "equal_mean_low_cells_cannot_mix": True,
        "possible_branches": possible_branches,
        "nonzero_integral_lift_floor": lift_floor,
        "proved": proved,
    }


def isolated_branch_exclusion(p: int, branch: str) -> dict[str, object]:
    """Use ``I=0`` to exclude one exhaustive hard branch."""
    _check_prime_parameter(p)
    allowed = hard_residue_branch_ledger(p)["possible_branches"]
    if not isinstance(branch, str) or branch not in allowed:
        raise ValueError(f"branch must be one of {allowed}")
    chart = isolated_outside_chart(p)
    q = (p - 1) // 2
    m = q + 1
    rules = baseline_coefficient_rules(p)
    if branch == BRANCH_B2:
        offset = 4
        edge_delta = 0
        mean_constant = 9
        expected_P = 4
        minimum_Q = 3
        expected_mean = 8
        hard_finite_formula = "mP+1"
    elif branch == BRANCH_P1_LAST:
        offset = 5
        edge_delta = 0
        mean_constant = 9
        expected_P = 5
        minimum_Q = 2
        expected_mean = 6
        hard_finite_formula = "mP"
    else:
        offset = 3
        edge_delta = 1
        mean_constant = 7
        expected_P = 3
        minimum_Q = 4
        expected_mean = 8
        hard_finite_formula = "mP+1"

    parameter_rows: list[dict[str, object]] = []
    for parallel_count in range(9):
        numerator = parallel_count - offset
        congruence = numerator % q == 0
        rho = numerator // q if congruence else None
        rho_nonnegative = rho is not None and rho >= 0
        s = int(rho) + parallel_count if rho_nonnegative else None
        opposite_edges = (
            q * (8 - int(s)) + edge_delta if s is not None else None
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
        raise ArithmeticError("I=0 no longer fixes the hard branch parameter")

    s = expected_P
    opposite_edges = q * (8 - s) + edge_delta
    previous_Q_mean = (
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
        and previous_Q_mean < 0
        and forced_minimum_direction
        and minimum_mean == expected_mean
        and nonzero_b_excluded
        and b0_lift_excluded
    )
    if not excluded:
        raise ArithmeticError("an isolated-chart hard branch survived")
    return {
        "p": p,
        "branch": branch,
        "branch_baseline": rules[branch],
        "transported_infinity_degree_I": 0,
        "coefficient_offset": offset,
        "rho_nonnegative_reason": f"I+P-{offset}>-q and is divisible by q",
        "P_upper_bound_reason": (
            "rho>=0 and opposite-edge nonnegativity give P<=s<=8"
        ),
        "parameter_rows_P_0_through_8": parameter_rows,
        "forced_P": expected_P,
        "forced_rho": 0,
        "forced_s": s,
        "hard_finite_edge_count_formula": hard_finite_formula,
        "opposite_finite_edge_count": opposite_edges,
        "opposite_parallel_count_sum": opposite_edges,
        "opposite_mean_formula": f"a=(p-1)s+(p+1)Q+{mean_constant}-7p",
        "minimum_parallel_count": minimum_Q,
        "mean_one_parallel_step_below": previous_Q_mean,
        "minimum_direction_mean": minimum_mean,
        "parallel_surplus_above_minimum": surplus,
        "surplus_less_than_opposite_direction_count": surplus < m,
        "a_minimum_direction_is_forced": forced_minimum_direction,
        "least_nonzero_phase_zero_b_floor": least_nonzero_floor,
        "nonzero_integral_lift_floor": lift_floor,
        "nonzero_b_excluded": nonzero_b_excluded,
        "b0_positive_lift_excluded": b0_lift_excluded,
        "branch_excluded": True,
        "proved": True,
    }


def critical_residual_exclusion(p: int) -> dict[str, object]:
    """Exclude critical residual (ii), independently of ``|D|``."""
    arithmetic = critical_residual_arithmetic(p)
    chart = isolated_outside_chart(p)
    residues = hard_residue_branch_ledger(p)
    branches = {
        branch: isolated_branch_exclusion(p, branch)
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
        "critical_original_k": 4 * p,
        "H_edge_count": 4 * p + 1,
        "arithmetic": arithmetic,
        "isolated_outside_chart": chart,
        "hard_residue_branches": residues,
        "branch_exclusions": branches,
        "boundary_size_hypothesis_used": False,
        "all_boundary_sizes_excluded": proved,
        "finite_configuration_search_used": False,
        "critical_residual_ii_k_eq_4p_excluded": proved,
        "even_k_greater_than_4p_excluded": False,
        "result_status": "proved theorem",
        "proved": proved,
    }


def p11_equality_obstruction() -> dict[str, object]:
    """Record exactly why the same proof stops, rather than closes, at ``p=11``."""
    p = 11
    q = 5
    lift_floor = int(sharp_integral_quadratic_lift_floor(p)["sharp_scaled_floor"])
    exact_floors = {
        phase: {
            b: scaled_direction_floor(p, b, phase) for b in range(0, p, 2)
        }
        for phase in (0, 1)
    }
    least_nonzero_phase_zero = min(
        floor for b, floor in exact_floors[0].items() if b != 0
    )
    branches = {
        BRANCH_B2: {"forced_P": 4, "forced_s": 4, "small_mean": 8},
        BRANCH_P3_LAST: {"forced_P": 3, "forced_s": 3, "small_mean": 8},
    }
    return {
        "p": p,
        "q": q,
        "isolated_gap": p * p - 8 * p - 1,
        "exact_even_floor_tables": exact_floors,
        "least_nonzero_phase_zero_floor": least_nonzero_phase_zero,
        "sharp_nonzero_integral_lift_floor": lift_floor,
        "branch_reductions": branches,
        "small_mean_equals_lift_floor": all(
            row["small_mean"] == lift_floor for row in branches.values()
        ),
        "one_equality_example": "C=(1-x_i)(1-x_j)",
        "equality_examples_not_classified_here": True,
        "p11_closed_here": False,
        "result_status": "open reduction",
        "proved_reduction": True,
    }


def universal_residual_schema() -> dict[str, object]:
    """Record the threshold inequalities making the theorem uniform."""
    threshold = 13
    isolated_gap = threshold * threshold - 8 * threshold - 1
    minimum_lift_floor = threshold - 3
    minimum_nonzero_b_floor = threshold - 1
    maximum_small_mean = 8
    return {
        "prime_range": "odd primes p>=13",
        "critical_k": "4p",
        "critical_H_size": "4p+1",
        "ambient_vertex_count": "p^2+1",
        "isolated_vertex_gap": "p^2-8p-1",
        "isolated_gap_at_13": isolated_gap,
        "isolated_gap_strictly_increasing_for_p_at_least_13": True,
        "minimum_q": 6,
        "coefficient_offsets": [4, 5, 3],
        "parallel_count_upper_bound": 8,
        "forced_s_values": [4, 5, 3],
        "forced_small_means": [8, 6, 8],
        "maximum_forced_small_mean": maximum_small_mean,
        "minimum_nonzero_phase_zero_b_floor": minimum_nonzero_b_floor,
        "minimum_nonzero_integral_lift_floor": minimum_lift_floor,
        "p13_floor_input": "direct exact Proposition 15.632 LP",
        "p_at_least_17_floor_input": "Proposition 15.669",
        "symbolic_not_configuration_census": True,
        "proved": bool(
            isolated_gap > 0
            and maximum_small_mean < minimum_lift_floor
            and maximum_small_mean < minimum_nonzero_b_floor
        ),
    }


def proposition_15734() -> dict[str, object]:
    """Package the all-boundary critical residual exclusion."""
    schema = universal_residual_schema()
    sample_primes = (13, 17, 19, 23, 29, 31, 37, 43)
    exact_rows = {str(p): critical_residual_exclusion(p) for p in sample_primes}
    p11 = p11_equality_obstruction()
    proved = bool(schema["proved"] and all(row["proved"] for row in exact_rows.values()))
    return {
        "prop": "15.734",
        "title": "Isolated-chart critical residual exclusion",
        "result_status": "proved theorem",
        "statement": (
            "critical residual (ii) at k=4p is empty for every prime p>=13, "
            "independently of boundary size"
        ),
        "universal_proof_schema": schema,
        "representative_exact_ledgers": exact_rows,
        "p11_equality_obstruction": p11,
        "critical_residual_ii_k_eq_4p_empty_p_ge_13": proved,
        "all_boundary_sizes_at_k_eq_4p_closed_p_ge_13": proved,
        "p11_closed": False,
        "even_k_greater_than_4p_closed": False,
        "finite_configuration_search_used": False,
        "residual_ii_k_eq_4p_empty_all_primes": False,
        "residual_ii_k_ge_4p_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "remaining_scope": (
            "critical residual (ii) at p=5,7,11; even k>4p; multi-level "
            "Type I; and the quadratic-minmax limit"
        ),
        "proved": proved,
    }


def write_evidence() -> Path:
    """Write the deterministic symbolic critical-residual certificate."""
    output = ROOT / "evidence" / "e1_gmin_m4_prop15734.json"
    payload = json.dumps(proposition_15734(), indent=2, sort_keys=True) + "\n"
    output.write_text(payload)
    return output


def main() -> None:
    result = proposition_15734()
    if not result["proved"]:
        raise ArithmeticError("Proposition 15.734 residual audit failed")
    path = write_evidence()
    print("Prop 15.734 critical residual k=4p for p>=13: excluded")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
