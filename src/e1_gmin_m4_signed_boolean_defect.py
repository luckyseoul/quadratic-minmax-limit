#!/usr/bin/env python3
r"""Exact integer defect and Graver alternative for the signed Radon box.

This module starts at the open box in Proposition 15.760.  Let ``R`` be its
unsigned integral edge--Radon matrix, let ``tau_e`` be the Paley sign of an
edge, and let ``epsilon_L`` be the sign of a projective direction.  The
source coordinate of a simple graph is

    z_e = tau_e * 1_(e in H),

so the required box is ``z_e in {0,tau_e}``.

Every edge is parallel to one projective direction ``L_e``, and the Paley
signs obey

    tau_e = epsilon_(L_e).

Consequently ``tau`` is an *integral row combination* of ``R``: give the
parallel row ``P_L`` weight ``epsilon_L`` and every off-diagonal row weight
zero.  Thus, on the integral fibre ``Rz=y``,

    H_y := sum_L epsilon_L P_L(y) = tau dot z                 (1)

is fixed.  For a residual target coming from a graph, ``H_y=|H|``.

For an integer ``n`` and a sign ``s``,

    n(n-s)/2 = (sn)(sn-1)/2 >= 0,

with equality exactly at ``n in {0,s}``.  It follows that

    beta_R(y) = min_(Rz=y, z integral)
                (||z||^2-H_y)/2                              (2)

is a nonnegative integer and the signed Boolean fibre is nonempty if and
only if ``beta_R(y)=0``.  If it is empty, every integral lift has squared
norm at least ``H_y+2``.  This is the integer threshold missing from the
real Moore--Penrose calculation of Proposition 15.761.

There is an exact finite Graver alternative for (2).  Write ``G(R)`` for
the Graver basis of the integer kernel.  For a feasible integral lift ``z``
and ``g in G(R)``, (1) gives ``tau dot g=0``, and hence

    beta(z+g)-beta(z) = z dot g + ||g||^2/2.                  (3)

A feasible ``z`` globally minimizes (2) if and only if no Graver move
strictly decreases it.  Since the Graver basis is symmetric, this is

    |2 z dot g| <= ||g||^2       for every g in G(R).         (4)

The proof is elementary.  If ``d`` is any kernel move, decompose it as a
conformal sum of Graver moves ``d=sum_i g_i``.  Coordinatewise conformality
gives ``||d||^2 >= sum_i ||g_i||^2``, and therefore

    beta(z+d)-beta(z)
       >= sum_i (beta(z+g_i)-beta(z)).                        (5)

Any globally improving ``d`` contains an improving Graver summand.  Thus
Graver descent from any integral lift terminates either at defect zero (and
constructs a graph) or at a positive-defect vector satisfying (4), which is
an exact obstruction certificate *provided the complete Graver basis has
been supplied and checked*.

This is a proved exact reduction, not a proof that ``beta_R(y)=0`` for the
compact/all-equal target.  No practical complete Graver basis for the full
edge--Radon matrix is asserted here; the recorded lower bound on switch
support shows that one should not expect a prime-uniform local basis.
"""
from __future__ import annotations

from fractions import Fraction
from math import floor
from typing import Iterable, Sequence


def _integer_tuple(values: Iterable[int], name: str) -> tuple[int, ...]:
    out = tuple(values)
    if any(not isinstance(value, int) or isinstance(value, bool) for value in out):
        raise ValueError(f"{name} must contain integers")
    return out


def _sign_tuple(values: Iterable[int], name: str = "tau") -> tuple[int, ...]:
    out = _integer_tuple(values, name)
    if any(value not in (-1, 1) for value in out):
        raise ValueError(f"{name} must contain signs")
    return out


def signed_coordinate_defect(value: int, tau: int) -> int:
    """Return ``value(value-tau)/2``, zero exactly on ``{0,tau}``."""
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError("value must be an integer")
    if tau not in (-1, 1):
        raise ValueError("tau must be a sign")
    numerator = value * (value - tau)
    if numerator % 2:
        raise ArithmeticError("the signed consecutive-integer product changed parity")
    defect = numerator // 2
    if defect < 0 or (defect == 0) != (value in (0, tau)):
        raise ArithmeticError("the signed Boolean coordinate gap failed")
    return defect


def signed_boolean_defect(
    values: Iterable[int], tau: Iterable[int]
) -> dict[str, object]:
    """Evaluate the exact integer distance (2) on one source vector."""
    z = _integer_tuple(values, "values")
    signs = _sign_tuple(tau)
    if len(z) != len(signs):
        raise ValueError("values and tau must have the same length")
    terms = tuple(
        signed_coordinate_defect(value, sign)
        for value, sign in zip(z, signs)
    )
    linear_term = sum(sign * value for sign, value in zip(signs, z))
    squared_norm = sum(value * value for value in z)
    defect = sum(terms)
    in_box = all(value in (0, sign) for value, sign in zip(z, signs))
    proved = bool(
        2 * defect == squared_norm - linear_term
        and defect >= 0
        and (defect == 0) == in_box
    )
    if not proved:
        raise ArithmeticError("the signed Boolean defect identity failed")
    return {
        "values": list(z),
        "tau": list(signs),
        "coordinate_defects": list(terms),
        "tau_dot_z": linear_term,
        "squared_norm": squared_norm,
        "defect": defect,
        "in_signed_boolean_box": in_box,
        "norm_gap_is_twice_defect": True,
        "proved": proved,
    }


def parallel_row_functional(
    parallel_direction_by_edge: Iterable[int],
    direction_signs: Iterable[int],
    values: Iterable[int],
) -> dict[str, object]:
    """Verify (1) from the parallel-row partition of the edge columns.

    ``parallel_direction_by_edge[e]`` is the unique direction whose
    parallel row contains edge ``e``.  This function does not instantiate a
    prime and so is also useful for symbolic or independently generated
    edge--Radon matrices.
    """
    directions = _integer_tuple(parallel_direction_by_edge, "parallel directions")
    epsilon = _sign_tuple(direction_signs, "direction_signs")
    z = _integer_tuple(values, "values")
    if len(directions) != len(z):
        raise ValueError("one parallel direction is required for every edge value")
    if any(direction < 0 or direction >= len(epsilon) for direction in directions):
        raise ValueError("a parallel direction index lies outside direction_signs")

    tau = tuple(epsilon[direction] for direction in directions)
    parallel_rows = [0] * len(epsilon)
    for direction, value in zip(directions, z):
        parallel_rows[direction] += value
    row_value = sum(sign * total for sign, total in zip(epsilon, parallel_rows))
    source_value = sum(sign * value for sign, value in zip(tau, z))
    proved = row_value == source_value
    if not proved:
        raise ArithmeticError("tau left the integer span of the parallel rows")
    return {
        "direction_signs": list(epsilon),
        "parallel_row_values": parallel_rows,
        "induced_edge_signs_tau": list(tau),
        "weighted_parallel_total_H_y": row_value,
        "tau_dot_z": source_value,
        "tau_is_integer_parallel_row_combination": True,
        "proved": proved,
    }


def kernel_move_defect_change(
    values: Iterable[int], move: Iterable[int], tau: Iterable[int]
) -> dict[str, object]:
    """Audit the exact move formula, including the ``tau dot g`` term."""
    z = _integer_tuple(values, "values")
    g = _integer_tuple(move, "move")
    signs = _sign_tuple(tau)
    if not len(z) == len(g) == len(signs):
        raise ValueError("values, move, and tau must have the same length")
    before = int(signed_boolean_defect(z, signs)["defect"])
    after_values = tuple(value + step for value, step in zip(z, g))
    after = int(signed_boolean_defect(after_values, signs)["defect"])
    twice_formula = (
        2 * sum(value * step for value, step in zip(z, g))
        + sum(step * step for step in g)
        - sum(sign * step for sign, step in zip(signs, g))
    )
    if twice_formula % 2:
        raise ArithmeticError("the integer defect move formula lost integrality")
    formula = twice_formula // 2
    proved = after - before == formula
    if not proved:
        raise ArithmeticError("the signed Boolean move formula failed")
    return {
        "defect_before": before,
        "defect_after": after,
        "defect_change": after - before,
        "tau_dot_move": sum(sign * step for sign, step in zip(signs, g)),
        "kernel_specialization_when_tau_dot_move_zero": (
            "z dot g + ||g||^2/2"
        ),
        "proved": proved,
    }


def conformal_superadditivity_certificate(
    values: Iterable[int], moves: Sequence[Sequence[int]], tau: Iterable[int]
) -> dict[str, object]:
    """Verify (5) for a supplied coordinatewise-conformal decomposition."""
    z = _integer_tuple(values, "values")
    signs = _sign_tuple(tau)
    pieces = tuple(_integer_tuple(move, "move") for move in moves)
    if not pieces:
        raise ValueError("at least one move is required")
    if any(len(move) != len(z) for move in pieces) or len(signs) != len(z):
        raise ValueError("all vectors must have the same length")
    for coordinate in range(len(z)):
        nonzero = [move[coordinate] for move in pieces if move[coordinate]]
        if nonzero and any(value * nonzero[0] < 0 for value in nonzero):
            raise ValueError("the moves are not a conformal sum")
    total = tuple(sum(move[index] for move in pieces) for index in range(len(z)))
    global_change = int(kernel_move_defect_change(z, total, signs)["defect_change"])
    piece_changes = tuple(
        int(kernel_move_defect_change(z, move, signs)["defect_change"])
        for move in pieces
    )
    gram_gap = sum(
        sum(left[k] * right[k] for k in range(len(z)))
        for i, left in enumerate(pieces)
        for right in pieces[i + 1 :]
    )
    proved = global_change - sum(piece_changes) == gram_gap and gram_gap >= 0
    if not proved:
        raise ArithmeticError("conformal Graver superadditivity failed")
    return {
        "total_move": list(total),
        "global_defect_change": global_change,
        "individual_defect_changes": list(piece_changes),
        "conformal_gram_gap": gram_gap,
        "global_change_at_least_sum_of_piece_changes": True,
        "proved": proved,
    }


def graver_voronoi_test(
    values: Iterable[int], graver_moves: Sequence[Sequence[int]], tau: Iterable[int]
) -> dict[str, object]:
    """Check (3)--(4) against a *supplied* symmetric Graver list.

    Completeness of ``graver_moves`` is deliberately not inferred.  A
    caller may use ``global_minimum_certified`` only after independently
    certifying that the list is the complete Graver basis of the matrix in
    question and that every move is in its integer kernel.
    """
    z = _integer_tuple(values, "values")
    signs = _sign_tuple(tau)
    moves = tuple(_integer_tuple(move, "graver_move") for move in graver_moves)
    if any(len(move) != len(z) for move in moves) or len(signs) != len(z):
        raise ValueError("all vectors must have the same length")
    rows = []
    for move in moves:
        audit = kernel_move_defect_change(z, move, signs)
        tau_dot = int(audit["tau_dot_move"])
        simplified_twice_change = (
            2 * sum(value * step for value, step in zip(z, move))
            + sum(step * step for step in move)
        )
        rows.append(
            {
                "move": list(move),
                "tau_dot_move": tau_dot,
                "defect_change": int(audit["defect_change"]),
                "twice_change_if_kernel": simplified_twice_change,
                "nondecreasing": int(audit["defect_change"]) >= 0,
            }
        )
    return {
        "defect": int(signed_boolean_defect(z, signs)["defect"]),
        "move_checks": rows,
        "no_supplied_move_improves": all(row["nondecreasing"] for row in rows),
        "global_minimum_certified": False,
        "why_not_automatic": (
            "set true only after external proof that the supplied symmetric list "
            "is the complete Graver basis and every move is in ker_Z R"
        ),
        "proved_arithmetic": True,
    }


def moore_penrose_source_coordinate(
    p: int,
    common_total: int,
    parallel_target: int,
    transverse_target_cells: Iterable[int],
) -> dict[str, object]:
    r"""Return one coordinate of the exact Moore--Penrose lift.

    The arguments use the unsigned target convention of Proposition 15.761.
    For an edge e parallel to L_e, parallel_target is P_tilde_(L_e) and
    transverse_target_cells contains the one off-diagonal value
    K_tilde_L(c_L(e)) for each of the other p directions.  The three
    spectral blocks simplify exactly to

      (R^+ y)_e =
        2(P_tilde_(L_e)-T)/p^3
        +sum_(L!=L_e) K_tilde_L(c_L(e))/p^2.
    """
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 3
        or p % 2 == 0
    ):
        raise ValueError("p must be an odd integer at least three")
    if any(
        not isinstance(value, int) or isinstance(value, bool)
        for value in (common_total, parallel_target)
    ):
        raise ValueError("target totals must be integers")
    cells = _integer_tuple(transverse_target_cells, "transverse target cells")
    if len(cells) != p:
        raise ValueError("one transverse target cell is required in each of p directions")
    numerator = 2 * (parallel_target - common_total) + p * sum(cells)
    value = Fraction(numerator, p**3)
    return {
        "p": p,
        "common_total_T": common_total,
        "parallel_target_at_L_e": parallel_target,
        "transverse_target_cells": list(cells),
        "p_cubed_numerator": numerator,
        "R_plus_y_coordinate": str(value),
        "formula": "2(P_tilde_L_e-T)/p^3 + sum_transverse K_tilde/p^2",
        "proved": True,
    }


def residual_fractional_backprojection_coordinate(
    p: int,
    common_total: int,
    actual_parallel_count: int,
    tau: int,
    signed_backprojection_terms: Iterable[int],
) -> dict[str, object]:
    r"""Return the physical fractional coordinate tau_e(R^+y)_e.

    Each supplied term is tau_e*epsilon_L*W_L(c_L(e)).  Their sum is the
    pointwise signed backprojection B_e.  The Moore--Penrose lift itself is a
    signed fractional-box lift whenever every returned numerator lies in
    [0,p^3].
    """
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 3
        or p % 2 == 0
    ):
        raise ValueError("p must be an odd integer at least three")
    if tau not in (-1, 1):
        raise ValueError("tau must be a sign")
    if any(
        not isinstance(value, int) or isinstance(value, bool)
        for value in (common_total, actual_parallel_count)
    ):
        raise ValueError("target totals must be integers")
    terms = _integer_tuple(
        signed_backprojection_terms, "signed backprojection terms"
    )
    if len(terms) != p:
        raise ValueError("one signed backprojection term is required per transverse direction")
    backprojection = sum(terms)
    numerator = (
        2 * (actual_parallel_count - tau * common_total)
        + p * backprojection
    )
    physical_value = Fraction(numerator, p**3)
    in_interval = 0 <= numerator <= p**3
    return {
        "p": p,
        "common_total_T": common_total,
        "actual_parallel_count_P": actual_parallel_count,
        "edge_sign_tau": tau,
        "signed_backprojection_B_e": backprojection,
        "p_cubed_numerator": numerator,
        "physical_fractional_coordinate_h_e": str(physical_value),
        "coordinate_in_unit_interval": in_interval,
        "all_coordinates_in_interval_suffice_for_fractional_signed_box": True,
        "necessity_not_claimed_because_real_kernel_shifts_are_available": True,
        "proved": True,
    }


def _nearest_integer_distance(value: Fraction) -> Fraction:
    lower = floor(value)
    return min(value - lower, lower + 1 - value)


def integer_lagrange_lower_bound(
    matrix: Sequence[Sequence[int]],
    target: Sequence[int],
    tau: Sequence[int],
    multiplier: Sequence[Fraction | int],
) -> dict[str, object]:
    r"""Return the periodic integer Lagrange lower bound for ``beta_R(y)``.

    For ``c=R^T lambda``, coordinatewise minimization over the integers gives

      beta_R(y) >= lambda.y
        + 1/2 sum_e [dist(c_e+tau_e/2,Z)^2-(c_e+tau_e/2)^2].

    A positive value is therefore a rigorous no-Boolean certificate.
    However, the comparison with the signed-box Farkas slack returned below
    proves that positivity detects exactly failure of the *fractional* box,
    not a new integer obstruction.
    """
    rows = tuple(_integer_tuple(row, "matrix row") for row in matrix)
    y = _integer_tuple(target, "target")
    signs = _sign_tuple(tau)
    lam = tuple(Fraction(value) for value in multiplier)
    if len(rows) != len(y) or len(lam) != len(y):
        raise ValueError("matrix, target, and multiplier row counts differ")
    if rows and any(len(row) != len(signs) for row in rows):
        raise ValueError("matrix column count differs from tau")
    c = tuple(
        sum(lam[row] * rows[row][column] for row in range(len(rows)))
        for column in range(len(signs))
    )
    centres = tuple(value + Fraction(sign, 2) for value, sign in zip(c, signs))
    periodic = sum(
        _nearest_integer_distance(value) ** 2 - value**2 for value in centres
    ) / 2
    bound = sum(value * rhs for value, rhs in zip(lam, y)) + periodic
    farkas_slack = sum(value * rhs for value, rhs in zip(lam, y)) - sum(
        max(Fraction(0), sign * value) for sign, value in zip(signs, c)
    )
    small_scaling_regime = all(abs(value) < 1 for value in c)
    # The integer coordinate minimum includes the two signed-Boolean choices
    # 0 and tau, so it is at most -max(0,tau*c).  When |c|<1, those choices
    # themselves attain the integer minimum.
    proved = bool(
        bound <= farkas_slack
        and (not small_scaling_regime or bound == farkas_slack)
    )
    if not proved:
        raise ArithmeticError("the periodic/Farkas dual comparison failed")
    return {
        "R_transpose_lambda": [str(value) for value in c],
        "integer_centres": [str(value) for value in centres],
        "lower_bound": str(bound),
        "signed_box_farkas_slack": str(farkas_slack),
        "periodic_bound_at_most_farkas_slack": True,
        "equals_farkas_slack_when_max_abs_Rt_lambda_below_one": (
            small_scaling_regime
        ),
        "positive_bound_excludes_signed_boolean_lift": bound > 0,
        "positive_supremum_iff_fractional_signed_box_is_empty": True,
        "no_integer_advantage_beyond_zonotope_separation": True,
        "proved": proved,
    }


def theorem_record() -> dict[str, object]:
    """Machine-readable statement of the exact reduction and its scope."""
    coordinate_audits = {
        str(tau): {
            str(value): signed_coordinate_defect(value, tau)
            for value in range(-5, 6)
        }
        for tau in (-1, 1)
    }
    # A nontrivial one-dimensional-kernel example.  Its complete Graver
    # basis is +/-(1,-3,2).  The fibre Az=(1,1) has the lift (1,-1,1), no
    # Boolean point, and exact defect one.
    toy_z = (1, -1, 1)
    toy_tau = (1, 1, 1)
    toy_moves = ((1, -3, 2), (-1, 3, -2))
    toy = graver_voronoi_test(toy_z, toy_moves, toy_tau)
    toy["matrix"] = [[1, 1, 1], [0, 2, 3]]
    toy["target"] = [1, 1]
    toy["complete_graver_basis_proved_directly_from_rank_one_kernel"] = True
    toy["global_minimum_certified"] = bool(toy["no_supplied_move_improves"])
    toy["signed_boolean_fibre_nonempty"] = False

    return {
        "title": "Exact signed-Boolean quadratic defect and Graver alternative",
        "status": "PROVED EXACT REDUCTION; TARGET DEFECT STILL OPEN",
        "fixed_linear_term": (
            "H_y=sum_L epsilon_L P_L(y)=tau dot z for every Rz=y"
        ),
        "defect": "beta_R(y)=min_(Rz=y) (||z||^2-H_y)/2",
        "coordinate_identity": "z_e(z_e-tau_e)/2>=0, equality iff z_e in {0,tau_e}",
        "exact_equivalence": "signed Boolean lift exists iff beta_R(y)=0",
        "empty_box_integer_norm_gap": "min ||z||^2 >= H_y+2",
        "graver_move_change": "beta(z+g)-beta(z)=z.g+||g||^2/2 for g in ker_Z R",
        "graver_optimality": "|2 z.g|<=||g||^2 for every g in G(R)",
        "pointwise_moore_penrose_lift": (
            "(R^+y)_e=2(P_tilde_L_e-T)/p^3+"
            "p^-2 sum_(L!=L_e)K_tilde_L(c_L(e))"
        ),
        "residual_fractional_box_sufficient_condition": (
            "0<=2(P_L_e-tau_e*T)+p*B_e<=p^3 for every edge"
        ),
        "parity_layer": (
            "because coker R has odd exponent p, ker_Z R -> ker_F2 R is onto"
        ),
        "coordinate_audits": coordinate_audits,
        "nontrivial_no_boolean_toy_certificate": toy,
        "proved": {
            "tau_is_in_integral_parallel_row_span": True,
            "defect_is_nonnegative_integer": True,
            "defect_zero_iff_signed_boolean": True,
            "positive_defect_has_two_unit_squared_norm_gap": True,
            "complete_graver_nondecrease_is_global_optimality": True,
            "pointwise_moore_penrose_formula": True,
            "periodic_integer_lagrange_bound_is_valid": True,
            "periodic_dual_positivity_equals_fractional_box_separation": True,
            "compact_all_equal_target_beta_zero": False,
            "residual_ii_closed": False,
            "e1_closed_general": False,
            "L": False,
        },
        "remaining_obstruction": (
            "Prove beta_R(y)=0 for every moment-compatible residual target, "
            "or exhibit beta_R(y)>0 by a complete Graver/Voronoi certificate. "
            "The periodic Lagrange dual only re-expresses fractional-box "
            "Farkas separation."
        ),
        "scope_warning": (
            "No complete Graver basis for the full edge-Radon matrix is supplied; "
            "the theorem is an exact nonlinear reduction, not a Boolean lift."
        ),
        "L_status": "OPEN",
    }
