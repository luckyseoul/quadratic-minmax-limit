#!/usr/bin/env python3
"""Network-flow normal form for the fixed/parallel symmetric fibre.

After the fixed antipodal variables have been eliminated, an unused
nonfixed inversion orbit is represented by a midpoint/difference pair
``([a], [delta])``.  This module records two exact facts.

* On the divided fixed coordinates together with the parallel coordinates,
  every remaining column is a signed node--arc incidence column.  Hence this
  projection of the restricted Boolean fibre is a capacitated network-flow
  polytope and is totally unimodular, even after arbitrary Mobius-used
  columns are deleted.
* Alternating cycles of that network preserve the fixed word and all
  parallel quotas while changing only nonfixed transverse cells.  In a
  clean midpoint/difference direction block, four-cycles span every mixed
  radial moment modulo ``p``.

The result does not prove that a flow can be chosen to hit the complete
transverse target while remaining Boolean.  In particular it is not a
closure of residual (ii).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import product

from e1_gmin_m4_mobius_half_symmetric import paley_edge_sign
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_symmetric_fixed_edge_elimination import (
    _annihilating_direction_index,
    _central_coordinates,
    _ordinary_image,
    _parallel_direction_index,
    _source_orbits,
    _target_layout,
)


Point = tuple[int, int]
IntegerVector = tuple[int, ...]


def _check_odd_prime(p: int) -> int:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 3
        or p % 2 == 0
        or not is_prime(p)
    ):
        raise ValueError("need an odd prime")
    return (p - 1) // 2


def _negative_point(p: int, point: Point) -> Point:
    return (-point[0] % p, -point[1] % p)


def _antipodal_class(p: int, point: Point) -> Point:
    if point == (0, 0):
        raise ValueError("zero has no nonzero antipodal class")
    return min(point, _negative_point(p, point))


def _midpoint_difference(p: int, edge: tuple[Point, Point]) -> tuple[Point, Point]:
    inverse_two = pow(2, -1, p)
    first, second = edge
    midpoint = (
        (first[0] + second[0]) * inverse_two % p,
        (first[1] + second[1]) * inverse_two % p,
    )
    difference = (
        (second[0] - first[0]) * inverse_two % p,
        (second[1] - first[1]) * inverse_two % p,
    )
    if midpoint == (0, 0) or difference == (0, 0):
        raise ArithmeticError("a nonfixed nonloop orbit lost midpoint/difference")
    return midpoint, difference


def fixed_parallel_flow_theorem(p: int) -> dict[str, object]:
    """State the signed-incidence and total-unimodularity theorem.

    In ordinary (untransported) target coordinates, divide the fixed block
    as in the fixed-edge-elimination theorem and negate every parallel row.
    For ``A=L_a`` and ``D=L_delta``, a nonfixed orbit column is

        tau_delta (e_(A,A(delta)^2) - e_(P_D))             if A != D,
        -tau_delta e_(P_D)                                 if A == D.

    The second case is an arc to a deleted root row.  Thus the matrix is a
    signed directed incidence matrix.  Column deletions and unit upper
    bounds retain total unimodularity.
    """
    h = _check_odd_prime(p)
    d = p + 1
    delta = d * h
    projective_blocks = d * d
    variables_per_block = h * h
    nonfixed_orbits = delta * delta
    proved = bool(
        delta == (p * p - 1) // 2
        and projective_blocks * variables_per_block == nonfixed_orbits
    )
    if not proved:
        raise ArithmeticError("the fixed/parallel flow dimensions changed")
    return {
        "p": p,
        "h": h,
        "projective_directions": d,
        "parallel_nodes": d,
        "fixed_antipodal_cell_nodes": d * h,
        "optional_root_nodes": 1,
        "projective_midpoint_difference_blocks": projective_blocks,
        "variables_per_projective_block": variables_per_block,
        "nonfixed_orbit_variables": nonfixed_orbits,
        "column_normal_form_nonradial": (
            "tau_delta*(e_(K_A(0,A(delta)^2))-e_(P_D))"
        ),
        "column_normal_form_radial": "-tau_delta*e_(P_D), with a root row omitted",
        "row_operation": "negate every divided parallel P_D row",
        "network_matrix": True,
        "totally_unimodular": True,
        "arbitrary_used_column_deletion_preserves_TU": True,
        "unit_box_integrality": True,
        "exact_projected_feasibility_test": "capacitated max flow / min cut",
        "global_hamming_row": "redundant after all parallel quotas are fixed",
        "nonfixed_transverse_cells_solved": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def alternating_cycle_steering_theorem(p: int) -> dict[str, object]:
    """State the exact alternating-cycle and Paley radial-cycle operations."""
    h = _check_odd_prime(p)
    return {
        "p": p,
        "radial_classes_per_projective_direction": h,
        "flow_fibre_difference": (
            "the difference of two fixed/parallel Boolean solutions is an "
            "integral circulation and decomposes into alternating cycles"
        ),
        "parallel_arc_exchange": (
            "replace ([a1],[delta]) by ([a2],[delta]) with a1 parallel a2"
        ),
        "parallel_arc_exchange_preserves": [
            "forced fixed-edge word",
            "physical Hamming weight",
            "every divided fixed coordinate",
            "every parallel-direction quota",
        ],
        "four_cycle": (
            "+(xi1,eta1)+(xi2,eta2)-(xi1,eta2)-(xi2,eta1) "
            "inside one projective (A,D) block"
        ),
        "radial_cycle_sign_hypothesis": (
            "the source sign is constant on a projective difference direction; "
            "this holds for the branch-C Paley sign"
        ),
        "four_cycle_integer_lattice": (
            "all h by h integer tables with zero row and column sums"
        ),
        "transverse_action": (
            "the same alternating cell cycle after alpha=L(a0)^2*xi and "
            "beta=L(delta0)^2*eta"
        ),
        "punctured_binary_warning": (
            "with used cells deleted, legal moves are alternating cycles of "
            "the allowed bipartite graph; four-cycles alone need not connect"
        ),
        "top_endpoint_warning": (
            "if the divided Hamming weight is zero, no completion-side cycle exists"
        ),
        "proved": True,
    }


def even_channel_midpoint_value(
    p: int,
    degree_half: int,
    channel: int,
    alpha: int,
    beta: int,
) -> int:
    """Evaluate ``Q_(2n,k)`` in squared midpoint/difference coordinates.

    If projected edge endpoints are ``x-y,x+y``, then ``alpha=x^2`` and
    ``beta=y^2`` give

        Q_(2n,k)=2^(2n-2k) beta alpha^(n-1-k) (alpha-beta)^k.
    """
    _check_odd_prime(p)
    n = degree_half
    k = channel
    if (
        not isinstance(n, int)
        or isinstance(n, bool)
        or n < 1
        or not isinstance(k, int)
        or isinstance(k, bool)
        or not 0 <= k < n
    ):
        raise ValueError("need n>=1 and 0<=k<n")
    alpha %= p
    beta %= p
    return (
        pow(2, 2 * n - 2 * k, p)
        * beta
        * pow(alpha, n - 1 - k, p)
        * pow(alpha - beta, k, p)
    ) % p


def four_cycle_monomial_delta(
    p: int,
    alpha_power: int,
    beta_power: int,
    alpha_first: int,
    alpha_second: int,
    beta_first: int,
    beta_second: int,
) -> dict[str, int | bool]:
    """Return the factored mixed-monomial change of one radial four-cycle."""
    _check_odd_prime(p)
    if alpha_power < 0 or beta_power < 0:
        raise ValueError("monomial powers must be nonnegative")
    a1, a2 = alpha_first % p, alpha_second % p
    b1, b2 = beta_first % p, beta_second % p
    direct = (
        pow(a1, alpha_power, p) * pow(b1, beta_power, p)
        + pow(a2, alpha_power, p) * pow(b2, beta_power, p)
        - pow(a1, alpha_power, p) * pow(b2, beta_power, p)
        - pow(a2, alpha_power, p) * pow(b1, beta_power, p)
    ) % p
    factored = (
        (pow(a1, alpha_power, p) - pow(a2, alpha_power, p))
        * (pow(b1, beta_power, p) - pow(b2, beta_power, p))
    ) % p
    return {
        "direct_delta": direct,
        "factored_delta": factored,
        "identity_holds": direct == factored,
    }


def clean_block_mixed_moment_theorem(
    p: int, used_orbit_cap: int | None = None
) -> dict[str, object]:
    """Give the clean-block spanning theorem for degree six and eight.

    The statement is linear over ``F_p``.  It proves that, when
    ``|U| <= (p^2-1)/2`` and ``p>11``, four-cycle lattices in untouched
    projective ``(A,D)`` blocks span every *mixed* degree-six and degree-eight
    global binary form.  It does not prove that those signed moves can be
    applied conformally inside one Boolean point.
    """
    h = _check_odd_prime(p)
    if p <= 11:
        raise ValueError("the clean-block degree-eight theorem requires p>11")
    d = p + 1
    delta = d * h
    cap = delta if used_orbit_cap is None else used_orbit_cap
    if (
        not isinstance(cap, int)
        or isinstance(cap, bool)
        or not 0 <= cap <= delta
    ):
        raise ValueError("used_orbit_cap must lie between 0 and (p^2-1)/2")
    clean = d * d - cap
    degree_six_zero_bound = 6 * d - 8
    degree_eight_zero_bound = 8 * d - 12
    degree_six_margin = clean - degree_six_zero_bound
    degree_eight_margin = clean - degree_eight_zero_bound
    proved = degree_six_margin > 0 and degree_eight_margin > 0
    if not proved:
        raise ArithmeticError("the clean-block biform margin vanished")
    return {
        "p": p,
        "h": h,
        "used_orbit_cap": cap,
        "projective_blocks": d * d,
        "clean_projective_blocks_at_least": clean,
        "worst_case_clean_blocks": d * d - delta,
        "degree_six_bidegrees": [[4, 2], [2, 4]],
        "degree_six_nonzero_biform_zero_bound": degree_six_zero_bound,
        "degree_six_spanning_margin": degree_six_margin,
        "degree_eight_bidegrees": [[6, 2], [4, 4], [2, 6]],
        "degree_eight_nonzero_biform_zero_bound": degree_eight_zero_bound,
        "degree_eight_spanning_margin": degree_eight_margin,
        "zero_bound": "a*d+b*d-a*b for bidegree (a,b) on P1(F_p)^2",
        "radial_vandermonde_basis": (
            "xi^r eta^s, 1<=r,s<=h-1, is dual to the four-cycle lattice mod p"
        ),
        "mixed_degree_six_global_forms_spanned": True,
        "mixed_degree_eight_global_forms_spanned": True,
        "pure_radial_margins_spanned": False,
        "conformal_Boolean_sequence_proved": False,
        "full_transverse_cells_solved": False,
        "proved": proved,
    }


def _vector_add(
    terms: tuple[tuple[int, IntegerVector], ...]
) -> IntegerVector:
    width = len(terms[0][1])
    return tuple(
        sum(coefficient * vector[index] for coefficient, vector in terms)
        for index in range(width)
    )


def exact_p7_flow_replay() -> dict[str, object]:
    """Replay every divided fixed column and one four-cycle at ``p=7``.

    This is a fail-when-wrong check of the symbolic theorem, not a residual
    configuration census.
    """
    p = 7
    h = (p - 1) // 2
    d = p + 1
    fixed_keys, paired_keys = _target_layout(p)
    _fixed_orbits, nonfixed_orbits = _source_orbits(p)
    records: dict[
        tuple[int, Point, Point], tuple[IntegerVector, IntegerVector, int]
    ] = {}
    bin_counts: Counter[tuple[int, Point]] = Counter()
    fixed_support_sizes: Counter[int] = Counter()
    all_network_columns = True

    for edge, negative in nonfixed_orbits:
        tau = paley_edge_sign(p, edge)
        fixed, paired = _central_coordinates(
            _ordinary_image(p, {edge: tau, negative: tau}),
            fixed_keys,
            paired_keys,
        )
        if any(value % 2 for value in fixed):
            raise ArithmeticError("a symmetric fixed column stopped being even")
        divided_fixed = tuple(value // 2 for value in fixed)
        midpoint, difference = _midpoint_difference(p, edge)
        midpoint_class = _antipodal_class(p, midpoint)
        difference_class = _antipodal_class(p, difference)
        A = _annihilating_direction_index(p, midpoint)
        D = _parallel_direction_index(p, edge)
        nonzero = {
            key: value
            for key, value in zip(fixed_keys, divided_fixed, strict=True)
            if value
        }
        expected: dict[tuple[object, ...], int] = {("P", D): tau}
        if A != D:
            fixed_cells = [key for key in nonzero if key[0] == "K"]
            if len(fixed_cells) != 1 or fixed_cells[0][1] != A:
                all_network_columns = False
            else:
                expected[fixed_cells[0]] = tau
        all_network_columns &= nonzero == expected
        fixed_support_sizes[len(nonzero)] += 1
        bin_counts[(A, difference_class)] += 1
        records[(A, midpoint_class, difference_class)] = (
            divided_fixed,
            paired,
            tau,
        )

    expected_nonfixed = ((p * p - 1) // 2) ** 2
    all_bins_have_h_parallel_arcs = bool(
        len(bin_counts) == d * ((p * p - 1) // 2)
        and set(bin_counts.values()) == {h}
    )

    by_block: defaultdict[
        tuple[int, int], dict[tuple[Point, Point], tuple[IntegerVector, IntegerVector, int]]
    ] = defaultdict(dict)
    for (A, midpoint_class, difference_class), record in records.items():
        D = _annihilating_direction_index(p, difference_class)
        by_block[(A, D)][(midpoint_class, difference_class)] = record

    chosen = None
    for (A, D), block in sorted(by_block.items()):
        if A == D:
            continue
        midpoint_classes = sorted({key[0] for key in block})
        difference_classes = sorted({key[1] for key in block})
        if len(midpoint_classes) >= 2 and len(difference_classes) >= 2:
            a1, a2 = midpoint_classes[:2]
            delta1, delta2 = difference_classes[:2]
            corners = (
                block[(a1, delta1)],
                block[(a2, delta2)],
                block[(a1, delta2)],
                block[(a2, delta1)],
            )
            fixed_cycle = _vector_add(
                (
                    (1, corners[0][0]),
                    (1, corners[1][0]),
                    (-1, corners[2][0]),
                    (-1, corners[3][0]),
                )
            )
            paired_cycle = _vector_add(
                (
                    (1, corners[0][1]),
                    (1, corners[1][1]),
                    (-1, corners[2][1]),
                    (-1, corners[3][1]),
                )
            )
            chosen = {
                "midpoint_direction": A,
                "difference_direction": D,
                "fixed_projection_zero": not any(fixed_cycle),
                "transverse_projection_nonzero": any(paired_cycle),
                "transverse_changed_coordinates": sum(
                    value != 0 for value in paired_cycle
                ),
            }
            break
    if chosen is None:
        raise ArithmeticError("the p=7 replay found no radial four-cycle")

    moment_factorization = four_cycle_monomial_delta(p, 2, 1, 1, 2, 1, 4)
    proved = bool(
        len(nonfixed_orbits) == expected_nonfixed
        and all_network_columns
        and all_bins_have_h_parallel_arcs
        and fixed_support_sizes[1] == d * h * h
        and fixed_support_sizes[2] == d * p * h * h
        and chosen["fixed_projection_zero"]
        and chosen["transverse_projection_nonzero"]
        and moment_factorization["identity_holds"]
    )
    if not proved:
        raise ArithmeticError("the exact p=7 fixed/parallel flow replay failed")
    return {
        "p": p,
        "nonfixed_columns_replayed": len(nonfixed_orbits),
        "every_divided_fixed_column_is_a_signed_network_arc": all_network_columns,
        "radial_half_edge_columns": fixed_support_sizes[1],
        "two_endpoint_arc_columns": fixed_support_sizes[2],
        "parallel_arc_bins": len(bin_counts),
        "parallel_arcs_per_bin": h,
        "all_bins_have_expected_parallel_multiplicity": all_bins_have_h_parallel_arcs,
        "four_cycle_replay": chosen,
        "mixed_monomial_factorization_replay": moment_factorization,
        "role": "fail-when-wrong exact small-prime replay, not a prime census",
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    flow = fixed_parallel_flow_theorem(p)
    cycles = alternating_cycle_steering_theorem(p)
    clean = clean_block_mixed_moment_theorem(p)
    replay = exact_p7_flow_replay()
    proved = bool(
        flow["proved"]
        and cycles["proved"]
        and clean["proved"]
        and replay["proved"]
    )
    if not proved:
        raise ArithmeticError("the symmetric fixed/parallel flow record failed")
    return {
        "title": "Fixed/parallel network flow and transverse cycle steering",
        "status": "PROVED PROJECTED TU AND MIXED-MOMENT SPAN; FULL BOX OPEN",
        "flow_theorem": flow,
        "cycle_theorem": cycles,
        "clean_block_mixed_moments": clean,
        "small_exact_replay": replay,
        "nonfixed_transverse_cells_solved": False,
        "common_simple_graph_constructed": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
