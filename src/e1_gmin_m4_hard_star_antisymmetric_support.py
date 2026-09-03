#!/usr/bin/env python3
"""Exact support bounds for the fixed hard-star antisymmetric target.

This is the first global support consequence of hard compact-residual
centrality.  It separates a Boolean edge set into antipodal single, double,
and fixed orbits; proves the sharp one-row ternary-support bound; and checks
that total edge counts, fixed-edge capacity, parallel pair totals, and the
real least-norm inequality do not by themselves contradict the balanced
branch-C ray.

The cross-directional antisymmetric condition isolated here was subsequently
solved constructively by the direction-localized Mobius trades in
``e1_gmin_m4_inversion_antisymmetric_radon``.  The live gate is now the
coupled symmetric half.
"""

from __future__ import annotations

from fractions import Fraction

from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_prop15758 import p3_local_survivor


def _check_p(p: int) -> int:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 31
        or p % 4 != 3
        or not is_prime(p)
    ):
        raise ValueError("need a prime p=4r+3 with p>=31")
    return (p - 3) // 4


def hard_star_antisymmetric_target(
    p: int, active_hard_direction_count: int
) -> dict[str, object]:
    """Return the exact norm and support of the fixed star differences."""
    r = _check_p(p)
    hard_direction_count = 2 * r + 2
    active = active_hard_direction_count
    if (
        not isinstance(active, int)
        or isinstance(active, bool)
        or not 0 <= active <= hard_direction_count
    ):
        raise ValueError("active hard direction count is out of range")

    cell_support_per_active_row = 2 * (p - 2)
    target_norm_squared = active * cell_support_per_active_row
    real_least_norm_squared = Fraction(target_norm_squared, p * p)
    combinatorial_support_floor = 0 if active == 0 else 2 * (p - 2)
    return {
        "p": p,
        "r": r,
        "hard_direction_count": hard_direction_count,
        "active_nonzero_star_centres": active,
        "zero_star_centres_have_zero_antisymmetric_target": True,
        "one_active_row_target": "plus/minus (S_j-S_-j), j!=0",
        "nonzero_cells_per_active_row": cell_support_per_active_row,
        "target_antisymmetric_norm_squared": target_norm_squared,
        "antisymmetric_RRt_eigenvalue": p * p,
        "exact_real_least_norm_squared": str(real_least_norm_squared),
        "real_bound_below_one_even_when_every_hard_row_is_active": (
            Fraction(
                hard_direction_count * cell_support_per_active_row,
                p * p,
            )
            < 1
        ),
        "ternary_support_floor": combinatorial_support_floor,
        "single_antipodal_orbit_floor": combinatorial_support_floor // 2,
        "support_floor_reason": (
            "the 2(p-2) nonzero cells in one active directional partition "
            "are disjoint and each needs a nonzero source coordinate"
        ),
        "proved": True,
    }


def exceptional_direction_support_bound(
    p: int,
    active_hard_direction_count: int,
    exceptional_active_incidence: int,
) -> dict[str, object]:
    """Refine the single-orbit floor by midpoint/difference exceptions.

    A nonfixed antipodal edge orbit has midpoint ``a!=0`` and difference
    class ``[delta]``.  It contributes to a nonself transverse cell in an
    active direction unless that direction annihilates ``a`` or ``delta``.
    Count the *distinct* active exceptional directions for every chosen
    single orbit and sum them as ``E``.  Incidence counting gives

        A*(p-2) <= A*c-E.
    """
    r = _check_p(p)
    active = active_hard_direction_count
    exceptional = exceptional_active_incidence
    maximum_active = 2 * r + 2
    if (
        not isinstance(active, int)
        or isinstance(active, bool)
        or not 1 <= active <= maximum_active
    ):
        raise ValueError("need at least one active hard direction")
    if (
        not isinstance(exceptional, int)
        or isinstance(exceptional, bool)
        or exceptional < 0
    ):
        raise ValueError("exceptional incidence must be nonnegative")
    extra = (exceptional + active - 1) // active
    single_floor = p - 2 + extra
    return {
        "p": p,
        "r": r,
        "active_hard_direction_count": active,
        "exceptional_active_incidence": exceptional,
        "incidence_inequality": "A*(p-2)<=A*c-E",
        "single_antipodal_orbit_floor": single_floor,
        "ternary_support_floor": 2 * single_floor,
        "equality_at_c_equals_p_minus_2_forces_E_zero": True,
        "equality_rigidity": (
            "every selected single orbit has both its midpoint direction "
            "and difference direction outside the active hard set, and its "
            "images form a bijective star-cell transversal in every active row"
        ),
        "proved": True,
    }


def equality_case_projective_pencil_lemma(
    p: int, active_hard_direction_count: int
) -> dict[str, object]:
    """Classify equality when at least nine hard star rows are active.

    Suppose ``c=p-2`` nonfixed antipodal source-edge orbits attain the
    one-row support floor.  In every active row each selected edge is then
    a transversal of exactly one of the two exceptional affine lines.
    For two selected edges ``{u,v}`` and ``{u',v'}``, their endpoint-square
    sets have a common member in every active direction ``L``.  Hence

        product_(a in {u,v}, b in {u',v'})
            (L(a)^2-L(b)^2)

    vanishes there.  This is a product of eight projective linear factors.
    Nine directions force one factor to vanish identically, so the two
    edges share an endpoint modulo sign.

    The projective endpoint pairs are therefore pairwise intersecting.
    A non-pencil family of two-subsets is contained in a triangle, and
    each projective endpoint pair has only two lifts modulo simultaneous
    negation.  It consequently contains at most six source-edge orbits,
    whereas ``p-2>=29``.  Thus all selected orbits share a projective
    source vertex ``[P]``.  Row bijectivity then forces
    ``j_L^2=L(P)^2`` in every active row.
    """
    r = _check_p(p)
    active = active_hard_direction_count
    maximum_active = 2 * r + 2
    if (
        not isinstance(active, int)
        or isinstance(active, bool)
        or not 9 <= active <= maximum_active
    ):
        raise ValueError("need between 9 and (p+1)/2 active hard directions")
    single_count = p - 2
    nonpencil_orbit_cap = 3 * 2
    proved = bool(
        active > 8
        and single_count > nonpencil_orbit_cap
        and p >= 31
    )
    if not proved:
        raise ArithmeticError("the equality-case pencil inequalities changed")
    return {
        "p": p,
        "r": r,
        "active_hard_direction_count": active,
        "hypothesis": "c=p-2 and every active row attains its star-cell floor",
        "pair_resultant": (
            "prod_{a in {u,v},b in {u',v'}} "
            "(L(a)^2-L(b)^2)"
        ),
        "pair_resultant_degree": 8,
        "projective_root_count_threshold": 9,
        "conclusion_pairwise_projective_endpoint_intersection": True,
        "lifts_per_projective_endpoint_pair": 2,
        "nonpencil_triangle_orbit_cap": nonpencil_orbit_cap,
        "single_orbit_count": single_count,
        "conclusion_common_projective_source_vertex": True,
        "conclusion_active_phase_coherence": "j_L^2=L(P)^2",
        "rows_with_j_L_zero_are_not_active_and_supply_no_root": True,
        "does_not_handle_fewer_than_nine_active_rows": True,
        "proved": proved,
    }


def _point_neg(point: tuple[int, int], p: int) -> tuple[int, int]:
    return ((-point[0]) % p, (-point[1]) % p)


def _edge_key(
    left: tuple[int, int], right: tuple[int, int]
) -> tuple[tuple[int, int], tuple[int, int]]:
    return tuple(sorted((left, right)))  # type: ignore[return-value]


def _edge_orbit_key(
    left: tuple[int, int], right: tuple[int, int], p: int
) -> tuple[tuple[int, int], tuple[int, int]]:
    edge = _edge_key(left, right)
    negative = _edge_key(_point_neg(left, p), _point_neg(right, p))
    return min(edge, negative)


def _projected_signed_radon(
    p: int,
    centre: tuple[int, int],
    signed_translates: list[tuple[tuple[int, int], int]],
    coordinate: int,
) -> dict[tuple[int, int], int]:
    """Project the antisymmetric pencil chain to x or y cells."""
    coefficients: dict[tuple[int, int], int] = {}
    for translate, sign in signed_translates:
        neighbour = (
            (translate[0] - centre[0]) % p,
            (translate[1] - centre[1]) % p,
        )
        cell = tuple(sorted((centre[coordinate], neighbour[coordinate])))
        negative_cell = tuple(
            sorted(((-centre[coordinate]) % p, (-neighbour[coordinate]) % p))
        )
        coefficients[cell] = coefficients.get(cell, 0) + sign
        coefficients[negative_cell] = coefficients.get(negative_cell, 0) - sign
    return {cell: value for cell, value in coefficients.items() if value}


def two_zero_direction_pencil_counterexample(p: int) -> dict[str, object]:
    """Falsify the claim that two zero rows kill a full-pencil equality case.

    This is deliberately only a counterexample to that proposed *local
    implication*.  It does not meet the active hard-row bijectivity
    conditions and is not a residual-(ii) lift.

    In translated neighbour coordinates ``z=w+P``, use a signed alternating
    path of five cells from row zero to column zero, followed by ``r-1``
    disjoint signed 2-by-2 rectangles.  There are ``5+4(r-1)=p-2``
    distinct nonzero points.  Every nonzero x- and y-fibre has signed sum
    zero.  The zero fibre is invisible to the antisymmetric edge Radon map,
    because ``{P,w}`` and ``{-P,-w}`` project to the same self-antipodal
    cell when ``L(z)=0``.
    """
    r = _check_p(p)
    signed_translates: list[tuple[tuple[int, int], int]] = [
        ((0, 1), 1),
        ((1, 1), -1),
        ((1, 2), 1),
        ((2, 2), -1),
        ((2, 0), 1),
    ]
    for index in range(1, r):
        low = 2 * index + 1
        high = 2 * index + 2
        signed_translates.extend(
            [
                ((low, low), 1),
                ((low, high), -1),
                ((high, low), -1),
                ((high, high), 1),
            ]
        )

    centre = (2 * r + 1, 2 * r + 1)
    points = [point for point, _ in signed_translates]
    neighbours = [
        ((point[0] - centre[0]) % p, (point[1] - centre[1]) % p)
        for point in points
    ]
    orbit_keys = [
        _edge_orbit_key(centre, neighbour, p) for neighbour in neighbours
    ]
    x_fibre_sums: dict[int, int] = {}
    y_fibre_sums: dict[int, int] = {}
    for (x_value, y_value), sign in signed_translates:
        x_fibre_sums[x_value] = x_fibre_sums.get(x_value, 0) + sign
        y_fibre_sums[y_value] = y_fibre_sums.get(y_value, 0) + sign
    x_radon = _projected_signed_radon(
        p, centre, signed_translates, coordinate=0
    )
    y_radon = _projected_signed_radon(
        p, centre, signed_translates, coordinate=1
    )
    twice_centre = ((2 * centre[0]) % p, (2 * centre[1]) % p)
    proved = bool(
        len(points) == p - 2
        and len(set(points)) == p - 2
        and (0, 0) not in points
        and twice_centre not in points
        and centre[0] != 0
        and centre[1] != 0
        and all(neighbour != centre for neighbour in neighbours)
        and all(neighbour != _point_neg(centre, p) for neighbour in neighbours)
        and len(set(orbit_keys)) == p - 2
        and x_fibre_sums.get(0) == 1
        and y_fibre_sums.get(0) == 1
        and all(value == 0 for key, value in x_fibre_sums.items() if key)
        and all(value == 0 for key, value in y_fibre_sums.items() if key)
        and not x_radon
        and not y_radon
    )
    if not proved:
        raise ArithmeticError("the two-zero-row pencil counterexample changed")
    return {
        "p": p,
        "r": r,
        "projective_pencil_centre": list(centre),
        "signed_translated_neighbours": [
            {"z": list(point), "sign": sign}
            for point, sign in signed_translates
        ],
        "single_nonfixed_edge_orbits": len(orbit_keys),
        "expected_single_orbits": p - 2,
        "x_nonzero_fibre_sums_zero": True,
        "y_nonzero_fibre_sums_zero": True,
        "x_zero_fibre_sum": x_fibre_sums[0],
        "y_zero_fibre_sum": y_fibre_sums[0],
        "x_antisymmetric_edge_Radon_zero": True,
        "y_antisymmetric_edge_Radon_zero": True,
        "what_is_falsified": (
            "a full projective pencil with p-2 single edge orbits cannot "
            "have zero antisymmetric projections in two directions"
        ),
        "active_hard_row_bijections_satisfied": False,
        "residual_ii_counterexample": False,
        "proved": proved,
    }


def antipodal_pair_total_ledger(p: int, t: int) -> dict[str, object]:
    """Check total-count/fixed-capacity compatibility at the support floor."""
    r = _check_p(p)
    row = p3_local_survivor(p, t)
    edge_count = int(row["H_edge_count"])
    fixed_edge_count = (p * p - 1) // 2
    nonfixed_pair_orbit_count = (p * p - 1) ** 2 // 4
    single_count = p - 2
    selected_fixed_count = 0
    double_count = (edge_count - single_count) // 2
    ternary_support = 2 * single_count
    antisymmetric_half_norm = Fraction(single_count, 2)
    symmetric_half_norm = Fraction(edge_count) - antisymmetric_half_norm

    proved = bool(
        row["proved_local_aggregate"]
        and edge_count % 2 == 1
        and single_count % 2 == 1
        and edge_count == selected_fixed_count + single_count + 2 * double_count
        and 0 <= selected_fixed_count <= fixed_edge_count
        and 0 <= double_count <= nonfixed_pair_orbit_count - single_count
        and antisymmetric_half_norm + symmetric_half_norm == edge_count
    )
    if not proved:
        raise ArithmeticError("the antipodal pair-total ledger changed")
    return {
        "p": p,
        "r": r,
        "t": t,
        "edge_count": edge_count,
        "available_fixed_antipodal_edges": fixed_edge_count,
        "available_nonfixed_antipodal_pair_orbits": nonfixed_pair_orbit_count,
        "single_pair_orbits_c": single_count,
        "double_pair_orbits_d": double_count,
        "selected_fixed_edges_f": selected_fixed_count,
        "total_count_identity": "|H|=f+c+2d",
        "ternary_support_identity": "|{e:x_e!=x_-e}|=2c",
        "ternary_support": ternary_support,
        "antisymmetric_half_squared_norm": str(antisymmetric_half_norm),
        "symmetric_half_squared_norm": str(symmetric_half_norm),
        "orthogonal_norm_identity": "||x^-||^2=c/2, ||x^+||^2=|H|-c/2",
        "total_and_fixed_capacity_contradiction": False,
        "this_is_only_an_orbit_count_ledger": True,
        "proved": proved,
    }


def parallel_pair_total_ledger(p: int, t: int) -> dict[str, object]:
    """Construct a scalar parallel-count ledger with all singles opposite.

    This does not construct edge orbits.  It proves that the exact equations
    ``P_D=f_D+c_D+2d_D`` and the fixed-edge capacities alone do not force a
    contradiction or force a single orbit to have hard difference direction.
    """
    r = _check_p(p)
    row = p3_local_survivor(p, t)
    fixed_capacity_per_direction = (p - 1) // 2
    hard_parallel = [int(item["parallel_P"]) for item in row["hard_rows"]]
    opposite_parallel = [int(item["Q"]) for item in row["opposite_rows"]]
    needed = p - 2
    opposite_singles: list[int] = []
    for parallel in opposite_parallel:
        take = min(parallel, needed)
        opposite_singles.append(take)
        needed -= take
    if needed:
        raise ArithmeticError("opposite parallel capacity fell below p-2")
    single_by_direction = [0] * len(hard_parallel) + opposite_singles
    parallel_counts = hard_parallel + opposite_parallel
    fixed_by_direction = [
        (parallel - single) % 2
        for parallel, single in zip(parallel_counts, single_by_direction)
    ]
    double_by_direction = [
        (parallel - fixed - single) // 2
        for parallel, fixed, single in zip(
            parallel_counts, fixed_by_direction, single_by_direction
        )
    ]
    proved = bool(
        row["proved_local_aggregate"]
        and sum(single_by_direction) == p - 2
        and not any(single_by_direction[: len(hard_parallel)])
        and all(0 <= fixed <= fixed_capacity_per_direction for fixed in fixed_by_direction)
        and all(double >= 0 for double in double_by_direction)
        and all(
            parallel == fixed + single + 2 * double
            for parallel, fixed, single, double in zip(
                parallel_counts,
                fixed_by_direction,
                single_by_direction,
                double_by_direction,
            )
        )
        and sum(parallel_counts) == int(row["H_edge_count"])
    )
    if not proved:
        raise ArithmeticError("the parallel pair-total ledger changed")
    return {
        "p": p,
        "r": r,
        "t": t,
        "fixed_edge_capacity_per_direction": fixed_capacity_per_direction,
        "parallel_counts": parallel_counts,
        "single_pair_orbits_by_direction": single_by_direction,
        "fixed_selected_edges_by_direction": fixed_by_direction,
        "double_pair_orbits_by_direction": double_by_direction,
        "parallel_identity": "P_D=f_D+c_D+2d_D",
        "all_single_difference_directions_are_opposite_type": True,
        "scalar_parallel_pair_totals_contradict": False,
        "common_edge_orbits_or_midpoints_constructed": False,
        "proved": proved,
    }


def symmetric_half_norm_barrier(p: int) -> dict[str, object]:
    """Show that the old full-target norm bound fits below the symmetric budget."""
    r = _check_p(p)
    minimum_edge_count = 4 * r * r + 8 * r + 9
    support_floor_single_count = p - 2
    symmetric_norm_budget = Fraction(minimum_edge_count) - Fraction(
        support_floor_single_count, 2
    )
    old_full_target_upper = (
        Fraction(
            (2 * r + 2) * ((10 * r - 4) ** 2 + (6 * r - 3) ** 2),
            p * p,
        )
        + 2
    )
    gap = symmetric_norm_budget - old_full_target_upper
    gap_numerator = (
        128 * r**4
        - 160 * r**3
        + 488 * r**2
        + 784 * r
        + 17
    )
    u = r - 7
    shifted_numerator = (
        128 * u**4
        + 3424 * u**3
        + 34760 * u**2
        + 159712 * u
        + 281865
    )
    proved = bool(
        gap == Fraction(gap_numerator, 2 * p * p)
        and gap_numerator == shifted_numerator
        and u >= 0
        and gap > 0
    )
    if not proved:
        raise ArithmeticError("the symmetric-half norm margin changed")
    return {
        "p": p,
        "r": r,
        "minimum_branch_C_edge_count": minimum_edge_count,
        "single_orbit_floor": support_floor_single_count,
        "symmetric_half_norm_budget_at_support_floor": str(
            symmetric_norm_budget
        ),
        "old_uniform_full_target_least_norm_upper_bound": str(
            old_full_target_upper
        ),
        "budget_minus_upper_bound": str(gap),
        "twice_p_squared_times_margin": gap_numerator,
        "shift_r_equals_u_plus_7_coefficients_low_to_high": [
            281865,
            159712,
            34760,
            3424,
            128,
        ],
        "every_shifted_coefficient_positive": True,
        "real_norm_and_symmetric_half_contradiction": False,
        "proved": proved,
    }


def balanced_ray_antisymmetric_support_barrier(p: int) -> dict[str, object]:
    """Record the support inequality and the superseding Mobius result."""
    r = _check_p(p)
    lower = 2 * r * r - 4 * r - 2
    upper = 4 * r * r - 2 * r - 5
    lower_ledger = antipodal_pair_total_ledger(p, lower)
    upper_ledger = antipodal_pair_total_ledger(p, upper)
    lower_parallel = parallel_pair_total_ledger(p, lower)
    upper_parallel = parallel_pair_total_ledger(p, upper)
    norm = symmetric_half_norm_barrier(p)
    proved = bool(
        lower_ledger["proved"]
        and upper_ledger["proved"]
        and lower_parallel["proved"]
        and upper_parallel["proved"]
        and norm["proved"]
        and lower_ledger["double_pair_orbits_d"] == 2 * r * r + 2 * r + 4
        and upper_ledger["double_pair_orbits_d"] == 4 * r * r + 4 * r + 1
        and upper_ledger["edge_count"]
        == upper_ledger["available_fixed_antipodal_edges"] - 1
    )
    if not proved:
        raise ArithmeticError("the full-ray antisymmetric support audit changed")
    return {
        "p": p,
        "r": r,
        "branch_C_t_interval": [lower, upper],
        "conditional_on_at_least_one_nonzero_hard_star_centre": True,
        "ternary_support_floor": 2 * (p - 2),
        "single_orbit_floor": p - 2,
        "floor_count_ledger_double_range": [
            lower_ledger["double_pair_orbits_d"],
            upper_ledger["double_pair_orbits_d"],
        ],
        "upper_edge_count_equals_fixed_capacity_minus_one": True,
        "total_edge_count_contradiction": False,
        "fixed_antipodal_capacity_contradiction": False,
        "parallel_pair_total_contradiction": False,
        "real_norm_contradiction": False,
        "former_antisymmetric_gate": (
            "realize common oriented single edge orbits across active hard "
            "and opposite/zero rows"
        ),
        "former_antisymmetric_gate_superseded_by": (
            "NOTE_2026-09-03_INVERSION_ANTISYMMETRIC_RADON.md, whose "
            "direction-localized Mobius trades give a disjoint ternary lift"
        ),
        "antisymmetric_ternary_lift_now_proved": True,
        "remaining_exact_gate": (
            "solve the coupled symmetric half: pair total s_e=1 on every "
            "used Mobius orbit, s_e in {0,2} on unused nonfixed orbits, and "
            "independent binary choices on fixed antipodal edges"
        ),
        "one_common_simple_graph_constructed": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def theorem_record() -> dict[str, object]:
    samples = {
        str(p): balanced_ray_antisymmetric_support_barrier(p)
        for p in (31, 43, 47)
    }
    return {
        "title": "Hard-star antisymmetric ternary-support barrier",
        "status": (
            "PROVED SUPPORT/EQUALITY STRUCTURE; ANTISYMMETRIC GATE "
            "SUPERSEDED BY MOBIUS LIFT; SYMMETRIC HALF OPEN"
        ),
        "proved": {
            "conditional_ternary_support_at_least_2p_minus_4": True,
            "exceptional_direction_refinement": True,
            "total_fixed_parallel_scalar_ledgers_remain_feasible": True,
            "real_norm_closes_balanced_branch_C": False,
            "one_common_simple_graph_constructed": False,
            "equality_with_at_least_nine_active_rows_forces_projective_pencil": True,
            "two_zero_rows_alone_exclude_that_pencil": False,
            "antisymmetric_ternary_lift_proved_by_subsequent_mobius_trade": True,
            "coupled_symmetric_half_proved": False,
            "residual_ii_closed": False,
        },
        "equality_case": equality_case_projective_pencil_lemma(31, 9),
        "two_zero_row_counterexample": two_zero_direction_pencil_counterexample(31),
        "sample_certificates": samples,
        "duplicate_work_guard": (
            "Do not retry the scalar total/fixed/parallel or Euclidean norm "
            "bounds, and do not claim that two zero directions alone kill a "
            "full pencil.  The Mobius construction has already solved the "
            "entire antisymmetric ternary target; attack the coupled "
            "symmetric pair totals instead."
        ),
        "L_status": "OPEN",
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), sort_keys=True, indent=2))
