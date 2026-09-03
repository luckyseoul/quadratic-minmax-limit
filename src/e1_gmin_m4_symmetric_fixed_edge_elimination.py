#!/usr/bin/env python3
"""Parity-forced elimination of the fixed symmetric Radon variables.

After a ternary antisymmetric lift has been chosen, a central completion has
one binary variable on each fixed antipodal edge and one binary *double-orbit*
variable on each unused nonfixed inversion orbit.  This module records the
exact block-triangular reduction which removes all fixed-edge variables.

The proof is all-prime and symbolic.  ``exact_elimination_replay`` builds
small matrices only as a fail-when-wrong check; it is not evidence for the
theorem and performs no residual configuration census.
"""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations, product

from e1_gmin_m4_inversion_antisymmetric_radon import (
    Edge,
    Point,
    _edge,
    _functional_value,
    _negative_edge,
    projective_functionals,
)
from e1_gmin_m4_inversion_symmetric_lattice import (
    symmetric_dimensions,
    symmetric_mod2_decomposition,
)
from e1_gmin_m4_mobius_half_symmetric import (
    mobius_parameter_edges,
    paley_edge_sign,
)
from e1_gmin_m4_prop15721 import is_prime


TargetKey = tuple[object, ...]
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


def fixed_edge_elimination_theorem(p: int) -> dict[str, object]:
    """Return the exact all-prime fixed-edge elimination theorem.

    The source sign ``tau`` is allowed to be any sign invariant under central
    inversion.  This includes the Paley column sign used in branch C.  In a
    basis consisting of fixed edges and nonfixed orbit sums, and a target
    basis consisting of fixed cells and nonfixed cell-pairs, the plus-map is

        [[A, 2B],
         [0,  C]].

    Modulo two, ``A`` is an isomorphism onto the compatible fixed-cell
    residues.  Thus a central target uniquely determines the fixed binary
    edge vector; subtracting it and dividing the fixed block by two leaves an
    equivalent zero-one system solely on unused nonfixed orbits.
    """
    h = _check_odd_prime(p)
    d = p + 1
    dimensions = symmetric_dimensions(p)
    mod2 = symmetric_mod2_decomposition(p)

    fixed_source = d * h
    fixed_target_coordinates = d * (h + 1)
    fixed_compatibility_equations = d
    compatible_fixed_rank = fixed_target_coordinates - fixed_compatibility_equations
    paired_target_rank = d * h * h

    proved = bool(
        dimensions["source_fixed_antipodal_edge_rank"] == fixed_source
        and mod2["fixed_edge_map_rank"] == fixed_source
        and mod2["fixed_edge_map_injective"]
        and compatible_fixed_rank == fixed_source
        and mod2["paired_nonfixed_target_rank"] == paired_target_rank
        and mod2["symmetric_mod2_surjective"]
    )
    if not proved:
        raise ArithmeticError("the fixed-edge elimination dimensions changed")

    return {
        "p": p,
        "h": h,
        "d": d,
        "hypotheses": {
            "target": "integral compatible inversion-symmetric T_U",
            "source_sign": "tau_e in {+1,-1} with tau_Je=tau_e",
            "used_nonfixed_orbits": "frozen after the antisymmetric lift",
            "unused_nonfixed_orbits": "binary choices 0 or tau_O(e+Je)",
            "fixed_antipodal_edges": "binary choices 0 or tau_f f",
        },
        "block_form": "R_plus=[[A,2B],[0,C]]",
        "block_reasons": {
            "lower_left_zero": (
                "a fixed antipodal edge maps only to inversion-fixed target cells"
            ),
            "upper_right_even": (
                "e and Je hit the same fixed target cell with the same tau sign"
            ),
        },
        "fixed_antipodal_variables": fixed_source,
        "fixed_target_coordinates": fixed_target_coordinates,
        "fixed_compatibility_equations": fixed_compatibility_equations,
        "compatible_fixed_residue_rank": compatible_fixed_rank,
        "fixed_map_mod2_rank": mod2["fixed_edge_map_rank"],
        "fixed_map_mod2_isomorphism": True,
        "unique_fixed_binary_vector": (
            "a(T_U)=A_bar^{-1}((T_U)_fix mod 2)"
        ),
        "explicit_inverse": (
            "a_[v]=g_(L_v)(0)+sum_L g_L(L(v)^2) (mod 2), "
            "where L_v(v)=0"
        ),
        "divided_target": (
            "T_hat=(((T_U)_fix-A*a(T_U))/2,(T_U)_pair)"
        ),
        "divided_unused_column": "B_hat_O=((B_O)_fix/2,(B_O)_pair)",
        "exact_equivalence": (
            "restricted symmetric fibre nonempty iff there is b in "
            "{0,1}^{Omega\\U} with sum_O b_O B_hat_O=T_hat"
        ),
        "hamming_slice": "2*sum_O b_O=|H|-|U|-|a(T_U)|",
        "directionwise_parallel_slice": (
            "Paley p=3 mod 4 specialization: n_L=(P_L-u_L-f_L)/2, "
            "with 0<=n_L<=d*h^2-u_L for every projective row L"
        ),
        "fixed_variables_eliminated": fixed_source,
        "remaining_boolean_variables": "unused nonfixed inversion orbits only",
        "geometric_fixed_word": (
            "Phi(a,[delta])=0 if a is parallel to delta; otherwise it is "
            "the p antipodal classes [delta+c*a], c in F_p"
        ),
        "support_weight_objective": (
            "|U|+|a_Y+sum_(O in U) Phi(O)|; physical feasibility requires "
            "this to be at most |H|"
        ),
        "restricted_symmetric_fibre_closed": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def hamming_slice_identity(
    used_orbit_count: int,
    fixed_edge_bits: tuple[int, ...] | list[int],
    unused_double_orbit_bits: tuple[int, ...] | list[int],
) -> dict[str, int | bool | str]:
    """Return the exact physical-edge count of a restricted completion."""
    if (
        not isinstance(used_orbit_count, int)
        or isinstance(used_orbit_count, bool)
        or used_orbit_count < 0
    ):
        raise ValueError("used_orbit_count must be a nonnegative integer")
    if any(bit not in (0, 1) or isinstance(bit, bool) for bit in fixed_edge_bits):
        raise ValueError("fixed-edge coefficients must be binary")
    if any(
        bit not in (0, 1) or isinstance(bit, bool)
        for bit in unused_double_orbit_bits
    ):
        raise ValueError("double-orbit coefficients must be binary")

    fixed_weight = sum(fixed_edge_bits)
    double_weight = sum(unused_double_orbit_bits)
    physical_edges = used_orbit_count + fixed_weight + 2 * double_weight
    return {
        "used_single_orbit_edges": used_orbit_count,
        "fixed_antipodal_edge_weight": fixed_weight,
        "unused_double_orbit_weight": double_weight,
        "physical_edge_count": physical_edges,
        "slice_equation": "2*double_weight=|H|-used_weight-fixed_weight",
        "slice_left": 2 * double_weight,
        "slice_right": physical_edges - used_orbit_count - fixed_weight,
        "proved": 2 * double_weight == physical_edges - used_orbit_count - fixed_weight,
    }


def forced_fixed_word_parity_theorem(
    p: int,
    graph_edge_count: int,
    used_orbit_count: int,
    used_parallel_orbit_count: int,
) -> dict[str, object]:
    """Return the exact parity forced by the fixed-cell inverse.

    If ``g_L(0)`` is the parallel fixed-cell bit and ``g_L(beta)`` are
    the nonzero fixed-cell bits, summing the explicit inverse over all
    fixed source classes gives

        |a(g)| = sum_(L,beta!=0) g_L(beta)                 (mod 2).

    For the target of a graph with ``H`` edges, every row total is ``H``
    modulo two and the parallel totals sum to ``H``.  Since ``p+1`` is
    even, ``|a_Y|=H`` modulo two.  A used nonparallel orbit contributes a
    ``p``-point word ``Phi``, whereas a used parallel orbit has ``Phi=0``.
    Thus, writing ``u_np`` and ``u_parallel`` for the two counts,

        |a(T_U)| = H + u_np                              (mod 2),
        H-|U|-|a(T_U)| = u_parallel                      (mod 2).

    Hence an even number of used parallel orbits is necessary for the
    divided Hamming slice.  An omitted compact fixed-cell contribution
    cannot be replaced by zero merely from centrality.
    """
    _check_odd_prime(p)
    for name, value in (
        ("graph_edge_count", graph_edge_count),
        ("used_orbit_count", used_orbit_count),
        ("used_parallel_orbit_count", used_parallel_orbit_count),
    ):
        if (
            not isinstance(value, int)
            or isinstance(value, bool)
            or value < 0
        ):
            raise ValueError(f"{name} must be a nonnegative integer")
    if used_orbit_count > graph_edge_count:
        raise ValueError("used_orbit_count cannot exceed graph_edge_count")
    if used_parallel_orbit_count > used_orbit_count:
        raise ValueError(
            "used_parallel_orbit_count cannot exceed used_orbit_count"
        )

    target_parity = graph_edge_count % 2
    used_nonparallel_orbit_count = (
        used_orbit_count - used_parallel_orbit_count
    )
    remainder_parity = (
        graph_edge_count + used_nonparallel_orbit_count
    ) % 2
    hamming_numerator_parity = used_parallel_orbit_count % 2
    return {
        "p": p,
        "fixed_inverse_parity": (
            "|a(g)|=sum_(L,beta!=0) g_L(beta) mod 2"
        ),
        "graph_target_fixed_word_parity": target_parity,
        "used_nonparallel_orbit_count": used_nonparallel_orbit_count,
        "used_parallel_orbit_count": used_parallel_orbit_count,
        "central_remainder_fixed_word_parity": remainder_parity,
        "central_remainder_identity": (
            "|a(T_U)|=|H|+|U_nonparallel| mod 2"
        ),
        "hamming_slice_numerator_parity": hamming_numerator_parity,
        "hamming_slice_numerator_even_automatically": False,
        "even_used_parallel_orbit_count_necessary": True,
        "hamming_slice_parity_feasible": hamming_numerator_parity == 0,
        "compact_fixed_cell_word_may_be_omitted": False,
        "parity_excludes_symmetric_completion": hamming_numerator_parity == 1,
        "proved": True,
    }


def p31_mobius_cancellation_parity_ladder(
    t: int,
    cancellation_units: int,
    zero_phi_cancellation_units: int | None = None,
) -> dict[str, object]:
    """Return the cancellation-aware fixed-word parity ladder at ``p=31``.

    There are sixteen localized Mobius halves.  Each contributes thirty raw
    orbit occurrences: one has ``Phi=0`` because its midpoint is parallel to
    its difference, and the other twenty-nine have odd ``Phi`` words.  Here
    *zero-Phi parallel* always means that midpoint/difference condition.  It
    is unrelated to the directionwise parallel-edge counts in
    :func:`directionwise_parallel_slices`.

    At an orbit with ``n_O`` raw occurrences and final ternary coefficient
    ``c_O``, put ``kappa_O=(n_O-|c_O|)/2``.  Thus cancellations of either
    type remove occurrences two at a time, even when three or more halves
    meet at one orbit.  If ``kappa_0`` is supplied, the exact post-cancellation
    counts are

        u_0  = 16  - 2*kappa_0,
        u_np = 464 - 2*(kappa-kappa_0).

    Without that optional split, their even parity is still forced.  Along
    the branch-C ray ``68 <= t <= 177``, the target has ``H=125+2t`` edges
    and the used support has ``480-2*kappa`` orbits.  Consequently the fixed
    word has odd weight and the Hamming numerator has even parity for every
    ternary cancellation pattern.  This is a necessary ledger only; it does
    not construct a completion or close residual (ii).
    """
    for name, value in (
        ("t", t),
        ("cancellation_units", cancellation_units),
    ):
        if not isinstance(value, int) or isinstance(value, bool):
            raise ValueError(f"{name} must be an integer")
    if not 68 <= t <= 177:
        raise ValueError("the p=31 branch-C ray requires 68<=t<=177")
    if not 0 <= cancellation_units <= 240:
        raise ValueError("cancellation_units must lie between 0 and 240")

    if zero_phi_cancellation_units is not None:
        if (
            not isinstance(zero_phi_cancellation_units, int)
            or isinstance(zero_phi_cancellation_units, bool)
        ):
            raise ValueError("zero_phi_cancellation_units must be an integer")
        if not 0 <= zero_phi_cancellation_units <= 8:
            raise ValueError(
                "zero_phi_cancellation_units must lie between 0 and 8"
            )
        if zero_phi_cancellation_units > cancellation_units:
            raise ValueError(
                "zero-Phi cancellation units cannot exceed total cancellations"
            )
        if cancellation_units - zero_phi_cancellation_units > 232:
            raise ValueError(
                "nonzero-Phi cancellation units cannot exceed 232"
            )

    graph_edge_count = 125 + 2 * t
    used_orbit_count = 480 - 2 * cancellation_units
    minimum_cancellation_units = 178 - t
    remaining_edge_capacity = graph_edge_count - used_orbit_count
    support_size_feasible = remaining_edge_capacity >= 0
    at_minimum_cancellation = (
        cancellation_units == minimum_cancellation_units
    )

    used_zero_phi_count: int | None = None
    used_nonzero_phi_count: int | None = None
    nonzero_phi_cancellation_units: int | None = None
    if zero_phi_cancellation_units is not None:
        nonzero_phi_cancellation_units = (
            cancellation_units - zero_phi_cancellation_units
        )
        used_zero_phi_count = 16 - 2 * zero_phi_cancellation_units
        used_nonzero_phi_count = (
            464 - 2 * nonzero_phi_cancellation_units
        )
        if used_zero_phi_count + used_nonzero_phi_count != used_orbit_count:
            raise ArithmeticError("the split cancellation ledger changed")

    fixed_word_parity = 1
    hamming_numerator_parity = 0
    minimum_fixed_word_weight = 1 if at_minimum_cancellation else None
    minimum_unused_double_orbit_count = (
        0 if at_minimum_cancellation else None
    )
    return {
        "p": 31,
        "t": t,
        "half_count": 16,
        "raw_orbit_occurrences": 480,
        "raw_zero_phi_occurrences": 16,
        "raw_nonzero_phi_occurrences": 464,
        "cancellation_units": cancellation_units,
        "zero_phi_cancellation_units": zero_phi_cancellation_units,
        "nonzero_phi_cancellation_units": nonzero_phi_cancellation_units,
        "used_orbit_count": used_orbit_count,
        "used_zero_phi_orbit_count": used_zero_phi_count,
        "used_nonzero_phi_orbit_count": used_nonzero_phi_count,
        "used_zero_phi_count_formula": "16-2*kappa_0",
        "used_nonzero_phi_count_formula": "464-2*(kappa-kappa_0)",
        "used_zero_phi_orbit_count_parity": 0,
        "used_nonzero_phi_orbit_count_parity": 0,
        "terminology": (
            "zero-Phi parallel means midpoint parallel to difference; it is "
            "not a directionwise parallel-edge count"
        ),
        "graph_edge_count": graph_edge_count,
        "minimum_cancellation_units_for_support_size": (
            minimum_cancellation_units
        ),
        "support_size_feasible": support_size_feasible,
        "remaining_edge_capacity": remaining_edge_capacity,
        "remaining_edge_capacity_formula": (
            "1+2*(kappa-(178-t))"
        ),
        "fixed_word_weight_parity": fixed_word_parity,
        "fixed_word_weight_is_odd_for_every_ternary_cancellation_pattern": (
            True
        ),
        "hamming_slice_numerator_parity": hamming_numerator_parity,
        "hamming_slice_parity_is_automatic_for_sixteen_halves": True,
        "at_minimum_cancellation": at_minimum_cancellation,
        "conditional_fixed_word_weight_if_completion": (
            minimum_fixed_word_weight
        ),
        "conditional_unused_double_orbit_count_if_completion": (
            minimum_unused_double_orbit_count
        ),
        "parity_excludes_mobius_completion": False,
        "symmetric_completion_constructed": False,
        "residual_ii_closed": False,
        "proved": True,
    }


def directionwise_parallel_slices(
    p: int,
    target_parallel_counts: tuple[int, ...] | list[int],
    used_single_parallel_counts: tuple[int, ...] | list[int],
    forced_fixed_parallel_counts: tuple[int, ...] | list[int],
) -> dict[str, object]:
    """Solve the exact parallel-coordinate slices of the reduced box.

    The entries are physical parallel counts in the normalized Paley rows.
    An edge parallel to row ``L`` has column sign ``tau_e=epsilon_L``;
    multiplying the row by ``epsilon_L`` therefore makes its contribution
    ``+1``.  A used orbit contributes one edge, a forced fixed edge one edge,
    and a selected unused symmetric orbit two edges.
    """
    h = _check_odd_prime(p)
    if p % 4 != 3:
        raise ValueError(
            "positive normalized Paley parallel slices require p=3 mod 4"
        )
    d = p + 1
    vectors = (
        target_parallel_counts,
        used_single_parallel_counts,
        forced_fixed_parallel_counts,
    )
    if any(len(vector) != d for vector in vectors):
        raise ValueError("parallel-count vectors must have p+1 entries")
    if any(
        not isinstance(value, int) or isinstance(value, bool) or value < 0
        for vector in vectors
        for value in vector
    ):
        raise ValueError("parallel counts must be nonnegative integers")

    total_nonfixed_orbit_capacity_per_direction = d * h * h
    if any(value > total_nonfixed_orbit_capacity_per_direction for value in used_single_parallel_counts):
        raise ValueError("used-orbit count exceeds its direction capacity")
    if any(value > h for value in forced_fixed_parallel_counts):
        raise ValueError("fixed-edge count exceeds its direction capacity")

    numerators = tuple(
        target - used - fixed
        for target, used, fixed in zip(
            target_parallel_counts,
            used_single_parallel_counts,
            forced_fixed_parallel_counts,
            strict=True,
        )
    )
    integral = tuple(value % 2 == 0 for value in numerators)
    selected = tuple(
        value // 2 if is_integral else None
        for value, is_integral in zip(numerators, integral, strict=True)
    )
    unused_capacities = tuple(
        total_nonfixed_orbit_capacity_per_direction - used
        for used in used_single_parallel_counts
    )
    nonnegative = tuple(
        value is not None and value >= 0 for value in selected
    )
    within_capacity = tuple(
        value is not None and 0 <= value <= capacity
        for value, capacity in zip(selected, unused_capacities, strict=True)
    )
    feasible = all(integral) and all(nonnegative) and all(within_capacity)
    selected_total = sum(value for value in selected if value is not None)
    global_numerator = (
        sum(target_parallel_counts)
        - sum(used_single_parallel_counts)
        - sum(forced_fixed_parallel_counts)
    )
    global_slice_matches = bool(
        feasible and 2 * selected_total == global_numerator
    )
    return {
        "p": p,
        "parallel_direction_count": d,
        "nonfixed_orbit_capacity_per_direction": (
            total_nonfixed_orbit_capacity_per_direction
        ),
        "target_parallel_counts": list(target_parallel_counts),
        "used_single_parallel_counts": list(used_single_parallel_counts),
        "forced_fixed_parallel_counts": list(forced_fixed_parallel_counts),
        "slice_numerators": list(numerators),
        "slice_integral": list(integral),
        "selected_unused_double_orbits": list(selected),
        "unused_double_orbit_capacities": list(unused_capacities),
        "slice_nonnegative": list(nonnegative),
        "slice_within_capacity": list(within_capacity),
        "all_direction_slices_feasible": feasible,
        "selected_unused_double_orbit_total": selected_total,
        "global_hamming_slice_recovered_by_summing_directions": (
            global_slice_matches
        ),
        "normalization_reason": (
            "for an edge parallel to L, epsilon_L*tau_e=1"
        ),
        "proved": True,
    }


def _antipodal_vector_class(p: int, point: Point) -> Point:
    point = (point[0] % p, point[1] % p)
    if point == (0, 0):
        raise ValueError("zero has no nonzero antipodal vector class")
    negative = ((-point[0]) % p, (-point[1]) % p)
    return min(point, negative)


def orbit_fixed_word(p: int, edge: Edge) -> dict[str, object]:
    """Return the fixed-edge parity word toggled by one nonfixed orbit.

    Write the edge as ``{a-delta,a+delta}``.  Applying the explicit inverse
    of the fixed-edge map to this selected edge's fixed-cell incidence gives
    zero when ``a`` and ``delta`` are parallel.  Otherwise its support is the
    affine line ``{[delta+c*a]:c in F_p}`` in the antipodal vector quotient.
    The result is unchanged by replacing the edge by its central negative or
    by replacing ``delta`` with ``-delta``.
    """
    _check_odd_prime(p)
    if len(edge) != 2 or len(edge[0]) != 2 or len(edge[1]) != 2:
        raise ValueError("an edge must contain two planar points")
    normalized = _edge(
        (edge[0][0] % p, edge[0][1] % p),
        (edge[1][0] % p, edge[1][1] % p),
    )
    if normalized[0] == normalized[1]:
        raise ValueError("source edges cannot be loops")
    if _negative_edge(p, normalized) == normalized:
        raise ValueError("the orbit must be nonfixed under central inversion")

    inverse_two = pow(2, -1, p)
    midpoint = (
        (normalized[0][0] + normalized[1][0]) * inverse_two % p,
        (normalized[0][1] + normalized[1][1]) * inverse_two % p,
    )
    difference = (
        (normalized[1][0] - normalized[0][0]) * inverse_two % p,
        (normalized[1][1] - normalized[0][1]) * inverse_two % p,
    )
    determinant = (
        midpoint[0] * difference[1] - midpoint[1] * difference[0]
    ) % p
    if determinant == 0:
        support: tuple[Point, ...] = ()
    else:
        support = tuple(sorted({
            _antipodal_vector_class(
                p,
                (
                    difference[0] + scalar * midpoint[0],
                    difference[1] + scalar * midpoint[1],
                ),
            )
            for scalar in range(p)
        }))
        if len(support) != p:
            raise ArithmeticError("the affine fixed word lost an antipodal class")

    negative = _negative_edge(p, normalized)
    negative_record_midpoint = (
        (negative[0][0] + negative[1][0]) * inverse_two % p,
        (negative[0][1] + negative[1][1]) * inverse_two % p,
    )
    orientation_invariant = bool(
        negative_record_midpoint
        == ((-midpoint[0]) % p, (-midpoint[1]) % p)
    )
    return {
        "p": p,
        "midpoint": list(midpoint),
        "difference_representative": list(difference),
        "midpoint_parallel_to_difference": determinant == 0,
        "fixed_word_support": [list(point) for point in support],
        "fixed_word_weight": len(support),
        "formula": (
            "0 if a||delta; otherwise {[delta+c*a]:c in F_p} "
            "in (V minus {0})/{+1,-1}"
        ),
        "independent_of_selected_orbit_side": orientation_invariant,
        "independent_of_difference_sign": True,
        "proved": bool(
            orientation_invariant and len(support) in (0, p)
        ),
    }


def fixed_word_block_basis_theorem(p: int) -> dict[str, object]:
    """Return the exact paired-affine-line block design over F2.

    The nonzero words from :func:`orbit_fixed_word` are the incidence
    vectors of paired non-origin affine lines in ``V/{+1,-1}``.  There are
    exactly as many block types as antipodal nonzero point classes.  Their
    square incidence matrix is orthogonal over ``F_2``.
    """
    h = _check_odd_prime(p)
    d = p + 1
    antipodal_classes = d * h
    block_types = antipodal_classes
    blocks_per_point = p
    points_per_block = p
    noncollinear_common_blocks = 2
    collinear_distinct_common_blocks = 0

    total_nonfixed_source_orbits = antipodal_classes * antipodal_classes
    zero_word_orbits = antipodal_classes * h
    nonzero_word_orbits = total_nonfixed_source_orbits - zero_word_orbits
    multiplicity_per_nonzero_block = p * h

    proved = bool(
        block_types == antipodal_classes
        and blocks_per_point % 2 == 1
        and noncollinear_common_blocks % 2 == 0
        and collinear_distinct_common_blocks % 2 == 0
        and nonzero_word_orbits == block_types * multiplicity_per_nonzero_block
        and antipodal_classes - h == p * h
    )
    if not proved:
        raise ArithmeticError("the fixed-word block design counts changed")
    return {
        "p": p,
        "h": h,
        "d": d,
        "antipodal_nonzero_point_classes": antipodal_classes,
        "paired_affine_line_block_types": block_types,
        "points_per_block": points_per_block,
        "blocks_per_point": blocks_per_point,
        "common_blocks_for_distinct_noncollinear_classes": (
            noncollinear_common_blocks
        ),
        "common_blocks_for_distinct_collinear_classes": (
            collinear_distinct_common_blocks
        ),
        "binary_gram_identity": "M*M^T=I, hence M^T*M=I",
        "block_vectors_form_basis": True,
        "total_nonfixed_source_orbits": total_nonfixed_source_orbits,
        "zero_word_orbits": zero_word_orbits,
        "nonzero_word_orbits": nonzero_word_orbits,
        "multiplicity_per_nonzero_block_type": multiplicity_per_nonzero_block,
        "disjoint_C_kernel_lifts_per_block_type": h,
        "columns_per_C_kernel_lift": p,
        "C_kernel_lift_construction": (
            "fix one of the h midpoint antipodal classes [a] in the block "
            "direction A and range [delta] over the p classes of B(A,c)"
        ),
        "zero_word_condition": "midpoint a is parallel to difference delta",
        "line_code_coordinates": (
            "if c(U) records the parity of selected orbits of each nonzero "
            "block type, then Phi(U)=M*c(U) and c(U)=M^T*Phi(U)"
        ),
        "proved": proved,
    }


def mobius_midpoint_direction_theorem(p: int) -> dict[str, object]:
    """Return the exact midpoint-direction profile of one Mobius half.

    In normalized ``(L,M)`` coordinates, the midpoint of parameter ``t`` has
    finite slope ``1-1/(t+1)^2``.  As ``t+1`` ranges through ``F_p^*``, each
    square occurs twice.  Thus the ``p-1`` midpoint orbits occupy exactly
    ``h`` spatial directions, two in each.
    """
    h = _check_odd_prime(p)
    proved = bool(p - 1 == 2 * h)
    if not proved:
        raise ArithmeticError("the Mobius midpoint-direction count changed")
    return {
        "p": p,
        "parameter_domain": "t in F_p minus {-1}",
        "normalized_midpoint_L_coordinate": "(1+t)/2",
        "normalized_midpoint_M_coordinate": "t*(t+2)/(2*(t+1))",
        "midpoint_slope_M_over_L": "1-1/(t+1)^2",
        "vertical_midpoint_direction_occurs": False,
        "distinct_midpoint_directions": h,
        "parameters_per_midpoint_direction": 2,
        "one_half_hits_any_midpoint_direction_at_most": 2,
        "proved": proved,
    }


def exact_mobius_midpoint_replay(p: int) -> dict[str, object]:
    """Replay the midpoint slope and exact two-to-one fibres at p=3,7."""
    _check_odd_prime(p)
    if p not in (3, 7):
        raise ValueError("the exact Mobius midpoint replay is limited to p=3 or p=7")
    edges = mobius_parameter_edges(p, (1, 0), (0, 1), 1)
    inverse_two = pow(2, -1, p)
    direction_counts: defaultdict[int, int] = defaultdict(int)
    formula_holds = True
    for parameter, edge in edges.items():
        midpoint = (
            (edge[0][0] + edge[1][0]) * inverse_two % p,
            (edge[0][1] + edge[1][1]) * inverse_two % p,
        )
        if midpoint[0] == 0:
            formula_holds = False
            break
        slope = midpoint[1] * pow(midpoint[0], -1, p) % p
        expected = (
            1 - pow(parameter + 1, -2, p)
        ) % p
        if slope != expected:
            formula_holds = False
            break
        direction_counts[slope] += 1
    h = (p - 1) // 2
    proved = bool(
        formula_holds
        and len(edges) == p - 1
        and len(direction_counts) == h
        and set(direction_counts.values()) == {2}
    )
    if not proved:
        raise ArithmeticError("the exact Mobius midpoint replay changed")
    return {
        "p": p,
        "parameter_count": len(edges),
        "midpoint_slope_formula_holds": formula_holds,
        "distinct_midpoint_directions": len(direction_counts),
        "multiplicity_set": sorted(set(direction_counts.values())),
        "role": "fail-when-wrong small-prime replay, not theorem evidence",
        "proved": proved,
    }


def _target_layout(
    p: int,
) -> tuple[tuple[TargetKey, ...], tuple[tuple[TargetKey, TargetKey], ...]]:
    """Return fixed target cells and canonical nonfixed target-cell pairs."""
    fixed: list[TargetKey] = []
    paired: list[tuple[TargetKey, TargetKey]] = []
    for direction_index in range(p + 1):
        fixed.append(("P", direction_index))
        for left, right in combinations(range(p), 2):
            key: TargetKey = ("K", direction_index, left, right)
            negative_values = sorted(((-left) % p, (-right) % p))
            negative: TargetKey = (
                "K",
                direction_index,
                negative_values[0],
                negative_values[1],
            )
            if key == negative:
                fixed.append(key)
            elif key < negative:
                paired.append((key, negative))
    return tuple(fixed), tuple(paired)


def _source_orbits(p: int) -> tuple[tuple[Edge, ...], tuple[tuple[Edge, Edge], ...]]:
    points = tuple(product(range(p), repeat=2))
    edges = tuple(_edge(points[i], points[j]) for i, j in combinations(range(p * p), 2))
    visited: set[Edge] = set()
    fixed: list[Edge] = []
    paired: list[tuple[Edge, Edge]] = []
    for edge in edges:
        if edge in visited:
            continue
        negative = _negative_edge(p, edge)
        if negative == edge:
            fixed.append(edge)
            visited.add(edge)
        else:
            first, second = sorted((edge, negative))
            paired.append((first, second))
            visited.update((first, second))
    return tuple((edge,) for edge in fixed), tuple(paired)


def _fixed_target_key_for_value(
    p: int, direction_index: int, value: int
) -> TargetKey:
    value %= p
    if value == 0:
        return ("P", direction_index)
    left, right = sorted((value, (-value) % p))
    return ("K", direction_index, left, right)


def _annihilating_direction_index(p: int, point: Point) -> int:
    hits = tuple(
        index
        for index, functional in enumerate(projective_functionals(p))
        if _functional_value(p, functional, point) == 0
    )
    if len(hits) != 1:
        raise ArithmeticError("a nonzero point lost its unique annihilating row")
    return hits[0]


def _explicit_fixed_inverse(
    p: int,
    fixed_orbits: tuple[tuple[Edge, ...], ...],
    fixed_target_by_key: dict[TargetKey, int],
) -> tuple[int, ...]:
    """Apply ``a_[v]=g_Lv(0)+sum_L g_L(L(v)^2)`` over F2."""
    directions = projective_functionals(p)
    recovered: list[int] = []
    for (edge,) in fixed_orbits:
        point = edge[0]
        if point == (0, 0):
            point = edge[1]
        annihilator = _annihilating_direction_index(p, point)
        value = fixed_target_by_key.get(("P", annihilator), 0)
        for direction_index, functional in enumerate(directions):
            projected = _functional_value(p, functional, point)
            key = _fixed_target_key_for_value(p, direction_index, projected)
            value += fixed_target_by_key.get(key, 0)
        recovered.append(value & 1)
    return tuple(recovered)


def _parallel_direction_index(p: int, edge: Edge) -> int:
    directions = projective_functionals(p)
    hits = tuple(
        index
        for index, functional in enumerate(directions)
        if _functional_value(p, functional, edge[0])
        == _functional_value(p, functional, edge[1])
    )
    if len(hits) != 1:
        raise ArithmeticError("an edge lost its unique parallel direction")
    return hits[0]


def _ordinary_image(p: int, source: dict[Edge, int]) -> dict[TargetKey, int]:
    image: defaultdict[TargetKey, int] = defaultdict(int)
    directions = projective_functionals(p)
    for edge, coefficient in source.items():
        for direction_index, functional in enumerate(directions):
            left = _functional_value(p, functional, edge[0])
            right = _functional_value(p, functional, edge[1])
            if left == right:
                key: TargetKey = ("P", direction_index)
            else:
                first, second = sorted((left, right))
                key = ("K", direction_index, first, second)
            image[key] += coefficient
    return {key: value for key, value in image.items() if value}


def _central_coordinates(
    image: dict[TargetKey, int],
    fixed_keys: tuple[TargetKey, ...],
    paired_keys: tuple[tuple[TargetKey, TargetKey], ...],
) -> tuple[IntegerVector, IntegerVector]:
    fixed = tuple(image.get(key, 0) for key in fixed_keys)
    paired: list[int] = []
    for first, second in paired_keys:
        first_value = image.get(first, 0)
        second_value = image.get(second, 0)
        if first_value != second_value:
            raise ArithmeticError("a symmetric source produced a noncentral target")
        paired.append(first_value)
    return fixed, tuple(paired)


def _column_add(
    columns: tuple[tuple[IntegerVector, IntegerVector], ...],
    bits: tuple[int, ...],
) -> tuple[IntegerVector, IntegerVector]:
    if len(columns) != len(bits):
        raise ValueError("column and bit counts differ")
    if not columns:
        return (), ()
    fixed_size = len(columns[0][0])
    paired_size = len(columns[0][1])
    fixed = [0] * fixed_size
    paired = [0] * paired_size
    for bit, (fixed_column, paired_column) in zip(bits, columns, strict=True):
        if bit not in (0, 1):
            raise ValueError("column coefficients must be binary")
        if bit:
            for index, value in enumerate(fixed_column):
                fixed[index] += value
            for index, value in enumerate(paired_column):
                paired[index] += value
    return tuple(fixed), tuple(paired)


def _binary_vector(values: IntegerVector) -> int:
    out = 0
    for index, value in enumerate(values):
        if value & 1:
            out |= 1 << index
    return out


def _binary_rank(vectors: tuple[int, ...]) -> int:
    basis: dict[int, int] = {}
    for original in vectors:
        value = original
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                value ^= basis[pivot]
            else:
                basis[pivot] = value
                break
    return len(basis)


def _xor_all(vectors: list[int] | tuple[int, ...]) -> int:
    out = 0
    for vector in vectors:
        out ^= vector
    return out


def _solve_independent_binary_columns(columns: tuple[int, ...], target: int) -> int:
    """Solve a target in the span of independent binary columns."""
    basis: dict[int, tuple[int, int]] = {}
    for column_index, original in enumerate(columns):
        value = original
        combination = 1 << column_index
        while value:
            pivot = value.bit_length() - 1
            if pivot in basis:
                basis_value, basis_combination = basis[pivot]
                value ^= basis_value
                combination ^= basis_combination
            else:
                basis[pivot] = (value, combination)
                break
        if value == 0:
            raise ArithmeticError("the fixed columns lost independence")

    value = target
    combination = 0
    while value:
        pivot = value.bit_length() - 1
        if pivot not in basis:
            raise ArithmeticError("the target fixed residue left the fixed-edge image")
        basis_value, basis_combination = basis[pivot]
        value ^= basis_value
        combination ^= basis_combination
    return combination


def exact_fixed_word_design_replay(p: int) -> dict[str, object]:
    """Replay the orbit-word formula and block orthogonality for p=3,7."""
    _check_odd_prime(p)
    if p not in (3, 7):
        raise ValueError("the exact fixed-word replay is limited to p=3 or p=7")

    h = (p - 1) // 2
    d = p + 1
    class_count = d * h
    fixed_keys, paired_keys = _target_layout(p)
    fixed_orbits, nonfixed_orbits = _source_orbits(p)
    fixed_classes = tuple(
        _antipodal_vector_class(p, edge[0]) for (edge,) in fixed_orbits
    )
    if len(set(fixed_classes)) != class_count:
        raise ArithmeticError("the fixed edges lost an antipodal vector class")

    block_multiplicities: defaultdict[frozenset[Point], int] = defaultdict(int)
    kernel_lift_columns: defaultdict[
        tuple[frozenset[Point], Point], list[int]
    ] = defaultdict(list)
    zero_words = 0
    direct_formula_agree = True
    for edge, negative in nonfixed_orbits:
        formula_record = orbit_fixed_word(p, edge)
        formula_support = frozenset(
            tuple(point) for point in formula_record["fixed_word_support"]
        )

        image = _ordinary_image(p, {edge: 1})
        fixed_target_by_key = {key: image.get(key, 0) for key in fixed_keys}
        direct_bits = _explicit_fixed_inverse(
            p, fixed_orbits, fixed_target_by_key
        )
        direct_support = frozenset(
            point
            for point, bit in zip(fixed_classes, direct_bits, strict=True)
            if bit
        )
        if direct_support != formula_support:
            direct_formula_agree = False
            break
        if formula_support:
            block_multiplicities[formula_support] += 1
            midpoint_class = _antipodal_vector_class(
                p, tuple(formula_record["midpoint"])
            )
            symmetric_image = _ordinary_image(p, {edge: 1, negative: 1})
            _fixed_part, paired_part = _central_coordinates(
                symmetric_image, fixed_keys, paired_keys
            )
            kernel_lift_columns[(formula_support, midpoint_class)].append(
                _binary_vector(paired_part)
            )
        else:
            zero_words += 1

    blocks = tuple(
        sorted(block_multiplicities, key=lambda block: tuple(sorted(block)))
    )
    point_classes = tuple(sorted(set(fixed_classes)))
    point_rows = tuple(
        sum(
            (1 << block_index)
            for block_index, block in enumerate(blocks)
            if point in block
        )
        for point in point_classes
    )
    block_columns = tuple(
        sum(
            (1 << point_index)
            for point_index, point in enumerate(point_classes)
            if point in block
        )
        for block in blocks
    )

    exact_intersections = True
    binary_gram_identity = True
    for first_index, first in enumerate(point_classes):
        for second_index, second in enumerate(point_classes):
            common = (point_rows[first_index] & point_rows[second_index]).bit_count()
            if first_index == second_index:
                expected = p
            else:
                determinant = (
                    first[0] * second[1] - first[1] * second[0]
                ) % p
                expected = 0 if determinant == 0 else 2
            if common != expected:
                exact_intersections = False
            if common % 2 != int(first_index == second_index):
                binary_gram_identity = False

    expected_zero = class_count * h
    expected_multiplicity = p * h
    kernel_lifts_exact = bool(
        len(kernel_lift_columns) == class_count * h
        and all(len(columns) == p for columns in kernel_lift_columns.values())
        and all(
            _xor_all(columns) == 0
            for columns in kernel_lift_columns.values()
        )
        and all(
            sum(1 for block, _midpoint in kernel_lift_columns if block == target_block)
            == h
            for target_block in blocks
        )
    )
    proved = bool(
        direct_formula_agree
        and len(nonfixed_orbits) == class_count * class_count
        and zero_words == expected_zero
        and len(blocks) == class_count
        and set(block_multiplicities.values()) == {expected_multiplicity}
        and all(len(block) == p for block in blocks)
        and all(row.bit_count() == p for row in point_rows)
        and exact_intersections
        and binary_gram_identity
        and _binary_rank(block_columns) == class_count
        and kernel_lifts_exact
    )
    if not proved:
        raise ArithmeticError("the exact fixed-word block replay changed")
    return {
        "p": p,
        "antipodal_point_classes": len(point_classes),
        "nonfixed_source_orbits": len(nonfixed_orbits),
        "zero_word_orbits": zero_words,
        "nonzero_block_types": len(blocks),
        "multiplicity_per_nonzero_block_type": expected_multiplicity,
        "points_per_block": p,
        "blocks_per_point": p,
        "direct_inverse_matches_affine_line_formula": direct_formula_agree,
        "exact_pair_intersection_numbers": exact_intersections,
        "binary_gram_is_identity": binary_gram_identity,
        "binary_block_rank": _binary_rank(block_columns),
        "disjoint_C_kernel_lifts_per_block_type": h,
        "columns_per_C_kernel_lift": p,
        "exact_C_kernel_lifts": kernel_lifts_exact,
        "role": "fail-when-wrong small-prime replay, not theorem evidence",
        "proved": proved,
    }


def exact_elimination_replay(p: int) -> dict[str, object]:
    """Build one small exact plus-map replay; intentionally restricted to p<=7."""
    _check_odd_prime(p)
    if p % 4 != 3 or p > 7:
        raise ValueError("the Paley replay is intentionally limited to p=3 or p=7")

    h = (p - 1) // 2
    d = p + 1
    fixed_keys, paired_keys = _target_layout(p)
    fixed_orbits, nonfixed_orbits = _source_orbits(p)

    fixed_columns: list[tuple[IntegerVector, IntegerVector]] = []
    for (edge,) in fixed_orbits:
        tau = paley_edge_sign(p, edge)
        fixed_columns.append(
            _central_coordinates(
                _ordinary_image(p, {edge: tau}), fixed_keys, paired_keys
            )
        )

    nonfixed_columns: list[tuple[IntegerVector, IntegerVector]] = []
    for edge, negative in nonfixed_orbits:
        tau = paley_edge_sign(p, edge)
        if paley_edge_sign(p, negative) != tau:
            raise ArithmeticError("the Paley sign changed under inversion")
        nonfixed_columns.append(
            _central_coordinates(
                _ordinary_image(p, {edge: tau, negative: tau}),
                fixed_keys,
                paired_keys,
            )
        )

    fixed_columns_tuple = tuple(fixed_columns)
    nonfixed_columns_tuple = tuple(nonfixed_columns)
    lower_left_zero = all(not any(column[1]) for column in fixed_columns_tuple)
    upper_right_even = all(
        all(value % 2 == 0 for value in column[0])
        for column in nonfixed_columns_tuple
    )
    fixed_binary_columns = tuple(
        _binary_vector(column[0]) for column in fixed_columns_tuple
    )
    fixed_rank = _binary_rank(fixed_binary_columns)

    fixed_bits = tuple(int(index % 3 == 1) for index in range(len(fixed_columns_tuple)))
    used_count = min(p - 1, len(nonfixed_columns_tuple))
    unused_bits = tuple(
        0 if index < used_count else int(index % (p + 2) in (1, 3))
        for index in range(len(nonfixed_columns_tuple))
    )
    fixed_target = _column_add(fixed_columns_tuple, fixed_bits)
    nonfixed_target = _column_add(nonfixed_columns_tuple, unused_bits)
    target_fixed = tuple(
        first + second
        for first, second in zip(fixed_target[0], nonfixed_target[0], strict=True)
    )
    target_paired = tuple(
        first + second
        for first, second in zip(fixed_target[1], nonfixed_target[1], strict=True)
    )

    recovered_mask = _solve_independent_binary_columns(
        fixed_binary_columns, _binary_vector(target_fixed)
    )
    recovered_bits = tuple(
        (recovered_mask >> index) & 1 for index in range(len(fixed_columns_tuple))
    )
    fixed_target_by_key = {
        key: value for key, value in zip(fixed_keys, target_fixed, strict=True)
    }
    explicit_recovered_bits = _explicit_fixed_inverse(
        p, fixed_orbits, fixed_target_by_key
    )
    recovered_fixed_target = _column_add(fixed_columns_tuple, recovered_bits)
    fixed_remainder = tuple(
        target - recovered
        for target, recovered in zip(
            target_fixed, recovered_fixed_target[0], strict=True
        )
    )
    if any(value % 2 for value in fixed_remainder):
        raise ArithmeticError("fixed-edge subtraction did not leave an even block")
    divided_target = (
        tuple(value // 2 for value in fixed_remainder),
        target_paired,
    )
    divided_columns = tuple(
        (
            tuple(value // 2 for value in fixed_column),
            paired_column,
        )
        for fixed_column, paired_column in nonfixed_columns_tuple
    )
    divided_reconstruction = _column_add(divided_columns, unused_bits)
    hamming = hamming_slice_identity(used_count, recovered_bits, unused_bits)

    used_parallel = [0] * d
    fixed_parallel = [0] * d
    double_parallel = [0] * d
    for index, (edge, _negative) in enumerate(nonfixed_orbits):
        direction_index = _parallel_direction_index(p, edge)
        if index < used_count:
            used_parallel[direction_index] += 1
        elif unused_bits[index]:
            double_parallel[direction_index] += 1
    for bit, (edge,) in zip(recovered_bits, fixed_orbits, strict=True):
        if bit:
            fixed_parallel[_parallel_direction_index(p, edge)] += 1
    target_parallel = tuple(
        used + fixed + 2 * double
        for used, fixed, double in zip(
            used_parallel, fixed_parallel, double_parallel, strict=True
        )
    )
    parallel_slices = directionwise_parallel_slices(
        p,
        target_parallel,
        tuple(used_parallel),
        tuple(fixed_parallel),
    )

    expected_fixed = d * h
    expected_fixed_coordinates = d * (h + 1)
    expected_paired = d * h * h
    proved = bool(
        len(fixed_orbits) == expected_fixed
        and len(fixed_keys) == expected_fixed_coordinates
        and len(paired_keys) == expected_paired
        and lower_left_zero
        and upper_right_even
        and fixed_rank == expected_fixed
        and recovered_bits == fixed_bits
        and explicit_recovered_bits == fixed_bits
        and divided_target == divided_reconstruction
        and hamming["proved"]
        and parallel_slices["all_direction_slices_feasible"]
        and parallel_slices["selected_unused_double_orbits"] == double_parallel
        and parallel_slices[
            "global_hamming_slice_recovered_by_summing_directions"
        ]
    )
    if not proved:
        raise ArithmeticError("the exact fixed-edge elimination replay changed")

    return {
        "p": p,
        "fixed_source_columns": len(fixed_orbits),
        "nonfixed_pair_columns": len(nonfixed_orbits),
        "fixed_target_coordinates": len(fixed_keys),
        "paired_target_coordinates": len(paired_keys),
        "lower_left_block_zero": lower_left_zero,
        "upper_right_block_even": upper_right_even,
        "fixed_map_mod2_rank": fixed_rank,
        "fixed_binary_vector_recovered_uniquely": recovered_bits == fixed_bits,
        "explicit_inverse_formula_recovers_fixed_vector": (
            explicit_recovered_bits == fixed_bits
        ),
        "fixed_remainder_even": True,
        "divided_target_reconstructed_exactly": divided_target == divided_reconstruction,
        "used_nonfixed_orbits_frozen": used_count,
        "hamming_slice": hamming,
        "directionwise_parallel_slices": parallel_slices,
        "role": "fail-when-wrong small-prime replay, not theorem evidence",
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    theorem = fixed_edge_elimination_theorem(p)
    block_basis = fixed_word_block_basis_theorem(p)
    midpoint_directions = mobius_midpoint_direction_theorem(p)
    fixed_word_parity = forced_fixed_word_parity_theorem(p, 481, 480, 16)
    p31_cancellation_parity = (
        p31_mobius_cancellation_parity_ladder(177, 1, 1)
        if p == 31
        else None
    )
    replay = exact_elimination_replay(7)
    block_replay = exact_fixed_word_design_replay(7)
    midpoint_replay = exact_mobius_midpoint_replay(7)
    proved = bool(
        theorem["proved"]
        and block_basis["proved"]
        and midpoint_directions["proved"]
        and fixed_word_parity["proved"]
        and (
            p31_cancellation_parity is None
            or p31_cancellation_parity["proved"]
        )
        and replay["proved"]
        and block_replay["proved"]
        and midpoint_replay["proved"]
    )
    if not proved:
        raise ArithmeticError("the fixed-edge elimination theorem record changed")
    return {
        "title": "Parity-forced fixed-edge elimination in the symmetric fibre",
        "status": "PROVED EXACT REDUCTION; REDUCED BOOLEAN FIBRE OPEN",
        "theorem": theorem,
        "fixed_word_block_basis": block_basis,
        "mobius_midpoint_directions": midpoint_directions,
        "forced_fixed_word_parity": fixed_word_parity,
        "p31_mobius_cancellation_parity": p31_cancellation_parity,
        "small_exact_replay": replay,
        "small_exact_fixed_word_replay": block_replay,
        "small_exact_mobius_midpoint_replay": midpoint_replay,
        "proved": {
            "block_triangular_form": True,
            "unique_fixed_edge_parity_vector": True,
            "explicit_fixed_edge_inverse": True,
            "divided_target_equivalence": True,
            "exact_hamming_slice": True,
            "forced_fixed_word_parity_identity": True,
            "p31_mobius_cancellation_parity_ladder": (
                p31_cancellation_parity is not None
            ),
            "directionwise_parallel_slices": True,
            "per_orbit_affine_fixed_word": True,
            "fixed_word_blocks_form_binary_basis": True,
            "disjoint_C_kernel_lifts": True,
            "mobius_midpoint_direction_multiplicity_two": True,
            "support_weight_coset_reduction": True,
            "reduced_unused_orbit_fibre_nonempty": False,
            "residual_ii_closed": False,
        },
        "duplicate_work_guard": (
            "The symmetric-lattice theorem gives mod-two surjectivity and "
            "the box separately. This theorem combines them to eliminate "
            "the fixed variables and divide the remaining fixed-cell block; "
            "do not retry first-layer or unrestricted parity/Smith "
            "obstructions. The punctured halved joint code remains open."
        ),
        "next_linear_gate": (
            "image of the punctured halved joint code [B mod 2; C mod 2], "
            "equivalently fixed cells mod 4 together with paired cells mod 2"
        ),
        "proved_all": proved,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
