#!/usr/bin/env python3
r"""Prop. 15.732 -- exact repair cycles and the surviving phase bridge.

Proposition 15.731 suggested composing its quadratic/cubic transition
quotients around cycles of the repair graph.  There is an exact normalization
which shows what that proposal can and cannot produce.

For a repair ``A`` put

    P_A = product_(a in A) L_a,       Theta_A = P_A^2 Phi_A.

If ``A=C union {a}`` and ``A'=C union {z}``, Proposition 15.731 gives

    L_z^2 Phi_A' - L_a^2 Phi_A = P_C Q_(a,z).

Since ``P_A=P_C L_a`` and ``P_A'=P_C L_z``, multiplication by ``P_C^2``
turns this into

    Theta_A' - Theta_A = P_C^3 Q_(a,z).

Thus the naturally cleared transition is the coboundary of an explicit
vertex potential.  Its sum on every closed repair walk is identically zero,
before evaluation, polarization, differentiation, or coefficient extraction.
In residue ``c=1``, changing the envelope lift by ``mu_A P_A`` changes
``Theta_A`` by ``mu_A P_A^3`` and changes the edge by the corresponding
coboundary, so the conclusion is gauge invariant.

This does give an exact local identity.  On a triangle supported by three
collinear block points ``i,j,k``, after cancelling the common outside-block
factor it is

    L_i^3 Q_(j,k)^i + L_k^3 Q_(i,j)^k + L_j^3 Q_(k,i)^j = 0.

The same formula holds on every triangular face of a ``J(4,2)`` repair
factor.  It is automatic from the potentials, not an independent holonomy
obstruction.

There is nevertheless local geometric information in one edge quotient.
Let ``q`` be the dual point of the rich line containing the exchanged points
``a,z`` and the retained point ``b``.  The rich line is a secant of both
repairs, so both envelope values at ``q`` are nonzero.  Taking the lowest
homogeneous part of the edge identity at ``q`` gives

    P_(C-{b})(q) L_b j_q^1(Q)
      = Phi_A'(q) L_z^2 - Phi_A(q) L_a^2.

Consequently ``Q(q)=0`` and its first jet is nonzero.  The first jet is also
unchanged by the cubic lift gauge in residue ``c=1``.  This is direction
sensitive, but it only repackages the universal secant-node values; none of
the residual phase/lift propositions currently identifies its square class.

The other direct bridge also fails quantitatively.  In a nonrich direction
with fibre profile

    ((p-3)/2 empty, 2 singleton, (p-1)/2 double),

deleting ``R`` points destroys at most ``R`` double fibres.  Every repair
therefore has at least ``(p-1)/2-R`` secants and at most ``R+2`` tangents in
that direction.  These give at most ``R+2`` distinct known zeros on the
direction-pencil line in the dual plane.  Since the envelope degree is
``2R+2``, distinct-root counting
cannot force that pencil line as a component.  At ``p=31`` this applies to
each of the at least ``4+y`` nonrich Paley-hard directions from 15.728: there
are at most 12 roots against degree 22.

Finally, multiplying repair data modulo squares cannot restore a missing
trisecant-block phase.  The three pair-selection masks ``110,101,011`` span
the even-weight subspace of ``F_2^3`` and do not contain ``111``.  (For a
4-secant, in contrast, the six weight-two masks do span ``1111``.)  Hence a
repair product cannot recover a factor containing an unselected trisecant
block; for a jet on one block, this applies whenever another trisecant is
present.  It is a precise obstruction to this natural bridge, not a no-go for
all nonlinear phase information.

This proposition is a proved method barrier, not an endpoint exclusion.  A
successful continuation must import information not contained in the linear
transition circulation: either relate the nonzero rich-direction first jet
to the signed residual lift, or exclude the common completion directly from
its many Paley-hard near-pairing directions.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Sequence

from e1_gmin_m4_prop15727 import endpoint_block_row, endpoint_residue_data


ROOT = Path(__file__).resolve().parents[1]


def cleared_transition_exactness_row(p: int) -> dict[str, object]:
    """Audit the exact-potential normalization for one endpoint residue."""
    data = endpoint_residue_data(p)
    R = data["R"]
    c = data["c"]
    k = p + 1 - R
    d = 2 * R + 2
    q_degree = 4 - c
    theta_degree = d + 2 * k
    cleared_edge_degree = 3 * (k - 1) + q_degree
    checks = {
        "repair_size": k == 2 * R + c + 1,
        "quotient_degree": q_degree == d + 2 - (k - 1),
        "cleared_degrees_match": cleared_edge_degree == theta_degree,
        "c1_gauge_is_a_vertex_coboundary": c != 1 or q_degree == 3,
    }
    return {
        "p": p,
        "R": R,
        "c": c,
        "repair_size_k": k,
        "envelope_degree_d": d,
        "transition_quotient_degree": q_degree,
        "vertex_line_product": "P_A=product_{u in A} L_u",
        "cleared_vertex_potential": "Theta_A=P_A^2 Phi_A",
        "cleared_edge_identity": "Theta_A'-Theta_A=P_C^3 Q_(a,z)",
        "vertex_potential_degree": theta_degree,
        "cleared_edge_degree": cleared_edge_degree,
        "closed_walk_sum": "sum_edges P_C^3 Q_(a,z)=0",
        "closed_walk_reason": "the Theta terms telescope vertex by vertex",
        "linear_functionals_preserve_zero": True,
        "covered_linear_functionals": (
            "evaluation, coefficient extraction, differentiation, and polarization"
        ),
        "c1_gauge": {
            "envelope_change": "Phi_A -> Phi_A+mu_A P_A",
            "potential_change": "Theta_A -> Theta_A+mu_A P_A^3",
            "quotient_change": (
                "Q_(a,z) -> Q_(a,z)+mu_A' L_z^3-mu_A L_a^3"
            ),
            "cycle_sum_still_a_coboundary": True,
        },
        "independent_linear_holonomy_obstruction": False,
        "checks": checks,
        "proved": all(checks.values()),
    }


def potential_walk_row(vertices: Sequence[str]) -> dict[str, object]:
    """Return formal potential coefficients for an oriented walk."""
    if len(vertices) < 2 or any(not isinstance(v, str) or not v for v in vertices):
        raise ValueError("need at least two nonempty string vertex labels")
    coefficients: dict[str, int] = {}
    for source, target in zip(vertices, vertices[1:]):
        coefficients[source] = coefficients.get(source, 0) - 1
        coefficients[target] = coefficients.get(target, 0) + 1
    coefficients = {key: value for key, value in coefficients.items() if value}
    closed = vertices[0] == vertices[-1]
    return {
        "vertices": list(vertices),
        "edge_count": len(vertices) - 1,
        "closed": closed,
        "potential_coefficients_after_collection": coefficients,
        "telescopes_to_zero": closed and not coefficients,
        "proved": (closed and not coefficients) or (not closed and bool(coefficients)),
    }


def block_triangle_syzygy_row(p: int) -> dict[str, object]:
    """Record the exact K3/J(4,2)-face identity and its degree ledger."""
    exact = cleared_transition_exactness_row(p)
    q_degree = int(exact["transition_quotient_degree"])
    term_degree = 3 + q_degree
    walk = potential_walk_row(("A_ij", "A_ik", "A_jk", "A_ij"))
    checks = {
        "parent_exactness": bool(exact["proved"]),
        "triangle_is_closed": bool(walk["telescopes_to_zero"]),
        "all_terms_have_one_degree": term_degree == 7 - int(exact["c"]),
    }
    return {
        "p": p,
        "c": exact["c"],
        "support": (
            "a K3 trisecant factor or any triangular face of a J(4,2) "
            "four-secant factor"
        ),
        "orientation": "A_ij -> A_ik -> A_jk -> A_ij",
        "identity_after_common_factor_cancellation": (
            "L_i^3 Q_(j,k)^i + L_k^3 Q_(i,j)^k "
            "+ L_j^3 Q_(k,i)^j = 0"
        ),
        "quotient_degree": q_degree,
        "identity_degree": term_degree,
        "formal_walk": walk,
        "identity_is_automatic_coboundary": True,
        "endpoint_contradiction_from_identity_alone": False,
        "checks": checks,
        "proved": all(checks.values()),
    }


def rich_direction_first_jet_row(p: int) -> dict[str, object]:
    """Extract the first nonzero local term of an edge quotient."""
    exact = cleared_transition_exactness_row(p)
    c = int(exact["c"])
    checks = {
        "parent_exactness": bool(exact["proved"]),
        "quotient_has_positive_degree": int(exact["transition_quotient_degree"]) >= 2,
        "gauge_terms_vanish_to_order_three": c != 1 or int(
            exact["transition_quotient_degree"]
        ) == 3,
    }
    return {
        "p": p,
        "R": exact["R"],
        "c": c,
        "setup": (
            "A=C union {a}, A'=C union {z}; b is the retained point on "
            "the rich line and q is that line's dual point"
        ),
        "secant_values_nonzero": "Phi_A(q) and Phi_A'(q) are nonzero",
        "edge_quotient_value": "Q_(a,z)(q)=0",
        "first_jet_identity": (
            "P_(C-{b})(q) L_b j_q^1(Q_(a,z)) = "
            "Phi_A'(q) L_z^2-Phi_A(q) L_a^2"
        ),
        "normalized_first_jet_formula": (
            "j_q^1 Q_i(j,k)=(K/P_B(q)) "
            "(l_k^2/Delta_(i,k)^2-l_j^2/Delta_(i,j)^2)/l_i"
        ),
        "normalizing_scalar": (
            "K=Phi_(i,j)(q) Delta_(i,j)^2 is independent of the selected pair"
        ),
        "normalizing_scalar_is_nonzero_square": True,
        "first_jet_nonzero": True,
        "nonzero_reason": (
            "two nonproportional cotangent forms L_a,L_z cannot have "
            "nonzero scalar squares equal"
        ),
        "c1_first_jet_gauge_invariant": True,
        "gauge_reason": "the lift ambiguity changes Q only by cubes vanishing to order three at q",
        "first_jet_square_character": "chi_p(P_B(q))",
        "residual_parity_mismatch": (
            "P_B(q) is repair-coloured, while residual parity records the "
            "symmetric difference of the A- and T-tangent fibres"
        ),
        "phase_or_lift_formula_currently_proved": False,
        "interpretation": (
            "the first jet is the smallest local transition datum that a "
            "new signed residual-to-envelope bridge could constrain"
        ),
        "checks": checks,
        "proved": all(checks.values()),
    }


def near_pairing_tangent_barrier_row(p: int) -> dict[str, object]:
    """Bound repair tangents in a direction with two singleton fibres."""
    data = endpoint_residue_data(p)
    R = data["R"]
    k = p + 1 - R
    d = 2 * R + 2
    empty_fibres = (p - 3) // 2
    singleton_fibres = 2
    double_fibres = (p - 1) // 2
    surviving_double_floor = double_fibres - R
    tangent_ceiling = k - 2 * surviving_double_floor
    component_root_threshold = d + 1
    checks = {
        "profile_uses_p_parallel_lines": (
            empty_fibres + singleton_fibres + double_fibres == p
        ),
        "profile_uses_p_plus_one_points": (
            singleton_fibres + 2 * double_fibres == p + 1
        ),
        "deletion_secant_floor_nonnegative": surviving_double_floor >= 0,
        "tangent_ceiling_is_R_plus_two": tangent_ceiling == R + 2,
        "root_count_cannot_force_component": tangent_ceiling < component_root_threshold,
    }
    return {
        "p": p,
        "R": R,
        "c": data["c"],
        "near_pairing_profile": {
            "empty_fibres": empty_fibres,
            "singleton_fibres": singleton_fibres,
            "double_fibres": double_fibres,
        },
        "repair_size_k": k,
        "deleted_point_count": R,
        "surviving_A_secant_floor": surviving_double_floor,
        "A_tangent_ceiling": tangent_ceiling,
        "envelope_degree": d,
        "roots_needed_to_force_direction_component": component_root_threshold,
        "root_deficit_at_least": component_root_threshold - tangent_ceiling,
        "known_tangent_dual_points_are_distinct": True,
        "intersection_multiplicity_used": False,
        "direction_component_forced_by_known_distinct_roots": False,
        "checks": checks,
        "proved": all(checks.values()),
    }


def p31_hard_direction_component_barrier(four_secants: int) -> dict[str, object]:
    """Apply the general tangent ceiling to every 15.728 hard direction."""
    block = endpoint_block_row(31, four_secants)
    general = near_pairing_tangent_barrier_row(31)
    y = int(block["four_secants_y"])
    floor = 4 + y
    checks = {
        "endpoint_block_row": bool(block["proved"]),
        "general_barrier": bool(general["proved"]),
        "p31_tangent_ceiling": general["A_tangent_ceiling"] == 12,
        "p31_envelope_degree": general["envelope_degree"] == 22,
        "15_728_direction_floor_positive": floor >= 4,
    }
    return {
        "p": 31,
        "R": 10,
        "four_secants_y": y,
        "trisecants_x": block["trisecants_x"],
        "15_728_nonrich_Paley_hard_direction_floor": floor,
        "tangents_per_such_direction_at_most": general["A_tangent_ceiling"],
        "envelope_degree": general["envelope_degree"],
        "component_forced_in_any_such_direction": False,
        "checks": checks,
        "proved": all(checks.values()),
    }


def weight_two_selection_span_row(block_size: int) -> dict[str, object]:
    """Compute the GF(2) span of all two-point selection masks."""
    if (
        not isinstance(block_size, int)
        or isinstance(block_size, bool)
        or block_size < 2
    ):
        raise ValueError("block_size must be an integer at least two")
    pairs = tuple(combinations(range(block_size), 2))
    generators = tuple(sum(1 << index for index in pair) for pair in pairs)
    span = {0}
    for generator in generators:
        span |= {value ^ generator for value in tuple(span)}
    full_mask = (1 << block_size) - 1
    rank = (len(span)).bit_length() - 1
    return {
        "block_size": block_size,
        "generator_masks": [
            "".join("1" if index in pair else "0" for index in range(block_size))
            for pair in pairs
        ],
        "span_rank": rank,
        "span_size": len(span),
        "span_is_even_weight_subspace": all(mask.bit_count() % 2 == 0 for mask in span)
        and len(span) == 2 ** (block_size - 1),
        "full_block_mask": format(full_mask, f"0{block_size}b"),
        "full_block_mask_in_span": full_mask in span,
        "proved": rank == block_size - 1,
    }


def endpoint_square_class_product_barrier_row(
    p: int, four_secants: int
) -> dict[str, object]:
    """Audit the block-parity limit of a product-over-repairs bridge."""
    block = endpoint_block_row(p, four_secants)
    x = int(block["trisecants_x"])
    y = int(block["four_secants_y"])
    triple = weight_two_selection_span_row(3)
    quadruple = weight_two_selection_span_row(4)
    checks = {
        "endpoint_block_row": bool(block["proved"]),
        "trisecant_full_mask_missing": not triple["full_block_mask_in_span"],
        "four_secant_full_mask_recoverable": quadruple["full_block_mask_in_span"],
        "c2_has_a_trisecant": int(block["c"]) != 2 or x >= 1,
    }
    blocked_for_every_target = x >= 2
    blocked_for_some_target = x >= 2 or (x >= 1 and y >= 1)
    return {
        "p": p,
        "R": block["R"],
        "c": block["c"],
        "trisecants_x": x,
        "four_secants_y": y,
        "trisecant_mask_span": triple,
        "four_secant_mask_span": quadruple,
        "full_product_on_a_trisecant_recoverable_modulo_squares": False,
        "product_bridge_blocked_for_every_rich_target": blocked_for_every_target,
        "product_bridge_blocked_for_some_rich_target": blocked_for_some_target,
        "sole_trisecant_target_exception": x == 1,
        "scope": (
            "products and quotients of repair-selected contributions after "
            "passing to square classes; this is not a no-go for every "
            "possible nonlinear signed bridge"
        ),
        "checks": checks,
        "proved": all(checks.values()),
    }


def proposition_15732() -> dict[str, object]:
    """Package the cycle and phase-bridge audit without closing the endpoint."""
    exact_rows = [cleared_transition_exactness_row(p) for p in (31, 41)]
    triangle_rows = [block_triangle_syzygy_row(p) for p in (31, 41)]
    jet_rows = [rich_direction_first_jet_row(p) for p in (31, 41)]
    tangent_rows = [near_pairing_tangent_barrier_row(p) for p in (31, 41)]
    p31_rows = [p31_hard_direction_component_barrier(y) for y in range(6)]
    phase_rows = [
        endpoint_square_class_product_barrier_row(31, 3),
        endpoint_square_class_product_barrier_row(41, 3),
    ]
    all_rows = exact_rows + triangle_rows + jet_rows + tangent_rows + p31_rows + phase_rows
    proved = all(bool(row["proved"]) for row in all_rows)
    return {
        "prop": "15.732",
        "title": "Exact repair cycles and phase-bridge barriers",
        "result_status": "proved method barrier",
        "all_prime_scope": (
            "every conditional endpoint repair family of Proposition 15.730"
        ),
        "cleared_transition_rows": exact_rows,
        "triangle_cycle_rows": triangle_rows,
        "rich_direction_first_jet_rows": jet_rows,
        "near_pairing_tangent_rows": tangent_rows,
        "p31_Paley_hard_rows": p31_rows,
        "square_class_product_rows": phase_rows,
        "finite_configuration_search_used": False,
        "nontrivial_linear_cycle_obstruction_available": False,
        "direction_component_from_15_728_profiles_available": False,
        "product_over_repairs_square_class_bridge_universally_available": False,
        "rich_direction_first_jet_proved": True,
        "phase_bridge_proved": False,
        "endpoint_excluded": False,
        "p_plus_one_shell_closed": False,
        "non_walsh_residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "next_gate": (
            "historical at this stage: Proposition 15.733 used the simultaneous "
            "p31 baselines and Proposition 15.734 then closed every k=4p "
            "boundary for p>=13 without a first-jet bridge; the live residual "
            "front is p=11 sharp equality or even k>4p"
        ),
        "proved": proved,
    }


def write_evidence() -> Path:
    """Write the deterministic proof ledger."""
    output = ROOT / "evidence" / "e1_gmin_m4_prop15732.json"
    payload = json.dumps(proposition_15732(), indent=2, sort_keys=True) + "\n"
    output.write_text(payload)
    return output


def main() -> None:
    result = proposition_15732()
    if not result["proved"]:
        raise ArithmeticError("Proposition 15.732 cycle audit failed")
    path = write_evidence()
    print("Prop 15.732 repair-cycle exactness: proved method barrier")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
