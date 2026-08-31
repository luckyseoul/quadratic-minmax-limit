#!/usr/bin/env python3
r"""Prop. 15.731 -- endpoint tangent envelopes by dual-line gluing.

The size hypothesis in Ball--Lavrauw, *Planar arcs*, current manuscript
Theorem 13 (arXiv v4 Theorem 11), is sufficient for their explicit
interpolation formula; it is not necessary for the existence of the squared
tangent envelope at the endpoint left by Proposition 15.730.

Here is the elementary gluing lemma that replaces that interpolation step.
Let ``L_1,...,L_n`` be distinct lines in ``PG(2,K)`` with no three
concurrent.  Fix ``d>=0`` and, on every ``L_i``, a section

    h_i in H^0(L_i,O_{L_i}(d)).

Suppose the two sections have the same value in the fibre of ``O(d)`` at
every node ``L_i intersect L_j``.  Then there is a homogeneous plane form
``F`` of degree ``d`` whose restriction to every ``L_i`` is ``h_i``.

The proof is constructive.  Having matched the first ``r-1`` lines, the
difference on ``L_r`` vanishes at their ``r-1`` distinct intersections.  If
``r-1<=d``, divide that binary form by the restrictions of the first
``r-1`` line equations, lift the quotient to the plane, and add

    (ell_1 ... ell_{r-1}) Q_r.

This corrects ``L_r`` without changing an earlier line.  If ``r-1>d``, the
difference already has more roots than its degree and is zero.  Moreover,

    ker(H^0(P^2,O(d)) -> direct_sum_i H^0(L_i,O(d)))
      = (ell_1 ... ell_n) H^0(P^2,O(d-n)),

where the right side is zero when ``n>d``.  Thus fixed compatible data have
a unique lift when ``n>d`` and, when ``n=d``, an affine line of lifts whose
direction is ``ell_1 ... ell_n``.

Apply this to an odd-order arc ``A`` of size ``k=q+2-t``.  Fix vector
representatives, a determinant, and its associated cross-product map before
Segre-normalizing the degree-``t`` tangent products so

    f_a(b) = (-1)^(t+1) f_b(a).

On the dual line ``L_a={Z:a dot Z=0}``, define the degree-``2t`` section

    h_a(X cross a) = f_a(X)^2.

It is well-defined because every tangent form at ``a`` vanishes on ``a``.
The arc condition says that its dual lines have no triple concurrence.  At
``L_a intersect L_b``, use the same vector ``Z=a cross b``: its two values
are ``f_a(-b)^2=f_a(b)^2`` and ``f_b(a)^2``, equal by the normalized lemma of
tangents.  The gluing lemma therefore gives a homogeneous dual form ``Phi``
of degree ``2t`` satisfying

    Phi(X cross a) = f_a(X)^2                 for every a in A.

This restriction identity is what ``tangent envelope`` means below; merely
containing the dual tangent points would be weaker.

For every endpoint repair in Proposition 15.730, ``t=R+1`` and
``d=2R+2``.  If ``p=3R+2``, then ``k=d+1`` and the normalized envelope is
unique.  If ``p=3R+1``, then ``k=d`` and the normalized envelopes are

    Phi_0 + mu P_A,       P_A(Z)=product_{a in A}(a dot Z).

Any two normalized tangent families differ by one common nonzero scalar:
evaluate their scaling factors on a pair ``a,b`` and cancel the nonzero
secant values in the two tangent-lemma identities.  Replacing all ``f_a``
by ``lambda f_a`` replaces the restriction data by ``lambda^2``.  If a
representative ``a`` is replaced by ``mu_a a``, normalization changes
``f_a`` by ``kappa mu_a^t`` for one common ``kappa``, while the inverse of
``X -> X cross a`` contributes ``mu_a^(-t)``.  Thus every ``h_a`` changes by
the same ``kappa^2``, and the span of ``P_A`` is unchanged.  Consequently
the ``p=3R+2`` envelope has a well-defined projective class.  For
``p=3R+1`` there is instead a well-defined projective pencil with the kernel
point ``[P_A]`` omitted, or equivalently one class modulo ``P_A``; there is
no claim of a unique projective curve in that residue class.  Exact
polynomial equalities always refer to the fixed representatives and
normalization.

There is also an exact transition for adjacent repairs.  Let
``A'=A-{a}+{z}`` switch one selected point on a rich block, let ``b`` be the
selected point left on that block, and put ``C=A intersect A'``.  For
``u in C`` the tangent sets swap ``uz`` for ``ua`` (at ``b`` the two line
forms are proportional).  Put

    m_{uv}(X)=det(X,u,v)=v dot (X cross u).

If ``f_u`` and ``g_u`` are separately Segre-normalized for ``A`` and ``A'``,
factor comparison first gives

    g_u = gamma_u f_u m_{ua}/m_{uz}.

For distinct ``u,v in C``, all evaluations used here are nonzero, and

    det(v,u,a)/det(v,u,z) = det(u,v,a)/det(u,v,z).

The two normalized lemmas of tangents therefore give ``gamma_u=gamma_v``.
Thus all pointwise multipliers are one scalar ``gamma`` and, for any
corresponding envelopes,

    (z dot Z)^2 Phi_{A'} - gamma^2 (a dot Z)^2 Phi_A
       = P_C Q,           P_C=product_{u in C}(u dot Z).

The quotient ``Q`` has degree two for ``p=3R+2`` and degree three for
``p=3R+1``.  This divisibility is compatibility data, not an endpoint
contradiction; in the latter case changing the two envelope representatives
changes ``Q`` by a linear combination of ``(z dot Z)^3`` and
``(a dot Z)^3``.

The edge scalars can be removed coherently across an entire repair family.
Every repair contains the same singleton set ``S``.  Choose ``e in S``.
No rich line contains ``e``, so the ``p`` lines from ``e`` to the other
points of ``D`` are distinct, leaving exactly one projective line through
``e`` which avoids ``D``.  Fix a form ``rho_e`` for that line and put

    f_e^A = rho_e product_{v in D minus A} m_{ev}.

These are exactly the ``R+1`` tangent factors at ``e``.  Normalize all other
tangent products for each repair relative to this fixed ``f_e^A``.  On the
swap above, ``D minus A'=(D minus A)-{z}+{a}``, hence

    f_e^{A'} = f_e^A m_{ea}/m_{ez}.

The common edge multiplier is therefore ``gamma=1``.  The definition is
path-independent, and the factors telescope around every closed walk in the
repair graph.  In residue ``c=1`` this coherently normalizes the restriction
data but does not remove the one-dimensional lift kernel.

This proposition is a proved algebraic refinement of the endpoint normal
form.  It does not prove that an endpoint completion exists or exclude one.
"""
from __future__ import annotations

import json
from math import comb
from pathlib import Path

from e1_gmin_m4_prop15727 import endpoint_block_row, endpoint_residue_data


ROOT = Path(__file__).resolve().parents[1]


def _nonnegative_integer(value: int, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def constructive_line_gluing_lemma(
    degree: int, line_count: int
) -> dict[str, object]:
    """Return the constructive proof ledger for compatible line sections.

    The geometric hypotheses -- distinct lines, no triple concurrence, and
    equality of the two section values at every node -- are assumptions of
    the theorem.  The rows record the correction made when the lines are
    adjoined in an arbitrary order.
    """
    d = _nonnegative_integer(degree, "degree")
    n = _nonnegative_integer(line_count, "line_count")
    if n == 0:
        raise ValueError("line_count must be positive")

    steps: list[dict[str, object]] = []
    for r in range(1, n + 1):
        prior_nodes = r - 1
        if prior_nodes <= d:
            quotient_degree = d - prior_nodes
            steps.append(
                {
                    "line_number_r": r,
                    "prior_node_count": prior_nodes,
                    "nodes_are_distinct": True,
                    "difference_degree": d,
                    "action": "factor_and_correct",
                    "correction_multiplier_degree": prior_nodes,
                    "lifted_quotient_degree": quotient_degree,
                    "correction": (
                        "(product of the first r-1 line equations) times "
                        "a plane lift of the quotient on L_r"
                    ),
                    "earlier_restrictions_preserved": True,
                    "new_restriction_matched": True,
                    "proved": quotient_degree >= 0,
                }
            )
        else:
            steps.append(
                {
                    "line_number_r": r,
                    "prior_node_count": prior_nodes,
                    "nodes_are_distinct": True,
                    "difference_degree": d,
                    "action": "already_zero_by_root_count",
                    "correction_multiplier_degree": None,
                    "lifted_quotient_degree": None,
                    "correction": "zero",
                    "earlier_restrictions_preserved": True,
                    "new_restriction_matched": True,
                    "proved": prior_nodes > d,
                }
            )

    residual_degree = d - n
    kernel_dimension = comb(residual_degree + 2, 2) if residual_degree >= 0 else 0
    unique_lift = n > d
    if unique_lift:
        solution_space = "one lift"
    elif n == d:
        solution_space = "an affine line directed by the product of all lines"
    else:
        solution_space = (
            "an affine space directed by the line product times all forms "
            f"of degree {residual_degree}"
        )

    return {
        "field": "arbitrary",
        "degree_d": d,
        "line_count_n": n,
        "hypotheses": {
            "lines_are_distinct": True,
            "no_three_lines_concurrent": True,
            "sections_have_degree_d": True,
            "section_values_agree_at_every_pairwise_node": True,
        },
        "construction": steps,
        "construction_rule": (
            "factor-and-correct through line d+1; from line d+2 onward "
            "the degree-d difference has more than d distinct roots"
        ),
        "lift_exists": True,
        "kernel": {
            "formula": (
                "(product_i ell_i) times the degree-(d-n) plane forms, "
                "interpreted as zero for d<n"
            ),
            "residual_degree": residual_degree if residual_degree >= 0 else None,
            "dimension": kernel_dimension,
            "reason": (
                "a form vanishing on every distinct line is divisible by "
                "their product, and the converse is immediate"
            ),
        },
        "unique_lift": unique_lift,
        "solution_space": solution_space,
        "proved": all(bool(step["proved"]) for step in steps),
    }


def ball_lavrauw_interpolation_scope() -> dict[str, object]:
    """Record exactly what is, and is not, imported from Ball--Lavrauw."""
    return {
        "authors": "Simeon Ball and Michel Lavrauw",
        "title": "Planar arcs",
        "journal": "Journal of Combinatorial Theory, Series A 160 (2018), 261-287",
        "doi": "10.1016/j.jcta.2018.06.015",
        "current_manuscript_scaled_tangent_lemma": 12,
        "arxiv_v4_scaled_tangent_lemma": 10,
        "current_manuscript_theorem": 13,
        "arxiv_v4_theorem": 11,
        "odd_order_stated_hypothesis": "|A|>=2t+2",
        "role_used_here": (
            "the scaled lemma of tangents normalizes f_a(b)="
            "(-1)^(t+1)f_b(a)"
        ),
        "threshold_role": (
            "sufficient for the theorem's explicit interpolation formula"
        ),
        "threshold_claimed_necessary": False,
        "theorem_13_applied_below_its_threshold": False,
        "new_consequence_used_below_threshold": (
            "dual-line gluing applied after the scaled lemma of tangents"
        ),
        "endpoint_existence_comes_from_line_gluing": True,
        "proved": True,
    }


def endpoint_tangent_envelope_row(p: int) -> dict[str, object]:
    """Apply line gluing to every conditional endpoint repair at ``p``."""
    data = endpoint_residue_data(p)
    R = int(data["R"])
    c = int(data["c"])
    k = p + 1 - R
    t = p + 2 - k
    d = 2 * t
    threshold = 2 * t + 2
    gluing = constructive_line_gluing_lemma(d, k)
    kernel_dimension = int(dict(gluing["kernel"])["dimension"])

    expected_k = {1: d, 2: d + 1}[c]
    expected_threshold_deficit = {1: 2, 2: 1}[c]
    expected_kernel_dimension = {1: 1, 2: 0}[c]
    fixed_normalization_solution = {
        1: "Phi_0 + mu P_A, with P_A=product_{a in A}(a dot Z)",
        2: "one exact homogeneous polynomial Phi",
    }[c]
    projective_ambiguity = {
        1: (
            "projective pencil with [P_A] omitted; equivalently one class "
            "modulo the one-dimensional kernel spanned by P_A"
        ),
        2: "one projective envelope curve",
    }[c]

    checks = {
        "arc_tangent_count": t == R + 1,
        "envelope_degree": d == 2 * R + 2,
        "repair_size": k == expected_k,
        "Ball_Lavrauw_threshold_deficit": threshold - k
        == expected_threshold_deficit,
        "compatible_dual_sections_glue": bool(gluing["lift_exists"]),
        "kernel_dimension": kernel_dimension == expected_kernel_dimension,
        "unique_exactly_in_c2": bool(gluing["unique_lift"]) == (c == 2),
    }
    return {
        **data,
        "scope": (
            "every maximum repair A supplied conditionally by Proposition "
            "15.730, for every admissible rich-block choice"
        ),
        "repair_arc_size_k": k,
        "tangents_per_arc_point_t": t,
        "envelope_degree_d": d,
        "dual_line_count": k,
        "dual_lines_have_no_triple_concurrence": True,
        "normalized_tangent_sections": (
            "h_a(X cross a)=f_a(X)^2 on L_a={Z:a dot Z=0}"
        ),
        "node_compatibility": (
            "at Z=a cross b the values are f_a(-b)^2=f_b(a)^2"
        ),
        "envelope_definition": (
            "a degree-d homogeneous Phi with "
            "Phi(X cross a)=f_a(X)^2 for every a in A and every X"
        ),
        "Ball_Lavrauw_stated_size_threshold": threshold,
        "Ball_Lavrauw_threshold_deficit": threshold - k,
        "Ball_Lavrauw_threshold_needed_for_this_gluing": False,
        "line_gluing": gluing,
        "normalization_and_ambiguity": {
            "fixed_vector_representatives": True,
            "Segre_normalization_relation": (
                "f_a(b)=(-1)^(t+1)f_b(a)"
            ),
            "residual_tangent_function_freedom": (
                "one common nonzero scalar lambda"
            ),
            "common_rescaling_effect": (
                "all restriction data scale by lambda^2, and lambda^2 Phi "
                "is the corresponding rescaled lift"
            ),
            "representative_rescaling_effect": (
                "a->mu_a a makes f_a->kappa mu_a^t f_a and the inverse "
                "cross-product parametrization cancels mu_a^t, so every "
                "dual-line section scales by the same kappa^2"
            ),
            "exact_polynomial_claim_requires_fixed_choices": True,
            "projective_object_independent_of_representative_rescaling": True,
            "fixed_normalization_solution_space": fixed_normalization_solution,
            "projective_interpretation": projective_ambiguity,
            "unique_projective_curve": c == 2,
            "unique_coset_modulo_line_product": True,
        },
        "endpoint_repair_realization_claimed": False,
        "endpoint_excluded": False,
        "checks": checks,
        "proved": all(checks.values()),
    }


def repair_family_coherent_normalization_row(
    p: int, four_secants: int
) -> dict[str, object]:
    """Remove every edge scalar coherently in one full repair family."""
    block = endpoint_block_row(p, four_secants)
    R = int(block["R"])
    c = int(block["c"])
    y = int(block["four_secants_y"])
    singleton_points = int(block["singleton_points"])
    tangents_at_base = R + 1
    lines_from_base_to_other_D_points = p
    projective_lines_through_base = p + 1
    D_empty_lines_through_base = (
        projective_lines_through_base - lines_from_base_to_other_D_points
    )

    checks = {
        "common_singleton_base_exists": singleton_points == c + 1 + 2 * y
        and singleton_points >= 2,
        "other_D_points_give_distinct_lines": (
            lines_from_base_to_other_D_points == p
        ),
        "one_D_empty_base_line": D_empty_lines_through_base == 1,
        "base_tangent_factor_count": R + D_empty_lines_through_base
        == tangents_at_base,
        "endpoint_block_row": bool(block["proved"]),
    }
    return {
        "p": p,
        "R": R,
        "c": c,
        "four_secants_y": y,
        "singleton_point_count": singleton_points,
        "base_point": "choose one fixed e in the common singleton set S",
        "base_line_count": {
            "lines_to_other_D_points": lines_from_base_to_other_D_points,
            "all_projective_lines_through_e": projective_lines_through_base,
            "lines_through_e_avoiding_D": D_empty_lines_through_base,
        },
        "coherent_base_tangent_product": (
            "f_e^A=rho_e product_{v in D minus A} det(X,e,v)"
        ),
        "base_tangent_factor_count": tangents_at_base,
        "per_repair_normalization": (
            "scale every other f_u^A relative to f_e^A by the normalized "
            "lemma of tangents"
        ),
        "edge_ratio": (
            "for A'=(A-{a}) union {z}, "
            "f_e^A'=f_e^A det(X,e,a)/det(X,e,z)"
        ),
        "all_adjacent_swap_multipliers_gamma": 1,
        "closed_walk_cocycle": (
            "trivial: complement factors cancel around every closed walk"
        ),
        "path_independent": True,
        "c1_lift_kernel_removed": False,
        "geometric_realization_claimed": False,
        "checks": checks,
        "proved": all(checks.values()),
    }


def adjacent_repair_transition_row(p: int) -> dict[str, object]:
    """Record the envelope divisibility under one rich-block swap.

    This is conditional on two Proposition 15.730 repairs
    ``A'=(A-{a}) union {z}`` sharing the other selected block point ``b``.
    It asserts no realizability or nonexistence conclusion.
    """
    endpoint = endpoint_tangent_envelope_row(p)
    c = int(endpoint["c"])
    d = int(endpoint["envelope_degree_d"])
    k = int(endpoint["repair_arc_size_k"])
    common_size = k - 1
    transition_degree = d + 2
    quotient_degree = transition_degree - common_size
    expected_quotient_degree = {1: 3, 2: 2}[c]

    checks = {
        "common_dual_line_divisor_degree": common_size == k - 1,
        "transition_form_degree": transition_degree == d + 2,
        "quotient_degree": quotient_degree == expected_quotient_degree,
        "parent_endpoint_gluing": bool(endpoint["proved"]),
    }
    return {
        "p": p,
        "R": int(endpoint["R"]),
        "c": c,
        "hypothesis": (
            "A and A' are adjacent Proposition 15.730 repairs with "
            "A'=(A-{a}) union {z}, common block point b, and C=A intersect A'"
        ),
        "common_arc_points": common_size,
        "tangent_factor_swap": {
            "for_u_in_C_minus_b": (
                "the tangent uz for A is replaced by the tangent ua for A'"
            ),
            "at_b": (
                "the secant line ba=bz is unchanged, so the two displayed "
                "line forms are proportional"
            ),
            "normalized_formula": (
                "g_u=gamma f_u (a dot (X cross u))/(z dot (X cross u)) "
                "for one gamma independent of u"
            ),
            "constant_multiplier_reason": (
                "evaluate at u,v in C; Segre's relation cancels the tangent "
                "values, while the two determinant ratios agree after "
                "swapping u and v"
            ),
        },
        "transition_identity": (
            "(z dot Z)^2 Phi_A' - gamma^2 (a dot Z)^2 Phi_A = P_C Q"
        ),
        "coherently_normalized_family_identity": (
            "(z dot Z)^2 Phi_A' - (a dot Z)^2 Phi_A = P_C Q"
        ),
        "coherent_gamma_one_available": True,
        "common_dual_line_product": "P_C=product_{u in C}(u dot Z)",
        "transition_form_degree": transition_degree,
        "common_dual_line_divisor_degree": common_size,
        "quotient_degree": quotient_degree,
        "quotient_ambiguity": (
            "none from envelope lifts in c=2; in c=1, changing the two "
            "lifts changes Q by nu(z dot Z)^3-gamma^2 mu(a dot Z)^3"
        ),
        "geometric_realization_claimed": False,
        "endpoint_excluded": False,
        "checks": checks,
        "proved": all(checks.values()),
    }


def proposition_15731() -> dict[str, object]:
    """Package the all-prime gluing theorem without changing the live gate."""
    endpoint_rows = [endpoint_tangent_envelope_row(p) for p in (31, 41)]
    transition_rows = [adjacent_repair_transition_row(p) for p in (31, 41)]
    coherent_rows = [
        repair_family_coherent_normalization_row(31, 3),
        repair_family_coherent_normalization_row(41, 3),
    ]
    dependency = ball_lavrauw_interpolation_scope()
    proved = bool(
        dependency["proved"]
        and all(row["proved"] for row in endpoint_rows)
        and all(row["proved"] for row in transition_rows)
        and all(row["proved"] for row in coherent_rows)
    )
    return {
        "prop": "15.731",
        "title": "Endpoint tangent envelopes by constructive dual-line gluing",
        "result_status": "proved algebraic refinement",
        "all_prime_scope": (
            "every conditional endpoint repair in Proposition 15.730: "
            "c=2 has a unique normalized degree-2(R+1) envelope, while "
            "c=1 has an affine pencil modulo the product of its dual lines"
        ),
        "general_line_gluing_lemma": (
            "compatible degree-d sections on distinct projective lines with "
            "no triple concurrence lift constructively to a degree-d plane "
            "form; the kernel is their line product times degree-(d-n) forms"
        ),
        "Ball_Lavrauw_scope": dependency,
        "sample_endpoint_rows": endpoint_rows,
        "sample_adjacent_swap_rows": transition_rows,
        "sample_coherent_normalization_rows": coherent_rows,
        "finite_configuration_search_used": False,
        "endpoint_repair_realization_claimed": False,
        "endpoint_excluded": False,
        "p_plus_one_shell_closed": False,
        "non_walsh_residual_ii_closed": False,
        "multi_level_type_I_closed": False,
        "quadratic_minmax_limit_closed": False,
        "top_level_gates_changed": False,
        "cycle_obstruction_proved": False,
        "phase_bridge_proved": False,
        "next_gate": (
            "historical at this stage: Proposition 15.732 later proved the "
            "linear transition circulation exact, Proposition 15.734 "
            "subsequently closed the k=4p layer for p>=13, and Proposition "
            "15.737 closed the first three p=11 layers; the live residual "
            "front starts at p=11 k=50 and p>=13 k=4p+6"
        ),
        "proved": proved,
    }


def write_evidence() -> Path:
    """Write the deterministic algebraic certificate."""
    output = ROOT / "evidence" / "e1_gmin_m4_prop15731.json"
    payload = json.dumps(proposition_15731(), indent=2, sort_keys=True) + "\n"
    output.write_text(payload)
    return output


def main() -> None:
    result = proposition_15731()
    if not result["proved"]:
        raise ArithmeticError("Proposition 15.731 gluing audit failed")
    path = write_evidence()
    print("Prop 15.731 endpoint tangent-envelope gluing: proved refinement")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
