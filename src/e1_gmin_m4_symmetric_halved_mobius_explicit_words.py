#!/usr/bin/env python3
"""Exclude the explicit weight-|Delta| dual words from Mobius supports.

For one nontrivial direction-localized Mobius half, the difference class
determines the Mobius parameter.  The midpoint and difference are proportional
only at the identity parameter.  Consequently a union of ``(p+1)/2`` halves
cannot contain a vertical fibre or a scalar graph in ``Delta x Delta``.

The proof recorded by :func:`explicit_word_intersection_theorem` is symbolic
for every odd prime whenever the displayed Mobius parameterization is used;
the Paley application has ``p=3 mod 4``.  The function
:func:`exact_small_explicit_word_replay` only checks the formulas at ``p=3,7``
and is not theorem evidence.
"""

from __future__ import annotations

import json

from e1_gmin_m4_inversion_antisymmetric_radon import (
    _functional_value,
    projective_functionals,
)
from e1_gmin_m4_mobius_half_symmetric import mobius_parameter_edges
from e1_gmin_m4_prop15721 import is_prime
from e1_gmin_m4_symmetric_halved_mobius_cover import (
    _antipodal_classes,
    _midpoint_difference,
)


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


def explicit_word_intersection_theorem(p: int) -> dict[str, object]:
    """Return the symbolic vertical-fibre and scalar-graph exclusion.

    In coordinates dual to independent functionals ``L,M``, put
    ``z=t+1``.  Before forgetting the sign of the half-difference, the
    parameterized half has

        a_z=(j/2)(z, z-z^-1),
        delta_z=(j/2)(z-2, z-2+z^-1).

    The theorem is uniform in ``L,M`` and every nonzero center ``j``.
    """
    h = _check_odd_prime(p)
    d = p + 1
    half_count = d // 2
    delta_size = d * h
    proved = bool(
        half_count == h + 1
        and delta_size == half_count * (p - 1)
        and half_count < delta_size
    )
    if not proved:
        raise ArithmeticError("the explicit-word intersection counts changed")
    return {
        "p": p,
        "h": h,
        "d": d,
        "delta_size": delta_size,
        "Mobius_half_count": half_count,
        "hypotheses": (
            "L,M are independent functionals, j is nonzero, and "
            "z=t+1 lies in F_p^*"
        ),
        "oriented_parameterization": {
            "midpoint": "a_z=(j/2)*(z, z-z^(-1)) in (L,M)-coordinates",
            "half_difference": (
                "delta_z=(j/2)*(z-2, z-2+z^(-1)) in "
                "(L,M)-coordinates"
            ),
        },
        "vertical_fibre": {
            "support": "Delta times {[delta_0]}",
            "size": delta_size,
            "one_half_intersection_at_most": 1,
            "injectivity_proof": (
                "delta_z'=epsilon*delta_z gives z'=z for epsilon=1; "
                "for epsilon=-1 the L and M-L coordinates give "
                "z'=4-z=-z, hence 4=0, impossible for odd p"
            ),
            "exact_hit_criterion": (
                "for a representative d of [delta_0], exactly one possible "
                "epsilon in {+1,-1} satisfies "
                "4*(M(d)-L(d))*(L(d)+epsilon*j)=j^2; then "
                "z=2+2*epsilon*L(d)/j"
            ),
            "union_intersection_at_most": half_count,
            "union_equality_criterion": (
                "every half hits [delta_0] and the resulting midpoint "
                "classes are pairwise distinct"
            ),
            "cannot_be_contained": True,
        },
        "scalar_graphs": {
            "support": "G_[c]={( [a],[delta] ):[a]=[c*delta]}",
            "size": delta_size,
            "determinant_formula": (
                "det_(L,M)(a_z,delta_z)=j^2*(z-1)/(2*z)"
            ),
            "identity_graph_one_half_intersection": 1,
            "identity_parameter": "z=1 (equivalently t=0), where a_1=-delta_1",
            "nonidentity_graph_one_half_intersection": 0,
            "identity_graph_union_intersection_at_most": half_count,
            "union_equality_criterion": (
                "the z=1 diagonal columns of all halves are pairwise distinct"
            ),
            "all_scalar_graphs_cannot_be_contained": True,
        },
        "actual_support_reason": (
            "post-cancellation support is contained in the raw union, so "
            "the same exclusions hold a fortiori"
        ),
        "directions_or_centers_coherence_used": False,
        "all_odd_prime_symbolic_proof": True,
        "Paley_application_scope": "p=3 mod 4",
        "structured_punctured_map_surjective_proved": False,
        "minimum_row_code_distance_proved": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def _determinant(p: int, first: tuple[int, int], second: tuple[int, int]) -> int:
    return (first[0] * second[1] - first[1] * second[0]) % p


def exact_small_explicit_word_replay(p: int) -> dict[str, object]:
    """Replay the formulas for every localized half at ``p=3,7`` only."""
    _check_odd_prime(p)
    if p not in (3, 7):
        raise ValueError("the exact replay is limited to p=3,7")

    directions = projective_functionals(p)
    classes = _antipodal_classes(p)
    half_count = 0
    all_difference_maps_injective = True
    all_hit_criteria_exact = True
    all_scalar_intersections_exact = True

    for direction in directions:
        for auxiliary in directions:
            if _determinant(p, direction, auxiliary) == 0:
                continue
            for center in range(1, p):
                parameter_edges = mobius_parameter_edges(
                    p, direction, auxiliary, center
                )
                columns = {
                    parameter: _midpoint_difference(p, edge)
                    for parameter, edge in parameter_edges.items()
                }
                differences = {
                    difference for _midpoint, difference in columns.values()
                }
                all_difference_maps_injective &= len(differences) == p - 1

                criterion_classes = set()
                for difference in classes:
                    l_value = _functional_value(p, direction, difference)
                    m_value = _functional_value(p, auxiliary, difference)
                    signs = tuple(
                        epsilon
                        for epsilon in (-1, 1)
                        if (
                            4
                            * (m_value - l_value)
                            * (l_value + epsilon * center)
                            - center * center
                        )
                        % p
                        == 0
                    )
                    if signs:
                        criterion_classes.add(difference)
                    all_hit_criteria_exact &= len(signs) <= 1
                all_hit_criteria_exact &= criterion_classes == differences

                proportional = tuple(
                    (parameter, midpoint, difference)
                    for parameter, (midpoint, difference) in columns.items()
                    if _determinant(p, midpoint, difference) == 0
                )
                all_scalar_intersections_exact &= bool(
                    len(proportional) == 1
                    and proportional[0][0] == 0
                    and proportional[0][1] == proportional[0][2]
                )
                half_count += 1

    expected_half_count = (p + 1) * p * (p - 1)
    proved = bool(
        half_count == expected_half_count
        and all_difference_maps_injective
        and all_hit_criteria_exact
        and all_scalar_intersections_exact
    )
    if not proved:
        raise ArithmeticError("the tiny explicit-word replay changed")
    return {
        "p": p,
        "localized_halves_checked": half_count,
        "expected_localized_halves": expected_half_count,
        "difference_class_map_injective": all_difference_maps_injective,
        "vertical_hit_criterion_exact": all_hit_criteria_exact,
        "one_identity_scalar_intersection_only": all_scalar_intersections_exact,
        "role": "fail-when-wrong p=3,7 replay, not theorem evidence",
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    """Return the symbolic theorem and its two tiny formula replays."""
    theorem = explicit_word_intersection_theorem(p)
    replay = {
        small_p: exact_small_explicit_word_replay(small_p)
        for small_p in (3, 7)
    }
    proved = bool(
        theorem["proved"] and all(record["proved"] for record in replay.values())
    )
    if not proved:
        raise ArithmeticError("the explicit-word theorem record changed")
    return {
        "title": "Mobius exclusion of explicit weight-|Delta| row-code words",
        "status": (
            "VERTICAL FIBRES AND SCALAR GRAPHS EXCLUDED; "
            "GENERAL LOW-WEIGHT SUPPORTS OPEN"
        ),
        "symbolic_theorem": theorem,
        "small_exact_replay": replay,
        "proved": {
            "vertical_fibre_containment_excluded": True,
            "scalar_graph_containment_excluded": True,
            "all_low_weight_row_words_excluded": False,
            "structured_punctured_surjectivity": False,
            "residual_ii_closed": False,
        },
        "proved_all_claimed_statements": proved,
    }


if __name__ == "__main__":
    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
