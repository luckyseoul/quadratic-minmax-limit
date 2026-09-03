#!/usr/bin/env python3
"""Exact exclusion of the all-active ``c=p-2`` pencil equality case.

The proof is symbolic.  It combines the equality pencil theorem from
``e1_gmin_m4_hard_star_antisymmetric_support`` with a two-point completion
and the prime-order Rédei--Megyesi direction theorem.  It performs no
finite configuration census and does not close the remaining Boolean lift.
"""
from __future__ import annotations

from e1_gmin_m4_prop15721 import is_prime


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


def all_active_pencil_equality_exclusion(p: int) -> dict[str, object]:
    """Prove that all active hard stars force ``c>=p-1``.

    The input ``c>=p-2`` and its equality-pencil classification are the
    exact theorems in ``e1_gmin_m4_hard_star_antisymmetric_support``.
    Under equality, orient the common projective endpoint as ``P`` and call
    the other endpoints ``Q``.  Hard-row bijectivity gives

        L(Q) = F_p minus {L(P),-L(P)}.

    The sign ``j_L=+/-L(P)``, Paley column signs, and physical orbit
    orientations do not change this support identity.  Hence adjoining
    ``P,-P`` produces a p-set with every hard direction undetermined.  The
    Rédei--Megyesi bound makes that p-set a line.  Its
    p-2 pencil edges are then parallel in one opposite row, exceeding that
    row's balanced quota.
    """
    r = _check_p(p)
    hard_count = 2 * r + 2
    opposite_count = hard_count
    equality_single_orbits = p - 2
    completed_point_count = p
    maximum_completed_directions = p + 1 - hard_count
    redei_noncollinear_minimum = (p + 3) // 2

    upper_t = 4 * r * r - 2 * r - 5
    upper_opposite_total = 10 * r + 6 + upper_t
    opposite_quota_maximum = (
        upper_opposite_total + opposite_count - 1
    ) // opposite_count

    proved = bool(
        hard_count >= 9
        and maximum_completed_directions == (p + 1) // 2
        and redei_noncollinear_minimum == (p + 3) // 2
        and maximum_completed_directions < redei_noncollinear_minimum
        and opposite_quota_maximum == 2 * r + 2
        and equality_single_orbits > opposite_quota_maximum
    )
    if not proved:
        raise ArithmeticError("the all-active pencil exclusion changed")
    return {
        "p": p,
        "r": r,
        "hypotheses": {
            "zero_odd_global_forms": True,
            "balanced_branch_C": True,
            "all_hard_star_centers_nonzero": True,
            "assumed_single_orbit_equality_c": equality_single_orbits,
        },
        "dependency": {
            "general_single_orbit_floor": "c>=p-2",
            "equality_pencil_threshold_rows": 9,
            "available_active_hard_rows": hard_count,
            "common_projective_endpoint": "[P]",
            "hard_phase_coherence": "j_L^2=L(P)^2",
        },
        "oriented_other_endpoint_set_size": equality_single_orbits,
        "hard_projection_identity": "L(Q)=F_p minus {L(P),-L(P)}",
        "hard_phase_sign_changes_support": False,
        "physical_orbit_orientation_does_not_change_support": True,
        "completion": "S=Q union {P,-P}",
        "completed_point_count": completed_point_count,
        "hard_undetermined_direction_count": hard_count,
        "maximum_completed_direction_count": maximum_completed_directions,
        "redei_megyesi_noncollinear_direction_minimum": (
            redei_noncollinear_minimum
        ),
        "completion_forced_collinear": True,
        "line_is_linear_because_it_contains_P_and_minus_P": True,
        "pencil_edges_parallel_to_line": equality_single_orbits,
        "line_annihilator_is_opposite": True,
        "full_ray_opposite_parallel_quota_maximum": opposite_quota_maximum,
        "parallel_quota_deficit": (
            equality_single_orbits - opposite_quota_maximum
        ),
        "equality_c_eq_p_minus_2_excluded": True,
        "conclusion": "c>=p-1 when every hard center is nonzero",
        "signed_boolean_lift_proved": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    return all_active_pencil_equality_exclusion(p)


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
