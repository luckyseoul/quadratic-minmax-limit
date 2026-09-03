#!/usr/bin/env python3
"""Odd-Radon centrality of the compact residual in branch-C hard rows.

For ``p=4r+3>=31``, a balanced branch-C hard row is a fixed unit star
plus ``e`` compact atoms, where ``e<=2r-2``.  The unit star is invisible
to every moment below degree ``p-1``.  This module records the exact
Couvreur/line/conic argument proving that, when all odd global forms
vanish, the compact residual has an integer centrally symmetric edge
chain.

The theorem concerns the compact residual, not the whole hard row.  It
does not address nonzero odd forms, even moments, a common edge lift, or
residual-(ii).
"""

from __future__ import annotations

from functools import lru_cache

from e1_gmin_m4_compact_ray_moment_gate import (
    p3_low_weight_line_peeling_certificate,
)
from e1_gmin_m4_conic_odd_radon import (
    conic_reduction_constants,
    theorem_record as conic_theorem_record,
)
from e1_gmin_m4_equianharmonic_component_packing import (
    equianharmonic_component_packing_certificate,
)
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


def _check_e(p: int, compact_count: int) -> tuple[int, int]:
    r = _check_p(p)
    if (
        not isinstance(compact_count, int)
        or isinstance(compact_count, bool)
        or not 0 <= compact_count <= 2 * r - 2
    ):
        raise ValueError("need 0<=compact_count<=2r-2")
    return r, compact_count


def p3_hard_compact_count_bound(p: int) -> dict[str, object]:
    """Certify ``e<=2r-2`` on the full balanced branch-C hard ray."""
    r = _check_p(p)
    direction_count = 2 * r + 2
    lower = 2 * r * r - 4 * r - 2
    upper = 4 * r * r - 2 * r - 5
    upper_numerator = upper + 1
    quotient, remainder = divmod(upper_numerator, direction_count)
    maximum_compact_count = quotient + bool(remainder)

    upper_profile = p3_local_survivor(p, upper)
    upper_counts = tuple(int(row["e"]) for row in upper_profile["hard_rows"])
    lower_profile = p3_local_survivor(p, lower)
    lower_counts = tuple(int(row["e"]) for row in lower_profile["hard_rows"])

    proved = bool(
        upper_profile["proved_local_aggregate"]
        and lower_profile["proved_local_aggregate"]
        and quotient == 2 * r - 3
        and remainder == 2
        and maximum_compact_count == 2 * r - 2
        and upper_counts.count(2 * r - 2) == 2
        and upper_counts.count(2 * r - 3) == 2 * r
        and max(lower_counts) <= maximum_compact_count
    )
    if not proved:
        raise ArithmeticError("the balanced hard-row compact bound changed")
    return {
        "p": p,
        "r": r,
        "branch_C_t_interval": [lower, upper],
        "balanced_hard_direction_count": direction_count,
        "hard_compact_count_formula": "e_L in balanced(t+1,2r+2)",
        "upper_endpoint_division": {
            "numerator": upper_numerator,
            "denominator": direction_count,
            "quotient": quotient,
            "remainder": remainder,
        },
        "upper_endpoint_compact_count_multiset": {
            str(2 * r - 3): 2 * r,
            str(2 * r - 2): 2,
        },
        "maximum_hard_compact_count": maximum_compact_count,
        "full_ray_bound": "0<=e<=2r-2",
        "proved": proved,
    }


def unit_star_odd_blind_certificate(p: int) -> dict[str, object]:
    """Record the field-power-sum proof that a unit star is odd-blind."""
    r = _check_p(p)
    maximum_degree = p - 2
    # Q_{d,k}(j,t) has t-degree d-k, hence at most d.  Every monomial
    # t^a with 0<=a<=p-2 has zero sum over F_p (including a=0).
    proved = maximum_degree < p - 1
    return {
        "p": p,
        "r": r,
        "moment_degree_interval": [2, maximum_degree],
        "summand_t_degree_upper_bound": "d-k<=d<=p-2",
        "field_power_sums_used": "sum_(t in F_p)t^a=0 for 0<=a<=p-2",
        "all_odd_degrees_through_p_minus_2_vanish": True,
        "indeed_all_degrees_through_p_minus_2_vanish": True,
        "proved": proved,
    }


def _three_to_one_l1_floor(size: int) -> int:
    quotient, remainder = divmod(size, 3)
    return 3 * quotient * (quotient + 1) // 2 + remainder * (quotient + 1)


@lru_cache(maxsize=1)
def _equianharmonic_dependency() -> dict[str, object]:
    """Cache the characteristic-uniform symbolic eight-case certificate."""
    return equianharmonic_component_packing_certificate()


def p3_hard_compact_line_exclusion(
    p: int, compact_count: int
) -> dict[str, object]:
    """Exclude every one- or two-maximal-line compact-only word."""
    r, e = _check_e(p, compact_count)
    h = 2 * r + 1
    m = h - 2
    occurrence_budget = 3 * e
    full_hard_budget = 6 * r - 6
    peeling = p3_low_weight_line_peeling_certificate(p)

    one_horizontal_diagonal_floor = h * (h + 1) // 2
    vertical_canonical_absolute_bound = occurrence_budget // h
    vertical_alternative_lift_floor = (h - 1) + (p - 1)
    double_vertical_demand = 2 * h
    compact_aligned_capacity = 2 * e

    different_family_floor = h * (h - 1) // 2
    same_vertical_injective_floor = (r + 1) ** 2
    same_vertical_one_zero_l1_floor = r * (r + 1)
    two_vertical_canonical_sum_bound = occurrence_budget // h
    two_vertical_alternative_lift_floor = 4 * h - 1
    three_to_one_floor = _three_to_one_l1_floor(h)
    reduced_three_to_one_floor = _three_to_one_l1_floor(h - 1)

    proved = bool(
        peeling["proved"]
        and occurrence_budget <= full_hard_budget
        and full_hard_budget == 3 * m - 3
        and one_horizontal_diagonal_floor > full_hard_budget
        and vertical_canonical_absolute_bound <= 2
        and vertical_alternative_lift_floor > full_hard_budget
        and compact_aligned_capacity < double_vertical_demand
        and different_family_floor > full_hard_budget
        and same_vertical_injective_floor > full_hard_budget
        and same_vertical_one_zero_l1_floor > full_hard_budget
        and two_vertical_canonical_sum_bound <= 2
        and two_vertical_alternative_lift_floor > full_hard_budget
        and three_to_one_floor > full_hard_budget
        and reduced_three_to_one_floor > full_hard_budget
    )
    if not proved:
        raise ArithmeticError("the hard compact line exclusion changed")
    return {
        "p": p,
        "r": r,
        "compact_atom_count": e,
        "signed_edge_occurrence_bound": occurrence_budget,
        "full_hard_signed_occurrence_bound": full_hard_budget,
        "dual_polynomial_total_degree": m,
        "line_peeling_dependency_proved": peeling["proved"],
        "one_horizontal_diagonal_l1_floor": one_horizontal_diagonal_floor,
        "one_horizontal_diagonal_margin": (
            one_horizontal_diagonal_floor - occurrence_budget
        ),
        "vertical_canonical_absolute_bound": vertical_canonical_absolute_bound,
        "vertical_alternative_lift_l1_floor": vertical_alternative_lift_floor,
        "unit_vertical_excluded_by_quotient_Euler_parity": True,
        "double_vertical_aligned_occurrence_demand": double_vertical_demand,
        "compact_aligned_occurrence_capacity": compact_aligned_capacity,
        "double_vertical_aligned_deficit": (
            double_vertical_demand - compact_aligned_capacity
        ),
        "different_family_two_line_l1_floor": different_family_floor,
        "same_vertical_injective_l1_floor": same_vertical_injective_floor,
        "same_vertical_h_minus_1_distinct_l1_floor": (
            same_vertical_one_zero_l1_floor
        ),
        "two_vertical_canonical_absolute_sum_bound": (
            two_vertical_canonical_sum_bound
        ),
        "two_vertical_alternative_lift_l1_floor": (
            two_vertical_alternative_lift_floor
        ),
        "same_vertical_constant_case_excluded_by_quotient_parity": True,
        "horizontal_diagonal_projective_fibre_bound": 3,
        "three_to_one_projective_l1_floor": three_to_one_floor,
        "h_minus_1_three_to_one_projective_l1_floor": (
            reduced_three_to_one_floor
        ),
        "all_one_maximal_line_supports_excluded": True,
        "all_two_maximal_line_supports_excluded": True,
        "reducible_conic_two_line_supports_excluded": True,
        "atom_profile_is_compact_only": True,
        "proved": proved,
    }


def p3_hard_compact_conic_exclusion(
    p: int, compact_count: int
) -> dict[str, object]:
    """Exclude the star, nonequianharmonic, and equianharmonic conics."""
    r, e = _check_e(p, compact_count)
    h = 2 * r + 1
    m = h - 2
    occurrence_budget = 3 * e
    target_size = p - 2
    conic_circuit_size = 2 * m + 2

    # The geometry and dual-weight classification are independent of the
    # opposite-row atom counts used by the original conic certificate.  At
    # b=r that certificate uses its largest allowed budget 3m, so every
    # geometric and l1 conclusion applies a fortiori here.
    geometry = conic_reduction_constants(p, r)
    conic_record = conic_theorem_record()
    equianharmonic_exists = p % 12 == 7
    packing = (
        _equianharmonic_dependency()
        if equianharmonic_exists
        else None
    )

    outside_support_upper_bound = occurrence_budget - conic_circuit_size
    nonconstant_l1_floor = h * (h - 1)
    nonunit_constant_l1_floor = 2 * target_size
    nonequianharmonic_score_bound = 2 * e + 2
    equianharmonic_score_bound = 2 * e
    packing_score_three_excluded = bool(
        packing is None
        or (
            packing["proved"]
            and packing["all_exceptional_characteristics_below_31"]
            and packing["score_three_compact_no_go"]["proved"]
            and not packing["score_three_compact_no_go"][
                "distinct_label_score_three_compact_exists"
            ]
        )
    )

    proved = bool(
        geometry["proved"]
        and geometry["forced_normal_form"]
        == "U=u*z^2, D=d*(z-1)^2 with u,d squares"
        and conic_record["proved"]["nonequianharmonic_constant_branch_is_excluded"]
        and occurrence_budget <= 3 * m - 3
        and outside_support_upper_bound <= m - 5
        and nonconstant_l1_floor > occurrence_budget
        and nonunit_constant_l1_floor > occurrence_budget
        and nonequianharmonic_score_bound < target_size
        and packing_score_three_excluded
        and (not equianharmonic_exists or equianharmonic_score_bound < target_size)
    )
    if not proved:
        raise ArithmeticError("the hard compact conic exclusion changed")
    return {
        "p": p,
        "r": r,
        "compact_atom_count": e,
        "signed_edge_occurrence_bound": occurrence_budget,
        "conic_circuit_size": conic_circuit_size,
        "outside_support_upper_bound": outside_support_upper_bound,
        "outside_nonzero_dual_minimum_support": m,
        "conic_peeling_closes": outside_support_upper_bound < m,
        "forced_normal_form": geometry["forced_normal_form"],
        "nonconstant_orbit_difference_l1_floor": nonconstant_l1_floor,
        "nonunit_constant_l1_floor": nonunit_constant_l1_floor,
        "constant_word_target_size": target_size,
        "star_constant_conic_excluded_by_quotient_Euler_parity": True,
        "nonequianharmonic_compact_score_bound": (
            nonequianharmonic_score_bound
        ),
        "nonequianharmonic_score_margin": (
            target_size - nonequianharmonic_score_bound
        ),
        "nonequianharmonic_dependency": (
            "NOTE_2026-09-03_CONIC_ODD_RADON_DICHOTOMY.md, Section 4"
        ),
        "equianharmonic_branch_exists": equianharmonic_exists,
        "equianharmonic_compact_atom_score_upper_bound": 2,
        "equianharmonic_compact_score_bound": equianharmonic_score_bound,
        "equianharmonic_score_margin": (
            target_size - equianharmonic_score_bound
        ),
        "equianharmonic_score_three_dependency_proved": (
            packing_score_three_excluded
        ),
        "all_irreducible_conic_supports_excluded": True,
        "proved": proved,
    }


def p3_hard_compact_odd_radon_centrality(
    p: int, compact_count: int
) -> dict[str, object]:
    """Prove integer centrality for one zero-odd hard compact residual."""
    r, e = _check_e(p, compact_count)
    h = 2 * r + 1
    m = h - 2
    occurrence_budget = 3 * e
    line = p3_hard_compact_line_exclusion(p, e)
    conic = p3_hard_compact_conic_exclusion(p, e)
    star = unit_star_odd_blind_certificate(p)
    individual_orbit_difference_bound = 2 * e

    proved = bool(
        star["proved"]
        and line["proved"]
        and conic["proved"]
        and occurrence_budget <= 3 * m - 3
        and occurrence_budget < 3 * m
        and individual_orbit_difference_bound < p
    )
    if not proved:
        raise ArithmeticError("the hard compact odd-Radon theorem changed")
    return {
        "p": p,
        "r": r,
        "compact_atom_count": e,
        "hard_row_decomposition": "one fixed unit star plus e compact atoms",
        "unit_star_odd_blind": star["proved"],
        "odd_degree_hypothesis": "all global odd forms d=3,5,...,p-2 vanish",
        "dual_polynomial_total_degree": m,
        "signed_edge_occurrence_bound": occurrence_budget,
        "couvreur_alternatives": [
            "m+2 collinear points",
            "2m+2 points on a conic",
            "3m-point cubic/degree-m complete intersection",
        ],
        "cubic_alternative_support_floor": 3 * m,
        "cubic_alternative_excluded_by_strict_support": True,
        "line_and_two_line_branches_excluded": line["proved"],
        "all_conic_branches_excluded": conic["proved"],
        "mod_p_compact_residual_word_is_zero": True,
        "individual_integer_orbit_difference_bound": (
            individual_orbit_difference_bound
        ),
        "integer_lift_margin": p - individual_orbit_difference_bound,
        "compact_residual_signed_edge_chain_is_centrally_symmetric": True,
        "whole_hard_row_is_centrally_symmetric": False,
        "whole_hard_row_equals": "fixed unit-star chain plus central residual",
        "assumes_balanced_branch_C_hard_template": True,
        "assumes_zero_odd_global_forms": True,
        "nonzero_odd_global_forms_ruled_out": False,
        "joint_degree_six_eight_ruled_out": False,
        "common_Fp_edge_lift_constructed": False,
        "Boolean_lift_constructed": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def p3_balanced_hard_rows_odd_radon_centrality(p: int) -> dict[str, object]:
    """Apply the compact-residual theorem to every balanced branch-C hard row."""
    count_bound = p3_hard_compact_count_bound(p)
    r = int(count_bound["r"])
    maximum = int(count_bound["maximum_hard_compact_count"])
    boundary = p3_hard_compact_odd_radon_centrality(p, maximum)
    proved = bool(
        count_bound["proved"]
        and boundary["proved"]
        and maximum == 2 * r - 2
        and boundary["integer_lift_margin"] > 0
    )
    if not proved:
        raise ArithmeticError("the balanced hard-row centrality theorem changed")
    return {
        "p": p,
        "r": r,
        "branch_C_t_interval": count_bound["branch_C_t_interval"],
        "hard_compact_count_interval": [0, maximum],
        "worst_case_certificate": boundary,
        "all_balanced_hard_compact_residuals_central_when_odd_forms_zero": True,
        "whole_hard_rows_retain_their_fixed_unit_stars": True,
        "assumes_zero_odd_global_forms": True,
        "nonzero_odd_global_forms_ruled_out": False,
        "joint_degree_six_eight_ruled_out": False,
        "unbalanced_hard_allocations_ruled_out": False,
        "common_Fp_edge_lift_constructed": False,
        "Boolean_lift_constructed": False,
        "residual_ii_closed": False,
        "proved": proved,
    }


def theorem_record() -> dict[str, object]:
    samples = {
        str(p): p3_balanced_hard_rows_odd_radon_centrality(p)
        for p in (31, 43, 47)
    }
    return {
        "title": "Branch-C hard-row compact-residual odd-Radon centrality",
        "status": "PROVED STRUCTURAL REDUCTION",
        "proved": {
            "unit_star_is_odd_blind_through_degree_p_minus_2": True,
            "every_balanced_hard_compact_residual_is_central_under_zero_odd_forms": True,
            "whole_hard_row_is_central": False,
            "nonzero_odd_global_forms_ruled_out": False,
            "joint_degree_six_eight_ruled_out": False,
            "common_Fp_edge_lift_constructed": False,
            "Boolean_lift_constructed": False,
            "residual_ii_closed": False,
        },
        "sample_certificates": samples,
        "remaining_gate": (
            "coordinate the fixed hard-row unit-star syndromes and the central "
            "compact residuals with the opposite rows in the even hierarchy, "
            "then solve the common integral/Boolean lift"
        ),
        "L_status": "OPEN",
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), sort_keys=True, indent=2))
