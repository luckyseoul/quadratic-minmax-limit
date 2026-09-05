"""Affine-parity cube mean floor and a closed infinite local branch.

If an integral nonnegative cube quadratic has affine parity on d variables,
its mean is at least ceil(d/2)/2. Averaging paired cubes through a minimum
of A on J(p,(p+1)/2), with an EVEN parity-support representative B, gives

    |B|*(p+1-|B|) <= 4*(min(A) + p*E[A]).

Consequently masses 2p*E[A] in {2p+4, 2p+6} exclude every even boundary
6 <= |B| <= p-5 for odd p >= 29. This is a local infinite-family theorem,
not a closure of the residual layer or of a global acceptance predicate.
No Boolean tables, graphs, primes, or equality families are enumerated.
"""
from __future__ import annotations

from fractions import Fraction


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def _integer(value: object, name: str) -> int:
    if type(value) is not int:
        raise ValueError(f"{name} must be an integer, not a Boolean")
    return value


def _odd_order(p: int) -> None:
    _integer(p, "p")
    if p < 5 or p % 2 == 0:
        raise ValueError("need odd p>=5; primality is not required locally")


def cube_affine_parity_mean_floor(d: int) -> Fraction:
    """Exact necessary mean floor, including the odd anchor degree."""
    _integer(d, "active parity count")
    if d < 0:
        raise ValueError("active parity count must be nonnegative")
    return Fraction(d + d % 2, 4)


def cube_affine_parity_certificate(d: int) -> dict[str, object]:
    """Record the weighted-degree proof, not a finite cube calculation."""
    floor = cube_affine_parity_mean_floor(d)
    odd_vertices = d + d % 2
    _require(odd_vertices % 2 == 0 and floor == Fraction(odd_vertices, 4),
             "homogenizing anchor parity was lost")
    return {
        "active_parity_variables": d,
        "odd_weighted_degree_vertices": odd_vertices,
        "anchor_has_odd_degree": bool(d % 2),
        "minimum_weighted_quadratic_upper_bound": -(odd_vertices // 2),
        "mean_floor": str(floor),
        "integral_multilinear_coefficients": True,
        "quadratic_multilinear_coefficients_are_even": True,
        "homogenization": "f=c+(1/2)*sum_(0<=i<j) w_ij*z_i*z_j, w_ij integral",
        "local_fields_at_a_global_minimum_are_nonpositive": True,
        "odd_local_fields_at_a_global_minimum_are_at_most_minus_one": True,
        "phase_independent": True,
        "dimension_free": True,
        "proved": True,
    }


def even_boundary_representative(p: int, b: int, phase: int) -> dict[str, int]:
    """Complement an odd support and adjust its phase on the middle slice."""
    _odd_order(p)
    _integer(b, "boundary size")
    _integer(phase, "phase")
    if not 0 <= b <= p or phase not in (0, 1):
        raise ValueError("need 0<=b<=p and phase in {0,1}")
    if b % 2:
        return {"boundary_size": p - b, "phase": (phase + (p + 1) // 2) % 2}
    return {"boundary_size": b, "phase": phase}


def paired_cube_operator_certificate(p: int) -> dict[str, object]:
    """Verify the averaged-cube identity on all monomial-position types."""
    _odd_order(p)
    m = (p + 1) // 2
    slice_single_mean = Fraction(m, p)
    slice_pair_mean = Fraction(m * (m - 1), p * (p - 1))
    single_inside = Fraction(m + 1, 2 * m)
    single_outside = Fraction(1, 2)
    pair_inside = Fraction(m + 2, 4 * m)
    pair_other = Fraction(1, 4)
    checks = [
        single_inside == (1 + p * slice_single_mean) / (p + 1),
        single_outside == p * slice_single_mean / (p + 1),
        pair_inside == (1 + p * slice_pair_mean) / (p + 1),
        pair_other == p * slice_pair_mean / (p + 1),
        Fraction(2, m) * Fraction(1, 2)
        + Fraction(m - 2, m) * Fraction(1, 4) == pair_inside,
        Fraction(1, m) * Fraction(1, 2)
        + Fraction(m - 2, m) * Fraction(1, 4) == pair_other,
    ]
    _require(all(checks), "paired-cube monomial identity failed")
    return {
        "p": p,
        "slice_size": m,
        "cube_dimension": m - 1,
        "identity": "T A(X)=(A(X)+p*E[A])/(p+1)",
        "monomial_position_types_verified": True,
        "includes_constant_linear_and_quadratic_terms": True,
        "proved": True,
    }


def paired_cube_parity_statistics(p: int, b: int, a: int) -> dict[str, object]:
    """Exact expectations at a fixed X with |X intersect B|=a."""
    _odd_order(p)
    _integer(b, "boundary size")
    _integer(a, "intersection size")
    m, q = (p + 1) // 2, (p - 1) // 2
    if b % 2 or not 0 <= b <= p - 1:
        raise ValueError("B must be the even parity-support representative")
    if not max(0, b - q) <= a <= min(b, m):
        raise ValueError("intersection size is infeasible on the middle slice")
    active_mean = Fraction(b * m - a * (2 * b + 1 - 2 * a), m)
    odd_probability = Fraction(a, m)
    floor = (active_mean + odd_probability) / 4
    pointwise_formula = Fraction(b * m - 2 * a * (b - a), 4 * m)
    universal_floor = Fraction(b * (p + 1 - b), 4 * (p + 1))
    _require(floor == pointwise_formula and floor >= universal_floor
             and floor - universal_floor == Fraction((2 * a - b) ** 2, 4 * (p + 1)),
             "even-boundary paired-cube floor failed")
    return {
        "p": p, "boundary_size": b, "intersection_size": a,
        "expected_active_parity_variables": str(active_mean),
        "probability_active_count_is_odd": str(odd_probability),
        "averaged_cube_mean_floor": str(floor),
        "intersection_free_cube_mean_floor": str(universal_floor),
        "exact_intersection_excess": str(floor - universal_floor),
        "unmatched_vertex_odd_degree_term_retained": True,
        "even_boundary_representative_required": True,
        "proved": True,
    }


def slice_affine_parity_minimum_budget(p: int, b: int) -> dict[str, object]:
    """Return the necessary bound b(p+1-b)<=4(min A+p E[A])."""
    _odd_order(p)
    _integer(b, "boundary size")
    if b % 2 or not 0 <= b <= p - 1:
        raise ValueError("B must be the even parity-support representative")
    cube = cube_affine_parity_certificate(3)
    operator = paired_cube_operator_certificate(p)
    balanced = paired_cube_parity_statistics(p, b, b // 2)
    _require(cube["proved"] and cube["anchor_has_odd_degree"]
             and cube["mean_floor"] == "1"
             and operator["proved"] and operator["monomial_position_types_verified"]
             and balanced["proved"]
             and balanced["unmatched_vertex_odd_degree_term_retained"],
             "cube floor or averaged-operator dependency failed")
    numerator = b * (p + 1 - b)
    return {
        "p": p, "boundary_size": b,
        "minimum_plus_p_mean_floor": str(Fraction(numerator, 4)),
        "inequality": "b*(p+1-b) <= 4*(min(A)+p*E[A])",
        "even_boundary_representative_required": True,
        "all_affine_parity_phases": True,
        "cube_floor_dependency": cube,
        "paired_operator_dependency": operator,
        "proved": True,
    }


def middle_boundary_mass_exclusion(p: int, excess: int) -> dict[str, object]:
    """Exclude 6<=b<=p-5 at mass 2p+excess when the strict margin is positive."""
    _odd_order(p)
    _integer(excess, "mass excess")
    if excess not in (4, 6):
        raise ValueError("this application treats exactly mass excess four or six")
    if p < 11:
        raise ValueError("need p>=11 so the middle-boundary interval is nonempty")
    budget = slice_affine_parity_minimum_budget(p, 6)
    mean = Fraction(2 * p + excess, 2 * p)
    minimum_upper_bound = mean.numerator // mean.denominator
    lower_numerator = 6 * (p - 5)
    upper_numerator = 4 * minimum_upper_bound + 2 * (2 * p + excess)
    margin = lower_numerator - upper_numerator
    _require(budget["proved"] and budget["even_boundary_representative_required"]
             and minimum_upper_bound == 1
             and margin == 2 * (p - 17 - excess),
             "local middle-boundary mass accounting failed")
    excluded = margin > 0
    return {
        "p": p, "scaled_mass": 2 * p + excess, "mass_excess": excess,
        "mean": str(mean), "minimum_A_at_most": minimum_upper_bound,
        "even_boundary_interval": [6, p - 5],
        "minimum_boundary_numerator": lower_numerator,
        "maximum_allowed_numerator": upper_numerator,
        "strict_margin": margin,
        "exact_strict_threshold": f"p>{17 + excess}",
        "all_affine_parity_phases": True,
        "all_middle_boundaries_excluded": excluded,
        "proved": excluded,
        "no_claim_at_zero_or_negative_margin": not excluded,
        "local_branch_only": True,
        "entire_residual_layer_closed": False,
        "residual_ii_closed_general": False,
        "limit_closed": False,
        "new_census_used": False,
    }


def first_uncovered_middle_boundary_closure(p: int) -> dict[str, object]:
    """Both low masses have no middle-boundary cell for every odd p>=29."""
    _odd_order(p)
    if p < 29:
        raise ValueError("the joint infinite-family application starts at odd p>=29")
    rows = [middle_boundary_mass_exclusion(p, excess) for excess in (4, 6)]
    _require(all(row["proved"] and row["all_middle_boundaries_excluded"]
                 and row["all_affine_parity_phases"] for row in rows),
             "one of the two low-mass branch closures failed")
    return {
        "p": p,
        "scaled_masses": [2 * p + 4, 2 * p + 6],
        "excluded_even_boundary_interval": [6, p - 5],
        "remaining_even_boundary_candidates": [0, 2, 4, p - 3, p - 1],
        "mass_records": rows,
        "physical_parallel_count_assumed": None,
        "signed_total_assumed": None,
        "both_prime_congruence_classes": True,
        "local_branch_closed": True,
        "entire_residual_layer_closed": False,
        "global_closure_claimed": False,
        "proved": True,
    }
