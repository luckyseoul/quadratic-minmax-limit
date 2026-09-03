#!/usr/bin/env python3
"""All-prime barrier to a standalone global Mobius-incidence closure.

For the ``m=(p+1)/2`` distinct hard target directions, choose every
auxiliary functional independently and uniformly among the ``p(p-1)``
functionals independent of its target.  Exact counting plus bounded
differences proves that some choice has pairwise-disjoint physical Mobius
supports while its nonzero fixed-word blocks have collision surplus at
least ``2m``.

Consequently distinct target directions, conic/Bezout incidence, arbitrary
triple block intersections, and physical ternarity alone cannot give an
upper bound contradicting the endpoint demand ``sigma >= kappa_z+m+q``.
Any successful use of that demand must additionally couple the block
incidences to the *required* physical cancellation count, or to the target
atom/even-moment equations.  This is a route barrier, not residual-(ii)
closure and not a construction of the endpoint target.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb

from e1_gmin_m4_prop15721 import is_prime


def _check_branch_prime(p: int) -> tuple[int, int]:
    if (
        not isinstance(p, int)
        or isinstance(p, bool)
        or p < 31
        or p % 4 != 3
        or not is_prime(p)
    ):
        raise ValueError("need a prime p=4r+3 with p>=31")
    return (p - 1) // 2, (p + 1) // 2


def _odd_binomial_probability(trials: int, probability: Fraction) -> Fraction:
    """Return ``P(Bin(trials, probability) is odd)`` exactly."""
    return (1 - (1 - 2 * probability) ** trials) / 2


def _expected_floor_half(trials: int, probability: Fraction) -> Fraction:
    """Return ``E floor(Bin(trials, probability)/2)`` exactly."""
    return (
        trials * probability
        - _odd_binomial_probability(trials, probability)
    ) / 2


def global_mobius_incidence_barrier(p: int) -> dict[str, object]:
    """Prove existence of a ternary all-target family with ``sigma>=2m``.

    The target directions may be any ``m`` distinct projective directions
    and the prescribed centers may be any nonzero values.  The branch-C hard
    directions are the intended specialization.

    A block is a nonzero functional modulo sign.  For one fixed target and
    one block not parallel to it, solving the normal form for the auxiliary
    gives ``p-2`` auxiliaries for ``N`` and ``p-2`` disjoint auxiliaries for
    ``-N``.  Thus its hit probability is

        rho = 2(p-2)/(p(p-1)).

    Let ``d_b`` be the number of selected half-conics through block ``b``
    and ``S=sum_b floor(d_b/2)``.  The exact expectation below separates
    blocks in the ``m`` target directions (one ineligible target) from the
    other ``m`` directions.

    The exact two-half physical-overlap count is obtained from the four
    endpoint matchings in normalized ``(q,r,A,B)`` coordinates.  Per sign,
    the direct zero branch contributes ``(p-1)^2`` choices, while the direct
    nonzero and swapped branches contribute ``(p-2)(p-3)`` each.
    Markov and McDiarmid then have a strictly positive common good event.
    """
    h, m = _check_branch_prime(p)
    auxiliary_count = p * (p - 1)
    hit_count = 2 * (p - 2)
    rho = Fraction(hit_count, auxiliary_count)
    a = 1 - 2 * rho

    floor_m = _expected_floor_half(m, rho)
    floor_m_minus_one = _expected_floor_half(m - 1, rho)
    expected_collision_surplus = h * m * (floor_m + floor_m_minus_one)

    # For each of the two physical signs: direct-zero, direct-nonzero,
    # swapped.  The candidate parametrization is bijective in each target's
    # p(p-1) possible auxiliaries.
    direct_zero = (p - 1) ** 2
    direct_nonzero = (p - 2) * (p - 3)
    swapped = direct_nonzero
    physical_overlap_pair_numerator = 2 * (
        direct_zero + direct_nonzero + swapped
    )
    physical_overlap_pair_probability = Fraction(
        physical_overlap_pair_numerator,
        auxiliary_count**2,
    )
    expected_total_pair_overlap = (
        comb(m, 2) * physical_overlap_pair_probability
    )

    # Replacing one half moves p-2 distinct block incidences.  Removed and
    # added blocks have equal cardinality, so the change in S is at most
    # p-2, rather than the cruder 2(p-2).
    bounded_difference = p - 2
    bad_integer_threshold = 2 * m - 1
    deviation = expected_collision_surplus - bad_integer_threshold
    mcdiarmid_exponent = Fraction(
        2 * deviation * deviation,
        m * bounded_difference**2,
    )

    # Uniform elementary bounds used to avoid a numerical probability
    # assertion.  Here a=1-2*rho.  For m>=16, a>=7/8 and
    # a^(m-1)>1/10.  The latter follows from
    # a>1-2/(m-1), monotonicity of (1-2/n)^n, and the n=15 value.
    a_at_least_seven_eighths = a >= Fraction(7, 8)
    a_power_above_one_tenth = a ** (m - 1) > Fraction(1, 10)
    expectation_lower_bound = Fraction(m * (35 * m - 67), 64)
    expectation_bound_holds = (
        expected_collision_surplus > expectation_lower_bound
    )
    overlap_markov_below_three_quarters = (
        expected_total_pair_overlap < Fraction(3, 4)
    )

    if m == 16:
        # At the least branch order, a^15>(7/8)^15>1/8 gives the sharper
        # exponent 4/3.  The exponential comparison is rational:
        # exp(4/3)>1+4/3+(4/3)^2/2+(4/3)^3/6=293/81, and
        # 81/293 < 1-E[T].
        boundary_power = a**15 > Fraction(1, 8)
        boundary_exponent = mcdiarmid_exponent > Fraction(4, 3)
        exponential_tail_rational_upper = Fraction(81, 293)
        probability_union_strictly_below_one = bool(
            boundary_power
            and boundary_exponent
            and exponential_tail_rational_upper
            < 1 - expected_total_pair_overlap
        )
        probability_proof = (
            "m=16: exponent>4/3, exp(-4/3)<81/293, and "
            "E[T]=10096/14415"
        )
    else:
        # Put n=m-18.  After the expectation lower bound is substituted,
        # exponent>3/2 is equivalent to positivity of
        # 1225*n^4+62262*n^3+1060117*n^2+6169740*n+2097892.
        n = m - 18
        exponent_polynomial = (
            1225 * n**4
            + 62262 * n**3
            + 1060117 * n**2
            + 6169740 * n
            + 2097892
        )
        boundary_power = True
        boundary_exponent = bool(
            m >= 18 and exponent_polynomial > 0
            and mcdiarmid_exponent > Fraction(3, 2)
        )
        exponential_tail_rational_upper = Fraction(1, 4)
        probability_union_strictly_below_one = bool(
            boundary_exponent
            and overlap_markov_below_three_quarters
        )
        probability_proof = (
            "m>=18: exponent>3/2, exp(-3/2)<1/4, and E[T]<3/4"
        )

    proved = bool(
        h == m - 1
        and auxiliary_count == p * (p - 1)
        and hit_count == 2 * (p - 2)
        and physical_overlap_pair_numerator == 6 * p * p - 24 * p + 26
        and a_at_least_seven_eighths
        and a_power_above_one_tenth
        and expectation_bound_holds
        and expected_collision_surplus > 2 * m
        and expected_total_pair_overlap < 1
        and overlap_markov_below_three_quarters
        and deviation > 0
        and probability_union_strictly_below_one
    )
    if not proved:
        raise ArithmeticError("the global Mobius-incidence barrier changed")

    return {
        "p": p,
        "h": h,
        "m": m,
        "target_quantifier": (
            "any m distinct projective target directions and any nonzero centers"
        ),
        "auxiliaries_per_target": auxiliary_count,
        "nonzero_blocks_per_half": p - 2,
        "signed_block_classes_total": (p * p - 1) // 2,
        "signed_block_classes_per_projective_direction": h,
        "fixed_eligible_block_hit_count": hit_count,
        "fixed_eligible_block_hit_probability": str(rho),
        "expected_raw_block_collision_surplus": str(
            expected_collision_surplus
        ),
        "expected_raw_block_collision_surplus_floor": (
            expected_collision_surplus.numerator
            // expected_collision_surplus.denominator
        ),
        "collision_surplus_definition": "S=sum_b floor(d_b/2)",
        "physical_pair_overlap_count_breakdown_per_sign": {
            "direct_zero": direct_zero,
            "direct_nonzero": direct_nonzero,
            "swapped": swapped,
        },
        "physical_pair_overlap_candidate_count": (
            physical_overlap_pair_numerator
        ),
        "physical_pair_overlap_probability": str(
            physical_overlap_pair_probability
        ),
        "expected_total_pair_overlap": str(expected_total_pair_overlap),
        "expected_total_pair_overlap_below_three_quarters": True,
        "bounded_difference_for_S": bounded_difference,
        "mcdiarmid_bad_event": "S<=2m-1",
        "mcdiarmid_exponent": str(mcdiarmid_exponent),
        "probability_proof": probability_proof,
        "bad_event_union_probability_strictly_below_one": True,
        "existence_conclusion": (
            "some auxiliary choice has no common physical inversion orbit "
            "between halves and S>=2m"
        ),
        "physical_supports_pairwise_disjoint": True,
        "sum_of_halves_is_ternary": True,
        "no_physical_cancellations_in_witness": True,
        "sigma_equals_S_for_the_witness": True,
        "sigma_lower_bound_in_witness": 2 * m,
        "maximum_kappa_zero_plus_m_plus_q_when_kappa_zero_is_zero": 2 * m,
        "standalone_global_incidence_bound_can_contradict_demand": False,
        "missing_coupling": (
            "the required physical cancellation count and its locations, "
            "or the target atom/even-moment equations"
        ),
        "endpoint_target_constructed": False,
        "residual_ii_closed": False,
        "proved": True,
    }


def theorem_record(p: int = 31) -> dict[str, object]:
    result = global_mobius_incidence_barrier(p)
    return {
        "title": "Global Mobius-incidence standalone-route barrier",
        "status": "PROVED ALL-PRIME NO-GO FOR INCIDENCE-ONLY CLOSURE",
        "result": result,
        "residual_ii_closed": False,
        "proved_all_claimed_statements": bool(result["proved"]),
    }


if __name__ == "__main__":
    from pprint import pprint

    pprint(theorem_record(), sort_dicts=True)
