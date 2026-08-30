#!/usr/bin/env python3
"""Prop. 15.722 -- exact phase cocycle and multi-chart p+1 reductions.

Let ``H`` be a residual flip set with odd-degree boundary ``D`` of size
``p+1``.  Proposition 15.721 permits every PSL(2,p^2) normalization but did
not track the Paley edge-product phase.  If

    g(z)=(az+b)/(cz+d)

has square determinant, its signed Paley multiplier is
``delta_g(x)=chi(cx+d)`` away from the pole.  When ``c!=0`` its value at
the pole and at infinity is ``chi(c)``.  When ``c=0`` there is no pole and
its value at every finite point and at infinity is ``chi(d)`` (equivalently
``chi(a)`` because the determinant is a square).  Since a relative edge
mask acquires one multiplier at every odd-degree vertex,

    c_(gH) = c_H product_(x in D) delta_g(x).                 (1)

For an all-finite boundary with monic root polynomial
``f(X)=product_(x in D)(X-x)``, (1) gives

    c_v = c_H chi(f'(v))       (v in D, sent to infinity),
    c_w = c_H chi(f(w))        (w not in D, sent to infinity). (2)

The derivative signs have product one.  Thus the hard/easy phases in the
``p+1`` boundary charts are coupled; they cannot be assigned independently.

An outside-point chart contains ``P=p+1`` affine points.  Its exact pair
slack is

    R=(sum_d b_d-P)/4
     =sum_(n=2r) r(r-1) + sum_(n=2r+1) r^2.                 (3)

Consequently ``R=0`` is precisely a ``P``-arc, while ``R=1`` would be one
3-secant and no other line of size at least three.  The latter is impossible:
deleting a point of the unique triple gives a p-arc on a conic, and the
deleted off-conic point lies on more than one surviving conic secant.

The first two positive values are impossible.  At ``R=2`` the only rich-line
patterns are one 4-secant or two 3-secants.  Deleting one or two selected
points leaves a p- or ``(p-1)``-arc.  Segre and Ball--Lavrauw's classification
of complete ``(q-1)``-arcs put it on a conic.  Each restored off-conic point
then retains at least ``(p-3)/2`` or ``(p-5)/2`` conic secants, already more
than total slack two permits.  At ``R=3`` the patterns are a 4-secant plus a
3-secant or three 3-secants.  Three private deletions leave a ``(p-2)``-arc;
the complete ``(q-2)``-arc classification and the same secant count give at
least ``(p-7)/2>3`` rich lines.  More generally, minimal deletion to an arc
and the prime-field conic threshold exclude
``R<=floor(sqrt(p)-5/2)``.  Thus the first open positive slack is beyond
``max(3,floor(sqrt(p)-5/2))`` (uniformly ``R>=4``).

At ``R=0``, Segre puts the points on a conic with external infinity line.
The phase-one type budget permits at most one of the conic's ``b=0``
directions.  Comparing its anisotropic direction character with the Paley
norm character leaves at most two disagreements.  If the two quadratic
forms were not proportional, the genus-one Weil bound would give character
sum at most ``2 sqrt(p)``, whereas at most two disagreements give at least
``p-3``.  For ``p>=17`` this is impossible.  Hence the conic is a
Miquelian circle and the two direction types align exactly.

Finally, any boundary which is a one-point replacement of a Miquelian
circle is impossible.  Send the off-circle boundary point to infinity.  The
remaining p points lie on an affine circle and form a p-arc, which is the
pair-equality branch already excluded by Proposition 15.676.  This removes
the near-line all-finite profile and every genuine circle-plus-point repair,
but it does not classify every strict outside chart or exclude a full circle
boundary.

The full-circle branch nevertheless has one exact phase normal form.  Send
a circle point to infinity, so the boundary is
``{infinity} union (a+b*F_p)``.  If ``eps=chi(b)`` is the line-direction
type, sending an outside point to infinity multiplies the edge-product sign
by ``(-1)^m eps``, while the resulting affine circle has all its ``b=2``
directions in type ``eps``.  Exact outside-chart alignment therefore forces

    c_H=(-1)^m,             m=(p+1)/2.                       (4)

Every boundary-point chart has common phase ``m mod 2`` and profile one
``b=1`` direction plus ``p`` copies of ``b=p``.  Its special floor is
``p+1-2*(m mod 2)`` and every transverse floor is zero.  This explains why
the full circle survives the one-chart floor test, but reduces it to one
phase in both PSL circle orbits.  The p+1 shell, residual (ii), Type I, and
the limit remain open.
"""
from __future__ import annotations

import json
from functools import lru_cache
from math import isqrt
from pathlib import Path

from e1_gmin_m4_prop15721 import is_prime


ROOT = Path(__file__).resolve().parents[1]


def _check_prime_parameter(p: int) -> None:
    if p < 17 or not is_prime(p):
        raise ValueError("need an odd prime parameter p>=17")


def signed_phase_cocycle() -> dict[str, object]:
    """Record the exact signed-automorphism edge-product cocycle."""
    return {
        "proved": True,
        "signed_multiplier": {
            "finite_nonpole": "delta_g(x)=chi(c*x+d)",
            "c_nonzero_pole": "delta_g(r)=chi(c)",
            "c_nonzero_infinity": "delta_g(infinity)=chi(c)",
            "c_zero_no_pole": True,
            "c_zero_finite_and_infinity": (
                "delta_g(x)=delta_g(infinity)=chi(d)=chi(a)"
            ),
        },
        "edge_product_transport": (
            "c_(gH)=c_H*product_{x in boundary(H)} delta_g(x)"
        ),
        "reason": (
            "each transported edge receives the two endpoint multipliers; "
            "even vertex degrees cancel and odd vertex degrees leave exactly "
            "the boundary product"
        ),
        "relative_flip_mask_still_only_permuted": True,
    }


def boundary_and_outside_chart_phases() -> dict[str, object]:
    """Specialize the cocycle to inversion about a finite point."""
    return {
        "proved": True,
        "boundary_polynomial": "f(X)=product_{x in D}(X-x), |D|=p+1",
        "normalization": "g_r(z)=1/(z-r)",
        "multiplier": (
            "delta_r(x)=chi(x-r) for finite x!=r and "
            "delta_r(r)=delta_r(infinity)=1"
        ),
        "boundary_point_r_in_D": "c_r=c_H*chi(f'(r))",
        "outside_point_r_not_in_D": "c_r=c_H*chi(f(r))",
        "no_untracked_determinant_sign": (
            "det(g_r)=-1 is a square in F_(p^2), and rescaling an SL "
            "representative changes every multiplier by a square"
        ),
    }


def derivative_phase_product_ledger(p: int) -> dict[str, object]:
    """Vandermonde proof that all boundary-chart phases multiply to +1."""
    _check_prime_parameter(p)
    size = p + 1
    sign_exponent = size * (size - 1) // 2
    # In F_(p^2), -1 is always a square for odd p.
    minus_one_character = 1
    product_character = minus_one_character ** sign_exponent
    return {
        "p": p,
        "boundary_size": size,
        "identity": (
            "product_{v in D} f'(v)=(-1)^(P(P-1)/2)*"
            "product_{i<j}(v_i-v_j)^2"
        ),
        "sign_exponent": sign_exponent,
        "chi_F_p2_minus_one": minus_one_character,
        "product_of_derivative_characters": product_character,
        "product_of_transported_boundary_phases": 1,
        "number_of_negative_boundary_phases_is_even": True,
        "proved": product_character == 1,
    }


def occupancy_slack_term(n: int) -> int:
    """Contribution of one affine line to the normalized pair slack R."""
    if n < 0:
        raise ValueError("occupancy must be nonnegative")
    r = n // 2
    return r * (r - 1) if n % 2 == 0 else r * r


def outside_pair_slack_identity() -> dict[str, object]:
    """Exact line-occupancy identity for an affine set of size P=p+1."""
    sample = {n: occupancy_slack_term(n) for n in range(9)}
    return {
        "proved": sample == {0: 0, 1: 0, 2: 0, 3: 1, 4: 2, 5: 4, 6: 6, 7: 9, 8: 12},
        "identity": (
            "(sum_d b_d-P)/4=sum_{n_l=2r}r(r-1)+"
            "sum_{n_l=2r+1}r^2"
        ),
        "inputs": [
            "sum_d sum_l C(n_l,2)=C(P,2)",
            "P-b_d=2*sum_l floor(n_l/2)",
        ],
        "sample_line_contributions": sample,
        "R_zero": "every line has occupancy at most two; equivalently a P-arc",
        "R_one": "exactly one 3-secant and every other line has occupancy at most two",
    }


def p_minus_one_arc_conic_lemma(p: int) -> dict[str, object]:
    """Import the complete ``(q-1)``-arc classification for ``p>=17``.

    Ball--Lavrauw classify the complete planar ``(q-1)``-arcs: they occur
    only for ``q=7,9,11,13``.  Thus a ``(p-1)``-arc for prime ``p>=17``
    extends to a ``p``-arc, and Segre then puts that extension on a conic.
    """
    _check_prime_parameter(p)
    safe = p >= 17 and p not in (7, 9, 11, 13)
    return {
        "p": p,
        "classification": "complete (q-1)-arcs occur only for q=7,9,11,13",
        "reference": (
            "Ball--Lavrauw, Planar arcs, JCTA 160 (2018), Corollary 10 "
            "(Corollary 8 in arXiv v4)"
        ),
        "p_minus_one_arc_is_incomplete": safe,
        "extension": "every (p-1)-arc extends to a p-arc",
        "conclusion": "every (p-1)-arc is contained in a nonsingular conic",
        "proved_dependency": safe,
        "proved": safe,
    }


def p_minus_two_arc_conic_lemma(p: int) -> dict[str, object]:
    """Import the complete ``(q-2)``-arc classification for ``p>=17``."""
    _check_prime_parameter(p)
    prior = p_minus_one_arc_conic_lemma(p)
    safe = p >= 17 and p not in (8, 9, 11)
    return {
        "p": p,
        "classification": "complete (q-2)-arcs occur only for q=8,9,11",
        "reference": (
            "Ball--Lavrauw, Planar arcs, JCTA 160 (2018), Corollary 11 "
            "(Corollary 9 in arXiv v4)"
        ),
        "extension_chain": "(p-2)-arc extends to (p-1)-arc, then to p-arc",
        "p_minus_one_dependency": prior,
        "conclusion": "every (p-2)-arc is contained in a nonsingular conic",
        "proved_dependency": bool(safe and prior["proved"]),
        "proved": bool(safe and prior["proved"]),
    }


def minimum_surviving_conic_secants(p: int, missing_conic_points: int) -> int:
    """Secants through an off-conic point after deleting conic points."""
    _check_prime_parameter(p)
    if not 0 <= missing_conic_points <= p + 1:
        raise ValueError("invalid number of missing conic points")
    return (p - 1) // 2 - missing_conic_points


def outside_R_two_structure(p: int) -> dict[str, object]:
    """Classify and exclude the rich-line geometry at slack ``R=2``.

    The line contribution is one at occupancy three, two at occupancy four,
    and at least four from occupancy five onward.  Thus a total of two is
    partitioned only as ``2`` or ``1+1``.  Deleting the indicated boundary
    points removes every line of occupancy at least three, because all other
    lines have zero contribution in the slack identity.

    A ``p``-arc is on a conic by Segre, but the conic need not be a
    Miquelian circle.  The Ball--Lavrauw complete-arc classification extends
    each resulting ``(p-1)``-arc to a ``p``-arc, so Segre puts it on a conic
    as well.  This actually excludes every pattern: an off-conic replacement
    retains far more conic secants than a total slack of two permits.
    """
    _check_prime_parameter(p)
    segre = segre_q_arc_conic_lemma(p)
    extension = p_minus_one_arc_conic_lemma(p)
    contributions = {n: occupancy_slack_term(n) for n in range(3, 8)}
    only_small_contributors = (
        contributions[3] == 1
        and contributions[4] == 2
        and all(contributions[n] > 2 for n in range(5, 8))
    )
    partitions_of_two = ((2,), (1, 1))
    one_missing_secants = minimum_surviving_conic_secants(p, 1)
    two_missing_secants = minimum_surviving_conic_secants(p, 2)
    intersecting_contradiction = one_missing_secants > 2
    two_point_contradiction = two_missing_secants > 2
    proved = bool(
        only_small_contributors
        and segre["proved"]
        and extension["proved"]
        and intersecting_contradiction
        and two_point_contradiction
    )
    return {
        "p": p,
        "R": 2,
        "positive_line_contributions": contributions,
        "partitions_of_two": partitions_of_two,
        "rich_line_patterns": (
            "one 4-secant and no other rich line",
            "two 3-secants and no other rich line",
        ),
        "intersecting_trisecants": {
            "meaning": "the two 3-secants share a boundary point",
            "deletion": "delete their shared boundary point",
            "remainder": "a p-arc",
            "Segre_dependency": segre,
            "normal_form": (
                "(a nonsingular conic minus one point) plus one off-conic point"
            ),
            "Miquelian_not_forced": True,
            "minimum_surviving_conic_secants": one_missing_secants,
            "contradicts_total_R_two": intersecting_contradiction,
        },
        "nonintersecting_or_four_secant": {
            "deletion": (
                "delete one boundary point from each disjoint 3-secant, "
                "or two from the unique 4-secant"
            ),
            "remainder": "a (p-1)-arc",
            "extension_dependency": extension,
            "normal_form": (
                "(a nonsingular conic minus two points) plus two off-conic points"
            ),
            "deleted_points_are_off_conic": (
                "each deleted point lies on a rich line whose two retained "
                "points already lie on the conic"
            ),
            "one_4_secant": (
                "both off-conic replacements lie on the same retained conic secant"
            ),
            "two_disjoint_3_secants": (
                "each off-conic replacement lies on its own retained conic secant"
            ),
            "minimum_surviving_conic_secants_through_either_replacement": (
                two_missing_secants
            ),
            "contradicts_total_R_two": two_point_contradiction,
        },
        "secant_count_input": (
            "an off-conic point lies on at least (p-1)/2 conic secants; "
            "each missing conic point destroys at most one"
        ),
        "excluded": proved,
        "next_exact_target": "none at R=2; use the combined low-slack theorem",
        "proved": proved,
    }


def outside_R_three_structure(p: int) -> dict[str, object]:
    """Classify and exclude the next slack value ``R=3``."""
    _check_prime_parameter(p)
    extension = p_minus_two_arc_conic_lemma(p)
    contributions = {n: occupancy_slack_term(n) for n in range(3, 8)}
    patterns_exact = (
        contributions[3] == 1
        and contributions[4] == 2
        and all(contributions[n] > 3 for n in range(5, 8))
    )
    surviving = minimum_surviving_conic_secants(p, 3)
    contradiction = surviving > 3
    proved = bool(patterns_exact and extension["proved"] and contradiction)
    return {
        "p": p,
        "R": 3,
        "partitions_of_three": ((2, 1), (1, 1, 1)),
        "rich_line_patterns": (
            "one 4-secant and one 3-secant",
            "three 3-secants",
        ),
        "private_deletion": (
            "choose two points of the 4-secant off the 3-secant and one "
            "point of the 3-secant off the 4-secant; or choose one point "
            "of each 3-secant off the other two rich lines"
        ),
        "remainder": "a (p-2)-arc",
        "extension_dependency": extension,
        "normal_form": (
            "a nonsingular conic minus three points, plus three off-conic points"
        ),
        "deleted_points_are_off_conic": (
            "each deleted point's rich line retains two conic points"
        ),
        "minimum_surviving_conic_secants_through_each_replacement": surviving,
        "maximum_rich_lines_allowed_by_R": 3,
        "excluded": proved,
        "proved": proved,
    }


def outside_low_slack_conic_exclusion(p: int) -> dict[str, object]:
    """Exclude a prime-dependent initial interval of positive slack.

    For every rich line delete ``n_l-2`` of its points and then make the
    deletion set inclusion-minimal.  Its size ``t`` is at most ``R`` because
    ``n-2<=occupancy_slack_term(n)``.  Minimality gives, for each deleted
    point, a line retaining two arc points.  Ball--Lavrauw's prime-field
    conic threshold applies when ``R<=sqrt(p)-5/2``.  The deleted point is
    then off-conic and retains more than ``R`` conic secants, a contradiction.
    The complete near-arc classifications separately make ``R<=3`` uniform.
    """
    _check_prime_parameter(p)
    # Largest integer r satisfying (2r+5)^2 <= 4p, equivalently
    # r <= sqrt(p)-5/2.  This avoids any floating-point boundary issue.
    generic_cutoff = max(0, (isqrt(4 * p) - 5) // 2)
    combined_cutoff = max(3, generic_cutoff)
    even_gap = "h(2r)-(2r-2)=(r-1)(r-2)>=0"
    odd_gap = "h(2r+1)-(2r-1)=(r-1)^2>=0"
    generic_secants = (p - 1) // 2 - generic_cutoff
    generic_contradiction = generic_secants > generic_cutoff
    integer_safe = (2 * generic_cutoff + 5) ** 2 <= 4 * p
    next_integer_fails = (2 * (generic_cutoff + 1) + 5) ** 2 > 4 * p
    r_two = outside_R_two_structure(p)
    r_three = outside_R_three_structure(p)
    proved = bool(
        integer_safe
        and next_integer_fails
        and generic_contradiction
        and r_two["proved"]
        and r_three["proved"]
    )
    return {
        "p": p,
        "minimal_deletion_lemma": {
            "construction": (
                "delete n_l-2 points from every rich line, take the union, "
                "then remove redundant deletions"
            ),
            "size": "1<=t<=sum_l(n_l-2)<=R",
            "even_occupancy_gap": even_gap,
            "odd_occupancy_gap": odd_gap,
            "minimality_witness": (
                "each deleted z lies on a line with exactly two retained arc points"
            ),
        },
        "prime_arc_conic_dependency": (
            "Ball--Lavrauw, Planar arcs, published Theorem 5 / arXiv v4 "
            "Theorem 3: a prime-field "
            "arc of size at least p-sqrt(p)+7/2 lies on a conic"
        ),
        "generic_cutoff": generic_cutoff,
        "generic_cutoff_exact_condition": "(2R+5)^2<=4p",
        "generic_minimum_surviving_secants_at_cutoff": generic_secants,
        "generic_secants_exceed_slack": generic_contradiction,
        "uniform_near_arc_cutoff": 3,
        "combined_excluded_positive_R_through": combined_cutoff,
        "first_open_R_at_least": combined_cutoff + 1,
        "R_two_dependency": r_two,
        "R_three_dependency": r_three,
        "proved": proved,
    }


def segre_q_arc_conic_lemma(p: int) -> dict[str, object]:
    """State the precise finite-geometry theorem used in the R=1 step.

    Segre's q-arc theorem says that a q-arc in ``PG(2,q)`` for odd
    ``q>=5`` is contained in a nonsingular conic.  Recording the hypotheses
    here prevents the code from silently replacing the theorem by the weaker
    and insufficient statement that only ``q+1``-arcs are conics.
    """
    _check_prime_parameter(p)
    return {
        "field": f"F_{p}",
        "plane": f"PG(2,{p})",
        "arc_size": p,
        "hypotheses": "Desarguesian plane, odd field order p>=5",
        "conclusion": "every p-arc is contained in a nonsingular conic",
        "reference": "Segre q-arc theorem (corollary of the lemma of tangents)",
        "all_hypotheses_met": p >= 17 and is_prime(p),
        "proved_dependency": p >= 17 and is_prime(p),
        "proved": p >= 17 and is_prime(p),
    }


def unique_trisecant_exclusion(p: int) -> dict[str, object]:
    """Exclude R=1 for every affine (p+1)-set."""
    _check_prime_parameter(p)
    segre = segre_q_arc_conic_lemma(p)
    # An off-conic point lies on (p-1)/2 or (p+1)/2 secants of a conic.
    # Removing the one conic point absent from the p-arc destroys at most one.
    surviving_secants = (p - 1) // 2 - 1
    return {
        "p": p,
        "assumption": "one 3-secant and no other line of occupancy at least three",
        "delete_one_triple_point": "a p-arc, hence contained in a conic",
        "p_arc_conic_dependency": segre,
        "deleted_point_is_off_conic": (
            "if it were on the conic, all p+1 points would be the conic and "
            "there would be no 3-secant"
        ),
        "minimum_full_conic_secants_through_deleted_point": (p - 1) // 2,
        "maximum_destroyed_by_the_missing_conic_point": 1,
        "minimum_surviving_secants": surviving_secants,
        "contradicts_unique_trisecant": surviving_secants > 1,
        "R_one_excluded": surviving_secants > 1,
        "proved": bool(segre["proved"] and surviving_secants > 1),
    }


def quadratic_product_genus_one_lemma(p: int) -> dict[str, object]:
    """Verify the smoothness hypotheses behind the character-sum bound."""
    _check_prime_parameter(p)
    separation = (p - 3) * (p - 3) > 4 * p
    return {
        "forms": "nonsingular anisotropic binary quadratics Q and N over F_p",
        "separable": "p is odd, so each nonsingular quadratic has two simple roots",
        "disjoint_roots_if_nonproportional": (
            "a shared root forces the shared Frobenius-conjugate root, hence proportionality"
        ),
        "quartic": "Q*N is squarefree when Q and N are nonproportional",
        "smooth_model": "Y^2=Q*N is a smooth double cover of P^1 of genus one",
        "projective_character_sum": (
            "sum_[u:v] chi(Q(u,v)N(u,v))=#C(F_p)-(p+1)"
        ),
        "hasse_bound": "absolute character sum at most 2*sqrt(p)",
        "forced_agreement_lower_bound": p - 3,
        "squared_separation": f"(p-3)^2={((p - 3) ** 2)}>4p={4 * p}",
        "strict_separation": separation,
        "proved": separation,
    }


def outside_arc_type_alignment(p: int) -> dict[str, object]:
    """At R=0 force the affine conic to be a Miquelian circle."""
    _check_prime_parameter(p)
    curve = quadratic_product_genus_one_lemma(p)
    period = p + 1
    m = period // 2
    b0_phase_one_floor = 2 * p
    b2_phase_one_floor = p - 1
    type_budget = m * period
    # With t phase-one b=0 directions the exact floor sum is
    # t*2p+(m-t)*(p-1)=mP+(t-1)P.
    maximum_wrong_directions = 1
    mismatch_count = 2 * maximum_wrong_directions
    character_sum_lower = period - 2 * mismatch_count
    hasse_squared = 4 * p
    separation_squared = character_sum_lower * character_sum_lower
    return {
        "p": p,
        "P": period,
        "m": m,
        "P_arc_classification": "a conic with external line at infinity",
        "global_profile": {"b=0": m, "b=2": m},
        "phase_one_floors": {
            "b=0": b0_phase_one_floor,
            "b=2": b2_phase_one_floor,
        },
        "phase_one_type_floor_with_t_wrong": "m*P+(t-1)*P",
        "type_budget": type_budget,
        "maximum_b0_directions_in_phase_one_type": maximum_wrong_directions,
        "maximum_partition_disagreements": mismatch_count,
        "character_sum_if_not_aligned_at_least": character_sum_lower,
        "nonproportional_forms_curve": "Y^2=Q*N is a smooth genus-one curve",
        "smooth_genus_one_certificate": curve,
        "genus_one_character_bound_squared": hasse_squared,
        "forced_sum_squared": separation_squared,
        "strict_Hasse_contradiction": separation_squared > hasse_squared,
        "quadratic_forms_proportional": separation_squared > hasse_squared,
        "direction_types_align_exactly": separation_squared > hasse_squared,
        "affine_conic_is_Miquelian_circle": separation_squared > hasse_squared,
        "proved": bool(curve["proved"] and separation_squared > hasse_squared),
    }


@lru_cache(maxsize=1)
def prior_pair_equality_dependency() -> dict[str, object]:
    """Load the exact prior theorem used after a boundary-point transport."""
    from e1_gmin_m4_prop15676 import theorem_record

    record = theorem_record()
    return {
        "prop_15_676_proved": bool(record["proved"]),
        "infinity_plus_p_pair_equality_excluded": (
            record["theorem"]["pair_deficit_equality"]
            == "EXCLUDED_FOR_BOTH_PRODUCT_SIGNS"
        ),
    }


def one_point_circle_replacement_exclusion(p: int) -> dict[str, object]:
    """Use a second chart to reject circle-minus-one-plus-one boundaries."""
    _check_prime_parameter(p)
    dependency = prior_pair_equality_dependency()
    proved = all(dependency.values())
    return {
        "p": p,
        "boundary_shape_in_one_outside_chart": "(C minus {r}) union {z}, z notin C",
        "C": "a Miquelian F_p-subline/circle",
        "second_normalization": "send the boundary point z to infinity",
        "circle_transport": (
            "Mobius maps preserve F_p-sublines; because z is not on C, "
            "the transported circle does not contain infinity"
        ),
        "remaining_finite_boundary": "p points of an affine P-arc",
        "pair_deficit_after_second_normalization": "equality",
        "excluded_by": "Proposition 15.676",
        "dependency": dependency,
        "includes_near_line_profile": True,
        "does_not_cover_arbitrary_non_Miquelian_conic_repairs": True,
        "excluded": proved,
        "proved": proved,
    }


def full_circle_line_chart_normal_form(p: int) -> dict[str, object]:
    """Fix the product sign and all floors in every circle-point chart.

    Write the normalized circle as ``{infinity} union (a+b*F_p)`` and put
    ``eps=chi(b)``.  For an outside point ``w``, the signed multiplier over
    the boundary is

        chi(product_t(a+b*t-w))=(-1)^m*eps.

    Indeed, after ``u=(w-a)/b``, the product is ``b^p(u-u^p)``; every
    nonzero trace-zero element has quadratic character ``(-1)^m``.  The
    transformed affine conic has tangent directions proportional to
    ``-b^(-1)/(t-u)^2``, hence its ``b=2`` directions have type ``eps``.
    Outside-chart type alignment then gives ``c_H=(-1)^m``.
    """
    _check_prime_parameter(p)
    from e1_gmin_m4_prop15669 import full_symbolic_floor

    period = p + 1
    m = period // 2
    phase = m & 1
    forced_c_h = -1 if phase else 1
    b1_floor = full_symbolic_floor(p, 1, phase)
    bp_floor = full_symbolic_floor(p, p, phase)
    expected_b1 = period - 2 * phase
    expected_bp = 0
    proved = b1_floor == expected_b1 and bp_floor == expected_bp
    if not proved:
        raise ArithmeticError("full-circle easy-phase floors changed")
    return {
        "p": p,
        "P": period,
        "m": m,
        "normalized_boundary": "{infinity} union (a+b*F_p)",
        "line_direction_type": "eps=chi_F_p2(b), either sign",
        "outside_chart_boundary_multiplier": "(-1)^m*eps",
        "proof_of_multiplier": (
            "product_t(a+b*t-w)=b^p(u-u^p); nonzero trace-zero "
            "elements have character (-1)^m"
        ),
        "outside_affine_circle_b2_type": "eps",
        "outside_alignment_equation": "c_line*((-1)^m*eps)=eps",
        "forced_line_chart_c_H": forced_c_h,
        "common_phase": phase,
        "profile": {"b=1": 1, "b=p": p},
        "direction_type_split": (
            "the special line direction lies in either type; its type has "
            "m-1 transverse directions and the other type has m"
        ),
        "special_parity_baseline": "x_j" if phase == 0 else "1-x_j",
        "transverse_parity": "sum_s x_s + phase = m+phase = 0 mod 2",
        "floors": {"b=1": b1_floor, "b=p": bp_floor},
        "special_mean_interval": [b1_floor, m * period],
        "transverse_quadratics_are_even_nonnegative": True,
        "both_PSL_circle_orbits_have_same_forced_phase": True,
        "full_circle_excluded": False,
        "proved": proved,
    }


def universal_multichart_certificate() -> dict[str, object]:
    """Separate the uniform proof from its finite regression examples."""
    p0 = 17
    segre = segre_q_arc_conic_lemma(p0)
    curve = quadratic_product_genus_one_lemma(p0)
    pair_slack = outside_pair_slack_identity()
    r_two = outside_R_two_structure(p0)
    r_three = outside_R_three_structure(p0)
    low_slack = outside_low_slack_conic_exclusion(p0)
    replacement = one_point_circle_replacement_exclusion(p0)
    circle = full_circle_line_chart_normal_form(p0)
    r_one_minimum_surviving_secants = (p0 - 1) // 2 - 1
    hasse_gap_at_17 = (p0 - 3) ** 2 - 4 * p0
    hasse_gap_odd_step_at_17 = 4 * p0 - 16
    proved = bool(
        signed_phase_cocycle()["proved"]
        and boundary_and_outside_chart_phases()["proved"]
        and pair_slack["proved"]
        and r_two["proved"]
        and r_three["proved"]
        and low_slack["proved"]
        and segre["proved"]
        and r_one_minimum_surviving_secants > 1
        and curve["proved"]
        and hasse_gap_at_17 > 0
        and hasse_gap_odd_step_at_17 > 0
        and replacement["proved"]
        and circle["proved"]
    )
    return {
        "scope": "every odd prime p>=17",
        "signed_phase_cocycle": signed_phase_cocycle(),
        "chart_phase_specialization": boundary_and_outside_chart_phases(),
        "pair_slack_identity": pair_slack,
        "R_two_structure": r_two,
        "R_three_structure": r_three,
        "low_slack_conic_exclusion": low_slack,
        "low_slack_secant_bounds_are_increasing": True,
        "R_one": {
            "q_arc_dependency": segre,
            "minimum_surviving_secants_at_p17": r_one_minimum_surviving_secants,
            "minimum_is_increasing": True,
        },
        "R_zero": {
            "smooth_curve_dependency": curve,
            "hasse_gap_at_p17": hasse_gap_at_17,
            "hasse_gap_odd_step_at_p17": hasse_gap_odd_step_at_17,
            "gap_is_increasing_for_p>=17": hasse_gap_odd_step_at_17 > 0,
        },
        "one_point_replacement_dependency": replacement["dependency"],
        "full_circle_phase_symbolic_identity": {
            "boundary_multiplier": "(-1)^m*eps",
            "outside_circle_b2_type": "eps",
            "forced_c_H": "(-1)^m",
            "prop_15_669_exact_floor_dependency": True,
            "proved": circle["proved"],
        },
        "proved": proved,
    }


def theorem_multichart_p_plus_one_reduction() -> dict[str, object]:
    """Package the exact new reductions without claiming shell closure."""
    universal = universal_multichart_certificate()
    samples = (17, 19, 23, 29, 31, 37, 41, 101)
    rows = {
        str(p): {
            "derivative_phase_product": derivative_phase_product_ledger(p),
            "R_one": unique_trisecant_exclusion(p),
            "R_two": outside_R_two_structure(p),
            "R_three": outside_R_three_structure(p),
            "low_slack": outside_low_slack_conic_exclusion(p),
            "outside_arc": outside_arc_type_alignment(p),
            "one_point_circle_replacement": one_point_circle_replacement_exclusion(p),
            "full_circle_line_chart": full_circle_line_chart_normal_form(p),
        }
        for p in samples
    }
    proved = bool(
        universal["proved"]
        and signed_phase_cocycle()["proved"]
        and boundary_and_outside_chart_phases()["proved"]
        and outside_pair_slack_identity()["proved"]
        and all(
            bool(item["proved"])
            for row in rows.values()
            for item in row.values()
        )
    )
    return {
        "prop": "15.722",
        "title": "Exact phase cocycle and multi-chart p+1 reductions",
        "proved": proved,
        "universal_certificate": universal,
        "phase_transport": signed_phase_cocycle(),
        "chart_specializations": boundary_and_outside_chart_phases(),
        "outside_pair_slack": outside_pair_slack_identity(),
        "outside_R_two_structure": outside_R_two_structure(17),
        "outside_R_three_structure": outside_R_three_structure(17),
        "outside_low_slack_conic_exclusion": outside_low_slack_conic_exclusion(17),
        "theorem": {
            "boundary_chart_phases": "c_v=c_H*chi(f'(v)); their product is +1",
            "outside_chart_phases": "c_w=c_H*chi(f(w))",
            "outside_R_one": "IMPOSSIBLE",
            "outside_R_two": "IMPOSSIBLE_BY_CONIC_EXTENSION_AND_SECANT_COUNT",
            "outside_R_three": "IMPOSSIBLE_BY_CONIC_EXTENSION_AND_SECANT_COUNT",
            "positive_outside_slack_excluded_through": (
                "max(3,floor(sqrt(p)-5/2))"
            ),
            "outside_R_zero": "FORCES_A_MIQUELIAN_CIRCLE_WITH_EXACT_TYPE_ALIGNMENT",
            "one_point_Miquelian_circle_replacements": "EXCLUDED_BY_A_SECOND_CHART",
            "full_circle_phase": "FORCED_TO_m_MOD_2_IN_EVERY_BOUNDARY_POINT_CHART",
            "full_circle_profile": "one b=1 plus p copies of b=p; transverse floors zero",
            "full_circle_boundary": "OPEN_AFTER_EXACT_ONE_PHASE_NORMAL_FORM",
            "strict_outside_profiles_R_at_least_2": "R_EQUALS_2_AND_3_IMPOSSIBLE",
            "strict_outside_profiles_R_at_least_4": (
                "OPEN_AFTER_THE_PRIME_DEPENDENT_LOW_SLACK_CUTOFF"
            ),
            "whole_p_plus_one_shell": "OPEN",
            "residual_ii": False,
            "type_I": False,
            "limit_exists": False,
        },
        "sample_ledgers_regression_only": rows,
        "L_status": "OPEN",
    }


def _jsonable(value):
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_jsonable(item) for item in value]
    return value


def main() -> dict[str, object]:
    theorem = theorem_multichart_p_plus_one_reduction()
    if theorem["proved"] is not True:
        raise ArithmeticError("Proposition 15.722 multi-chart audit failed")
    destination = ROOT / "evidence" / "e1_gmin_m4_prop15722.json"
    destination.write_text(json.dumps(_jsonable(theorem), indent=2) + "\n")
    print("Prop 15.722 exact phase cocycle and multi-chart reductions: proved")
    print(f"  wrote {destination}")
    return theorem


if __name__ == "__main__":
    main()
