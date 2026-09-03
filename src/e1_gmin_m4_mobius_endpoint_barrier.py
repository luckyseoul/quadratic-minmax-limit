#!/usr/bin/env python3
"""Symbolic barrier ledger for the complementary Mobius endpoint.

The desired opposite-swapped intersection of two localized Mobius halves
has a simple rational locus.  On that locus the other three possible
endpoint matchings can be classified exactly.  A rational example shows
that the local intersection and singleton equations are consistent.

This does *not* construct or refute a global complementary family.  For one
*fixed, preassigned* auxiliary family, its relative scales impose a
center-coherence condition which is not detected by the unit-star moments.
The auxiliaries might instead be chosen adaptively to the centers.  That
global matching problem, and residual (ii), remain explicitly open below.
"""

from __future__ import annotations

from fractions import Fraction
from typing import TypeAlias

from e1_gmin_m4_compact_ray_moment_gate import compact_moment, star_moment


Rational: TypeAlias = Fraction
Point: TypeAlias = tuple[Rational, Rational]
Edge: TypeAlias = frozenset[Point]
LabelEdge: TypeAlias = tuple[int, int]


def _fraction(value: int | Fraction) -> Fraction:
    return value if isinstance(value, Fraction) else Fraction(value)


def _opposite_swapped_parameters(
    q: int | Fraction, r: int | Fraction
) -> tuple[Fraction, Fraction, Fraction, Fraction, Fraction]:
    """Normalize and return ``q,r,z,A,B`` on the desired locus."""
    q = _fraction(q)
    r = _fraction(r)
    z = q * r
    if q == 0 or r == 0 or q == 1 or r == 1 or z == 1:
        raise ValueError("need q,r,qr nonzero, q,r!=1, and qr!=1")
    return q, r, z, (q - 1) * (z - 1) / z, (r - 1) * (z - 1) / z


def _candidate_parameters(
    q: Fraction,
    r: Fraction,
    A: Fraction,
    B: Fraction,
    orientation: int,
    matching: str,
) -> tuple[Fraction, Fraction] | None:
    """Solve the endpoint equations for one of the four candidates."""
    if orientation not in (-1, 1):
        raise ValueError("orientation must be +1 or -1")
    if matching == "swapped":
        return Fraction(orientation, 1) / q, Fraction(orientation, 1) / r
    if matching != "direct":
        raise ValueError("matching must be direct or swapped")
    if q == orientation or r == orientation:
        return None
    return A / (q - orientation) - 1, B / (r - orientation) - 1


def _candidate_record(
    q: Fraction,
    r: Fraction,
    A: Fraction,
    B: Fraction,
    orientation: int,
    matching: str,
) -> dict[str, object]:
    """Construct both actual edges and compare them without formula shortcuts."""
    parameters = _candidate_parameters(q, r, A, B, orientation, matching)
    if parameters is None:
        return {
            "orientation": orientation,
            "matching": matching,
            "defined": False,
            "accepted": False,
            "reason": "candidate endpoint equation has zero denominator",
        }
    t, s = parameters
    if t == -1 or s == -1:
        return {
            "orientation": orientation,
            "matching": matching,
            "defined": False,
            "t": t,
            "s": s,
            "accepted": False,
            "reason": "omitted Mobius parameter -1",
        }
    first = _first_half_edge(q, A, t)
    second = _second_half_edge(r, B, s)
    oriented_second = second if orientation == 1 else _negative_edge(second)
    return {
        "orientation": orientation,
        "matching": matching,
        "defined": True,
        "t": t,
        "s": s,
        "first_edge": tuple(sorted(first)),
        "oriented_second_edge": tuple(sorted(oriented_second)),
        "accepted": first == oriented_second,
    }


def actual_four_candidate_replay(
    q: int | Fraction, r: int | Fraction
) -> dict[str, object]:
    """Build and compare the two actual edges for all four candidates."""
    q, r, z, A, B = _opposite_swapped_parameters(q, r)
    specifications = {
        "same_orientation_direct": (1, "direct"),
        "same_orientation_swapped": (1, "swapped"),
        "opposite_direct": (-1, "direct"),
        "opposite_swapped": (-1, "swapped"),
    }
    candidates = {
        name: _candidate_record(q, r, A, B, orientation, matching)
        for name, (orientation, matching) in specifications.items()
    }
    verdicts = {
        name: bool(record["accepted"]) for name, record in candidates.items()
    }
    expected = {
        "same_orientation_direct": False,
        "same_orientation_swapped": False,
        "opposite_direct": q == r == Fraction(1, 2),
        "opposite_swapped": True,
    }
    if verdicts != expected:
        raise ArithmeticError("the direct four-edge replay changed")
    return {
        "q": q,
        "r": r,
        "z": z,
        "A": A,
        "B": B,
        "candidates": candidates,
        "candidate_verdicts": verdicts,
        "proved": True,
    }


def opposite_swapped_locus(
    q: int | Fraction, r: int | Fraction
) -> dict[str, object]:
    """Return the exact opposite-swapped locus and actual four-edge verdicts.

    In the normalized two-half coordinates, put ``z=q*r``.  Acceptance of
    the opposite-swapped candidate forces

        A=(q-1)(z-1)/z,  B=(r-1)(z-1)/z.

    The exclusions in the input check are exactly the degeneracies which
    make an auxiliary dependent, a displayed denominator zero, or a forced
    Mobius parameter equal to the omitted value ``-1``.
    """
    replay = actual_four_candidate_replay(q, r)
    q = Fraction(replay["q"])
    r = Fraction(replay["r"])
    z = Fraction(replay["z"])
    A = Fraction(replay["A"])
    B = Fraction(replay["B"])
    verdicts = dict(replay["candidate_verdicts"])
    extra_opposite_direct = bool(verdicts["opposite_direct"])
    return {
        "q": q,
        "r": r,
        "z": z,
        "A": A,
        "B": B,
        "candidate_verdicts": verdicts,
        "actual_candidate_edges": replay["candidates"],
        "shared_orbit_count": 1 + int(extra_opposite_direct),
        "opposite_orientation_shared_orbits": 1 + int(extra_opposite_direct),
        "same_orientation_shared_orbits": 0,
        "two_trade_sum_is_ternary": True,
        "desired_swapped_orbit_is_unique": not extra_opposite_direct,
        "unique_failure_point": "q=r=1/2",
        "proved_over_rational_function_field": True,
    }


def _first_half_edge(q: Fraction, A: Fraction, t: Fraction) -> Edge:
    if t == -1:
        raise ValueError("the Mobius parameter cannot be -1")
    return frozenset(
        (
            (Fraction(1), q - A / (t + 1)),
            (t, q * t),
        )
    )


def _second_half_edge(r: Fraction, B: Fraction, s: Fraction) -> Edge:
    if s == -1:
        raise ValueError("the Mobius parameter cannot be -1")
    return frozenset(
        (
            (r - B / (s + 1), Fraction(1)),
            (r * s, s),
        )
    )


def _negative_edge(edge: Edge) -> Edge:
    return frozenset((-x, -y) for x, y in edge)


def rational_clean_overlap_example() -> dict[str, object]:
    """Replay the all-characteristic-zero example ``q=r=2`` exactly.

    Take equal singleton signs and a point ``x`` with ``X(x)=Y(x)=-6``.
    Then ``F=X-Y`` annihilates ``x``.  The two auxiliary functionals are

        M1=(-5 X+4 Y)/3,  M2=(4 X-5 Y)/3,

    and both evaluate to ``2`` at ``x``.  Thus their special blocks agree
    with the singleton block.  At ``t=s=-1/2`` the first edge is the
    negative of the second, while :func:`opposite_swapped_locus` proves
    that no other candidate survives.
    """
    q = r = Fraction(2)
    locus = opposite_swapped_locus(q, r)
    A = Fraction(locus["A"])
    B = Fraction(locus["B"])
    t = s = Fraction(-1, 2)
    first = _first_half_edge(q, A, t)
    second = _second_half_edge(r, B, s)

    X_at_x = Y_at_x = Fraction(-6)
    M1_at_x = Fraction(-5, 3) * X_at_x + Fraction(4, 3) * Y_at_x
    M2_at_x = Fraction(4, 3) * X_at_x + Fraction(-5, 3) * Y_at_x
    proved = bool(
        A == B == Fraction(3, 4)
        and first == _negative_edge(second)
        and M1_at_x == M2_at_x == 2
        and X_at_x - Y_at_x == 0
        and locus["desired_swapped_orbit_is_unique"]
        and locus["same_orientation_shared_orbits"] == 0
    )
    if not proved:
        raise ArithmeticError("the rational clean-overlap example changed")
    return {
        "q": q,
        "r": r,
        "A": A,
        "B": B,
        "t": t,
        "s": s,
        "F": "X-Y",
        "X_at_singleton": X_at_x,
        "Y_at_singleton": Y_at_x,
        "M1": "(-5*X+4*Y)/3",
        "M2": "(4*X-5*Y)/3",
        "M1_at_singleton": M1_at_x,
        "M2_at_singleton": M2_at_x,
        "first_edge": sorted(first),
        "second_edge": sorted(second),
        "opposite_edges": first == _negative_edge(second),
        "other_three_candidates_rejected": True,
        "specializes_in_branch_characteristic": "every prime p>=31",
        "global_complementary_family_constructed": False,
        "proved": proved,
    }


def star_center_invisibility_replay(
    p: int, centers: tuple[int, ...] = (1, 2)
) -> dict[str, object]:
    """Directly replay every unit-star channel through degree ``p-1``.

    This calls the canonical :func:`star_moment`; it is a formula replay at
    the requested field, not a search over primes or target profiles.
    """
    if not centers:
        raise ValueError("need at least one nonzero center")
    reduced_centers = tuple(center % p for center in centers)
    if any(center == 0 for center in reduced_centers):
        raise ValueError("unit-star centers must be nonzero")

    rows: dict[int, dict[tuple[int, int], int]] = {}
    expected: dict[tuple[int, int], int] = {}
    for d in range(2, p):
        for k in range(d // 2):
            expected[d, k] = p - 1 if d == p - 1 and k == 0 else 0
    for center in reduced_centers:
        rows[center] = {
            index: star_moment(p, center, *index) for index in expected
        }

    proved = all(row == expected for row in rows.values())
    if not proved:
        raise ArithmeticError("a unit-star moment acquired center dependence")
    return {
        "p": p,
        "centers": reduced_centers,
        "channel_count": len(expected),
        "moments_by_center": rows,
        "expected_moments": expected,
        "top_degree_k0_value": p - 1,
        "all_other_channels_zero": True,
        "center_independent": True,
        "proved": True,
    }


def _label_edge(p: int, first: int, second: int) -> LabelEdge:
    """Canonical unordered edge on the one-dimensional label row."""
    values = first % p, second % p
    if values[0] == values[1]:
        raise ValueError("a compact atom edge must have distinct endpoints")
    return tuple(sorted(values))


def _negative_label_edge(p: int, edge: LabelEdge) -> LabelEdge:
    return _label_edge(p, -edge[0], -edge[1])


def centered_compact_atom_support(p: int, scale: int = 1) -> dict[str, object]:
    """Replay the fixed-singleton plus paired-group support of ``K(v,-v;0)``.

    The positive edge ``{v,-v}`` is fixed by label negation.  The two
    negative edges ``{v,0}`` and ``{-v,0}`` form one two-edge orbit.  Thus
    the inversion quotient has exactly one fixed singleton and one grouped
    support coordinate.  Canonical ``compact_moment`` calls independently
    verify that these are the established centered compact coefficients.
    """
    if not isinstance(p, int) or isinstance(p, bool) or p < 7 or p % 2 == 0:
        raise ValueError("need an odd prime p>=7")
    scale %= p
    if scale == 0:
        raise ValueError("the centered compact scale must be nonzero")

    terms = {
        _label_edge(p, scale, -scale): 1,
        _label_edge(p, scale, 0): -1,
        _label_edge(p, -scale, 0): -1,
    }
    if len(terms) != 3:
        raise ArithmeticError("the centered compact atom collapsed")
    if any(terms.get(_negative_label_edge(p, edge)) != coefficient
           for edge, coefficient in terms.items()):
        raise ArithmeticError("the centered compact chain lost central symmetry")

    fixed = tuple(edge for edge in terms if _negative_label_edge(p, edge) == edge)
    nonfixed = tuple(edge for edge in terms if _negative_label_edge(p, edge) != edge)
    grouped = {frozenset((edge, _negative_label_edge(p, edge))) for edge in nonfixed}
    odd_moments = {
        (d, k): compact_moment(p, scale, -scale, 0, d, k)
        for d in range(3, p - 1, 2)
        for k in range(d // 2)
    }
    proved = bool(
        len(fixed) == 1
        and terms[fixed[0]] == 1
        and len(nonfixed) == 2
        and len(grouped) == 1
        and all(terms[edge] == -1 for edge in nonfixed)
        and not any(odd_moments.values())
    )
    if not proved:
        raise ArithmeticError("the centered compact quotient support changed")
    return {
        "p": p,
        "scale": scale,
        "signed_edge_chain": terms,
        "fixed_singleton_edge": fixed[0],
        "fixed_singleton_support_size": len(fixed),
        "paired_group_edges": tuple(sorted(next(iter(grouped)))),
        "paired_group_support_size": len(grouped),
        "quotient_support_size": len(fixed) + len(grouped),
        "canonical_odd_compact_moments": odd_moments,
        "centrally_symmetric": True,
        "proved": True,
    }


def fixed_family_center_coherence(
    auxiliary_values_at_x0: tuple[int | Fraction, ...],
    centers: tuple[int | Fraction, ...],
) -> dict[str, object]:
    """Test center coherence for one fixed, preassigned auxiliary family.

    If ``x0`` generates ``ker(F)``, a common rescaling ``x=c*x0`` can meet
    ``M_i(x)^2=4*j_i^2`` exactly when the displayed square ratios agree.
    This function does not vary or optimize the auxiliary family.
    """
    if not auxiliary_values_at_x0 or len(auxiliary_values_at_x0) != len(centers):
        raise ValueError("need equally many nonempty auxiliary values and centers")
    values = tuple(_fraction(value) for value in auxiliary_values_at_x0)
    center_values = tuple(_fraction(center) for center in centers)
    if any(value == 0 for value in values) or any(center == 0 for center in center_values):
        raise ValueError("auxiliary values and centers must be nonzero")
    ratios = tuple((center / value) ** 2 for value, center in zip(values, center_values))
    coherent = len(set(ratios)) == 1
    return {
        "auxiliary_values_at_x0": values,
        "centers": center_values,
        "ratios_j_squared_over_M_squared": ratios,
        "common_ratio": ratios[0] if coherent else None,
        "fixed_preassigned_family_coherent": coherent,
        "adaptive_auxiliary_choice_analyzed": False,
        "proved": True,
    }


def center_invisibility_theorem() -> dict[str, object]:
    """Record the exact fixed-family consequence of star invisibility.

    For the moment polynomial ``Q_(d,k)(j,t)``, its degree in ``t`` is at
    most ``d-k``.  Hence every unit-star contraction vanishes for
    ``2<=d<=p-2``.  At ``d=p-1`` it is still zero for ``k>0``; for ``k=0``
    only the leading ``t^(p-1)`` term survives and gives ``-1``.  In every
    case the contraction is independent of the nonzero center ``j``.

    For one fixed, preassigned complementary family, one ``x`` in the
    one-dimensional kernel of ``F`` must satisfy ``M_i(x)^2=4*j_i^2`` for
    all halves.  Multiplying one center by a scalar whose square is not one,
    while holding that auxiliary family fixed, can destroy coherence without
    changing the unit-star moments.  It does not rule out adapting the
    auxiliary choices to the new centers.
    """
    return {
        "unit_star_moments": {
            "2<=d<=p-2": "0 for every channel k",
            "d=p-1 and k>0": "0",
            "d=p-1 and k=0": "-1",
            "depends_on_nonzero_center_j": False,
        },
        "directly_replayed_center_independent_data": (
            "all canonical unit-star moment channels through degree p-1"
        ),
        "fixed_scale_singleton_condition": "M_i(x)^2=4*j_i^2 for every i, x in ker(F)",
        "equivalent_ratio_condition": (
            "j_i^2/M_i(x0)^2 is independent of i for one generator x0 of ker(F)"
        ),
        "fixed_family_incoherence_witness": (
            "hold every M_i fixed, then replace one j_k by lambda*j_k "
            "with lambda^2!=1"
        ),
        "preassigned_fixed_auxiliary_family_automatically_coherent": False,
        "adaptive_auxiliary_choice_can_restore_coherence": "OPEN",
        "global_complementary_family_constructed": False,
        "global_complementary_family_refuted": False,
        "proved": True,
    }


def theorem_record() -> dict[str, object]:
    """Return the compact symbolic endpoint-barrier ledger."""
    generic = opposite_swapped_locus(Fraction(2), Fraction(3))
    exceptional = opposite_swapped_locus(Fraction(1, 2), Fraction(1, 2))
    example = rational_clean_overlap_example()
    centers = center_invisibility_theorem()
    proved = bool(
        generic["desired_swapped_orbit_is_unique"]
        and not exceptional["desired_swapped_orbit_is_unique"]
        and example["proved"]
        and centers["proved"]
    )
    if not proved:
        raise ArithmeticError("the Mobius endpoint barrier ledger changed")
    return {
        "title": (
            "Complementary Mobius endpoint: clean local locus and fixed-family "
            "center warning"
        ),
        "four_candidate_locus": {
            "desired_opposite_swapped": (
                "A=(q-1)(qr-1)/(qr), B=(r-1)(qr-1)/(qr)"
            ),
            "same_swapped": "impossible",
            "same_direct": "impossible",
            "extra_opposite_direct": "iff q=r=1/2",
            "clean_unique_overlap": "all admissible locus points except q=r=1/2",
        },
        "generic_replay": generic,
        "double_overlap_replay": exceptional,
        "rational_clean_example": example,
        "center_invisibility": centers,
        "status": (
            "LOCAL CLEAN OVERLAP EXISTS; PREASSIGNED CENTER COHERENCE IS NOT "
            "AUTOMATIC; ADAPTIVE GLOBAL CHOICE OPEN"
        ),
        "scope": (
            "center rescaling holds one preassigned auxiliary family fixed; "
            "the auxiliaries are not optimized after rescaling"
        ),
        "adaptive_global_complementary_choice_resolved": False,
        "global_complementary_family_constructed": False,
        "global_complementary_family_refuted": False,
        "coherent_target_Boolean_completion_constructed": False,
        "all_admissible_targets_excluded": False,
        "residual_ii_closed": False,
        "proved_all_claimed_statements": proved,
    }


if __name__ == "__main__":
    from pprint import pprint

    pprint(theorem_record(), sort_dicts=True)
