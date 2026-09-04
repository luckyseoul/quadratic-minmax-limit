#!/usr/bin/env python3
"""Exact parity/support split for the p=31 top ``j=1, f=3, d=0`` ledger.

This module deliberately stops before physical synchronization.  It combines
the top parallel parity word, three fixed antipodal edges, two cancellation
units, and the adaptive kernel-selector product.  The result classifies every
possible fixed-direction parity support and the forced auxiliary
multiplicities.  It does not assert that any surviving case has compatible
Mobius scales, centres, or compact atoms.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations, combinations_with_replacement
from typing import Iterable, Sequence

from e1_gmin_m4_inversion_antisymmetric_radon import projective_functionals
from e1_gmin_m4_mobius_half_symmetric import paley_direction_sign


P = 31
M = (P + 1) // 2
FIXED_EDGE_COUNT = 3
CANCELLATION_UNITS = 2

# One labelled representative of the top quota assignment.  Its multisets
# are hard 14^14 15^2 and opposite 15^3 16^13.
TARGET_PROFILE = (
    15, 14, 14, 15, 16, 15, 16, 14,
    16, 14, 14, 16, 16, 15, 15, 14,
    14, 16, 16, 16, 16, 14, 14, 16,
    14, 16, 16, 16, 14, 14, 14, 14,
)

DIRECTIONS = tuple(projective_functionals(P))
DIRECTION_SIGNS = tuple(paley_direction_sign(P, row) for row in DIRECTIONS)
HARD = frozenset(i for i, sign in enumerate(DIRECTION_SIGNS) if sign == 1)
OPPOSITE = frozenset(range(P + 1)) - HARD

# The established parallel-slice parity word is v=P_D+1_hard (mod 2).
V_SUPPORT = frozenset(
    i
    for i, value in enumerate(TARGET_PROFILE)
    if (value + int(i in HARD)) % 2
)
V_HARD = V_SUPPORT & HARD
V_OPPOSITE = V_SUPPORT & OPPOSITE


def fixed_parity_support(fixed_counts: Sequence[int]) -> frozenset[int]:
    """Return ``b_D=f_D mod 2`` for an exact three-fixed-edge allocation."""
    if len(fixed_counts) != P + 1:
        raise ValueError("fixed_counts must have 32 entries")
    if any(not isinstance(value, int) or value < 0 for value in fixed_counts):
        raise ValueError("fixed counts must be nonnegative integers")
    if sum(fixed_counts) != FIXED_EDGE_COUNT:
        raise ValueError("the j=1,f=3 branch has exactly three fixed edges")
    return frozenset(i for i, value in enumerate(fixed_counts) if value % 2)


def classify_support(b_support: Iterable[int]) -> dict[str, object]:
    """Classify one possible fixed-direction parity support.

    Let ``a=v+b`` be the parity vector of the sixteen auxiliary direction
    occurrences.  Necessarily ``wt(a)<=16``.  Weight sixteen forces an SDR;
    weight fourteen forces exactly one direction ``A`` to receive the only
    extra pair of occurrences.  The adaptive-selector product gives

        z = 1 + |b intersect hard|  (mod 2),

    where ``z`` is the parity of origin-containing cancellation units.
    An SDR permits none.  In the weight-fourteen case the unique repeated
    auxiliary group permits at most one, so this parity determines its exact
    value.
    """
    b = frozenset(int(value) for value in b_support)
    if any(value not in range(P + 1) for value in b):
        raise ValueError("a direction index is out of range")
    if len(b) not in (1, 3):
        raise ValueError("three fixed edges give parity support of weight 1 or 3")

    a = V_SUPPORT ^ b
    k = len(V_SUPPORT & b)
    hard_b = len(HARD & b)
    hard_auxiliary_parity = len(HARD & a) % 2
    required_origin_parity = (hard_auxiliary_parity + FIXED_EDGE_COUNT) % 2
    if required_origin_parity != (1 + hard_b) % 2:
        raise ArithmeticError("the selector parity simplification changed")

    record: dict[str, object] = {
        "b_support": tuple(sorted(b)),
        "b_weight": len(b),
        "b_intersection_v": k,
        "b_hard_count": hard_b,
        "auxiliary_parity_support": tuple(sorted(a)),
        "auxiliary_parity_weight": len(a),
        "required_origin_cancellation_parity": required_origin_parity,
        "possible": False,
    }

    if len(a) > M:
        record.update(
            classification="capacity_excluded",
            reason="auxiliary parity support exceeds sixteen occurrences",
        )
        return record
    if len(a) not in (M, M - 2):
        raise ArithmeticError("an unexpected auxiliary parity weight survived")

    if len(a) == M:
        record["auxiliary_multiplicity_form"] = "each support direction once (SDR)"
        record["origin_cancellation_units"] = 0
        if required_origin_parity:
            record.update(
                classification="selector_excluded_sdr",
                reason="the selector product requires an origin correction but the SDR forbids one",
            )
            return record
        record.update(
            possible=True,
            nonorigin_cancellation_units=2,
            correction_signature_possible_weights=(0, 2, 4),
            classification=(
                "singleton_sdr_two_nonorigin"
                if len(b) == 1
                else "triple_k2_sdr_two_nonorigin"
            ),
        )
        return record

    # With fourteen odd directions and sixteen occurrences, n=a+2r and
    # sum(r)=1.  Thus one direction A is uniquely repeated: multiplicity 3
    # if A is in supp(a), and multiplicity 2 otherwise.
    origin_units = required_origin_parity
    record.update(
        possible=True,
        auxiliary_multiplicity_form=(
            "choose one A: multiplicity 3 at A and 1 on the other 13 support directions "
            "if A lies in the support; otherwise multiplicity 2 at A and 1 on all 14 support directions"
        ),
        unique_repeated_auxiliary_direction=True,
        origin_cancellation_units=origin_units,
        nonorigin_cancellation_units=CANCELLATION_UNITS - origin_units,
        correction_signature_possible_weights=(1, 3) if origin_units else (0, 2, 4),
        classification=(
            "triple_k3_repeated_one_origin_one_nonorigin"
            if origin_units
            else "triple_k3_repeated_two_nonorigin"
        ),
    )
    return record


def classify_fixed_counts(fixed_counts: Sequence[int]) -> dict[str, object]:
    """Add the actual fixed multiplicities and raw-profile equation."""
    result = classify_support(fixed_parity_support(fixed_counts))
    result["fixed_counts"] = tuple(fixed_counts)
    result["raw_profile_equation"] = "R_D = P_D + 2*kappa_D - f_D, sum(kappa)=2"
    return result


def required_raw_profile(
    fixed_counts: Sequence[int], cancellation_counts: Sequence[int]
) -> tuple[int, ...]:
    """Return the raw sixteen-half profile forced by fixed/cancellation data."""
    fixed_parity_support(fixed_counts)
    if len(cancellation_counts) != P + 1:
        raise ValueError("cancellation_counts must have 32 entries")
    if any(not isinstance(value, int) or value < 0 for value in cancellation_counts):
        raise ValueError("cancellation counts must be nonnegative integers")
    if sum(cancellation_counts) != CANCELLATION_UNITS:
        raise ValueError("the j=1,f=3 branch has exactly two cancellation units")
    raw = tuple(
        TARGET_PROFILE[i] + 2 * cancellation_counts[i] - fixed_counts[i]
        for i in range(P + 1)
    )
    if sum(raw) != M * (P - 1):
        raise ArithmeticError("the raw profile lost its 480-edge total")
    return raw


def theorem_record() -> dict[str, object]:
    """Exhaust all supports and all 5,984 allocations of three fixed edges."""
    support_classes: Counter[str] = Counter()
    support_examples: dict[str, tuple[int, ...]] = {}
    for weight in (1, 3):
        for support in combinations(range(P + 1), weight):
            row = classify_support(support)
            name = str(row["classification"])
            support_classes[name] += 1
            support_examples.setdefault(name, support)

    allocation_classes: Counter[str] = Counter()
    for directions in combinations_with_replacement(range(P + 1), 3):
        fixed = [0] * (P + 1)
        for direction in directions:
            fixed[direction] += 1
        row = classify_fixed_counts(fixed)
        allocation_classes[str(row["classification"])] += 1

    possible_supports = sum(
        count
        for name, count in support_classes.items()
        if not name.startswith("capacity_") and not name.startswith("selector_")
    )
    possible_allocations = sum(
        count
        for name, count in allocation_classes.items()
        if not name.startswith("capacity_") and not name.startswith("selector_")
    )
    proved = bool(
        len(HARD) == len(OPPOSITE) == M
        and len(V_HARD) == 14
        and len(V_OPPOSITE) == 3
        and len(V_SUPPORT) == 17
        and sum(support_classes.values()) == 32 + 4960
        and support_classes["singleton_sdr_two_nonorigin"] == 14
        and support_classes["triple_k2_sdr_two_nonorigin"] == 734
        and support_classes["triple_k3_repeated_one_origin_one_nonorigin"] == 274
        and support_classes["triple_k3_repeated_two_nonorigin"] == 406
        and sum(allocation_classes.values()) == 5984
        and possible_supports == 1428
        and possible_allocations == 1862
    )
    if not proved:
        raise ArithmeticError("the p31 top j=1,f=3 case split changed")
    return {
        "p": P,
        "m": M,
        "fixed_edges": FIXED_EDGE_COUNT,
        "cancellation_units": CANCELLATION_UNITS,
        "target_profile": TARGET_PROFILE,
        "hard_directions": tuple(sorted(HARD)),
        "opposite_directions": tuple(sorted(OPPOSITE)),
        "v_hard_support": tuple(sorted(V_HARD)),
        "v_opposite_support": tuple(sorted(V_OPPOSITE)),
        "support_class_counts": dict(sorted(support_classes.items())),
        "support_class_examples": dict(sorted(support_examples.items())),
        "possible_support_count": possible_supports,
        "fixed_allocation_class_counts": dict(sorted(allocation_classes.items())),
        "fixed_allocation_count": sum(allocation_classes.values()),
        "possible_fixed_allocation_count": possible_allocations,
        "nonorigin_signature": "e_U+e_V with U!=V and spatial direction D not in {U,V}",
        "origin_signature": "e_A at the repeated auxiliary direction A",
        "scope": (
            "complete parity/support classification conditional on the p31 top localized-Mobius "
            "j=1,f=3,d=0 ledger; scale, centre, collision-locus, atom, and common-graph synchronization remain open"
        ),
        "residual_ii_closed": False,
        "proved": proved,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(theorem_record(), indent=2, sort_keys=True))
