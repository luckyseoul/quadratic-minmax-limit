#!/usr/bin/env python3
"""Exact cross-ratio reduction for the open ``p=31``, weight-eight case.

This module does not run the proposed large meet-in-the-middle search.  It
proves the finite symmetry reduction which cuts that search from all possible
fourth silent directions to four generic anharmonic-orbit representatives,
and records the exact coverage conditions for a later ``4+4`` certificate.
"""

from __future__ import annotations

from math import comb


P = 31
H = 15
DIRECTION_COUNT = 32
POINT_CLASS_COUNT = 480
FIXED_SILENT_TRIPLE = (0, 1, 31)


def _inverse(value: int) -> int:
    value %= P
    if value == 0:
        raise ValueError("zero has no inverse")
    return pow(value, -1, P)


def anharmonic_orbit(value: int) -> tuple[int, ...]:
    """Return the orbit of a fourth point under the stabilizer of ``{0,1,inf}``.

    The stabilizer induces the six standard cross-ratio transforms.  The
    input is therefore a finite direction different from zero and one.
    """
    value %= P
    if value in (0, 1):
        raise ValueError("the fourth direction must differ from 0, 1, infinity")
    orbit = {
        value,
        (1 - value) % P,
        _inverse(value),
        _inverse(1 - value),
        value * _inverse(value - 1) % P,
        (value - 1) * _inverse(value) % P,
    }
    return tuple(sorted(orbit))


def p31_anharmonic_partition() -> tuple[tuple[int, ...], ...]:
    """Partition the 29 possible fourth directions into anharmonic orbits."""
    remaining = set(range(2, P))
    parts: list[tuple[int, ...]] = []
    while remaining:
        orbit = anharmonic_orbit(min(remaining))
        if not set(orbit) <= remaining:
            raise ArithmeticError("anharmonic orbits failed to partition")
        parts.append(orbit)
        remaining.difference_update(orbit)
    return tuple(parts)


def fourth_silent_direction_reduction() -> dict[str, object]:
    """Prove that four generic fourth-direction cases cover every counterexample.

    A weight-eight counterexample has at least nine silent directions.  After
    moving any three of them to ``{0,1,infinity}``, six silent directions
    remain.  The harmonic and equianharmonic orbits together contain only
    five directions, so one remaining direction belongs to a generic orbit.
    """
    partition = p31_anharmonic_partition()
    harmonic = anharmonic_orbit(2)
    equianharmonic = anharmonic_orbit(6)
    exceptional = tuple(sorted(set(harmonic) | set(equianharmonic)))
    generic = tuple(orbit for orbit in partition if orbit not in (harmonic, equianharmonic))
    generic_representatives = tuple(orbit[0] for orbit in generic)
    extra_silent_count = 9 - len(FIXED_SILENT_TRIPLE)
    proved = bool(
        tuple(map(len, partition)) == (3, 6, 6, 6, 2, 6)
        and harmonic == (2, 16, 30)
        and equianharmonic == (6, 26)
        and len(exceptional) == 5
        and extra_silent_count == 6
        and generic_representatives == (3, 4, 5, 12)
        and set().union(*map(set, partition)) == set(range(2, P))
    )
    if not proved:
        raise ArithmeticError("the p=31 cross-ratio reduction changed")
    return {
        "p": P,
        "point_weight": 8,
        "counterexample_requires_silent_directions": 9,
        "fixed_silent_triple": list(FIXED_SILENT_TRIPLE),
        "remaining_silent_directions": extra_silent_count,
        "anharmonic_orbits": [list(orbit) for orbit in partition],
        "harmonic_orbit": list(harmonic),
        "equianharmonic_orbit": list(equianharmonic),
        "exceptional_union_size": len(exceptional),
        "generic_orbit_representatives": list(generic_representatives),
        "coverage_reason": (
            "six remaining silent directions cannot all lie in the five-point "
            "union of the harmonic and equianharmonic orbits"
        ),
        "large_search_run": False,
        "proved": proved,
    }


def quadruple_mitm_certificate_plan() -> dict[str, object]:
    """Record the exact ``4+4`` coverage and storage counts, without running it."""
    quadruples = comb(POINT_CLASS_COUNT, 4)
    unordered_partitions = comb(8, 4) // 2
    signature_bits = 4 * H
    record_bytes_aligned = 16
    aligned_gib = quadruples * record_bytes_aligned / 2**30
    parallel_array_bytes = quadruples * (8 + 4)
    proved = bool(
        quadruples == 2_184_297_480
        and unordered_partitions == 35
        and signature_bits == 60
        and quadruples.bit_length() == 32
        and aligned_gib < 33
    )
    if not proved:
        raise ArithmeticError("the p=31 weight-eight MITM counts changed")
    return {
        "p": P,
        "point_weight": 8,
        "fourth_direction_cases": [3, 4, 5, 12],
        "signature": (
            "the 60 parity bits of the fifteen nonzero squared-projection "
            "blocks in directions 0,1,infinity,lambda"
        ),
        "quadruples": quadruples,
        "quadruple_index_bits": quadruples.bit_length(),
        "unordered_4_plus_4_partitions_per_support": unordered_partitions,
        "coverage": (
            "a weight-eight support is silent in the four fixed directions "
            "iff the two halves in every 4+4 partition have equal signatures"
        ),
        "zero_fibre_reason": (
            "even total weight and even nonzero fibres force the omitted zero "
            "fibre parity to be even"
        ),
        "aligned_16_byte_record_gib": aligned_gib,
        "parallel_key_and_rank_array_gib": parallel_array_bytes / 2**30,
        "required_final_replay": (
            "recompute all 32 direction syndromes and reject unless at least "
            "nine direction groups are silent"
        ),
        "large_search_run": False,
        "proved": proved,
    }


def theorem_record() -> dict[str, object]:
    reduction = fourth_silent_direction_reduction()
    plan = quadruple_mitm_certificate_plan()
    return {
        "p": P,
        "proved": {
            "four_generic_cross_ratio_cases_are_complete": reduction["proved"],
            "four_plus_four_signature_coverage": plan["proved"],
            "weight_eight_counterexample_excluded": False,
            "group_support_lemma": False,
            "row_code_minimum_distance": False,
            "residual_ii_closed": False,
        },
        "large_search_run": False,
    }
