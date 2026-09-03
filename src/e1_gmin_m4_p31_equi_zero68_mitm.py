#!/usr/bin/env python3
"""Cheap structural replay for the p=31 equianharmonic zero-form certificate.

The exhaustive generator is ``evidence/p31_equianharmonic_zero68_mitm.cpp``.
This module rederives the alignment census and records the disjoint exhaustive
partition totals.  Its conclusion is finite and deliberately does not close
residual (ii) or the global limit problem.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations
from math import comb

from e1_gmin_m4_conic_odd_radon import tangent_conic_target


P = 31
K = 11


def _edge(a: int, b: int) -> tuple[int, int]:
    a %= P
    b %= P
    if a == b:
        raise ValueError("loop is not an atom edge")
    return (a, b) if a < b else (b, a)


def _orbit(edge: tuple[int, int]) -> tuple[tuple[int, int] | None, int]:
    negative = _edge(-edge[0], -edge[1])
    if negative == edge:
        return None, 0
    representative = min(edge, negative)
    return representative, 1 if representative == edge else -1


def _orbit_chain(
    occurrences: tuple[tuple[int, int, int], ...],
) -> dict[tuple[int, int], int]:
    chain: dict[tuple[int, int], int] = defaultdict(int)
    for a, b, coefficient in occurrences:
        representative, sign = _orbit(_edge(a, b))
        if representative is not None:
            chain[representative] += coefficient * sign
    return {edge: value for edge, value in chain.items() if value}


def _alignment_score(
    occurrences: tuple[tuple[int, int, int], ...],
    target: dict[tuple[int, int], int],
) -> int:
    score = 0
    for a, b, coefficient in occurrences:
        representative, sign = _orbit(_edge(a, b))
        if representative is not None:
            score += coefficient * sign * target.get(representative, 0)
    return score


def _ae_occurrences(triple: tuple[int, int, int]):
    a, b, c = triple
    return ((a, b, 1), (a, c, 1), (b, c, 1))


def _compact_occurrences(
    triple: tuple[int, int, int], distinguished: int
):
    a, b = (value for value in triple if value != distinguished)
    return (
        (a, b, 1),
        (a, distinguished, -1),
        (b, distinguished, -1),
    )


def p31_equianharmonic_alignment_census() -> dict[str, object]:
    """Recompute the complete 4,495/13,485 atom alignment census."""
    target = tangent_conic_target(P, K)
    ae_scores: Counter[int] = Counter()
    compact_scores: Counter[int] = Counter()
    maximal_ae_edges = []
    maximal_compact_chains = []

    for triple in combinations(range(P), 3):
        ae_occurrences = _ae_occurrences(triple)
        ae_score = _alignment_score(ae_occurrences, target)
        ae_scores[ae_score] += 1
        if ae_score == 3:
            chain = _orbit_chain(ae_occurrences)
            maximal_ae_edges.append(frozenset(chain))

        for distinguished in triple:
            occurrences = _compact_occurrences(triple, distinguished)
            compact_score = _alignment_score(occurrences, target)
            compact_scores[compact_score] += 1
            if compact_score == 2:
                maximal_compact_chains.append(_orbit_chain(occurrences))

    target_edges_covered = set().union(*maximal_ae_edges)
    maximal_compact_off = []
    maximal_compact_antipodal = 0
    for chain in maximal_compact_chains:
        off = {edge: value for edge, value in chain.items() if edge not in target}
        if off:
            maximal_compact_off.append(next(iter(off)))
        else:
            maximal_compact_antipodal += 1

    return {
        "p": P,
        "b": 7,
        "k": K,
        "target_support": len(target),
        "ae_atom_count": sum(ae_scores.values()),
        "compact_atom_count": sum(compact_scores.values()),
        "ae_score_counts": dict(sorted(ae_scores.items())),
        "compact_score_counts": dict(sorted(compact_scores.items())),
        "maximum_ae_score": max(ae_scores),
        "maximum_compact_score": max(compact_scores),
        "maximal_ae_cycle_count": len(maximal_ae_edges),
        "maximal_ae_cycles_are_edge_disjoint": sum(
            map(len, maximal_ae_edges)
        ) == len(target_edges_covered),
        "maximal_ae_target_edges_covered": len(target_edges_covered),
        "broken_target_edges": len(target) - len(target_edges_covered),
        "maximal_compact_count": len(maximal_compact_chains),
        "maximal_compact_antipodal_off_edge_count": (
            maximal_compact_antipodal
        ),
        "maximal_compact_supported_off_orbit_count": len(
            set(maximal_compact_off)
        ),
    }


def p31_equianharmonic_zero68_mitm_certificate() -> dict[str, object]:
    """Return the exact exhaustive totals and the resulting finite theorem."""
    census = p31_equianharmonic_alignment_census()
    maximum_total_score = 6 * census["maximum_ae_score"] + 7 * census[
        "maximum_compact_score"
    ]
    total_deficit = maximum_total_score - census["target_support"]

    d3_candidates = census["ae_score_counts"][0] + census[
        "compact_score_counts"
    ][-1]
    d2_candidates = census["ae_score_counts"][1] + census[
        "compact_score_counts"
    ][0]
    d1_candidates = census["ae_score_counts"][2] + census[
        "compact_score_counts"
    ][1]
    compact_d1 = census["compact_score_counts"][1]
    d1_with_ae = 1 + compact_d1 + comb(compact_d1 + 1, 2)
    d1_all_compact = comb(compact_d1 + 2, 3)

    partitions = {
        "deficit_3": {
            "exceptional_candidates": d3_candidates,
            "completions": 13_528_344,
            "edge_hits": 60,
        },
        "deficit_2_plus_1": {
            "exceptional_pairs": d2_candidates * d1_candidates,
            "off_compatible_pairs": 79_918,
            "completions": 87_840_508,
            "edge_hits": 2_160,
        },
        "three_deficit_1_with_ae": {
            "exceptional_multisets": d1_with_ae,
            "off_compatible_multisets": 24_828,
            "completions": 20_465_801,
            "edge_hits": 392,
        },
        "three_compact_deficit_1": {
            "all_unordered_multisets": d1_all_compact,
            "unsupported_projection_pair_keys": 1_031_232,
            "unsupported_projection_3sum_hits": 2_027_542,
            "supported_off_completable_multisets": 1_089_526,
            "completions": 108_480_057,
            "edge_hits": 14_464,
        },
    }
    completions = sum(row["completions"] for row in partitions.values())
    edge_hits = sum(row["edge_hits"] for row in partitions.values())

    return {
        "scope": "p=31,b=7,k=11 constant tangent-conic fiber only",
        "alignment_census": census,
        "maximum_total_score": maximum_total_score,
        "required_target_score": census["target_support"],
        "total_deficit": total_deficit,
        "deficit_partitions": ((3,), (2, 1), (1, 1, 1)),
        "partitions": partitions,
        "maximal_completions_tested": completions,
        "edge_hits": edge_hits,
        "zero_degree_six_and_eight_hits": 0,
        "p31_b7_k11_constant_conic_zero68_fiber": "UNSAT",
        "proved": completions == 230_314_710 and edge_hits == 17_076,
        "residual_ii_closed": False,
        "L_status": "OPEN",
    }
