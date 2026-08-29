#!/usr/bin/env python3
"""Prop. 15.689 -- p=19 low-slack conic reduction.

Proposition 15.688's complete residue-zero census has 143 phase-labelled
profiles.  This proposition excludes every profile of pair slack at most
12, leaving exactly 14 profiles with slack 16 through 32.

The external finite-geometric input is the complete-arc spectrum of
``PG(2,19)``: complete arcs have sizes 10,11,12,13,14,20 and no others.
Consequently every arc of size at least 15 is contained in the unique
20-arc, a nondegenerate conic.

For slack zero the boundary is a 16-arc.  Three undetermined infinity
points give two overlapping 18-arc extensions whose conics coincide and
would contain three collinear points.  With one or two undetermined
directions, the containing conic is respectively tangent or secant to the
line at infinity.  Every other direction then has at least six affine
secants, hence odd-fibre count at most four; the exact profiles violate
this bound.

Slack four repairs to a 15-arc by one deletion.  Its deleted point is off
the containing conic and lies on at least four retained conic secants, so
the pair slack is at least 16, a contradiction.

At slack eight or twelve, repair deletes at most two or three points.  The
profiles have at least two undetermined directions; adjoining those points
at infinity produces an arc of size at least 16 or 15 and hence a conic.
If ``j`` deleted points are off that conic, retained-secant counting gives
pair slack at least ``4*j*(5-j)``, which is 16,24,24 for ``j=1,2,3``.
"""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from e1_gmin_m4_prop15688 import p19_residue_zero_profiles


ROOT = Path(__file__).resolve().parents[1]
P = 19
BOUNDARY_SIZE = 16
EXCLUDED_SLACKS = (0, 4, 8, 12)


def p19_complete_arc_spectrum() -> dict[str, object]:
    """The published complete-arc spectrum used by the reduction."""
    return {
        "external_dependency": True,
        "source": (
            "G. Faina, S. Marcugini, A. Milani, and F. Pambianco, The "
            "spectrum of values k for complete k-arcs in PG(2,q) for "
            "q<=23, Ars Combinatoria 47 (1997), 3-11"
        ),
        "independent_full_classification": (
            "H. Sticker, Classification of Arcs in Desarguesian Projective "
            "Planes, PhD thesis, Ghent University, Table for PG(2,19)"
        ),
        "complete_arc_sizes": [10, 11, 12, 13, 14, 20],
        "no_complete_arc_sizes": [15, 16, 17, 18, 19],
        "unique_size_twenty_arc": "nondegenerate conic",
        "consequence": "every arc of size at least 15 is conic-contained",
        "proved_conditional_on_external_classification": True,
    }


def p19_low_slack_profile_ledger() -> dict[str, object]:
    """Extract and audit the 129 profiles handled geometrically."""
    census = p19_residue_zero_profiles()
    profiles = list(census["profiles"])
    counts = Counter(
        int(row["pair_slack"])
        for row in profiles
        if int(row["pair_slack"]) in EXCLUDED_SLACKS
    )
    expected = {0: 54, 4: 37, 8: 25, 12: 13}
    if dict(sorted(counts.items())) != expected:
        raise ArithmeticError("p=19 low-slack block changed")

    arc_rows = [row for row in profiles if int(row["pair_slack"]) == 0]
    small_t_arc_rows = [
        row for row in arc_rows if int(row["undetermined_directions"]) <= 2
    ]
    small_t_high_b = []
    for row in small_t_arc_rows:
        high = max(
            int(b)
            for b, count in row["global_b_profile"].items()
            if int(count) and int(b) != BOUNDARY_SIZE
        )
        small_t_high_b.append(high)
    if len(small_t_arc_rows) != 29 or min(small_t_high_b) <= 4:
        raise ArithmeticError("p=19 arc profile obstruction changed")

    minimum_t = {
        slack: min(
            int(row["undetermined_directions"])
            for row in profiles
            if int(row["pair_slack"]) == slack
        )
        for slack in (8, 12)
    }
    if minimum_t != {8: 2, 12: 2}:
        raise ArithmeticError("p=19 repaired profiles lost infinity points")
    return {
        "profile_counts_by_slack": expected,
        "excluded_profile_count": sum(expected.values()),
        "slack_zero_profiles_with_at_least_three_undetermined": (
            len(arc_rows) - len(small_t_arc_rows)
        ),
        "slack_zero_profiles_with_one_or_two_undetermined": len(
            small_t_arc_rows
        ),
        "small_t_arc_minimum_high_nonundetermined_b": min(small_t_high_b),
        "minimum_undetermined_directions": minimum_t,
        "proved": True,
    }


def p19_low_slack_geometric_exclusion() -> dict[str, object]:
    """Proposition 15.689."""
    classification = p19_complete_arc_spectrum()
    ledger = p19_low_slack_profile_ledger()
    repair_rows = {
        4: {
            "repair_deletion_bound": 1,
            "repaired_arc_minimum_size": 15,
            "off_conic_points": [1],
            "pair_slack_floors": [16],
        },
        8: {
            "repair_deletion_bound": 2,
            "adjoined_undetermined_infinity_points": 2,
            "extended_arc_minimum_size": 16,
            "off_conic_points": [1, 2],
            "pair_slack_floors": [16, 24],
        },
        12: {
            "repair_deletion_bound": 3,
            "adjoined_undetermined_infinity_points": 2,
            "extended_arc_minimum_size": 15,
            "off_conic_points": [1, 2, 3],
            "pair_slack_floors": [16, 24, 24],
        },
    }
    for slack, row in repair_rows.items():
        if min(row["pair_slack_floors"]) <= slack:
            raise ArithmeticError("p=19 conic-core floor stopped contradicting slack")

    census = p19_residue_zero_profiles()
    remaining_histogram = dict(census["pair_slack_histogram"])
    for slack in EXCLUDED_SLACKS:
        remaining_histogram.pop(slack)
    remaining_count = sum(remaining_histogram.values())
    if remaining_count != 14 or remaining_histogram != {
        16: 7,
        20: 4,
        24: 1,
        28: 1,
        32: 1,
    }:
        raise ArithmeticError("p=19 post-conic remainder changed")
    return {
        "proposition": "15.689",
        "p": P,
        "boundary_size": BOUNDARY_SIZE,
        "profile_count_before": int(census["phase_labelled_profile_count"]),
        "excluded_pair_slacks": list(EXCLUDED_SLACKS),
        "profile_count_excluded": int(ledger["excluded_profile_count"]),
        "profile_count_after": remaining_count,
        "remaining_pair_slack_histogram": remaining_histogram,
        "slack_zero_exclusion": {
            "at_least_three_undetermined": (
                "two overlapping 18-arc extensions force one conic through "
                "three collinear infinity points"
            ),
            "one_or_two_undetermined": (
                "the containing tangent/secant conic forces every other "
                "direction to have b<=4; every exact profile has b>=6"
            ),
        },
        "repair_conic_core_rows": repair_rows,
        "retained_secant_inequality": (
            "with two adjoined infinity points and j off-conic boundary "
            "points, slack >= 4*j*(5-j)"
        ),
        "classification": classification,
        "profile_ledger": ledger,
        "p19_second_all_finite_endpoint_closed": False,
        "remaining_same_boundary_primes": [17, 19, 23],
        "top_level_gates_changed": False,
        "proved_conditional_on_external_classification": True,
    }


def main() -> None:
    theorem = p19_low_slack_geometric_exclusion()
    target = ROOT / "evidence" / "e1_gmin_m4_prop15689.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.689: p=19 residue-zero profiles "
        f"{theorem['profile_count_before']} -> {theorem['profile_count_after']}"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
