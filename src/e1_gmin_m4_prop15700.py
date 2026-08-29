#!/usr/bin/env python3
"""Prop. 15.700 -- reduce the p=17 second boundary, slack zero to two rows.

The second even all-finite size at p=17 is s=16.  Exact same-type
quotient arithmetic and Proposition 15.688's sharp integral lift floor leave
phase-zero residues 0,7,8 and phase-one residues 0,8.  Completion-bounded
profile enumeration, the pair budget, and slack divisibility give 1,575
phase-labelled rows, including 247 rows of pair slack zero.

Slack zero makes the sixteen affine boundary points a projective 16-arc.
Sticker's exhaustive PG(2,17) classification has one 16-arc class, represented
by a conic with two points deleted.  Fixing one conic and enumerating every
line at infinity and every eligible deleted pair gives 21,267 affine cases.
Their exact Paley-phase directional census has 53 labelled profiles (including
the nonsquare phase swap).  Only two of the 247 arithmetic profiles occur;
both are tangent-at-infinity conic-minus-two cases:

    phase 0 {0:1,2:7,16:1}, phase 1 {2:9};
    phase 0 {0:1,2:8},      phase 1 {2:8,16:1}.

Thus 245 profiles are excluded and the exact p=17 remainder drops from 1,575
to 1,330.  This is a strict reduction, not endpoint closure: the two conic
profiles and every positive-slack profile remain.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from pathlib import Path

from e1_gmin_m4_prop15632 import field_direction_data, projective_directions
from e1_gmin_m4_prop15669 import full_symbolic_floor
from e1_gmin_m4_prop15678 import p17_arc_classification_ledger
from e1_gmin_m4_prop15688 import sharp_integral_quadratic_lift_floor


ROOT = Path(__file__).resolve().parents[1]
P = 17
M = 9
S = 16
PERIOD = 18
PAIR_DEFICIT_BUDGET = S * (S - 1)
Point = tuple[int, int, int]
ProfileKey = tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]


def _histogram(values: tuple[int, ...]) -> dict[int, int]:
    return dict(sorted(Counter(values).items()))


@lru_cache(maxsize=None)
def _profile_rows(
    phase: int, u: int, deficit_cap: int
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Enumerate all exact p=17 profiles within one type deficit cap."""
    if phase not in (0, 1) or not 0 <= u < M:
        raise ValueError("phase must be 0/1 and 0<=u<9")
    target = M - u
    options: list[tuple[int, int, int]] = []
    for b in range(0, S + 1, 2):
        floor = full_symbolic_floor(P, b, phase)
        for quotient in range(target + 1):
            excess = 2 * u + PERIOD * quotient - floor
            if excess >= 0 and excess != 2:
                options.append((quotient, S - b, b))

    infinity = deficit_cap + S * M + 1
    completion = [[infinity] * (target + 1) for _ in range(M + 1)]
    completion[0][0] = 0
    for count in range(1, M + 1):
        for quotient_sum in range(target + 1):
            completion[count][quotient_sum] = min(
                (
                    added + completion[count - 1][quotient_sum - quotient]
                    for quotient, added, _b in options
                    if quotient <= quotient_sum
                ),
                default=infinity,
            )

    states: set[tuple[int, int, tuple[int, ...]]] = {(0, 0, ())}
    for count in range(M):
        next_states: set[tuple[int, int, tuple[int, ...]]] = set()
        for used, deficit, profile in states:
            for quotient, added, b in options:
                new_used = used + quotient
                new_deficit = deficit + added
                remaining_count = M - count - 1
                remaining_sum = target - new_used
                if (
                    new_used <= target
                    and new_deficit <= deficit_cap
                    and new_deficit
                    + completion[remaining_count][remaining_sum]
                    <= deficit_cap
                ):
                    next_states.add(
                        (
                            new_used,
                            new_deficit,
                            tuple(sorted(profile + (b,))),
                        )
                    )
        states = next_states
    return tuple(
        sorted(
            (deficit, profile)
            for used, deficit, profile in states
            if used == target
        )
    )


def residue_and_lift_ledger() -> dict[str, object]:
    """Record all relaxed minima and the sharp-lift residue reduction."""
    minima: dict[str, dict[int, int]] = {"0": {}, "1": {}}
    row_counts: dict[str, dict[int, int]] = {"0": {}, "1": {}}
    for phase in (0, 1):
        for u in range(M):
            rows = _profile_rows(phase, u, PAIR_DEFICIT_BUDGET)
            if rows:
                minima[str(phase)][u] = int(rows[0][0])
                row_counts[str(phase)][u] = len(rows)
    expected_minima = {
        "0": {0: 78, 2: 32, 3: 48, 4: 64, 5: 80, 6: 96, 7: 112, 8: 0},
        "1": {0: 0, 8: 112},
    }
    expected_counts = {
        "0": {0: 375, 2: 540, 3: 300, 4: 144, 5: 100, 6: 28, 7: 15, 8: 73},
        "1": {0: 1, 8: 9},
    }
    if minima != expected_minima or row_counts != expected_counts:
        raise ArithmeticError("p=17 residue ledger changed")

    lift = sharp_integral_quadratic_lift_floor(P)
    excluded_rows = []
    least_positive_floor = min(
        full_symbolic_floor(P, b, 0) for b in range(2, S + 1, 2)
    )
    for u in range(2, 7):
        scaled_mean = 2 * u
        excluded_rows.append(
            {
                "u": u,
                "quotient_sum": M - u,
                "quotient_zero_forced": M - u < M,
                "zero_quotient_scaled_mean": scaled_mean,
                "least_positive_b_floor": least_positive_floor,
                "therefore_b_zero": scaled_mean < least_positive_floor,
                "sharp_nonzero_lift_floor": int(lift["sharp_scaled_floor"]),
                "excluded": int(lift["sharp_scaled_floor"]) > scaled_mean,
            }
        )
    if least_positive_floor != 16 or not all(row["excluded"] for row in excluded_rows):
        raise ArithmeticError("p=17 sharp-lift reduction changed")
    return {
        "p": P,
        "s": S,
        "phase_minimum_deficits": minima,
        "phase_profile_row_counts_at_full_cap": row_counts,
        "phase_zero_positive_residues_excluded": excluded_rows,
        "remaining_phase_zero_residues": [0, 7, 8],
        "remaining_phase_one_residues": [0, 8],
        "sharp_lift": lift,
        "proved": True,
    }


@lru_cache(maxsize=1)
def p17_second_boundary_profile_census() -> dict[str, object]:
    """Enumerate the complete post-lift pair-budget profile ledger."""
    ledger = residue_and_lift_ledger()
    minima = ledger["phase_minimum_deficits"]
    candidates = []
    for u0 in ledger["remaining_phase_zero_residues"]:
        for u1 in ledger["remaining_phase_one_residues"]:
            phase_zero = _profile_rows(
                0, int(u0), PAIR_DEFICIT_BUDGET - int(minima["1"][u1])
            )
            phase_one = _profile_rows(
                1, int(u1), PAIR_DEFICIT_BUDGET - int(minima["0"][u0])
            )
            for deficit_zero, profile_zero in phase_zero:
                for deficit_one, profile_one in phase_one:
                    slack = PAIR_DEFICIT_BUDGET - deficit_zero - deficit_one
                    if slack < 0 or slack % 4:
                        continue
                    candidates.append(
                        {
                            "u0": int(u0),
                            "u1": int(u1),
                            "phase_deficits": {
                                "0": deficit_zero,
                                "1": deficit_one,
                            },
                            "phase_profiles_b": {
                                "0": _histogram(profile_zero),
                                "1": _histogram(profile_one),
                            },
                            "pair_slack": slack,
                        }
                    )
    candidates.sort(
        key=lambda row: (
            int(row["pair_slack"]),
            int(row["u0"]),
            int(row["u1"]),
            tuple(row["phase_profiles_b"]["0"].items()),
            tuple(row["phase_profiles_b"]["1"].items()),
        )
    )
    slack_histogram = dict(
        sorted(Counter(int(row["pair_slack"]) for row in candidates).items())
    )
    residue_histogram = dict(
        sorted(Counter((int(row["u0"]), int(row["u1"])) for row in candidates).items())
    )
    slack_zero = [row for row in candidates if int(row["pair_slack"]) == 0]
    slack_zero_residues = dict(
        sorted(Counter((int(row["u0"]), int(row["u1"])) for row in slack_zero).items())
    )
    if (
        len(candidates) != 1575
        or len(slack_zero) != 247
        or residue_histogram
        != {(0, 0): 181, (0, 8): 1062, (7, 0): 9, (7, 8): 9, (8, 0): 37, (8, 8): 277}
        or slack_zero_residues != {(0, 8): 234, (7, 8): 4, (8, 8): 9}
    ):
        raise ArithmeticError("p=17 second-boundary profile census changed")
    canonical = json.dumps(candidates, sort_keys=True, separators=(",", ":"))
    return {
        "p": P,
        "boundary_size": S,
        "pair_deficit_budget": PAIR_DEFICIT_BUDGET,
        "phase_labelled_profile_count": len(candidates),
        "pair_slack_histogram": slack_histogram,
        "residue_pair_histogram": residue_histogram,
        "slack_zero_profile_count": len(slack_zero),
        "slack_zero_residue_pair_histogram": slack_zero_residues,
        "canonical_profile_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
        "profiles": candidates,
        "residue_and_lift_ledger": ledger,
        "proved": True,
    }


def _dot(left: Point, right: Point) -> int:
    return sum(a * b for a, b in zip(left, right)) % P


def _determinant(rows: tuple[Point, Point, Point]) -> int:
    a, b, c = rows
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    ) % P


def _projective_points_or_lines() -> tuple[Point, ...]:
    return tuple(
        [(1, y, z) for y in range(P) for z in range(P)]
        + [(0, 1, z) for z in range(P)]
        + [(0, 0, 1)]
    )


def _chart(line: Point) -> tuple[Point, Point, Point]:
    basis: tuple[Point, ...] = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    for first in basis:
        for second in basis:
            rows = (first, second, line)
            if first != second and _determinant(rows):
                return rows
    raise ArithmeticError("failed to construct affine chart")


def _affine_point(matrix: tuple[Point, Point, Point], point: Point) -> tuple[int, int]:
    image = tuple(_dot(row, point) for row in matrix)
    inverse = pow(image[2], -1, P)
    return image[0] * inverse % P, image[1] * inverse % P


def _profile_key(points: tuple[tuple[int, int], ...]) -> ProfileKey:
    phases = {0: Counter(), 1: Counter()}
    for direction in projective_directions(P):
        eps, labels = field_direction_data(P, direction)
        phase = 0 if eps == 1 else 1
        counts = [0] * P
        for x, y in points:
            counts[labels[y * P + x]] += 1
        phases[phase][sum(value & 1 for value in counts)] += 1
    return tuple(sorted(phases[0].items())), tuple(sorted(phases[1].items()))


@lru_cache(maxsize=1)
def conic_minus_two_affine_profile_census() -> dict[str, object]:
    """Enumerate every affine chart/deleted pair of the unique 16-arc."""
    conic: tuple[Point, ...] = tuple(
        [(t * t % P, t, 1) for t in range(P)] + [(1, 0, 0)]
    )
    lines = _projective_points_or_lines()
    profiles: set[ProfileKey] = set()
    raw_by_intersection = Counter()
    examples: dict[ProfileKey, dict[str, object]] = {}
    raw_case_count = 0
    for line in lines:
        on_line = {
            index for index, point in enumerate(conic) if _dot(line, point) == 0
        }
        matrix = _chart(line)
        for removed in combinations(range(len(conic)), 2):
            if not on_line <= set(removed):
                continue
            points = tuple(
                _affine_point(matrix, point)
                for index, point in enumerate(conic)
                if index not in removed
            )
            if len(points) != S or len(set(points)) != S:
                raise ArithmeticError("conic-minus-two affine image changed")
            key = _profile_key(points)
            swapped = key[1], key[0]
            profiles.update((key, swapped))
            raw_case_count += 1
            raw_by_intersection[len(on_line)] += 1
            example = {
                "line_at_infinity": list(line),
                "removed_conic_indices": list(removed),
                "line_conic_intersection_size": len(on_line),
            }
            examples.setdefault(key, example)
            examples.setdefault(swapped, example)
    expected_raw = {0: 20808, 1: 306, 2: 153}
    if (
        len(lines) != 307
        or raw_case_count != 21267
        or dict(raw_by_intersection) != expected_raw
        or len(profiles) != 53
    ):
        raise ArithmeticError("p=17 conic profile census changed")
    return {
        "canonical_conic": [list(point) for point in conic],
        "projective_line_count": len(lines),
        "raw_affine_case_count": raw_case_count,
        "raw_case_count_by_line_intersection": dict(sorted(raw_by_intersection.items())),
        "phase_labelled_profile_count_including_swap": len(profiles),
        "profiles": profiles,
        "examples": examples,
        "proved": True,
    }


def _key_from_row(row: dict[str, object]) -> ProfileKey:
    profiles = row["phase_profiles_b"]
    return (
        tuple(sorted((int(b), int(count)) for b, count in profiles["0"].items())),
        tuple(sorted((int(b), int(count)) for b, count in profiles["1"].items())),
    )


def p17_slack_zero_conic_reduction() -> dict[str, object]:
    """Proposition 15.700."""
    arithmetic = p17_second_boundary_profile_census()
    geometry = conic_minus_two_affine_profile_census()
    classification = p17_arc_classification_ledger()
    if int(classification["pgl_classes_in_pg2_17"]["16"]) != 1:
        raise ArithmeticError("p=17 16-arc classification changed")
    slack_zero = [
        row for row in arithmetic["profiles"] if int(row["pair_slack"]) == 0
    ]
    geometric_keys = geometry["profiles"]
    survivors = [row for row in slack_zero if _key_from_row(row) in geometric_keys]
    expected = [
        ({0: 1, 2: 7, 16: 1}, {2: 9}),
        ({0: 1, 2: 8}, {2: 8, 16: 1}),
    ]
    observed = [
        (
            dict(row["phase_profiles_b"]["0"]),
            dict(row["phase_profiles_b"]["1"]),
        )
        for row in survivors
    ]
    if observed != expected or any((row["u0"], row["u1"]) != (0, 8) for row in survivors):
        raise ArithmeticError("p=17 conic/arithmetic intersection changed")
    survivor_records = []
    for row in survivors:
        key = _key_from_row(row)
        example = geometry["examples"][key]
        if int(example["line_conic_intersection_size"]) != 1:
            raise ArithmeticError("surviving p=17 conic profile is not tangent")
        survivor_records.append({**row, "conic_example": example})
    before = int(arithmetic["phase_labelled_profile_count"])
    excluded = len(slack_zero) - len(survivors)
    after = before - excluded
    if before != 1575 or excluded != 245 or after != 1330:
        raise ArithmeticError("p=17 reduction accounting changed")
    return {
        "proposition": "15.700",
        "p": P,
        "boundary_size": S,
        "profile_count_before": before,
        "slack_zero_profile_count_before": len(slack_zero),
        "slack_zero_profile_count_after": len(survivors),
        "profiles_excluded_here": excluded,
        "profile_count_after": after,
        "surviving_slack_zero_profiles": survivor_records,
        "survivor_geometry": "conic minus its tangent point at infinity and one affine point",
        "all_positive_slack_profiles_remain": True,
        "p17_second_all_finite_endpoint_closed": False,
        "profile_census": arithmetic,
        "conic_profile_census": {
            key: value
            for key, value in geometry.items()
            if key not in ("profiles", "examples")
        },
        "external_classification": classification,
        "proved_conditional_on_external_classification": True,
    }


def _jsonable(value):
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list, set)):
        return [_jsonable(item) for item in value]
    return value


def main() -> None:
    theorem = _jsonable(p17_slack_zero_conic_reduction())
    target = ROOT / "evidence" / "e1_gmin_m4_prop15700.json"
    target.write_text(json.dumps(theorem, indent=2, sort_keys=True) + "\n")
    print(
        "Prop. 15.700: p=17 second-boundary profiles "
        f"{theorem['profile_count_before']} -> {theorem['profile_count_after']}; "
        "slack zero 247 -> 2"
    )
    print(f"wrote {target}")


if __name__ == "__main__":
    main()
