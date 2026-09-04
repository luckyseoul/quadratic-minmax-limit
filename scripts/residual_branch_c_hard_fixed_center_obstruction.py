#!/usr/bin/env python3
"""Exact centre obstruction for one hard-fixed complement-pair witness.

The selected sixteen labelled halves have the exact hard-fixed parallel
profile and the required two-bit boundary correction.  This module checks
all nonzero centre pairs and proves that none of their shared inversion
orbits has the required spatial direction zero.  Since every repeated orbit
among two or more halves appears in a pairwise intersection, the fixed half
family cannot realize its required cancellation for any centre tuple.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from e1_gmin_m4_p31_top_mobius_boundary_parity import (  # noqa: E402
    half_kernel_sigma_formula,
)
from io_atomic import write_json_atomic  # noqa: E402
from scripts import residual_branch_c_auxiliary_transverse_gpu as auxiliary  # noqa: E402


P = 31
FIXED_DIRECTION_INDEX = 0
SOURCE_WITNESS_SHA256 = (
    "9e241aeefee9d91132833b83d0d2e6cf2e4f875b531ddd039551fafb0b5cea91"
)
OPTION_CATALOG_SHA256 = (
    "cec2298fc950ed6a97a60b203756c6caddcb15f239cc78e6ab2de79b6cdea5e0"
)
HALF_CHOICE_ROWS = (
    (0, 22, 12),
    (1, 27, 27),
    (2, 15, 11),
    (3, 28, 7),
    (7, 16, 9),
    (9, 4, 5),
    (10, 7, 20),
    (15, 26, 16),
    (16, 24, 28),
    (21, 29, 15),
    (22, 9, 19),
    (24, 10, 15),
    (28, 21, 21),
    (29, 31, 27),
    (30, 2, 1),
    (31, 3, 25),
)
EXPECTED_RAW_PROFILE = (
    15, 15, 14, 14, 15, 16, 16, 14,
    16, 14, 14, 16, 16, 16, 16, 14,
    14, 16, 16, 16, 16, 14, 14, 16,
    14, 16, 15, 15, 14, 14, 15, 14,
)
EXPECTED_AGGREGATE_SIGNATURE = 0x00800005
EXPECTED_CORRECTION_SIGNATURE = 0x00800004


def fixed_family_center_obstruction_certificate() -> dict[str, object]:
    design = tuple(auxiliary.HalfChoice(*row) for row in HALF_CHOICE_ROWS)
    target_indices = tuple(choice.target_index for choice in design)
    auxiliary_indices = tuple(choice.auxiliary_index for choice in design)
    raw_profile = auxiliary.raw_parallel_profile(design)

    aggregate_signature = 0
    for choice in design:
        for kernel_index, kernel in enumerate(auxiliary.DIRECTIONS):
            if half_kernel_sigma_formula(
                choice.target, choice.auxiliary, kernel
            ) == -1:
                aggregate_signature ^= 1 << kernel_index
    correction_signature = aggregate_signature ^ (1 << FIXED_DIRECTION_INDEX)

    cache = auxiliary._orbit_cache(design)
    intersection_size_histogram: Counter[int] = Counter()
    spatial_direction_histogram: Counter[int] = Counter()
    orientation_histogram: Counter[tuple[int, int, int]] = Counter()
    pair_center_cases = 0
    for first, second in combinations(range(len(design)), 2):
        for first_center in range(1, P):
            first_map = cache[first][first_center - 1]
            for second_center in range(1, P):
                second_map = cache[second][second_center - 1]
                pair_center_cases += 1
                common = set(first_map) & set(second_map)
                intersection_size_histogram[len(common)] += 1
                for orbit in common:
                    direction = auxiliary.frozen._spatial_direction_index(orbit)
                    spatial_direction_histogram[direction] += 1
                    orientation_histogram[
                        (direction, first_map[orbit], second_map[orbit])
                    ] += 1

    expected_pair_center_cases = len(tuple(combinations(range(16), 2))) * 30**2
    direction_zero_incidence = spatial_direction_histogram[FIXED_DIRECTION_INDEX]
    exact_input_replay = bool(
        target_indices == auxiliary.HARD
        and len(set(auxiliary_indices)) == 16
        and raw_profile == EXPECTED_RAW_PROFILE
        and aggregate_signature == EXPECTED_AGGREGATE_SIGNATURE
        and correction_signature == EXPECTED_CORRECTION_SIGNATURE
        and correction_signature.bit_count() == 2
        and not (correction_signature >> FIXED_DIRECTION_INDEX & 1)
    )
    exhaustive_pair_replay = bool(
        pair_center_cases == expected_pair_center_cases
        and sum(intersection_size_histogram.values()) == pair_center_cases
        and max(intersection_size_histogram, default=0) <= 1
    )
    physically_excluded = bool(
        exact_input_replay
        and exhaustive_pair_replay
        and direction_zero_incidence == 0
    )
    if not physically_excluded:
        raise ArithmeticError("the fixed-family centre obstruction changed")

    return {
        "schema": "resii_p31_hard_fixed_center_obstruction_v1",
        "classification": "rigorous obstruction for one fixed 16-half family",
        "p": P,
        "source_witness_sha256": SOURCE_WITNESS_SHA256,
        "option_catalog_sha256": OPTION_CATALOG_SHA256,
        "half_choices_target_auxiliary_scale": HALF_CHOICE_ROWS,
        "fixed_and_required_cancellation_direction_index": FIXED_DIRECTION_INDEX,
        "raw_parallel_profile": raw_profile,
        "exact_profile_replay": raw_profile == EXPECTED_RAW_PROFILE,
        "aggregate_signature_hex": f"{aggregate_signature:08x}",
        "correction_signature_hex": f"{correction_signature:08x}",
        "correction_signature_support": tuple(
            index
            for index in range(P + 1)
            if correction_signature >> index & 1
        ),
        "target_count": len(design),
        "nonzero_centers_per_half": P - 1,
        "half_pair_count": len(tuple(combinations(range(len(design)), 2))),
        "pair_center_cases_checked": pair_center_cases,
        "expected_pair_center_cases": expected_pair_center_cases,
        "intersection_size_histogram": dict(sorted(intersection_size_histogram.items())),
        "shared_orbit_incidence_count": sum(spatial_direction_histogram.values()),
        "shared_orbit_spatial_direction_histogram": dict(
            sorted(spatial_direction_histogram.items())
        ),
        "shared_orbit_orientation_histogram": {
            f"{direction}:{first}:{second}": count
            for (direction, first, second), count in sorted(
                orientation_histogram.items()
            )
        },
        "required_direction_shared_orbit_incidence_count": direction_zero_incidence,
        "pairwise_intersections_capture_every_multi_half_repeated_orbit": True,
        "two_half_cancellation_in_required_direction_exists": False,
        "three_or_more_half_cancellation_in_required_direction_exists": False,
        "fixed_16_half_family_physically_excluded": physically_excluded,
        "scope": (
            "excludes only the displayed sixteen labelled halves selected from "
            "the eight complement pairs; it does not exclude other complement-pair choices"
        ),
        "residual_ii_closed": False,
    }


def main() -> dict[str, object]:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = fixed_family_center_obstruction_certificate()
    if args.output is not None:
        write_json_atomic(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


if __name__ == "__main__":
    main()
