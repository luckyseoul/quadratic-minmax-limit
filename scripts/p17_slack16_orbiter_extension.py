#!/usr/bin/env python3
"""Exact PGL-orbit extension census for the p=17 slack-sixteen remainder."""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
from collections import Counter
from itertools import combinations, product
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from e1_gmin_m4_prop15700 import (  # noqa: E402
    _affine_point,
    _chart,
    _dot,
    _profile_key,
    _projective_points_or_lines,
    p17_second_boundary_profile_census,
)

from p17_orbiter_complete12_filter import (  # noqa: E402
    P,
    classify,
    determinant,
    orbiter_unrank,
)


N = P * P + P + 1


def atomic_write(path: Path, value: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def json_profile_key(key: tuple) -> list[list[list[int]]]:
    return [[list(pair) for pair in phase] for phase in key]


def main() -> None:
    source = REPO / "evidence" / "p17_arcs_d2_reps_lvl_12.csv"
    output = REPO / "evidence" / "p17_slack16_orbiter_extension.json"
    source_data = classify(source)
    records = (
        source_data["complete_representatives"]
        + source_data["extendible_representatives"]
    )
    records.sort(key=lambda row: int(row["orbit_index"]))
    if len(records) != 629:
        raise ArithmeticError("12-arc orbit census changed")

    points = tuple(orbiter_unrank(rank) for rank in range(N))
    lines = _projective_points_or_lines()
    if len(points) != N or len(lines) != N:
        raise ArithmeticError("PG(2,17) point/line count changed")

    target_rows = [
        row
        for row in p17_second_boundary_profile_census()["profiles"]
        if int(row["pair_slack"]) == 16
        and sum(
            int(row["phase_profiles_b"][phase].get(16, 0))
            for phase in ("0", "1")
        )
        == 0
    ]
    target_keys = {
        (
            tuple(sorted((int(b), int(n)) for b, n in row["phase_profiles_b"]["0"].items())),
            tuple(sorted((int(b), int(n)) for b, n in row["phase_profiles_b"]["1"].items())),
        ): index
        for index, row in enumerate(target_rows)
    }
    if len(target_rows) != 13 or len(target_keys) != 13:
        raise ArithmeticError("final p17 slack-sixteen target block changed")

    secant_index_histograms: Counter[tuple[tuple[int, int], ...]] = Counter()
    raw_quadruples = 0
    occupancy_valid = 0
    occupancy_patterns: Counter[tuple[int, int]] = Counter()
    unique_boundaries: set[tuple[int, ...]] = set()
    disjoint_chart_cases = 0
    profile_histogram: Counter[tuple] = Counter()
    hits: list[dict[str, object]] = []

    for ordinal, record in enumerate(records):
        core_ranks = tuple(int(rank) for rank in record["representative"])
        core = tuple(points[rank] for rank in core_ranks)
        core_set = set(core_ranks)
        candidates_by_index: list[list[int]] = [[] for _ in range(5)]
        for rank, point in enumerate(points):
            if rank in core_set:
                continue
            secant_index = sum(
                determinant(points[left], points[right], point) == 0
                for left, right in combinations(core_ranks, 2)
            )
            if secant_index <= 4:
                candidates_by_index[secant_index].append(rank)
        secant_index_histograms[
            tuple((index, len(bucket)) for index, bucket in enumerate(candidates_by_index))
        ] += 1

        deletion_iterators = []
        for n0 in range(5):
            for n1 in range(5 - n0):
                for n2 in range(5 - n0 - n1):
                    for n3 in range(5 - n0 - n1 - n2):
                        n4 = 4 - n0 - n1 - n2 - n3
                        counts = (n0, n1, n2, n3, n4)
                        if sum(index * count for index, count in enumerate(counts)) > 4:
                            continue
                        if any(
                            count > len(candidates_by_index[index])
                            for index, count in enumerate(counts)
                        ):
                            continue
                        deletion_iterators.append(
                            product(
                                *(
                                    combinations(candidates_by_index[index], count)
                                    for index, count in enumerate(counts)
                                )
                            )
                        )

        for choices in deletion_iterators:
            for parts in choices:
                deleted = tuple(rank for part in parts for rank in part)
                if len(deleted) != 4:
                    raise ArithmeticError("deletion composition changed")
                raw_quadruples += 1
                boundary_ranks = tuple(sorted(core_ranks + deleted))
                boundary = core + tuple(points[rank] for rank in deleted)
                occupancies = tuple(
                    sum(_dot(line, point) == 0 for point in boundary) for line in lines
                )
                if max(occupancies) > 4:
                    continue
                pattern = (occupancies.count(3), occupancies.count(4))
                if pattern not in ((0, 2), (2, 1), (4, 0)):
                    continue
                occupancy_valid += 1
                occupancy_patterns[pattern] += 1
                unique_boundaries.add(boundary_ranks)

                for line, occupancy in zip(lines, occupancies):
                    if occupancy:
                        continue
                    disjoint_chart_cases += 1
                    matrix = _chart(line)
                    affine = tuple(_affine_point(matrix, point) for point in boundary)
                    if len(set(affine)) != 16:
                        raise ArithmeticError("disjoint chart collapsed boundary points")
                    key = _profile_key(affine)
                    profile_histogram[key] += 1
                    labelled = ((key, False), ((key[1], key[0]), True))
                    for labelled_key, swapped in labelled:
                        if labelled_key not in target_keys:
                            continue
                        hits.append(
                            {
                                "target_index": target_keys[labelled_key],
                                "orbit_index": int(record["orbit_index"]),
                                "core": list(core_ranks),
                                "restored": list(deleted),
                                "line_at_infinity": list(line),
                                "phase_swapped": swapped,
                                "profile": json_profile_key(labelled_key),
                            }
                        )

        if (ordinal + 1) % 100 == 0:
            print(
                f"orbits={ordinal + 1}/629 raw={raw_quadruples} "
                f"valid={occupancy_valid} charts={disjoint_chart_cases}",
                flush=True,
            )

    expected = (97122, 47, 10, 6345, 317, 0)
    observed = (
        raw_quadruples,
        occupancy_valid,
        len(unique_boundaries),
        disjoint_chart_cases,
        len(profile_histogram),
        len(hits),
    )
    if observed != expected or dict(occupancy_patterns) != {(4, 0): 47}:
        raise ArithmeticError(f"p17 orbit-extension census changed: {observed}")
    result = {
        "experiment": "p17_slack16_orbiter_twelve_arc_extension_census",
        "p": P,
        "source": str(source.relative_to(REPO)),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "complete_twelve_arc_orbits": int(
            source_data["complete_twelve_arc_orbit_count"]
        ),
        "extendible_twelve_arc_orbits": int(
            source_data["extendible_twelve_arc_orbit_count"]
        ),
        "orbiter_build": int(source_data["orbiter_build"]),
        "twelve_arc_orbits": len(records),
        "target_profile_count": len(target_rows),
        "distinct_core_secant_index_histograms": len(secant_index_histograms),
        "raw_four_point_extensions_with_core_secant_charge_at_most_four": raw_quadruples,
        "occupancy_valid_extensions": occupancy_valid,
        "occupancy_pattern_histogram": {
            f"n3={key[0]},n4={key[1]}": value
            for key, value in sorted(occupancy_patterns.items())
        },
        "unique_boundary_rank_sets": len(unique_boundaries),
        "disjoint_line_affine_chart_cases": disjoint_chart_cases,
        "distinct_unlabelled_phase_profiles": len(profile_histogram),
        "phase_labelled_target_hits": len(hits),
        "hits": hits,
        "all_thirteen_profiles_excluded": not hits,
        "conditional_on_orbiter_exhaustive_pgl_orbit_census": True,
    }
    atomic_write(output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
