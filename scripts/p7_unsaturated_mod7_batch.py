#!/usr/bin/env python3
"""Exact mod-seven exhaustion of all unsaturated p=7 four-finite cases.

For every fixed boundary and elevated-direction choice, each complete slack
catalog fixes the 280 affine bad-edge counts.  The common linear system has
282 equations (edge count, distinguished edge, and those 280 counts) in
1225 edge indicators.  Its rank over F_7 is only 147, so 135 exact left-null
dependencies constrain the catalog right-hand side.

There are at most two non-singleton directional catalogs in any remaining
case.  Their dependency syndromes are joined by exact hashing, avoiding a
Cartesian enumeration as large as 2233*1764.  Zero matching tuples is a
rigorous finite infeasibility certificate for the case.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from p7_no_infinity_unsaturated_cpsat import atomic_write  # noqa: E402
from p7_no_infinity_unsaturated_orbit_batch import elevation_cases  # noqa: E402
from p7_unsaturated_gf2_catalog_filter import direction_scope  # noqa: E402
from p7_unsaturated_modular_catalog_filter import (  # noqa: E402
    equation_matrix,
    left_dependencies,
)
from p7_unsaturated_slack_catalog import exact_slack_catalog_values  # noqa: E402


POINTS = tuple(itertools.combinations(range(7), 4))
POINT_INDEX = {point: index for index, point in enumerate(POINTS)}


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mapped_slack_catalog(
    odd_fibres: int,
    phase: int,
    scaled_mean: int,
    B: set[int],
) -> np.ndarray:
    """Return complete slack values with canonical fibres relabelled to B."""
    if len(B) != odd_fibres:
        raise ValueError("odd-fibre set has the wrong size")
    canonical = np.asarray(
        exact_slack_catalog_values(odd_fibres, phase, scaled_mean),
        dtype=np.int16,
    )
    actual_B = sorted(B)
    actual_complement = sorted(set(range(7)) - B)
    permutation = dict(zip(range(odd_fibres), actual_B)) | dict(
        zip(range(odd_fibres, 7), actual_complement)
    )
    inverse = {target: source for source, target in permutation.items()}
    source_columns = [
        POINT_INDEX[tuple(sorted(inverse[value] for value in point))]
        for point in POINTS
    ]
    mapped = canonical[:, source_columns]
    if len({bytes(row) for row in mapped.astype(np.int8)}) != len(mapped):
        raise AssertionError("fibre relabeling collapsed a catalog")
    return mapped


def contribution_matrix(
    dependencies: np.ndarray,
    direction_index: int,
    odd_fibres: int,
    phase: int,
    scaled_mean: int,
    B: set[int],
) -> np.ndarray:
    """Dependency syndromes, one column per complete catalog row."""
    values = mapped_slack_catalog(
        odd_fibres, phase, scaled_mean, B
    ).astype(np.int64)
    # epsilon*S=3+2A and 29 selected edges imply bad=13-A exactly.
    bad_counts = 13 - values
    if np.any(bad_counts < 0):
        raise AssertionError("catalog score exceeds the 29-edge range")
    block = dependencies[
        :, 2 + 35 * direction_index : 2 + 35 * (direction_index + 1)
    ]
    return (block @ bad_counts.T % 7).astype(np.uint8)


def count_zero_sum_tuples(
    base: np.ndarray,
    contributions: list[np.ndarray],
) -> int:
    """Count exact catalog tuples whose total dependency syndrome is zero."""
    target = (-base.astype(np.int16)) % 7
    if len(contributions) == 1:
        matrix = contributions[0]
        return sum(
            bytes(matrix[:, index]) == bytes(target.astype(np.uint8))
            for index in range(matrix.shape[1])
        )
    if len(contributions) != 2:
        raise AssertionError("remaining cases must have one or two catalogs")
    first, second = contributions
    first_counts = Counter(
        bytes(first[:, index]) for index in range(first.shape[1])
    )
    total = 0
    for index in range(second.shape[1]):
        needed = (target - second[:, index].astype(np.int16)) % 7
        total += first_counts.get(bytes(needed.astype(np.uint8)), 0)
    return total


def run_batch(source: Path, c_h: int) -> dict:
    started = time.time()
    payload = json.loads(source.read_text())
    if int(payload["p"]) != 7 or int(payload["c_H"]) != c_h:
        raise ValueError("orbit source does not match p=7 and c_H")

    matrix = equation_matrix()
    rank, dependencies = left_dependencies(matrix, 7)
    if rank != 147 or dependencies.shape != (135, 282):
        raise AssertionError(
            f"unexpected mod-seven linear dimensions rank={rank}, "
            f"dependencies={dependencies.shape}"
        )
    if np.any(dependencies @ (matrix % 7) % 7):
        raise AssertionError("left-null dependency audit failed")
    # Calibration: every right side manufactured as A*x must pass all cuts.
    calibration = np.arange(matrix.shape[1], dtype=np.int64) % 2
    if np.any(dependencies @ ((matrix @ calibration) % 7) % 7):
        raise AssertionError("manufactured consistent right side was rejected")

    base_edge_syndrome = (
        dependencies[:, :2] @ np.asarray([29, 1], dtype=np.int64) % 7
    ).astype(np.uint8)
    contribution_cache: dict[tuple, np.ndarray] = {}
    rows = []
    catalog_pattern_counts: Counter[tuple[int, ...]] = Counter()

    unsaturated = [
        (orbit_index, orbit)
        for orbit_index, orbit in enumerate(payload["orbits"])
        if any(int(value) != 32 for value in orbit["type_costs"].values())
    ]
    for orbit_index, orbit in unsaturated:
        boundary = tuple(int(value) for value in orbit["representative_vertices"])
        for elevated in elevation_cases(orbit):
            direction_rows, type_floors = direction_scope(c_h, boundary, elevated)
            elevated_set = set(elevated)
            base = base_edge_syndrome.astype(np.int16).copy()
            variable_contributions = []
            variable_rows = []
            for direction_index, direction_row in enumerate(direction_rows):
                eps = int(direction_row["eps"])
                odd_fibres = len(direction_row["B"])
                phase = int(direction_row["phase"])
                scaled_mean = int(direction_row["floor"]) + (
                    8 if direction_index in elevated_set else 0
                )
                key = (
                    direction_index,
                    odd_fibres,
                    phase,
                    scaled_mean,
                    tuple(sorted(direction_row["B"])),
                )
                if key not in contribution_cache:
                    contribution_cache[key] = contribution_matrix(
                        dependencies,
                        direction_index,
                        odd_fibres,
                        phase,
                        scaled_mean,
                        set(direction_row["B"]),
                    )
                contribution = contribution_cache[key]
                if contribution.shape[1] == 1:
                    base = (base + contribution[:, 0]) % 7
                else:
                    variable_contributions.append(contribution)
                    variable_rows.append(
                        {
                            "direction_index": direction_index,
                            "eps": eps,
                            "b": odd_fibres,
                            "phase": phase,
                            "scaled_mean": scaled_mean,
                            "catalog_rows": int(contribution.shape[1]),
                        }
                    )
            pattern = tuple(
                sorted(
                    (int(matrix.shape[1]) for matrix in variable_contributions),
                    reverse=True,
                )
            )
            catalog_pattern_counts[pattern] += 1
            consistent = count_zero_sum_tuples(
                base.astype(np.uint8), variable_contributions
            )
            rows.append(
                {
                    "orbit_index": orbit_index,
                    "orbit_size": int(orbit["size"]),
                    "representative_vertices": list(boundary),
                    "type_floor_sums": {
                        str(key): value for key, value in type_floors.items()
                    },
                    "elevated_directions": list(elevated),
                    "non_singleton_catalogs": variable_rows,
                    "catalog_tuple_count": int(np.prod(pattern, dtype=object)),
                    "mod7_consistent_catalog_tuples": int(consistent),
                    "mod7_infeasible": consistent == 0,
                }
            )

    surviving = [row for row in rows if not row["mod7_infeasible"]]
    return {
        "experiment": "p7_unsaturated_mod7_batch",
        "status": "complete_exact_mod_seven_catalog_syndrome_exhaustion",
        "p": 7,
        "c_H": c_h,
        "source": str(source),
        "source_sha256": source_hash(source),
        "linear_system": {
            "equations": int(matrix.shape[0]),
            "edge_variables": int(matrix.shape[1]),
            "modulus": 7,
            "rank": rank,
            "left_dependency_dimension": int(len(dependencies)),
            "left_null_audit": True,
            "manufactured_rhs_calibration": True,
        },
        "unsaturated_orbits": len(unsaturated),
        "unsaturated_boundary_size_sum": sum(
            int(orbit["size"]) for _index, orbit in unsaturated
        ),
        "fixed_elevation_cases": len(rows),
        "catalog_pattern_counts": {
            "x".join(map(str, pattern)): count
            for pattern, count in sorted(catalog_pattern_counts.items())
        },
        "syndrome_contribution_cache_entries": len(contribution_cache),
        "mod7_infeasible_cases": len(rows) - len(surviving),
        "surviving_cases": len(surviving),
        "all_cases_mod7_infeasible": not surviving,
        "survivors": surviving,
        "elapsed_seconds": time.time() - started,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--c-h", type=int, choices=(-1, 1), default=-1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = run_batch(args.source, args.c_h)
    if args.output is not None:
        atomic_write(args.output, out)
    summary = {key: value for key, value in out.items() if key != "rows"}
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
