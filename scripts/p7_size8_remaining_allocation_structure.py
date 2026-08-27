#!/usr/bin/env python3
"""Classify the exact post-15.664 p=7 size-eight allocation leaves.

This is a structural prepass for the next GPU sieve.  It reads the complete
floor-profile census, removes the already closed conic, forced-floor, and
four-allocation strata, and records how many directional catalogs are raised
in each of the 11-, 16-, 24-, and 44-allocation leaves.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import itertools
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from p7_fixed_boundary_mean_allocation_batch import allocations  # noqa: E402
from p7_fixed_boundary_modular_cpsat import linear_data  # noqa: E402
from p7_exceptional_omit_high_catalogs import modular_rank  # noqa: E402
from p7_size8_one_elevation_gpu import abstract_rows, type_costs  # noqa: E402
from p7_size8_floor_profile_gpu import direction_tables  # noqa: E402


EXPECTED_BOUNDARIES = {11: 154_056, 16: 1_194_816, 24: 1_176, 44: 69_384}
EXPECTED_PROFILES = {11: 248, 16: 516, 24: 8, 44: 110}


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def counter_json(counter: Counter) -> dict[str, int]:
    return {str(key): int(value) for key, value in sorted(counter.items(), key=lambda row: str(row[0]))}


def run(source_path: Path) -> dict:
    source = json.loads(source_path.read_text())
    if (
        source.get("experiment") != "p7_size8_floor_profile_gpu"
        or source.get("status") != "complete_exact_floor_profile_census"
        or int(source.get("p", 0)) != 7
        or int(source.get("c_H", 0)) != -1
    ):
        raise ValueError("source is not the complete c_H=-1 p=7 floor census")
    _labels, epsilon_array = direction_tables()
    epsilons = tuple(int(value) for value in epsilon_array)

    records: dict[int, dict[str, Counter | int]] = {
        count: {
            "profiles": 0,
            "boundaries": 0,
            "leaves": 0,
            "floor_pairs_by_boundary": Counter(),
            "odd_secants_by_boundary": Counter(),
            "raised_support_by_leaf": Counter(),
            "raised_amounts_by_leaf": Counter(),
            "raised_type_support_by_leaf": Counter(),
            "remaining_36row_floor_catalogs_by_leaf": Counter(),
            "profile_signature_sets": Counter(),
        }
        for count in EXPECTED_BOUNDARIES
    }
    for source_row in source["survivor_ordered_profiles"]:
        profile = tuple(int(value) for value in source_row["b_by_direction"])
        boundary_count = int(source_row["count"])
        floors = type_costs(profile, epsilons)
        if sum(profile) == 8 or floors == (32, 32):
            continue
        rows = abstract_rows(profile, epsilons)
        leaves = allocations(rows)
        if len(leaves) == 4:
            continue
        if len(leaves) not in records:
            raise AssertionError(f"unexpected remaining allocation count {len(leaves)}")
        record = records[len(leaves)]
        record["profiles"] = int(record["profiles"]) + 1
        record["boundaries"] = int(record["boundaries"]) + boundary_count
        record["leaves"] = int(record["leaves"]) + boundary_count * len(leaves)
        record["floor_pairs_by_boundary"][floors] += boundary_count
        record["odd_secants_by_boundary"][sum(profile)] += boundary_count

        floor_vector = tuple(row["floor"] for row in rows)
        profile_signatures: Counter = Counter()
        for leaf in leaves:
            differences = tuple(value - floor for value, floor in zip(leaf, floor_vector))
            support = tuple(index for index, value in enumerate(differences) if value)
            amounts = tuple(sorted(differences[index] for index in support))
            type_support = (
                sum(epsilons[index] == -1 for index in support),
                sum(epsilons[index] == 1 for index in support),
            )
            record["raised_support_by_leaf"][len(support)] += boundary_count
            record["raised_amounts_by_leaf"][amounts] += boundary_count
            record["raised_type_support_by_leaf"][type_support] += boundary_count
            remaining_variable = sum(
                epsilons[index] == -1 and profile[index] == 4 and index not in support
                for index in range(8)
            )
            record["remaining_36row_floor_catalogs_by_leaf"][remaining_variable] += boundary_count
            profile_signatures[(len(support), amounts, type_support)] += 1
        record["profile_signature_sets"][tuple(sorted(profile_signatures.items()))] += 1

    output_records = {}
    for allocation_count, record in records.items():
        if int(record["boundaries"]) != EXPECTED_BOUNDARIES[allocation_count]:
            raise AssertionError(f"boundary count changed for {allocation_count}")
        if int(record["profiles"]) != EXPECTED_PROFILES[allocation_count]:
            raise AssertionError(f"profile count changed for {allocation_count}")
        expected_leaves = allocation_count * EXPECTED_BOUNDARIES[allocation_count]
        if int(record["leaves"]) != expected_leaves:
            raise AssertionError(f"leaf count changed for {allocation_count}")
        output_records[str(allocation_count)] = {
            "ordered_profiles": int(record["profiles"]),
            "boundaries_per_sign": int(record["boundaries"]),
            "allocation_leaves_per_sign": int(record["leaves"]),
            "floor_pairs_by_boundary": counter_json(record["floor_pairs_by_boundary"]),
            "odd_secants_by_boundary": counter_json(record["odd_secants_by_boundary"]),
            "raised_support_by_leaf": counter_json(record["raised_support_by_leaf"]),
            "raised_amounts_by_leaf": counter_json(record["raised_amounts_by_leaf"]),
            "raised_type_support_by_leaf": counter_json(record["raised_type_support_by_leaf"]),
            "remaining_36row_floor_catalogs_by_leaf": counter_json(
                record["remaining_36row_floor_catalogs_by_leaf"]
            ),
            "distinct_profile_signature_sets": len(record["profile_signature_sets"]),
            "profile_signature_set_histogram": counter_json(record["profile_signature_sets"]),
        }

    _matrix, dependencies, _linear_rows = linear_data((7,))
    dependency = dependencies[7]
    conditioned_dimensions: dict[str, dict[str, dict[str, int]]] = {}
    for size in range(1, 6):
        by_type: dict[str, Counter] = {}
        for subset in itertools.combinations(range(8), size):
            columns = [
                column
                for direction in subset
                for column in range(2 + 35 * direction, 2 + 35 * (direction + 1))
            ]
            dimension = 135 - modular_rank(dependency[:, columns])
            type_support = (
                sum(epsilons[index] == -1 for index in subset),
                sum(epsilons[index] == 1 for index in subset),
            )
            key = str(type_support)
            by_type.setdefault(key, Counter())[dimension] += 1
        conditioned_dimensions[str(size)] = {
            key: counter_json(histogram) for key, histogram in sorted(by_type.items())
        }

    return {
        "experiment": "p7_size8_remaining_allocation_structure",
        "status": "complete_exact_post_15664_structure",
        "source": str(source_path),
        "source_sha256": sha256(source_path),
        "p": 7,
        "c_H": -1,
        "remaining_boundaries_per_sign": sum(EXPECTED_BOUNDARIES.values()),
        "remaining_allocation_leaves_per_sign": sum(
            count * boundaries for count, boundaries in EXPECTED_BOUNDARIES.items()
        ),
        "conditioned_mod7_dependency_dimensions_by_omitted_support": conditioned_dimensions,
        "strata": output_records,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = run(args.source)
    atomic_json(args.output, output)
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
