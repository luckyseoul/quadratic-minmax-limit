#!/usr/bin/env python3
"""Independent coverage audit for the complete p=7 six-finite exclusion."""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15654 import p7_nonsquare_signed_permutation  # noqa: E402
from p7_no_infinity_unsaturated_cpsat import atomic_write  # noqa: E402
from p7_size6_finite_deep_allocation_batch import exact_mean_allocations  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_hash(value) -> str:
    return hashlib.sha256(json.dumps(value, separators=(",", ":")).encode()).hexdigest()


def load_many(paths: tuple[Path, ...]) -> list[dict]:
    return [json.loads(path.read_text()) for path in paths]


def audit(
    source: Path,
    nuka_summary: Path,
    ordinary_shards: tuple[Path, ...],
    deep_shards: tuple[Path, ...],
    allocation_shards: tuple[Path, ...],
    join_shards: tuple[Path, ...],
) -> dict:
    payload = json.loads(source.read_text())
    digest = sha256(source)
    if (
        int(payload["p"]) != 7
        or int(payload["c_H"]) != -1
        or int(payload["boundary_size"]) != 6
        or int(payload["infinity_value"]) != 0
        or int(payload["orbit_count"]) != 80_704
        or int(payload["orbit_size_sum"]) != 3_856_300
    ):
        raise ValueError("orbit source has the wrong or incomplete scope")
    ordinary_indices = []
    deep_indices = []
    for index, orbit in enumerate(payload["orbits"]):
        costs = {int(key): int(value) for key, value in orbit["type_costs"].items()}
        (ordinary_indices if all(value in (24, 32) for value in costs.values()) else deep_indices).append(index)
    if len(ordinary_indices) != 80_519 or len(deep_indices) != 185:
        raise AssertionError("ordinary/deep orbit partition mismatch")

    nuka = json.loads(nuka_summary.read_text())
    nuka_ok = bool(
        nuka["host"].upper() == "NUKA"
        and nuka["floor"]["backend"] == "numpy"
        and int(nuka["floor"]["checked_boundaries"]) == 13_983_816
        and int(nuka["floor"]["floor_surviving_boundaries"]) == 3_856_300
        and nuka["floor"]["survivor_sha256"]
        == payload["upstream_floor_sieve"]["survivor_sha256"]
        and int(nuka["orbits"]["orbit_count"]) == 80_704
        and int(nuka["orbits"]["orbit_size_sum"]) == 3_856_300
        and nuka["orbits"]["ordered_orbits_canonical_sha256"]
        == canonical_hash(payload["orbits"])
        and nuka["orbits"]["profile_histogram_canonical_sha256"]
        == canonical_hash(payload["profile_histogram"])
    )

    ordinary = load_many(ordinary_shards)
    ordinary_shard_indices = sorted(int(row["shard_index"]) for row in ordinary)
    ordinary_floor_pairs: Counter[str] = Counter()
    ordinary_patterns: Counter[str] = Counter()
    ordinary_processed = ordinary_cases = ordinary_rejected = ordinary_survivors = 0
    ordinary_scope_ok = len(ordinary) == 16 and ordinary_shard_indices == list(range(16))
    for recording in ordinary:
        shard_index = int(recording["shard_index"])
        shard_count = int(recording["shard_count"])
        expected_indices = [index for index in ordinary_indices if index % shard_count == shard_index]
        expected_cases = sum(
            4 ** sum(int(value) == 24 for value in payload["orbits"][index]["type_costs"].values())
            for index in expected_indices
        )
        ordinary_scope_ok &= bool(
            recording["source_sha256"] == digest
            and recording["moduli" if "moduli" in recording else "linear_system"] is not None
            and int(recording["processed_ordinary_orbits"]) == len(expected_indices)
            and int(recording["elevation_cases"]) == expected_cases
            and int(recording["modular_infeasible_cases"]) == expected_cases
            and int(recording["surviving_cases"]) == 0
            and recording["survivors"] == []
        )
        ordinary_processed += int(recording["processed_ordinary_orbits"])
        ordinary_cases += int(recording["elevation_cases"])
        ordinary_rejected += int(recording["modular_infeasible_cases"])
        ordinary_survivors += int(recording["surviving_cases"])
        ordinary_floor_pairs.update({key: int(value) for key, value in recording["floor_pair_counts"].items()})
        ordinary_patterns.update({key: int(value) for key, value in recording["catalog_pattern_counts"].items()})
    ordinary_ok = bool(
        ordinary_scope_ok
        and ordinary_processed == 80_519
        and ordinary_cases == ordinary_rejected == 160_745
        and ordinary_survivors == 0
        and dict(ordinary_floor_pairs) == nuka["ordinary"]["floor_pair_counts"]
        and dict(ordinary_patterns) == nuka["ordinary"]["catalog_pattern_counts"]
        and int(nuka["ordinary"]["processed_ordinary_orbits"]) == 80_519
        and int(nuka["ordinary"]["elevation_cases"]) == 160_745
        and int(nuka["ordinary"]["surviving_cases"]) == 0
    )

    deep = load_many(deep_shards)
    deep_rows = [row for recording in deep for row in recording["rows"]]
    deep_by_index = {int(row["orbit_index"]): row for row in deep_rows}
    deep_scope_ok = bool(
        len(deep) == 16
        and sorted(int(row["shard_index"]) for row in deep) == list(range(16))
        and len(deep_rows) == len(deep_by_index) == 185
        and set(deep_by_index) == set(deep_indices)
        and all(recording["source_sha256"] == digest for recording in deep)
    )
    for index, row in deep_by_index.items():
        orbit = payload["orbits"][index]
        deep_scope_ok &= bool(
            tuple(row["fixed_boundary"]) == tuple(orbit["representative_vertices"])
            and row["type_floor_sums"] == orbit["type_costs"]
            and row["solver_status"] in {"INFEASIBLE", "UNKNOWN"}
            and not row["feasible"]
        )
    initial_infeasible = {index for index, row in deep_by_index.items() if row["solver_status"] == "INFEASIBLE"}
    initial_unknown = {index for index, row in deep_by_index.items() if row["solver_status"] == "UNKNOWN"}

    allocations = load_many(allocation_shards)
    allocation_rows = [row for recording in allocations for row in recording["rows"]]
    allocation_by_index = {int(row["orbit_index"]): row for row in allocation_rows}
    allocation_scope_ok = bool(
        len(allocations) == 16
        and sorted(int(row["shard_index"]) for row in allocations) == list(range(16))
        and len(allocation_rows) == len(allocation_by_index) == len(initial_unknown) == 93
        and set(allocation_by_index) == initial_unknown
        and all(recording["source_sha256"] == digest for recording in allocations)
    )
    unknown_leaf_keys = set()
    allocation_infeasible = 0
    allocation_unknown = 0
    for index, recording in allocation_by_index.items():
        expected = {
            tuple(sorted(row.items()))
            for row in exact_mean_allocations(deep_by_index[index]["direction_rows"])
        }
        leaves = recording["leaves"]
        observed = {
            tuple(sorted((int(key), int(value)) for key, value in leaf["fixed_scaled_means"].items()))
            for leaf in leaves
        }
        allocation_scope_ok &= bool(
            int(recording["allocation_count"]) == len(expected) == 10
            and observed == expected
            and all(leaf["solver_status"] in {"INFEASIBLE", "UNKNOWN"} and not leaf["feasible"] for leaf in leaves)
        )
        for leaf in leaves:
            key = (
                index,
                tuple(sorted((int(k), int(v)) for k, v in leaf["fixed_scaled_means"].items())),
            )
            if leaf["solver_status"] == "INFEASIBLE":
                allocation_infeasible += 1
            else:
                allocation_unknown += 1
                unknown_leaf_keys.add(key)

    joins = load_many(join_shards)
    join_rows = [row for recording in joins for row in recording["rows"]]
    join_by_key = {
        (
            int(row["orbit_index"]),
            tuple(sorted((int(k), int(v)) for k, v in row["fixed_scaled_means"].items())),
        ): row
        for row in join_rows
    }
    join_ok = bool(
        len(joins) == 8
        and sorted(int(row["shard_index"]) for row in joins) == list(range(8))
        and len(join_rows) == len(join_by_key) == len(unknown_leaf_keys) == 120
        and set(join_by_key) == unknown_leaf_keys
        and all(recording["source_sha256"] == digest for recording in joins)
        and all(
            row["modularly_infeasible"]
            and int(row["modular_consistent_catalog_tuples"]) == 0
            for row in join_rows
        )
    )

    anti = p7_nonsquare_signed_permutation()
    sign_transfer_ok = bool(
        anti["fixes_infinity"]
        and anti["fixes_finite_zero"]
        and anti["fixes_distinguished_edge"]
        and anti["signed_conference_anti_isometry"]
    )
    deep_ok = bool(
        deep_scope_ok
        and len(initial_infeasible) == 92
        and len(initial_unknown) == 93
        and allocation_scope_ok
        and allocation_infeasible == 810
        and allocation_unknown == 120
        and join_ok
    )
    proved = bool(nuka_ok and ordinary_ok and deep_ok and sign_transfer_ok)
    return {
        "experiment": "p7_size6_finite_global_audit",
        "status": "complete_independent_coverage_audit" if proved else "incomplete_audit",
        "proved": proved,
        # Keep the compact certificate independent of the machine-specific
        # location from which the permanent raw source was replayed.
        "source": source.name,
        "source_sha256": digest,
        "source_floor_survivors": int(payload["candidate_boundaries"]),
        "source_orbits": int(payload["orbit_count"]),
        "nuka_independent_reproduction": {"proved": nuka_ok, "summary_sha256": sha256(nuka_summary)},
        "ordinary": {
            "proved": ordinary_ok,
            "orbits": ordinary_processed,
            "elevation_cases": ordinary_cases,
            "rejected_cases": ordinary_rejected,
            "surviving_cases": ordinary_survivors,
            "floor_pair_counts": dict(ordinary_floor_pairs),
            "catalog_pattern_counts": dict(ordinary_patterns),
        },
        "deep": {
            "proved": deep_ok,
            "orbits": len(deep_by_index),
            "initial_infeasible_orbits": len(initial_infeasible),
            "initial_unknown_orbits": len(initial_unknown),
            "allocation_leaves": allocation_infeasible + allocation_unknown,
            "allocation_infeasible_leaves": allocation_infeasible,
            "allocation_unknown_leaves": allocation_unknown,
            "low_catalog_join_leaves": len(join_rows),
            "low_catalog_join_survivors": sum(not row["modularly_infeasible"] for row in join_rows),
        },
        "sign_transfer": {
            "proved": sign_transfer_ok,
            "nonsquare_multiplier": anti["nonsquare_multiplier"],
            "fixes_distinguished_edge": anti["fixes_distinguished_edge"],
            "signed_conference_anti_isometry": anti["signed_conference_anti_isometry"],
        },
        "p7_six_finite_both_product_signs_closed": proved,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--nuka-summary", type=Path, required=True)
    parser.add_argument("--ordinary-shards", type=Path, nargs="+", required=True)
    parser.add_argument("--deep-shards", type=Path, nargs="+", required=True)
    parser.add_argument("--allocation-shards", type=Path, nargs="+", required=True)
    parser.add_argument("--join-shards", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = audit(
        args.source,
        args.nuka_summary,
        tuple(args.ordinary_shards),
        tuple(args.deep_shards),
        tuple(args.allocation_shards),
        tuple(args.join_shards),
    )
    atomic_write(args.output, out)
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
