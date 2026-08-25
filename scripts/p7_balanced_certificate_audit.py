#!/usr/bin/env python3
"""Audit the complete balanced-p7 negative two-point finite certificate."""
from __future__ import annotations

import argparse
import glob
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np


POINTS = tuple(itertools.combinations(range(7), 4))
PAIRS = tuple(itertools.combinations(range(7), 2))
EXPECTED_LIFT_COUNTS = {
    (2, 2, 2, 2, 2): 56,
    (1, 1, 1, 1, 1, 1, 2, 2): 280,
    (1, 1, 1, 1, 1, 1, 1, 1, 2): 420,
    (1, 1, 1, 1, 1, 1, 1, 1, 1, 1): 1008,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def vector_is_quadratic(vector: tuple[int, ...]) -> bool:
    vertex_mass = [
        sum(vector[a] for a, X in enumerate(POINTS) if i in X) for i in range(7)
    ]
    pair_mass = {
        pair: sum(
            vector[a]
            for a, X in enumerate(POINTS)
            if pair[0] in X and pair[1] in X
        )
        for pair in PAIRS
    }
    return all(
        6 * vector[a]
        == 2 * sum(pair_mass[pair] for pair in itertools.combinations(X, 2))
        - 3 * sum(vertex_mass[i] for i in X)
        + 36
        for a, X in enumerate(POINTS)
    )


def audit(
    star_files: dict[tuple[int, int], Path],
    shard_patterns: dict[tuple[int, int], str],
    retry_path: Path,
    lift_jsons: list[Path],
    lift_npzs: list[Path],
    combined_table: Path,
) -> dict:
    files = [*star_files.values(), retry_path, *lift_jsons, *lift_npzs, combined_table]
    pairs = {}
    all_unknown_keys = set()
    total_main_infeasible = 0
    total_main_unknown = 0
    for pair in sorted(star_files):
        source = json.loads(star_files[pair].read_text())
        assert tuple(source["exception_pair"]) == pair
        assert source["all_star_count"] == 18424
        assert source["survivor_count"] == 18214
        assert source["stabilizer_size"] == 6
        assert source["orbit_count"] == 3038
        assert sum(row["size"] for row in source["orbits"]) == 18214

        shard_paths = [Path(path) for path in sorted(glob.glob(shard_patterns[pair]))]
        files.extend(shard_paths)
        assert len(shard_paths) == 32
        rows = []
        for path in shard_paths:
            payload = json.loads(path.read_text())
            assert tuple(payload["exception_pair"]) == pair
            assert payload["completed_count"] == payload["assigned_count"]
            assert payload["feasible_count"] == 0
            rows.extend(payload["rows"])
        assert len(rows) == 3038
        by_index = {row["orbit_index"]: row for row in rows}
        assert len(by_index) == 3038
        assert set(by_index) == set(range(3038))
        for index, orbit in enumerate(source["orbits"]):
            row = by_index[index]
            assert row["representative"] == orbit["representative"]
            assert row["orbit_size"] == orbit["size"]
            assert not row["feasible"]
            assert row["solver_status"] in ("INFEASIBLE", "UNKNOWN")
        statuses = Counter(row["solver_status"] for row in rows)
        unknown_keys = {
            (pair, row["orbit_index"])
            for row in rows
            if row["solver_status"] == "UNKNOWN"
        }
        all_unknown_keys |= unknown_keys
        total_main_infeasible += statuses["INFEASIBLE"]
        total_main_unknown += statuses["UNKNOWN"]
        pairs[f"{pair[0]}_{pair[1]}"] = {
            "all_stars": source["all_star_count"],
            "l1_filtered": source["all_star_count"] - source["survivor_count"],
            "surviving_stars": source["survivor_count"],
            "orbit_count": source["orbit_count"],
            "main_status_counts": dict(statuses),
        }

    retry = json.loads(retry_path.read_text())
    assert retry["status"] == "complete_all_infeasible"
    assert retry["assigned_count"] == retry["completed_count"]
    assert retry["unknown_count"] == retry["feasible_count"] == 0
    retry_keys = {
        (tuple(row["exception_pair"]), row["orbit_index"]) for row in retry["rows"]
    }
    assert retry_keys == all_unknown_keys
    assert all(row["solver_status"] == "INFEASIBLE" for row in retry["rows"])

    enumerated_vectors = []
    lift_counts = {}
    for json_path, npz_path in zip(lift_jsons, lift_npzs):
        metadata = json.loads(json_path.read_text())
        assert metadata["complete"] and metadata["solver_status"] == "OPTIMAL"
        vectors = np.asarray(np.load(npz_path)["values"], dtype=np.int8)
        histogram = tuple(metadata["positive_values"])
        expected = EXPECTED_LIFT_COUNTS[histogram]
        assert metadata["solution_count"] == len(vectors) == expected
        assert len(np.unique(vectors, axis=0)) == expected
        assert all(tuple(sorted(row[row > 0].tolist())) == histogram for row in vectors)
        enumerated_vectors.append(vectors)
        lift_counts[",".join(map(str, histogram))] = expected
    union = np.unique(np.concatenate(enumerated_vectors), axis=0)
    table = np.asarray(np.load(combined_table)["values"], dtype=np.int8)
    assert len(union) == len(table) == sum(EXPECTED_LIFT_COUNTS.values()) == 1764
    assert np.array_equal(union, table)
    assert all(vector_is_quadratic(tuple(map(int, row))) for row in table)

    hashes = {path.name: sha256(path) for path in sorted(set(files))}
    manifest_text = "".join(f"{name} {digest}\n" for name, digest in sorted(hashes.items()))
    return {
        "experiment": "p7_balanced_certificate_audit",
        "proved": True,
        "exception_pair_orbit_representatives": [[0, 1], [0, 3]],
        "pairs": pairs,
        "main_infeasible": total_main_infeasible,
        "main_unknown_retried": total_main_unknown,
        "retry_infeasible": len(retry["rows"]),
        "final_unknown": 0,
        "final_feasible": 0,
        "lift_histogram_counts": lift_counts,
        "total_lift_vectors": len(table),
        "file_sha256": hashes,
        "manifest_sha256": hashlib.sha256(manifest_text.encode()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--star-e01", type=Path, required=True)
    parser.add_argument("--star-e03", type=Path, required=True)
    parser.add_argument("--shards-e01", required=True)
    parser.add_argument("--shards-e03", required=True)
    parser.add_argument("--retry", type=Path, required=True)
    parser.add_argument("--lift-json", type=Path, action="append", required=True)
    parser.add_argument("--lift-npz", type=Path, action="append", required=True)
    parser.add_argument("--combined-table", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if len(args.lift_json) != 4 or len(args.lift_npz) != 4:
        raise ValueError("supply four lift JSON/NPZ pairs in matching order")
    result = audit(
        {(0, 1): args.star_e01, (0, 3): args.star_e03},
        {(0, 1): args.shards_e01, (0, 3): args.shards_e03},
        args.retry,
        args.lift_json,
        args.lift_npz,
        args.combined_table,
    )
    rendered = json.dumps(result, indent=2)
    print(rendered, flush=True)
    args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
