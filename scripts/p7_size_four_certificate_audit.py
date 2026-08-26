#!/usr/bin/env python3
"""Independent audit of the p=7 infinity-plus-three orbit certificates."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p7_size_four_slack_classify import classify_three_odd_fibres  # noqa: E402
from residual_size_four_boundary_orbits import classify  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(
        f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp"
    )
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def audit(orbit_path: Path, result_paths: list[Path]) -> dict:
    source = json.loads(orbit_path.read_text())
    if (
        int(source["p"]) != 7
        or int(source["c_H"]) != 1
        or int(source["infinity_value"]) != 1
    ):
        raise ValueError("unexpected orbit source scope")
    orbits = source["orbits"]
    expected = {
        index: tuple(int(v) for v in row["representative_vertices"])
        for index, row in enumerate(orbits)
    }
    boundary_to_index = {boundary: index for index, boundary in expected.items()}
    valid_infeasible: set[int] = set()
    malformed = []
    statuses: dict[int, list[str]] = {index: [] for index in expected}
    rows_seen = 0
    for result_path in result_paths:
        payload = json.loads(result_path.read_text())
        payload_rows = payload.get("rows")
        if payload_rows is None:
            payload_rows = [payload]
        for row in payload_rows:
            rows_seen += 1
            boundary = tuple(int(v) for v in row["fixed_boundary"])
            index = int(row.get("orbit_index", boundary_to_index.get(boundary, -1)))
            status = str(row["solver_status"])
            if index not in expected:
                malformed.append({"path": str(result_path), "orbit_index": index})
                continue
            statuses[index].append(status)
            structurally_valid = bool(
                boundary == expected[index]
                and int(row["p"]) == 7
                and int(row["c_H"]) == 1
                and int(row["boundary_size"]) == 4
                and int(row["infinity_value"]) == 1
                and row["distinguished_edge"] == [0, 1]
                and row["shell_mode"] == "affine"
                and int(row["p7_saturated_slack_equalities"]) == 177
                and int(row["n_score_constraints"]) == 0
            )
            if not structurally_valid:
                malformed.append(
                    {"path": str(result_path), "orbit_index": index, "row": row}
                )
            elif status == "INFEASIBLE" and bool(
                row["finite_infeasibility_certificate"]
            ):
                valid_infeasible.add(index)

    missing = sorted(set(expected) - valid_infeasible)
    slack = classify_three_odd_fibres()
    fresh_orbits = classify(7, 1, 1)
    orbit_reclassification_matches = bool(
        int(fresh_orbits["candidate_boundaries"]) == 18424
        and int(fresh_orbits["stabilizer_size"]) == 48
        and int(fresh_orbits["orbit_count"]) == 416
        and [
            (row["representative_vertices"], int(row["size"]))
            for row in fresh_orbits["orbits"]
        ]
        == [
            (row["representative_vertices"], int(row["size"]))
            for row in orbits
        ]
    )
    orbit_coverage = bool(
        int(source["candidate_boundaries"]) == 18424
        and int(source["stabilizer_size"]) == 48
        and int(source["orbit_count"]) == 416
        and int(source["orbit_size_sum"]) == 18424
        and sum(int(row["size"]) for row in orbits) == 18424
    )
    proved = bool(
        orbit_coverage
        and orbit_reclassification_matches
        and slack["proved"]
        and not malformed
        and not missing
        and len(valid_infeasible) == 416
    )
    return {
        "experiment": "p7_size_four_certificate_audit",
        "status": "independent_orbit_and_certificate_audit",
        "proved": proved,
        "orbit_source": str(orbit_path),
        "orbit_source_sha256": sha256(orbit_path),
        "result_files": [
            {"path": str(path), "sha256": sha256(path)} for path in result_paths
        ],
        "candidate_boundaries": int(source["candidate_boundaries"]),
        "stabilizer_size": int(source["stabilizer_size"]),
        "orbit_count": int(source["orbit_count"]),
        "orbit_size_sum": int(source["orbit_size_sum"]),
        "orbit_coverage_valid": orbit_coverage,
        "fresh_orbit_reclassification_matches": orbit_reclassification_matches,
        "slack_classification": {
            key: value for key, value in slack.items() if key != "survivors"
        },
        "certificate_rows_seen_including_duplicates": rows_seen,
        "valid_infeasible_orbits": len(valid_infeasible),
        "missing_orbits": missing,
        "malformed_rows": malformed,
        "statuses_by_missing_orbit": {
            str(index): statuses[index] for index in missing
        },
        "feasible_rows": sum(
            status in ("FEASIBLE", "OPTIMAL")
            for rows in statuses.values()
            for status in rows
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orbits", type=Path, required=True)
    parser.add_argument("--results", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = audit(args.orbits, args.results)
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        atomic_write(args.output, out)


if __name__ == "__main__":
    main()
