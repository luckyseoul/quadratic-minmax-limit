#!/usr/bin/env python3
"""Independent coverage audit for the unsaturated p=7 orbit certificates."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15654 import p7_nonsquare_signed_permutation  # noqa: E402
from p7_no_infinity_unsaturated_orbit_batch import elevation_cases  # noqa: E402
from p7_no_infinity_unsaturated_partition_retry import (  # noqa: E402
    initial_intervals,
    interval_leaves,
)
from p7_unsaturated_slack_catalog import exact_slack_catalog_values  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def catalog_audit() -> dict:
    points = tuple(itertools.combinations(range(7), 4))
    excess = set(exact_slack_catalog_values(0, 0, 8))
    translated_counts = {}
    translated_matches = True
    for phase, mean in ((0, 16), (1, 14)):
        parity = tuple(
            (sum(vertex in point for vertex in (0, 1)) + phase) & 1
            for point in points
        )
        translated = {
            tuple(row[index] - parity[index] for index in range(35))
            for row in exact_slack_catalog_values(2, phase, mean)
        }
        translated_counts[f"b2_phase{phase}_mean{mean}"] = len(translated)
        translated_matches &= translated == excess

    odd_set = {0, 1, 2, 3}
    minimum = tuple(
        (len(set(point) & odd_set) - 2) ** 2 for point in points
    )
    deficits = Counter()
    for row in exact_slack_catalog_values(4, 0, 16):
        negative = tuple(
            value
            for value in (
                row[index] - minimum[index] for index in range(35)
            )
            if value < 0
        )
        deficits[negative] += 1
    b4_split = {
        "regular": deficits[()],
        "single_minus_two": deficits[(-2,)],
        "single_minus_four": deficits[(-4,)],
    }
    proved = bool(
        len(excess) == 1764
        and translated_matches
        and translated_counts
        == {"b2_phase0_mean16": 1764, "b2_phase1_mean14": 1764}
        and b4_split
        == {"regular": 1764, "single_minus_two": 448, "single_minus_four": 21}
        and sum(b4_split.values()) == 2233
        and len(exact_slack_catalog_values(4, 1, 14)) == 36
    )
    return {
        "proved": proved,
        "universal_even_excess_count": len(excess),
        "b2_translated_counts": translated_counts,
        "b2_translations_match_universal_excess": translated_matches,
        "b4_mean16_split": b4_split,
        "b4_mean16_total": sum(b4_split.values()),
        "b4_phase1_mean14_count": len(exact_slack_catalog_values(4, 1, 14)),
    }


def partition_audit(
    source_payload: dict,
    source_digest: str,
    paths: list[Path],
) -> tuple[set[tuple[int, tuple[int, ...]]], list[dict], list[dict]]:
    closed: set[tuple[int, tuple[int, ...]]] = set()
    summaries = []
    malformed = []
    for path in paths:
        payload = json.loads(path.read_text())
        summary = {
            "path": str(path),
            "sha256": sha256(path),
            "orbit_index": payload.get("orbit_index"),
            "elevated_directions": payload.get("elevated_directions"),
            "partition_direction": payload.get("partition_direction"),
            "catalog_total": payload.get("catalog_total"),
            "attempted_intervals": payload.get("attempted_intervals"),
            "terminal_status_counts": payload.get("terminal_status_counts"),
        }
        summaries.append(summary)
        try:
            orbit_index = int(payload["orbit_index"])
            if not 0 <= orbit_index < len(source_payload["orbits"]):
                raise ValueError("orbit index")
            orbit = source_payload["orbits"][orbit_index]
            boundary = tuple(int(value) for value in orbit["representative_vertices"])
            elevated = tuple(int(value) for value in payload["elevated_directions"])
            partition_direction = int(payload["partition_direction"])
            catalog_total = int(payload["catalog_total"])
            initial_chunk_size = int(payload["initial_chunk_size"])
            min_chunk_size = int(payload["min_chunk_size"])
            if initial_chunk_size <= 0 or min_chunk_size <= 0:
                raise ValueError("partition chunk sizes")
            if (
                payload["source_sha256"] != source_digest
                or int(payload["c_H"]) != -1
                or payload.get("catalog_partition_basis")
                != "mapped_target_catalog_rows_lexicographic_v1"
                or tuple(payload["representative_vertices"]) != boundary
                or elevated not in elevation_cases(orbit)
                or partition_direction not in elevated
            ):
                raise ValueError("partition scope")
            direction = orbit["direction_rows"][partition_direction]
            phase = int(int(direction["eps"]) == -1)
            expected_total = len(
                exact_slack_catalog_values(
                    int(direction["b"]),
                    phase,
                    int(direction["floor"]) + 8,
                )
            )
            if catalog_total != expected_total:
                raise ValueError("catalog total")
            rows_by_interval = {}
            for row in payload["rows"]:
                interval = (int(row["catalog_start"]), int(row["catalog_stop"]))
                if interval in rows_by_interval:
                    raise ValueError("duplicate interval")
                start, stop = interval
                result = row["result"]
                if not 0 <= start < stop <= catalog_total:
                    raise ValueError("interval outside catalog")
                if (
                    int(row["catalog_rows"]) != stop - start
                    or tuple(result["fixed_boundary"]) != boundary
                    or int(result["c_H"]) != -1
                    or tuple(result["fixed_elevated_directions"]) != elevated
                    or result.get("catalog_ranges")
                    != {str(partition_direction): [start, stop]}
                    or result.get("catalog_partition_basis")
                    != "mapped_target_catalog_rows_lexicographic_v1"
                    or bool(result.get("direct_score_cuts"))
                    != bool(payload.get("direct_score_cuts"))
                    or bool(result.get("pointwise_score_equalities"))
                    != bool(payload.get("pointwise_score_equalities"))
                    or bool(result.get("pointwise_only"))
                    != bool(payload.get("pointwise_only"))
                    or int(
                        result["direction_rows"][partition_direction][
                            "target_option_total_count"
                        ]
                    )
                    != catalog_total
                    or int(
                        result["direction_rows"][partition_direction][
                            "target_option_count"
                        ]
                    )
                    != stop - start
                ):
                    raise ValueError("interval result scope")
                rows_by_interval[interval] = row
            roots = initial_intervals(catalog_total, initial_chunk_size)
            leaves = interval_leaves(roots, rows_by_interval, min_chunk_size)
            cursor = 0
            for start, stop, status in leaves:
                if start != cursor or stop <= start:
                    raise ValueError("terminal partition coverage")
                cursor = stop
                row = rows_by_interval.get((start, stop))
                if (
                    status != "INFEASIBLE"
                    or row is None
                    or row["result"]["finite_infeasibility_certificate"] is not True
                    or row["result"]["feasible"] is not False
                ):
                    raise ValueError("terminal interval is not infeasible")
            if cursor != catalog_total or not leaves:
                raise ValueError("terminal partition incomplete")
            if payload.get("proved") is not True:
                raise ValueError("payload did not mark the reconstructed proof")
            key = (orbit_index, elevated)
            if key in closed:
                raise ValueError("duplicate partition proof for case")
            closed.add(key)
            summary["independently_verified"] = True
        except (IndexError, KeyError, TypeError, ValueError) as exc:
            summary["independently_verified"] = False
            malformed.append({"path": str(path), "reason": str(exc)})
    return closed, summaries, malformed


def audit(source: Path, shards: list[Path], partitions: list[Path] | None = None) -> dict:
    source_payload = json.loads(source.read_text())
    source_digest = sha256(source)
    if int(source_payload["p"]) != 7 or int(source_payload["c_H"]) != -1:
        raise ValueError("source must be the p=7, c_H=-1 no-infinity orbit file")
    unsaturated = {
        index: orbit
        for index, orbit in enumerate(source_payload["orbits"])
        if any(int(value) != 32 for value in orbit["type_costs"].values())
    }
    expected = {
        (index, elevated)
        for index, orbit in unsaturated.items()
        for elevated in elevation_cases(orbit)
    }
    partition_closed, partition_summaries, partition_malformed = partition_audit(
        source_payload,
        source_digest,
        list(partitions or []),
    )

    attempts = {}
    shard_rows = []
    malformed = list(partition_malformed)
    for path in shards:
        payload = json.loads(path.read_text())
        shard_rows.append(
            {
                "path": str(path),
                "sha256": sha256(path),
                "shard_index": int(payload["shard_index"]),
                "shard_count": int(payload["shard_count"]),
                "completed_cases": int(payload["completed_cases"]),
                "status_counts": payload["status_counts"],
            }
        )
        if payload["source_sha256"] != source_digest:
            malformed.append({"path": str(path), "reason": "source hash"})
        for row in payload["rows"]:
            key = (
                int(row["orbit_index"]),
                tuple(int(value) for value in row["elevated_directions"]),
            )
            attempts.setdefault(key, []).append(row)

    final = {}
    conflicting = []
    for key, rows in attempts.items():
        statuses = {row["result"]["solver_status"] for row in rows}
        infeasible = [
            row for row in rows if row["result"]["solver_status"] == "INFEASIBLE"
        ]
        if infeasible:
            final[key] = infeasible[-1]
        else:
            final[key] = rows[-1]
        if "INFEASIBLE" in statuses and statuses & {"FEASIBLE", "OPTIMAL"}:
            conflicting.append({"key": key, "statuses": sorted(statuses)})

    for key, row in final.items():
        orbit_index, elevated = key
        if orbit_index not in unsaturated:
            malformed.append({"key": key, "reason": "orbit outside scope"})
            continue
        orbit = unsaturated[orbit_index]
        result = row["result"]
        if tuple(row["representative_vertices"]) != tuple(
            orbit["representative_vertices"]
        ):
            malformed.append({"key": key, "reason": "representative"})
        if elevated not in elevation_cases(orbit):
            malformed.append({"key": key, "reason": "elevation case"})
        if tuple(result["fixed_boundary"]) != tuple(
            orbit["representative_vertices"]
        ):
            malformed.append({"key": key, "reason": "result boundary"})
        if int(result["c_H"]) != -1:
            malformed.append({"key": key, "reason": "product sign"})
        if key not in partition_closed and not (
            result["solver_status"] == "INFEASIBLE"
            and result["finite_infeasibility_certificate"] is True
            and result["feasible"] is False
        ):
            malformed.append({"key": key, "reason": "not infeasible"})

    covered = set(final) | partition_closed
    missing = sorted(expected - covered)
    unexpected = sorted(set(final) - expected)
    unknown = sorted(
        key
        for key in expected & set(final)
        if key not in partition_closed
        and final[key]["result"]["solver_status"] == "UNKNOWN"
    )
    feasible = sorted(
        key
        for key in expected & set(final)
        if final[key]["result"]["solver_status"] in {"FEASIBLE", "OPTIMAL"}
    )
    catalog = catalog_audit()
    anti_isometry = p7_nonsquare_signed_permutation()
    proved = bool(
        len(unsaturated) == 518
        and sum(int(orbit["size"]) for orbit in unsaturated.values()) == 23520
        and len(expected) == 2408
        and len(covered) == len(expected)
        and not missing
        and not unexpected
        and not unknown
        and not feasible
        and not malformed
        and not conflicting
        and catalog["proved"]
        and anti_isometry["fixes_distinguished_edge"]
        and anti_isometry["signed_conference_anti_isometry"]
    )
    return {
        "experiment": "p7_no_infinity_unsaturated_audit",
        "status": "complete_independent_coverage_audit" if proved else "incomplete_audit",
        "proved": proved,
        "source": str(source),
        "source_sha256": source_digest,
        "shards": shard_rows,
        "partition_certificates": partition_summaries,
        "partition_closed_cases": len(partition_closed),
        "unsaturated_orbits": len(unsaturated),
        "unsaturated_boundary_size_sum": sum(
            int(orbit["size"]) for orbit in unsaturated.values()
        ),
        "expected_elevation_cases": len(expected),
        "final_case_rows": len(final),
        "missing": len(missing),
        "unexpected": len(unexpected),
        "unknown": len(unknown),
        "feasible": len(feasible),
        "malformed": malformed,
        "conflicting": conflicting,
        "catalog": catalog,
        "sign_transfer": {
            "fixes_distinguished_edge": anti_isometry["fixes_distinguished_edge"],
            "signed_conference_anti_isometry": anti_isometry[
                "signed_conference_anti_isometry"
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--shards", type=Path, nargs="+", required=True)
    parser.add_argument("--partitions", type=Path, nargs="*", default=[])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = audit(args.source, args.shards, args.partitions)
    rendered = json.dumps(out, indent=2)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
