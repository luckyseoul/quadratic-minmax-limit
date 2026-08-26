#!/usr/bin/env python3
"""Audit the generalized layered SCIP certificate for a p=5 circle class."""
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

from p5_full_shell_circle_crossing_resplit_scip import (  # noqa: E402
    quotient_patterns as second_row_quotient,
)
from p5_full_shell_circle_crossing_split_scip import (  # noqa: E402
    quotient_patterns as first_row_quotient,
)
from p5_full_shell_circle_degree_shards import degree_orbits  # noqa: E402
from p5_full_shell_circle_pattern_shards import internal_pattern_orbits  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def tuple_edges(rows: list[list[int]] | None) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(tuple(sorted(int(value) for value in edge)) for edge in (rows or [])))


class GeneralCircleAudit:
    def __init__(
        self,
        source: Path,
        orbit_index: int,
        workspace: Path,
        layer_manifest: Path | None = None,
    ):
        self.source_path = source
        self.source_sha = sha256(source)
        self.source = json.loads(source.read_text())
        self.orbit_index = orbit_index
        self.workspace = workspace
        orbit = self.source["orbits"][orbit_index]
        self.boundary = tuple(
            sorted(int(value) for value in orbit["representative_vertices"])
        )
        self.c_h = int(self.source["c_H"])
        self.edges = tuple(tuple(edge) for edge in geometry(5, "full")["edges"])
        self.edge_index = {edge: index for index, edge in enumerate(self.edges)}
        self.used: set[Path] = {source}
        self.layer_manifest_path = layer_manifest
        self.layer_manifest = None
        if layer_manifest is not None:
            self.layer_manifest = json.loads(layer_manifest.read_text())
            if self.layer_manifest.get("source_sha256") != self.source_sha:
                raise AssertionError("layer manifest source hash mismatch")
            if int(self.layer_manifest.get("orbit_index", -1)) != orbit_index:
                raise AssertionError("layer manifest orbit mismatch")
            self.used.add(layer_manifest)

    def load(self, path: Path, experiment: str | None = None) -> dict:
        payload = json.loads(path.read_text())
        if experiment is not None and payload.get("experiment") != experiment:
            raise AssertionError(f"wrong experiment in {path}")
        if payload.get("source_sha256") != self.source_sha:
            raise AssertionError(f"source hash mismatch in {path}")
        if int(payload.get("orbit_index", -1)) != self.orbit_index:
            raise AssertionError(f"orbit mismatch in {path}")
        if tuple(int(value) for value in payload.get("boundary", [])) != self.boundary:
            raise AssertionError(f"boundary mismatch in {path}")
        if int(payload.get("c_H", 0)) != self.c_h:
            raise AssertionError(f"c_H mismatch in {path}")
        self.used.add(path)
        return payload

    @staticmethod
    def solver_infeasible(row: dict) -> bool:
        return bool(
            row.get("solver_status") == "infeasible"
            and row.get("finite_infeasibility_certificate") is True
            and row.get("feasible") is False
        )

    @staticmethod
    def solver_unknown(row: dict) -> bool:
        return bool(
            row.get("solver_status") != "infeasible"
            and row.get("feasible") is not True
        )

    @staticmethod
    def assert_no_witness(rows: list[dict], label: str) -> None:
        if any(row.get("feasible") is True for row in rows):
            raise AssertionError(f"{label} contains a feasible witness")

    def direct_count(self, internal_edges: int) -> dict:
        path = self.workspace / f"p5_scip_screen_noinf_i{self.orbit_index}_a{internal_edges}.json"
        row = self.load(path, "p5_full_shell_fixed_boundary_scip")
        if int(row["boundary_internal_edges"]) != internal_edges:
            raise AssertionError("direct count has wrong internal-edge count")
        if row.get("boundary_cross_edges") is not None:
            raise AssertionError("direct count unexpectedly fixes crossing count")
        if row.get("feasible") is True:
            raise AssertionError("direct count found a witness")
        return row

    def cross_partition(self, internal_edges: int) -> tuple[list[int], list[dict]]:
        all_counts = list(range(0, 22 - internal_edges, 2))
        rows = []
        hard = []
        for cross_edges in all_counts:
            path = self.workspace / (
                f"p5_scip_crossscreen_noinf_i{self.orbit_index}_a{internal_edges}_b{cross_edges}.json"
            )
            row = self.load(path, "p5_full_shell_fixed_boundary_scip")
            if (
                int(row["boundary_internal_edges"]) != internal_edges
                or int(row["boundary_cross_edges"]) != cross_edges
                or row.get("fixed_internal_edges") is not None
                or row.get("boundary_cross_degrees") is not None
            ):
                raise AssertionError("bad direct crossing-count scope")
            if row.get("feasible") is True:
                raise AssertionError("crossing-count screen found a witness")
            if not self.solver_infeasible(row):
                hard.append(cross_edges)
            rows.append(row)
        return hard, rows

    def verify_pattern_layer(self, internal_edges: int, hard_counts: list[int]) -> dict:
        path = self.workspace / f"p5_scip_patterns_noinf_i{self.orbit_index}_a{internal_edges}.json"
        payload = self.load(path, "p5_full_shell_circle_pattern_scip")
        if int(payload["boundary_internal_edges"]) != internal_edges:
            raise AssertionError("pattern layer has wrong internal-edge count")
        if [int(value) for value in payload["cross_counts"]] != hard_counts:
            raise AssertionError("pattern layer does not match hard crossing counts")
        fresh = internal_pattern_orbits(self.boundary, internal_edges)
        if payload["internal_pattern_orbits"] != fresh:
            raise AssertionError("internal-pattern quotient mismatch")
        expected = {
            (pattern_index, cross_edges)
            for pattern_index in range(len(fresh))
            for cross_edges in hard_counts
        }
        rows = {
            (
                int(row["internal_pattern_orbit_index"]),
                int(row["boundary_cross_edges"]),
            ): row
            for row in payload["rows"]
        }
        if set(rows) != expected or int(payload["case_count"]) != len(expected):
            raise AssertionError("pattern row coverage mismatch")
        self.assert_no_witness(list(rows.values()), "pattern layer")
        unknown = {key for key, row in rows.items() if self.solver_unknown(row)}
        recorded_unknown = {
            tuple(int(value) for value in row) for row in payload["unknown_cases"]
        }
        if unknown != recorded_unknown:
            raise AssertionError("pattern unknown list mismatch")
        if any(
            not self.solver_infeasible(row) and key not in unknown
            for key, row in rows.items()
        ):
            raise AssertionError("pattern row has an unclassified status")
        return {
            "path": path,
            "payload": payload,
            "patterns": fresh,
            "rows": rows,
            "unknown": unknown,
        }

    def verify_merged_retry_shards(
        self, payload: dict, parent: dict, merged_path: Path
    ) -> dict[tuple[int, ...], dict]:
        records = payload.get("merged_shards")
        if not records:
            return {}
        shard_count = len(records)
        if {int(record["shard_index"]) for record in records} != set(
            range(shard_count)
        ):
            raise AssertionError(f"merged shard indices are incomplete in {merged_path}")
        ordered_parent = [
            tuple(int(value) for value in row)
            for row in parent["payload"]["unknown_branches"]
        ]
        if len(ordered_parent) != len(set(ordered_parent)):
            raise AssertionError("parent transition has duplicate unknown branch keys")
        parent_hash = sha256(parent["path"])
        covered: set[tuple[int, ...]] = set()
        shard_rows: dict[tuple[int, ...], dict] = {}
        for record in records:
            recorded_path = Path(record["path"])
            shard_path = (
                recorded_path
                if recorded_path.exists()
                else self.workspace / recorded_path.name
            )
            if sha256(shard_path) != record["sha256"]:
                raise AssertionError(f"merged shard hash mismatch for {shard_path}")
            shard = self.load(
                shard_path, "p5_full_shell_circle_branch_retry_scip"
            )
            shard_index = int(record["shard_index"])
            if (
                int(shard["shard_index"]) != shard_index
                or int(shard["shard_count"]) != shard_count
                or int(shard["all_input_unknown_branch_count"])
                != len(ordered_parent)
                or shard.get("split_source_sha256") != parent_hash
            ):
                raise AssertionError(f"bad merged shard scope in {shard_path}")
            expected = {
                key
                for position, key in enumerate(ordered_parent)
                if position % shard_count == shard_index
            }
            actual = {
                tuple(int(value) for value in row)
                for row in shard["input_unknown_branches"]
            }
            if actual != expected or len(actual) != len(shard["input_unknown_branches"]):
                raise AssertionError(f"merged shard partition mismatch in {shard_path}")
            rows = {
                tuple(int(value) for value in row["source_branch_key"]): row
                for row in shard["rows"]
            }
            if set(rows) != actual or len(rows) != len(shard["rows"]):
                raise AssertionError(f"merged shard row coverage mismatch in {shard_path}")
            self.assert_no_witness(list(rows.values()), f"merged shard {shard_path}")
            unknown = {key for key, row in rows.items() if self.solver_unknown(row)}
            recorded_unknown = {
                tuple(int(value) for value in row) for row in shard["unknown_branches"]
            }
            if unknown != recorded_unknown or any(
                not self.solver_infeasible(row) and key not in unknown
                for key, row in rows.items()
            ):
                raise AssertionError(f"merged shard status mismatch in {shard_path}")
            if (
                int(record["branch_count"]) != len(rows)
                or int(record["infeasible_count"])
                != sum(self.solver_infeasible(row) for row in rows.values())
                or int(record["unknown_count"]) != len(unknown)
                or int(record["feasible_count"]) != 0
            ):
                raise AssertionError(f"merged shard summary mismatch in {shard_path}")
            if covered & actual:
                raise AssertionError("merged retry shards overlap")
            covered |= actual
            shard_rows.update(rows)
        if covered != set(ordered_parent):
            raise AssertionError("merged retry shards do not exhaust the parent")
        return shard_rows

    def verify_generic_transition(self, parent: dict, path: Path) -> dict:
        payload = self.load(path)
        experiment = payload.get("experiment")
        supported = {
            "p5_full_shell_circle_crossing_resplit_scip",
            "p5_full_shell_circle_crossing_recursive_scip",
            "p5_full_shell_circle_branch_retry_scip",
        }
        if experiment not in supported:
            raise AssertionError(f"unsupported recursive transition in {path}")
        if experiment == "p5_full_shell_circle_crossing_resplit_scip":
            parent_hash = payload.get("parent_split_source_sha256")
        elif experiment == "p5_full_shell_circle_crossing_recursive_scip":
            parent_hash = payload.get("parent_source_sha256")
        else:
            parent_hash = payload.get("split_source_sha256")
        if parent_hash != sha256(parent["path"]):
            raise AssertionError(f"transition parent hash mismatch in {path}")
        input_unknown = {
            tuple(int(value) for value in row)
            for row in payload["input_unknown_branches"]
        }
        if input_unknown != parent["unknown"]:
            raise AssertionError(f"transition input mismatch in {path}")

        if experiment == "p5_full_shell_circle_branch_retry_scip":
            shard_rows = self.verify_merged_retry_shards(payload, parent, path)
            rows = {
                tuple(int(value) for value in row["source_branch_key"]): row
                for row in payload["rows"]
            }
            if set(rows) != input_unknown or len(rows) != len(payload["rows"]):
                raise AssertionError(f"retry row coverage mismatch in {path}")
            self.assert_no_witness(list(rows.values()), f"retry transition {path}")
            unknown = {key for key, row in rows.items() if self.solver_unknown(row)}
            recorded = {
                tuple(int(value) for value in row)
                for row in payload["unknown_branches"]
            }
            if unknown != recorded:
                raise AssertionError(f"retry unknown mismatch in {path}")
            if shard_rows and rows != shard_rows:
                raise AssertionError(f"merged rows differ from source shards in {path}")
            return {
                "path": path,
                "payload": payload,
                "flat_rows": rows,
                "unknown": unknown,
                "transition_type": "retry",
            }

        profiles = {
            tuple(int(value) for value in row["parent_branch_key"]): row
            for row in payload["profiles"]
        }
        if set(profiles) != input_unknown or len(profiles) != len(payload["profiles"]):
            raise AssertionError(f"partition profile coverage mismatch in {path}")
        rows_by_parent: dict[tuple[int, ...], dict[int, dict]] = {
            key: {} for key in input_unknown
        }
        flat_rows = {}
        for row in payload["rows"]:
            if experiment == "p5_full_shell_circle_crossing_resplit_scip":
                parent_key = (
                    *tuple(int(value) for value in row["source_profile_key"]),
                    int(row["parent_representative_index"]),
                )
                child_key = (*parent_key, int(row["representative_index"]))
            else:
                child_key = tuple(int(value) for value in row["source_branch_key"])
                parent_key = child_key[:-1]
            index = child_key[-1]
            if parent_key not in rows_by_parent or index in rows_by_parent[parent_key]:
                raise AssertionError(f"duplicate or stray partition row in {path}")
            rows_by_parent[parent_key][index] = row
            flat_rows[child_key] = row
        if len(flat_rows) != len(payload["rows"]):
            raise AssertionError(f"duplicate child key in {path}")

        unknown = set()
        for parent_key, metadata in profiles.items():
            parent_row = parent["flat_rows"][parent_key]
            fixed = self.fixed_internal_indices(parent_row)
            degrees = self.row_degrees(parent_row)
            inherited = self.required_indices(parent_row)
            vertex = int(metadata["enumerated_crossing_vertex"])
            representatives, weights, raw_total, stabilizer = second_row_quotient(
                self.boundary, fixed, degrees, inherited, vertex
            )
            if (
                int(metadata["crossing_raw_pattern_count"]) != raw_total
                or int(metadata["crossing_stabilizer_size"]) != stabilizer
                or int(metadata["crossing_representative_count"])
                != len(representatives)
                or set(rows_by_parent[parent_key]) != set(range(len(representatives)))
            ):
                raise AssertionError(f"recursive quotient mismatch in {path}")
            for index, representative in enumerate(representatives):
                row = rows_by_parent[parent_key][index]
                expected_required = tuple(sorted(set(inherited) | set(representative)))
                if (
                    int(row["orbit_weight"]) != weights[index]
                    or self.required_indices(row) != expected_required
                ):
                    raise AssertionError(f"recursive branch mismatch in {path}")
                if row.get("feasible") is True:
                    raise AssertionError(f"recursive transition found witness in {path}")
                if self.solver_unknown(row):
                    unknown.add((*parent_key, index))
        recorded = {
            tuple(int(value) for value in row) for row in payload["unknown_branches"]
        }
        if unknown != recorded:
            raise AssertionError(f"recursive unknown mismatch in {path}")
        return {
            "path": path,
            "payload": payload,
            "flat_rows": flat_rows,
            "unknown": unknown,
            "transition_type": "partition",
        }

    def verify_degree_layer(self, internal_edges: int, pattern: dict) -> dict | None:
        if not pattern["unknown"]:
            return None
        path = self.workspace / f"p5_scip_degrees_noinf_i{self.orbit_index}_a{internal_edges}.json"
        payload = self.load(path, "p5_full_shell_circle_degree_scip")
        if payload.get("pattern_source_sha256") != sha256(pattern["path"]):
            raise AssertionError("degree layer references wrong pattern file")
        input_unknown = {
            tuple(int(value) for value in row)
            for row in payload["input_unknown_pattern_cases"]
        }
        if input_unknown != pattern["unknown"]:
            raise AssertionError("degree layer input mismatch")
        cubes = {}
        expected = {}
        for cube in payload["cube_records"]:
            pattern_index = int(cube["internal_pattern_orbit_index"])
            cross_edges = int(cube["boundary_cross_edges"])
            indices = tuple(
                int(value)
                for value in pattern["patterns"][pattern_index]["representative_indices"]
            )
            fresh = degree_orbits(self.boundary, indices, cross_edges)
            if cube["degree_orbits"] != fresh:
                raise AssertionError("crossing-degree quotient mismatch")
            cubes[(pattern_index, cross_edges)] = cube
            for degree_index, degree_row in enumerate(fresh):
                expected[(pattern_index, cross_edges, degree_index)] = tuple(
                    int(value) for value in degree_row["representative"]
                )
        if set(cubes) != pattern["unknown"]:
            raise AssertionError("degree cubes do not cover pattern unknowns")
        rows = {
            (
                int(row["internal_pattern_orbit_index"]),
                int(row["boundary_cross_edges"]),
                int(row["cross_degree_orbit_index"]),
            ): row
            for row in payload["rows"]
        }
        if set(rows) != set(expected) or int(payload["case_count"]) != len(expected):
            raise AssertionError("degree row coverage mismatch")
        for key, row in rows.items():
            observed = tuple(
                int(row["boundary_cross_degrees"][str(vertex)])
                for vertex in self.boundary
            )
            if observed != expected[key]:
                raise AssertionError("degree row profile mismatch")
        self.assert_no_witness(list(rows.values()), "degree layer")
        unknown = {key for key, row in rows.items() if self.solver_unknown(row)}
        if unknown != {
            tuple(int(value) for value in row) for row in payload["unknown_cases"]
        }:
            raise AssertionError("degree unknown list mismatch")
        return {
            "path": path,
            "payload": payload,
            "cubes": cubes,
            "rows": rows,
            "unknown": unknown,
        }

    def audit_hard_count_manifest(
        self, internal_edges: int, hard_counts: list[int], transition_names: list[str]
    ) -> dict:
        pattern = self.verify_pattern_layer(internal_edges, hard_counts)
        degree = self.verify_degree_layer(internal_edges, pattern)
        transitions = []
        retry = first = None
        if degree is None:
            if transition_names:
                raise AssertionError("manifest has transitions after a closed pattern layer")
            unresolved = set()
        else:
            retry = self.verify_retry_layer(internal_edges, degree)
            unresolved = retry["unknown"] if retry is not None else degree["unknown"]
            first = self.verify_first_split(internal_edges, degree, unresolved, retry)
            if first is None:
                if transition_names:
                    raise AssertionError("manifest has transitions after a closed degree layer")
                unresolved = set()
            else:
                parent = first
                for name in transition_names:
                    parent = self.verify_generic_transition(
                        parent, self.workspace / str(name)
                    )
                    transitions.append(parent)
                unresolved = parent["unknown"]
        if unresolved:
            raise AssertionError(
                f"manifest chain for internal count {internal_edges} leaves "
                f"{len(unresolved)} unresolved branches"
            )
        return {
            "internal_edges": internal_edges,
            "hard_cross_counts": hard_counts,
            "internal_pattern_orbits": len(pattern["patterns"]),
            "pattern_cases": len(pattern["rows"]),
            "pattern_infeasible": sum(
                self.solver_infeasible(row) for row in pattern["rows"].values()
            ),
            "degree_profiles": 0 if degree is None else len(degree["rows"]),
            "initial_degree_infeasible": (
                0
                if degree is None
                else sum(
                    self.solver_infeasible(row) for row in degree["rows"].values()
                )
            ),
            "degree_retry_cases": 0 if retry is None else len(retry["rows"]),
            "first_split_branches": (
                0 if first is None else len(first["flat_rows"])
            ),
            "recursive_transitions": [
                {
                    "file": transition["path"].name,
                    "type": transition["transition_type"],
                    "rows": len(transition["flat_rows"]),
                    "unknown_after": len(transition["unknown"]),
                }
                for transition in transitions
            ],
        }

    def verify_retry_layer(self, internal_edges: int, degree: dict) -> dict | None:
        path = self.workspace / f"p5_scip_degree_retry_noinf_i{self.orbit_index}_a{internal_edges}.json"
        if not path.exists():
            return None
        payload = self.load(path, "p5_full_shell_circle_degree_retry_scip")
        if payload.get("degree_source_sha256") != sha256(degree["path"]):
            raise AssertionError("retry layer references wrong degree file")
        input_unknown = {
            tuple(int(value) for value in row) for row in payload["input_unknown_cases"]
        }
        if input_unknown != degree["unknown"]:
            raise AssertionError("retry input does not equal degree unknowns")
        rows = {
            (
                int(row["internal_pattern_orbit_index"]),
                int(row["boundary_cross_edges"]),
                int(row["cross_degree_orbit_index"]),
            ): row
            for row in payload["rows"]
        }
        if set(rows) != input_unknown:
            raise AssertionError("retry row coverage mismatch")
        self.assert_no_witness(list(rows.values()), "degree retry")
        unknown = {key for key, row in rows.items() if self.solver_unknown(row)}
        if unknown != {
            tuple(int(value) for value in row) for row in payload["unknown_cases"]
        }:
            raise AssertionError("retry unknown list mismatch")
        return {"path": path, "payload": payload, "rows": rows, "unknown": unknown}

    def fixed_internal_indices(self, row: dict) -> tuple[int, ...]:
        return tuple(sorted(self.edge_index[edge] for edge in tuple_edges(row["fixed_internal_edges"])))

    def row_degrees(self, row: dict) -> tuple[int, ...]:
        return tuple(
            int(row["boundary_cross_degrees"][str(vertex)])
            for vertex in self.boundary
        )

    def required_indices(self, row: dict) -> tuple[int, ...]:
        return tuple(sorted(self.edge_index[edge] for edge in tuple_edges(row["required_cross_edges"])))

    def verify_first_split(
        self,
        internal_edges: int,
        degree: dict,
        input_unknown: set[tuple[int, int, int]],
        retry: dict | None,
    ) -> dict | None:
        if not input_unknown:
            return None
        path = self.workspace / f"p5_scip_crosssplit_noinf_i{self.orbit_index}_a{internal_edges}.json"
        payload = self.load(path, "p5_full_shell_circle_crossing_split_scip")
        if payload.get("degree_source_sha256") != sha256(degree["path"]):
            raise AssertionError("first split references wrong degree file")
        if retry is None:
            if payload.get("retry_source") is not None:
                raise AssertionError("unexpected retry parent")
        elif payload.get("retry_source_sha256") != sha256(retry["path"]):
            raise AssertionError("first split references wrong retry file")
        if {
            tuple(int(value) for value in row)
            for row in payload["input_unknown_profiles"]
        } != input_unknown:
            raise AssertionError("first split input mismatch")
        profiles = {
            tuple(int(value) for value in row["profile_key"]): row
            for row in payload["profiles"]
        }
        if set(profiles) != input_unknown:
            raise AssertionError("first split profile metadata mismatch")
        rows_by_profile: dict[tuple[int, int, int], dict[int, dict]] = {
            key: {} for key in input_unknown
        }
        for row in payload["rows"]:
            key = tuple(int(value) for value in row["source_profile_key"])
            index = int(row["representative_index"])
            if index in rows_by_profile[key]:
                raise AssertionError("duplicate first-row representative")
            rows_by_profile[key][index] = row
        unknown = set()
        for key, metadata in profiles.items():
            base = degree["rows"][key]
            fixed = self.fixed_internal_indices(base)
            degrees = self.row_degrees(base)
            vertex = int(metadata["enumerated_crossing_vertex"])
            representatives, weights, raw_total, stabilizer = first_row_quotient(
                self.boundary, fixed, degrees, vertex
            )
            if (
                int(metadata["crossing_raw_pattern_count"]) != raw_total
                or int(metadata["crossing_stabilizer_size"]) != stabilizer
                or int(metadata["crossing_representative_count"]) != len(representatives)
                or set(rows_by_profile[key]) != set(range(len(representatives)))
            ):
                raise AssertionError("first-row quotient metadata mismatch")
            for index, representative in enumerate(representatives):
                row = rows_by_profile[key][index]
                if int(row["orbit_weight"]) != weights[index]:
                    raise AssertionError("first-row orbit weight mismatch")
                if self.required_indices(row) != tuple(sorted(representative)):
                    raise AssertionError("first-row required-edge mismatch")
                if row.get("feasible") is True:
                    raise AssertionError("first split found a witness")
                if self.solver_unknown(row):
                    unknown.add((*key, index))
        if unknown != {
            tuple(int(value) for value in row) for row in payload["unknown_branches"]
        }:
            raise AssertionError("first split unknown list mismatch")
        return {
            "path": path,
            "payload": payload,
            "profiles": profiles,
            "rows": rows_by_profile,
            "unknown": unknown,
            "flat_rows": {
                (*key, index): row
                for key, indexed in rows_by_profile.items()
                for index, row in indexed.items()
            },
        }

    def verify_second_split(
        self, internal_edges: int, degree: dict, first: dict
    ) -> dict | None:
        if not first["unknown"]:
            return None
        path = self.workspace / f"p5_scip_crossresplit_noinf_i{self.orbit_index}_a{internal_edges}.json"
        payload = self.load(path, "p5_full_shell_circle_crossing_resplit_scip")
        if payload.get("degree_source_sha256") != sha256(degree["path"]):
            raise AssertionError("second split references wrong degree file")
        if payload.get("parent_split_source_sha256") != sha256(first["path"]):
            raise AssertionError("second split references wrong first split")
        if {
            tuple(int(value) for value in row)
            for row in payload["input_unknown_branches"]
        } != first["unknown"]:
            raise AssertionError("second split input mismatch")
        profiles = {
            tuple(int(value) for value in row["parent_branch_key"]): row
            for row in payload["profiles"]
        }
        if set(profiles) != first["unknown"]:
            raise AssertionError("second split profile metadata mismatch")
        rows_by_parent: dict[tuple[int, ...], dict[int, dict]] = {
            key: {} for key in first["unknown"]
        }
        for row in payload["rows"]:
            key = (
                *tuple(int(value) for value in row["source_profile_key"]),
                int(row["parent_representative_index"]),
            )
            index = int(row["representative_index"])
            if index in rows_by_parent[key]:
                raise AssertionError("duplicate second-row representative")
            rows_by_parent[key][index] = row
        unknown = set()
        for parent_key, metadata in profiles.items():
            profile_key = parent_key[:3]
            base = degree["rows"][profile_key]
            parent = first["rows"][profile_key][parent_key[3]]
            fixed = self.fixed_internal_indices(base)
            degrees = self.row_degrees(base)
            inherited = self.required_indices(parent)
            vertex = int(metadata["enumerated_crossing_vertex"])
            representatives, weights, raw_total, stabilizer = second_row_quotient(
                self.boundary, fixed, degrees, inherited, vertex
            )
            if (
                int(metadata["crossing_raw_pattern_count"]) != raw_total
                or int(metadata["crossing_stabilizer_size"]) != stabilizer
                or int(metadata["crossing_representative_count"]) != len(representatives)
                or set(rows_by_parent[parent_key]) != set(range(len(representatives)))
            ):
                raise AssertionError("second-row quotient metadata mismatch")
            for index, representative in enumerate(representatives):
                row = rows_by_parent[parent_key][index]
                expected_required = tuple(sorted(set(inherited) | set(representative)))
                if (
                    int(row["orbit_weight"]) != weights[index]
                    or self.required_indices(row) != expected_required
                ):
                    raise AssertionError("second-row branch mismatch")
                if row.get("feasible") is True:
                    raise AssertionError("second split found a witness")
                if self.solver_unknown(row):
                    unknown.add((*parent_key, index))
        if unknown != {
            tuple(int(value) for value in row) for row in payload["unknown_branches"]
        }:
            raise AssertionError("second split unknown list mismatch")
        return {
            "path": path,
            "payload": payload,
            "profiles": profiles,
            "rows": rows_by_parent,
            "unknown": unknown,
            "flat_rows": {
                (*key, index): row
                for key, indexed in rows_by_parent.items()
                for index, row in indexed.items()
            },
        }

    def verify_branch_retry(
        self, internal_edges: int, second: dict
    ) -> dict | None:
        if not second["unknown"]:
            return None
        path = self.workspace / f"p5_scip_branch_retry_noinf_i{self.orbit_index}_a{internal_edges}.json"
        payload = self.load(path, "p5_full_shell_circle_branch_retry_scip")
        if payload.get("split_source_sha256") != sha256(second["path"]):
            raise AssertionError("branch retry references wrong second split")
        if {
            tuple(int(value) for value in row)
            for row in payload["input_unknown_branches"]
        } != second["unknown"]:
            raise AssertionError("branch retry input mismatch")
        rows = {
            tuple(int(value) for value in row["source_branch_key"]): row
            for row in payload["rows"]
        }
        if set(rows) != second["unknown"]:
            raise AssertionError("branch retry row coverage mismatch")
        self.assert_no_witness(list(rows.values()), "branch retry")
        unknown = {key for key, row in rows.items() if self.solver_unknown(row)}
        if unknown != {
            tuple(int(value) for value in row) for row in payload["unknown_branches"]
        }:
            raise AssertionError("branch retry unknown list mismatch")
        return {"path": path, "payload": payload, "rows": rows, "unknown": unknown}

    def verify_final_leaves(
        self, internal_edges: int, branch_retry: dict
    ) -> list[dict]:
        if not branch_retry["unknown"]:
            return []
        paths = sorted(
            self.workspace.glob(
                f"p5_scip_final_noinf_i{self.orbit_index}_a{internal_edges}_*.json"
            )
        )
        leaves = [self.load(path, "p5_full_shell_fixed_boundary_scip") for path in paths]
        unused = set(range(len(leaves)))
        matched = []
        for key in sorted(branch_retry["unknown"]):
            parent = branch_retry["rows"][key]
            required = tuple_edges(parent["required_cross_edges"])
            candidates = [
                index
                for index in unused
                if tuple_edges(leaves[index].get("required_cross_edges")) == required
                and int(leaves[index]["boundary_internal_edges"]) == internal_edges
                and int(leaves[index]["boundary_cross_edges"])
                == int(parent["boundary_cross_edges"])
                and leaves[index].get("boundary_cross_degrees")
                == parent.get("boundary_cross_degrees")
            ]
            if len(candidates) != 1:
                raise AssertionError("final leaf does not uniquely match its branch")
            index = candidates[0]
            if not self.solver_infeasible(leaves[index]):
                raise AssertionError("final leaf is not an infeasibility certificate")
            unused.remove(index)
            matched.append({"branch_key": list(key), "file": paths[index].name})
        if unused:
            raise AssertionError("unmatched final leaf file")
        return matched

    def audit_hard_count(self, internal_edges: int, hard_counts: list[int]) -> dict:
        pattern = self.verify_pattern_layer(internal_edges, hard_counts)
        degree = self.verify_degree_layer(internal_edges, pattern)
        if degree is None:
            unresolved = set()
            retry = first = second = branch_retry = None
            final = []
        else:
            retry = self.verify_retry_layer(internal_edges, degree)
            unresolved = retry["unknown"] if retry is not None else degree["unknown"]
            first = self.verify_first_split(
                internal_edges, degree, unresolved, retry
            )
            if first is None:
                second = branch_retry = None
                final = []
            else:
                second = self.verify_second_split(internal_edges, degree, first)
                if second is None:
                    branch_retry = None
                    final = []
                else:
                    branch_retry = self.verify_branch_retry(internal_edges, second)
                    final = (
                        self.verify_final_leaves(internal_edges, branch_retry)
                        if branch_retry is not None
                        else []
                    )
        initial_degree_infeasible = (
            0
            if degree is None
            else sum(self.solver_infeasible(row) for row in degree["rows"].values())
        )
        return {
            "internal_edges": internal_edges,
            "hard_cross_counts": hard_counts,
            "internal_pattern_orbits": len(pattern["patterns"]),
            "pattern_cases": len(pattern["rows"]),
            "pattern_infeasible": sum(
                self.solver_infeasible(row) for row in pattern["rows"].values()
            ),
            "degree_profiles": 0 if degree is None else len(degree["rows"]),
            "initial_degree_infeasible": initial_degree_infeasible,
            "degree_retry_cases": 0 if retry is None else len(retry["rows"]),
            "first_split_profiles": 0 if first is None else len(first["profiles"]),
            "first_split_branches": (
                0
                if first is None
                else sum(len(rows) for rows in first["rows"].values())
            ),
            "second_split_parent_branches": (
                0 if second is None else len(second["profiles"])
            ),
            "second_split_branches": (
                0
                if second is None
                else sum(len(rows) for rows in second["rows"].values())
            ),
            "branch_retry_cases": (
                0 if branch_retry is None else len(branch_retry["rows"])
            ),
            "final_direct_leaves": final,
        }

    def audit(self) -> dict:
        if int(self.source["p"]) != 5 or self.c_h != -1:
            raise AssertionError("unexpected source scope")
        direct = {}
        hard_internal = []
        for internal_edges in range(16):
            row = self.direct_count(internal_edges)
            direct[internal_edges] = row
            if self.solver_infeasible(row):
                continue
            if not self.solver_unknown(row):
                raise AssertionError("internal-count screen has unclassified status")
            hard_internal.append(internal_edges)
        levels = []
        manifest_chains = (
            None if self.layer_manifest is None else self.layer_manifest.get("chains", {})
        )
        for internal_edges in hard_internal:
            hard_cross, cross_rows = self.cross_partition(internal_edges)
            if not hard_cross:
                raise AssertionError("hard internal count had no hard crossing count")
            if manifest_chains is None:
                levels.append(self.audit_hard_count(internal_edges, hard_cross))
            else:
                key = str(internal_edges)
                if key not in manifest_chains:
                    raise AssertionError(f"manifest omits hard internal count {key}")
                levels.append(
                    self.audit_hard_count_manifest(
                        internal_edges, hard_cross, list(manifest_chains[key])
                    )
                )
            if any(
                not self.solver_infeasible(row)
                and int(row["boundary_cross_edges"]) not in hard_cross
                for row in cross_rows
            ):
                raise AssertionError("cross partition has an unclassified row")
        covered = sorted(
            [a for a, row in direct.items() if self.solver_infeasible(row)]
            + [row["internal_edges"] for row in levels]
        )
        proved = covered == list(range(16))
        artifacts = [
            {
                "file": path.name,
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in sorted(self.used, key=lambda value: value.name)
        ]
        relation = (
            "internal"
            if 0 in self.boundary and 1 in self.boundary
            else "crossing"
            if (0 in self.boundary) != (1 in self.boundary)
            else "outside"
        )
        return {
            "experiment": "p5_size6_general_circle_scip_audit",
            "status": "exact_generalized_layered_finite_coverage_audit",
            "proved": proved,
            "scope": {
                "p": 5,
                "c_H": self.c_h,
                "orbit_index": self.orbit_index,
                "boundary": list(self.boundary),
                "mandatory_edge_relation": relation,
                "source": str(self.source_path),
                "source_sha256": self.source_sha,
            },
            "directly_closed_internal_counts": sorted(
                a for a, row in direct.items() if self.solver_infeasible(row)
            ),
            "layered_internal_counts": hard_internal,
            "covered_internal_edge_counts": covered,
            "levels": levels,
            "used_artifact_count": len(artifacts),
            "used_artifacts": artifacts,
            "solver_claim_limit": (
                "The audit reconstructs every finite quotient and verifies every "
                "recorded SCIP infeasibility status; it is not a separately "
                "checkable SCIP proof-log certificate."
            ),
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--orbit-index", type=int, required=True)
    parser.add_argument("--workspace", type=Path, default=Path("/tmp"))
    parser.add_argument("--layer-manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = GeneralCircleAudit(
        args.source, args.orbit_index, args.workspace, args.layer_manifest
    ).audit()
    atomic_write(args.output, result)
    print(
        json.dumps(
            {
                "experiment": result["experiment"],
                "proved": result["proved"],
                "covered_internal_edge_counts": result[
                    "covered_internal_edge_counts"
                ],
                "used_artifact_count": result["used_artifact_count"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
