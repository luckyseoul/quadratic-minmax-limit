#!/usr/bin/env python3
"""Audit the layered SCIP certificate for one hard p=5 six-point circle.

The certificate partitions first by the number of boundary-internal edges,
then by the number of crossing edges, an internal-pattern orbit, and a
crossing-degree orbit.  Remaining SCIP timeouts may be partitioned further
by exact outside-neighbour sets of selected boundary vertices.  This audit
rebuilds the finite group quotients, verifies every shard range and branch,
and records hashes for every artifact actually used in the proof.

The current certificate layout is for the normalized circle whose mandatory
edge (0, 1) is internal.  The solver model itself is more general, but this
auditor deliberately rejects a different geometry instead of silently
claiming coverage it has not reconstructed.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import os
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p5_full_shell_circle_degree_shards import degree_orbits  # noqa: E402
from p5_full_shell_circle_pattern_shards import internal_pattern_orbits  # noqa: E402
from p5_full_shell_fixed_boundary_cpsat import boundary_edge_stabilizers  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402


EXPECTED_SOURCE_SHA256 = "a3ec787dbbb6ad08213573088274a744e086b64cb5998cea3a3ed77dd5040834"


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


def edge_tuple(edge: list[int] | tuple[int, int]) -> tuple[int, int]:
    a, b = (int(value) for value in edge)
    return tuple(sorted((a, b)))


def edge_set(rows: list[list[int]] | None) -> frozenset[tuple[int, int]]:
    return frozenset(edge_tuple(row) for row in (rows or []))


def degree_tuple(payload: dict, boundary: tuple[int, ...]) -> tuple[int, ...]:
    degrees = payload.get("boundary_cross_degrees")
    if not isinstance(degrees, dict):
        raise AssertionError("missing crossing-degree profile")
    return tuple(int(degrees[str(vertex)]) for vertex in boundary)


@dataclass(frozen=True, order=True)
class Profile:
    internal_edges: int
    fixed_internal_edges: tuple[tuple[int, int], ...]
    cross_edges: int
    crossing_degrees: tuple[int, ...]


@dataclass(frozen=True)
class Record:
    path: Path
    payload: dict


def profile_from_payload(payload: dict, boundary: tuple[int, ...]) -> Profile:
    fixed = tuple(sorted(edge_set(payload.get("fixed_internal_edges"))))
    return Profile(
        int(payload["boundary_internal_edges"]),
        fixed,
        int(payload["boundary_cross_edges"]),
        degree_tuple(payload, boundary),
    )


class CircleAudit:
    def __init__(self, source: Path, orbit_index: int, workspace: Path):
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
        self.data = geometry(5, "full")
        self.edges = tuple(edge_tuple(edge) for edge in self.data["edges"])
        self.edge_index = {edge: index for index, edge in enumerate(self.edges)}
        self.boundary_set = set(self.boundary)
        self.cross_indices = tuple(
            index
            for index, (a, b) in enumerate(self.edges)
            if (a in self.boundary_set) != (b in self.boundary_set)
        )
        self.cross_groups = {
            vertex: frozenset(
                index for index in self.cross_indices if vertex in self.edges[index]
            )
            for vertex in self.boundary
        }
        self.group_to_vertex = {
            group: vertex for vertex, group in self.cross_groups.items()
        }
        self.used_paths: set[Path] = {source}
        self.records_by_profile: dict[Profile, list[Record]] = defaultdict(list)
        self.branch_cache: dict[
            tuple[Profile, frozenset[tuple[int, int]]], tuple[bool, dict]
        ] = {}
        self.active_branches: set[
            tuple[Profile, frozenset[tuple[int, int]]]
        ] = set()

    def check_scope(self, payload: dict, path: Path, require_hash: bool = True) -> None:
        if int(payload.get("orbit_index", -1)) != self.orbit_index:
            raise AssertionError(f"orbit mismatch in {path}")
        if tuple(int(value) for value in payload.get("boundary", [])) != self.boundary:
            raise AssertionError(f"boundary mismatch in {path}")
        if int(payload.get("c_H", 0)) != self.c_h:
            raise AssertionError(f"c_H mismatch in {path}")
        if require_hash and payload.get("source_sha256") != self.source_sha:
            raise AssertionError(f"source hash mismatch in {path}")

    def load_scip_records(self, globs: list[str]) -> None:
        paths = sorted({path for pattern in globs for path in self.workspace.glob(pattern)})
        if not paths:
            raise AssertionError("no SCIP records matched")
        for path in paths:
            payload = json.loads(path.read_text())
            if payload.get("experiment") != "p5_full_shell_fixed_boundary_scip":
                raise AssertionError(f"wrong SCIP experiment in {path}")
            self.check_scope(payload, path)
            profile = profile_from_payload(payload, self.boundary)
            self.records_by_profile[profile].append(Record(path, payload))

    def validate_required_edges(
        self, profile: Profile, required: frozenset[tuple[int, int]]
    ) -> None:
        if any(self.edge_index[edge] not in self.cross_indices for edge in required):
            raise AssertionError("a branch required a non-crossing edge")
        for vertex, degree in zip(self.boundary, profile.crossing_degrees):
            incident = sum(vertex in edge for edge in required)
            if incident > degree:
                raise AssertionError("required branch edges exceed a boundary degree")

    @lru_cache(maxsize=None)
    def crossing_representatives(
        self,
        profile: Profile,
        required: frozenset[tuple[int, int]],
        enumerated_vertices: tuple[int, ...],
    ) -> tuple[tuple[tuple[int, ...], ...], tuple[int, ...], int, int]:
        """Rebuild the exact crossing-pattern quotient used by SCIP shards."""
        self.validate_required_edges(profile, required)
        if not enumerated_vertices or len(set(enumerated_vertices)) != len(
            enumerated_vertices
        ):
            raise AssertionError("bad enumerated boundary-vertex set")
        if not set(enumerated_vertices) <= self.boundary_set:
            raise AssertionError("enumerated vertex is outside the boundary")
        degree_by_vertex = dict(zip(self.boundary, profile.crossing_degrees))
        choices = [
            tuple(
                itertools.combinations(
                    tuple(self.cross_groups[vertex]), degree_by_vertex[vertex]
                )
            )
            for vertex in enumerated_vertices
        ]
        raw_total = math.prod(len(values) for values in choices)
        fixed_internal = {
            self.edge_index[edge] for edge in profile.fixed_internal_edges
        }
        required_indices = {self.edge_index[edge] for edge in required}
        enumerated_set = set(enumerated_vertices)
        stabilizers = []
        for permutation in boundary_edge_stabilizers(self.boundary):
            if {permutation[index] for index in fixed_internal} != fixed_internal:
                continue
            valid = True
            for vertex, group in self.cross_groups.items():
                image = frozenset(permutation[index] for index in group)
                target = self.group_to_vertex.get(image)
                if (
                    target is None
                    or degree_by_vertex[target] != degree_by_vertex[vertex]
                    or ((vertex in enumerated_set) != (target in enumerated_set))
                ):
                    valid = False
                    break
            if valid and required_indices and {
                permutation[index] for index in required_indices
            } != required_indices:
                valid = False
            if valid:
                stabilizers.append(permutation)
        if not stabilizers:
            raise AssertionError("crossing quotient lost the identity")

        weights: dict[tuple[int, ...], int] = {}
        for per_vertex in itertools.product(*choices):
            pattern = tuple(sorted(index for group in per_vertex for index in group))
            canonical = min(
                tuple(sorted(permutation[index] for index in pattern))
                for permutation in stabilizers
            )
            weights[canonical] = weights.get(canonical, 0) + 1
        if sum(weights.values()) != raw_total:
            raise AssertionError("crossing quotient lost raw patterns")
        representatives = tuple(sorted(weights))
        return (
            representatives,
            tuple(weights[representative] for representative in representatives),
            raw_total,
            len(stabilizers),
        )

    @staticmethod
    def record_required(record: Record) -> frozenset[tuple[int, int]]:
        return edge_set(record.payload.get("required_cross_edges"))

    def direct_leaf(
        self,
        records: list[Record],
        required: frozenset[tuple[int, int]],
    ) -> Record | None:
        candidates = []
        for record in records:
            payload = record.payload
            if payload.get("status") == "exact_scip_crossing_orbit_shard":
                continue
            if self.record_required(record) != required:
                continue
            if payload.get("fixed_cross_edges") not in (None, []):
                continue
            if (
                payload.get("solver_status") == "infeasible"
                and payload.get("finite_infeasibility_certificate") is True
                and payload.get("feasible") is False
            ):
                candidates.append(record)
        return min(candidates, key=lambda record: record.path.name) if candidates else None

    def enumeration_groups(
        self,
        records: list[Record],
        required: frozenset[tuple[int, int]],
    ) -> list[list[Record]]:
        groups: dict[tuple, list[Record]] = defaultdict(list)
        for record in records:
            payload = record.payload
            if payload.get("status") != "exact_scip_crossing_orbit_shard":
                continue
            if self.record_required(record) != required:
                continue
            # The first SCIP crossing enumerator predated the explicit
            # metadata field and always enumerated every boundary vertex.
            enumerated = payload.get("enumerated_crossing_vertices")
            if enumerated is None:
                enumerated = self.boundary
            key = (
                tuple(int(value) for value in enumerated),
                int(payload["crossing_representative_count"]),
                int(payload["crossing_raw_pattern_count"]),
                int(payload["crossing_stabilizer_size"]),
            )
            groups[key].append(record)
        return [groups[key] for key in sorted(groups)]

    def close_enumeration_group(
        self,
        profile: Profile,
        required: frozenset[tuple[int, int]],
        records: list[Record],
    ) -> tuple[bool, dict]:
        first = records[0].payload
        first_vertices = first.get("enumerated_crossing_vertices")
        if first_vertices is None:
            first_vertices = self.boundary
        vertices = tuple(int(value) for value in first_vertices)
        representatives, weights, raw_total, stabilizer_size = (
            self.crossing_representatives(profile, required, vertices)
        )
        representative_count = len(representatives)
        if (
            int(first["crossing_representative_count"]) != representative_count
            or int(first["crossing_raw_pattern_count"]) != raw_total
            or int(first["crossing_stabilizer_size"]) != stabilizer_size
        ):
            return False, {"reason": "crossing quotient metadata mismatch"}

        states: dict[int, list[tuple[str, frozenset[tuple[int, int]], Record]]] = (
            defaultdict(list)
        )
        for record in records:
            payload = record.payload
            payload_vertices = payload.get("enumerated_crossing_vertices")
            if payload_vertices is None:
                payload_vertices = self.boundary
            if (
                tuple(int(value) for value in payload_vertices)
                != vertices
                or int(payload["crossing_representative_count"])
                != representative_count
                or int(payload["crossing_raw_pattern_count"]) != raw_total
                or int(payload["crossing_stabilizer_size"]) != stabilizer_size
            ):
                return False, {"reason": "inconsistent crossing shard family"}
            start = int(payload["crossing_start"])
            stop = int(payload["crossing_stop"])
            unknown_rows = {
                int(row["representative_index"]): row
                for row in payload.get("unknown_cases", [])
            }
            if not 0 <= start < stop <= representative_count:
                return False, {"reason": "invalid or empty shard range"}
            if int(payload["crossing_attempted"]) != stop - start:
                return False, {"reason": "incomplete shard attempt count"}
            if int(payload["unknown_case_count"]) != len(unknown_rows):
                return False, {"reason": "unknown-count mismatch"}
            if payload.get("feasible") is True:
                return False, {"reason": "feasible crossing representative"}
            inferred_infeasible = 0
            for index in range(start, stop):
                if index in unknown_rows:
                    row = unknown_rows[index]
                    pattern_indices = representatives[index]
                    selected = edge_set(row["selected_crossing_edges"])
                    expected = frozenset(self.edges[value] for value in pattern_indices)
                    if selected != expected or int(row["orbit_weight"]) != weights[index]:
                        return False, {"reason": "unknown representative mismatch"}
                    for vertex in vertices:
                        expected_degree = profile.crossing_degrees[
                            self.boundary.index(vertex)
                        ]
                        if sum(vertex in edge for edge in selected) != expected_degree:
                            return False, {"reason": "enumerated row is not exact"}
                    states[index].append(("unknown", required | selected, record))
                else:
                    states[index].append(("infeasible", required, record))
                    inferred_infeasible += 1
            if inferred_infeasible != int(payload["infeasible_representatives"]):
                return False, {"reason": "infeasible representative mismatch"}

        child_summaries = []
        covered_files = set()
        for index in range(representative_count):
            options = states.get(index, [])
            exact = next((option for option in options if option[0] == "infeasible"), None)
            if exact is not None:
                covered_files.add(exact[2].path)
                continue
            unknown = [option for option in options if option[0] == "unknown"]
            if not unknown:
                return False, {"reason": f"uncovered representative {index}"}
            children = {option[1] for option in unknown}
            if len(children) != 1:
                return False, {"reason": f"inconsistent branch at representative {index}"}
            child_required = next(iter(children))
            closed, child = self.close_branch(profile, child_required)
            if not closed:
                return False, {
                    "reason": f"unclosed child representative {index}",
                    "child": child,
                }
            covered_files.update(option[2].path for option in unknown)
            child_summaries.append(
                {
                    "representative_index": index,
                    "required_cross_edges": [list(edge) for edge in sorted(child_required)],
                    "proof": child,
                }
            )
        self.used_paths.update(covered_files)
        return True, {
            "method": "exact_crossing_orbit_partition",
            "enumerated_crossing_vertices": list(vertices),
            "raw_patterns": raw_total,
            "stabilizer_size": stabilizer_size,
            "representatives": representative_count,
            "child_branches": child_summaries,
            "shard_files": sorted(path.name for path in covered_files),
        }

    def close_branch(
        self,
        profile: Profile,
        required: frozenset[tuple[int, int]],
    ) -> tuple[bool, dict]:
        key = (profile, required)
        if key in self.branch_cache:
            return self.branch_cache[key]
        if key in self.active_branches:
            return False, {"reason": "cyclic branch proof"}
        self.active_branches.add(key)
        self.validate_required_edges(profile, required)
        records = self.records_by_profile.get(profile, [])
        leaf = self.direct_leaf(records, required)
        if leaf is not None:
            self.used_paths.add(leaf.path)
            result = (
                True,
                {
                    "method": "direct_scip_infeasible",
                    "required_cross_edges": [list(edge) for edge in sorted(required)],
                    "file": leaf.path.name,
                },
            )
        else:
            result = (False, {"reason": "no complete crossing partition"})
            for group in self.enumeration_groups(records, required):
                candidate = self.close_enumeration_group(profile, required, group)
                if candidate[0]:
                    result = candidate
                    break
                result = candidate
        self.active_branches.remove(key)
        self.branch_cache[key] = result
        return result

    def validate_count_file(self, path: Path, expected_a: int) -> dict:
        payload = json.loads(path.read_text())
        self.check_scope(payload, path)
        valid = bool(
            payload.get("experiment") == "p5_full_shell_fixed_boundary_scip"
            and int(payload["boundary_internal_edges"]) == expected_a
            and payload.get("boundary_cross_edges") is None
            and payload.get("solver_status") == "infeasible"
            and payload.get("finite_infeasibility_certificate") is True
            and payload.get("feasible") is False
        )
        if not valid:
            raise AssertionError(f"invalid direct internal-count closure in {path}")
        self.used_paths.add(path)
        return {
            "internal_edges": expected_a,
            "file": path.name,
            "elapsed_seconds": payload.get("elapsed_seconds"),
        }

    def validate_cross_file(self, path: Path, expected_a: int, expected_b: int) -> dict:
        payload = json.loads(path.read_text())
        self.check_scope(payload, path)
        valid = bool(
            payload.get("experiment") == "p5_full_shell_fixed_boundary_scip"
            and int(payload["boundary_internal_edges"]) == expected_a
            and int(payload["boundary_cross_edges"]) == expected_b
            and payload.get("fixed_internal_edges") is None
            and payload.get("boundary_cross_degrees") is None
            and payload.get("solver_status") == "infeasible"
            and payload.get("finite_infeasibility_certificate") is True
            and payload.get("feasible") is False
        )
        if not valid:
            raise AssertionError(f"invalid crossing-count closure in {path}")
        self.used_paths.add(path)
        return {"cross_edges": expected_b, "file": path.name}

    def audit_profile_level(self, internal_edges: int) -> dict:
        pattern_path = self.workspace / f"p5_size6_patternshards_i{self.orbit_index}_a{internal_edges}.json"
        degree_path = self.workspace / f"p5_size6_degreeshards_i{self.orbit_index}_a{internal_edges}.json"
        pattern_data = json.loads(pattern_path.read_text())
        degree_data = json.loads(degree_path.read_text())
        self.check_scope(pattern_data, pattern_path, require_hash=False)
        self.check_scope(degree_data, degree_path, require_hash=False)
        if pattern_data.get("source") != str(self.source_path):
            raise AssertionError("pattern source path mismatch")
        if degree_data.get("source") != str(self.source_path):
            raise AssertionError("degree source path mismatch")

        fresh_patterns = internal_pattern_orbits(self.boundary, internal_edges)
        if pattern_data["internal_pattern_orbits"] != fresh_patterns:
            raise AssertionError("internal-pattern orbit reconstruction mismatch")
        if int(pattern_data["full_boundary_stabilizer_size"]) != len(
            boundary_edge_stabilizers(self.boundary)
        ):
            raise AssertionError("boundary stabilizer size mismatch")
        pattern_rows = {
            (
                int(row["internal_pattern_orbit_index"]),
                int(row["boundary_cross_edges"]),
            ): row
            for row in pattern_data["rows"]
        }
        expected_pattern_rows = {
            (pattern_index, int(cross_edges))
            for pattern_index in range(len(fresh_patterns))
            for cross_edges in pattern_data["cross_counts"]
        }
        if set(pattern_rows) != expected_pattern_rows:
            raise AssertionError("pattern-shard row coverage mismatch")
        pattern_unknown = {
            tuple(int(value) for value in row)
            for row in pattern_data["unknown_cases"]
        }
        if pattern_unknown != {
            key for key, row in pattern_rows.items() if row["solver_status"] == "UNKNOWN"
        }:
            raise AssertionError("pattern unknown-case list mismatch")
        if any(row["solver_status"] in {"OPTIMAL", "FEASIBLE"} for row in pattern_rows.values()):
            raise AssertionError("pattern layer contains a feasible case")

        cube_records = {}
        expected_profiles: dict[tuple[int, int, int], Profile] = {}
        for cube in degree_data["cube_records"]:
            pattern_index = int(cube["internal_pattern_orbit_index"])
            cross_edges = int(cube["boundary_cross_edges"])
            pattern = tuple(
                int(value)
                for value in fresh_patterns[pattern_index]["representative_indices"]
            )
            fresh_degrees = degree_orbits(self.boundary, pattern, cross_edges)
            if cube["degree_orbits"] != fresh_degrees:
                raise AssertionError("degree-orbit reconstruction mismatch")
            if int(cube["degree_orbit_count"]) != len(fresh_degrees):
                raise AssertionError("degree-orbit count mismatch")
            if int(cube["all_degree_vectors"]) != sum(
                int(row["orbit_size"]) for row in fresh_degrees
            ):
                raise AssertionError("raw degree-vector count mismatch")
            cube_records[(pattern_index, cross_edges)] = cube
            fixed_edges = tuple(
                sorted(
                    edge_tuple(edge)
                    for edge in fresh_patterns[pattern_index]["representative_edges"]
                )
            )
            for degree_index, degree_record in enumerate(fresh_degrees):
                expected_profiles[(pattern_index, cross_edges, degree_index)] = Profile(
                    internal_edges,
                    fixed_edges,
                    cross_edges,
                    tuple(int(value) for value in degree_record["representative"]),
                )
        if set(cube_records) != pattern_unknown:
            raise AssertionError("degree layer does not match pattern timeouts")

        degree_rows = {
            (
                int(row["internal_pattern_orbit_index"]),
                int(row["boundary_cross_edges"]),
                int(row["cross_degree_orbit_index"]),
            ): row
            for row in degree_data["rows"]
        }
        if set(degree_rows) != set(expected_profiles):
            raise AssertionError("degree-shard row coverage mismatch")
        if int(degree_data["case_count"]) != len(degree_rows):
            raise AssertionError("degree case count mismatch")
        initial_infeasible = 0
        scip_closures = []
        for key, profile in sorted(expected_profiles.items()):
            row = degree_rows[key]
            if profile_from_payload(row, self.boundary) != profile:
                raise AssertionError("degree row has the wrong profile")
            if row.get("source_sha256") != self.source_sha:
                raise AssertionError("degree row source hash mismatch")
            status = row["solver_status"]
            if status == "INFEASIBLE" and row.get("finite_infeasibility_certificate") is True:
                initial_infeasible += 1
                continue
            if status in {"OPTIMAL", "FEASIBLE"}:
                raise AssertionError("degree layer contains a feasible case")
            closed, proof = self.close_branch(profile, frozenset())
            if not closed:
                raise AssertionError(f"unclosed profile {profile}: {proof}")
            scip_closures.append(
                {
                    "pattern_index": key[0],
                    "cross_edges": key[1],
                    "degree_index": key[2],
                    "degree_vector": list(profile.crossing_degrees),
                    "proof": proof,
                }
            )

        selected_cross_counts = {int(value) for value in pattern_data["cross_counts"]}
        all_cross_counts = set(range(0, 22 - internal_edges, 2))
        omitted = sorted(all_cross_counts - selected_cross_counts)
        cross_closures = []
        for cross_edges in omitted:
            path = self.workspace / (
                f"p5_scip_crossclosure_i{self.orbit_index}_a{internal_edges}_b{cross_edges}.json"
            )
            cross_closures.append(
                self.validate_cross_file(path, internal_edges, cross_edges)
            )

        self.used_paths.update({pattern_path, degree_path})
        return {
            "internal_edges": internal_edges,
            "all_even_cross_counts": sorted(all_cross_counts),
            "profiled_cross_counts": sorted(selected_cross_counts),
            "direct_cross_count_closures": cross_closures,
            "internal_pattern_orbits": len(fresh_patterns),
            "degree_profiles": len(expected_profiles),
            "initial_cpsat_infeasible_profiles": initial_infeasible,
            "scip_closed_profiles": len(scip_closures),
            "scip_profile_proofs": scip_closures,
        }

    def audit(self) -> dict:
        if self.source_sha != EXPECTED_SOURCE_SHA256:
            raise AssertionError("unexpected normalized-circle source hash")
        if (int(self.source["p"]), self.c_h, self.boundary) != (
            5,
            -1,
            (0, 1, 2, 3, 4, 5),
        ):
            raise AssertionError("this audit only covers normalized circle index 0")
        if self.edges[0] != (0, 1) or (0, 1) not in {
            edge for edge in self.edges if set(edge) <= self.boundary_set
        }:
            raise AssertionError("mandatory edge is not internal to this circle")

        self.load_scip_records(
            [
                "p5_scip_direct_i0_a1_*.json",
                "p5_scip_enum_i0_a1*.json",
                "p5_scip_full_i0_a1_*.json",
                "p5_scip_hard_i0_a1_*.json",
                "p5_scip_retry_i0_a1_*.json",
                "p5_scip_retry2_i0_a1_*.json",
                "p5_scip_direct_i0_a2_*.json",
                "p5_scip_retry_i0_a2_*.json",
                "p5_scip_split_i0_a2_*.json",
                "p5_scip_split2_i0_a2_*.json",
                "p5_scip_direct_i0_a3_*.json",
                "p5_scip_retry_i0_a3_*.json",
            ]
        )

        direct_counts = []
        for internal_edges in (0, *range(4, 16)):
            path = self.workspace / (
                f"p5_scip_countclosure_i{self.orbit_index}_a{internal_edges}.json"
            )
            direct_counts.append(self.validate_count_file(path, internal_edges))
        profile_levels = [self.audit_profile_level(value) for value in (1, 2, 3)]
        covered_counts = sorted(
            [row["internal_edges"] for row in direct_counts]
            + [row["internal_edges"] for row in profile_levels]
        )
        proved = covered_counts == list(range(16))
        artifacts = [
            {
                "file": path.name,
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in sorted(self.used_paths, key=lambda value: value.name)
        ]
        return {
            "experiment": "p5_size6_circle_scip_audit",
            "status": "exact_layered_finite_coverage_audit",
            "proved": proved,
            "scope": {
                "p": 5,
                "c_H": self.c_h,
                "orbit_index": self.orbit_index,
                "boundary": list(self.boundary),
                "mandatory_edge": [0, 1],
                "source": str(self.source_path),
                "source_sha256": self.source_sha,
            },
            "covered_internal_edge_counts": covered_counts,
            "direct_internal_count_closures": direct_counts,
            "profile_levels": profile_levels,
            "used_artifact_count": len(artifacts),
            "used_artifacts": artifacts,
            "solver_claim_limit": (
                "The audit proves exhaustive finite partition coverage and checks each "
                "recorded SCIP infeasibility status. It is not a separately checkable "
                "SCIP proof-log certificate."
            ),
        }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--orbit-index", type=int, default=0)
    parser.add_argument("--workspace", type=Path, default=Path("/tmp"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = CircleAudit(args.source, args.orbit_index, args.workspace).audit()
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
