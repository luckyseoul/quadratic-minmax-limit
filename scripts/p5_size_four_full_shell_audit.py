#!/usr/bin/env python3
"""Structural audit of the complete p=5 four-point certificate.

This independently rebuilds the boundary orbits, full-shell matrices,
parity/lift masses, and nonsquare sign-transfer bijection.  It then checks
that every claimed finite-solver exclusion has the expected scope and that
the sole mod-five timeout is exactly the orbit closed by the separate
mod-seven certificate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15598 import field_ctx  # noqa: E402
from p7_unsaturated_modular_catalog_filter import left_dependencies  # noqa: E402
from residual_boundary_four_lift_cpsat import geometry  # noqa: E402
from residual_size_four_boundary_orbits import classify, stabilizer_permutations  # noqa: E402


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


def source_rows(path: Path, c_h: int, infinity_value: int) -> tuple[dict, list[dict], bool]:
    source = json.loads(path.read_text())
    if (int(source["p"]), int(source["c_H"]), int(source["infinity_value"])) != (
        5,
        c_h,
        infinity_value,
    ):
        raise ValueError(f"unexpected source scope in {path}")
    fresh = classify(5, c_h, infinity_value)
    rows = list(source["orbits"])
    matches = bool(
        int(source["candidate_boundaries"]) == int(fresh["candidate_boundaries"])
        and int(source["stabilizer_size"]) == int(fresh["stabilizer_size"]) == 24
        and int(source["orbit_count"]) == int(fresh["orbit_count"]) == len(rows)
        and int(source["orbit_size_sum"])
        == int(fresh["orbit_size_sum"])
        == sum(int(row["size"]) for row in rows)
        and [
            (row["representative_vertices"], int(row["size"])) for row in rows
        ]
        == [
            (row["representative_vertices"], int(row["size"]))
            for row in fresh["orbits"]
        ]
    )
    return source, rows, matches


def shell_representatives() -> tuple[dict, dict]:
    data = geometry(5, "full")
    edges = data["edges"]
    C = data["C"]
    left = np.asarray([a for a, _b in edges], dtype=np.int16)
    right = np.asarray([b for _a, b in edges], dtype=np.int16)
    edge_count = np.ones(len(edges), dtype=np.int16)
    fixed_edge = np.zeros(len(edges), dtype=np.int16)
    fixed_edge[edges.index((0, 1))] = 1
    reps = {}
    bad_rows = {}
    dimensions = {}
    combined_rows = [edge_count, fixed_edge]
    for eps in (-1, 1):
        Y = data["shells"][eps]
        features = (Y[:, left] * Y[:, right] * C[left, right]).astype(np.int8)
        unique, indices, counts = np.unique(
            np.ascontiguousarray(features),
            axis=0,
            return_index=True,
            return_counts=True,
        )
        normalized = eps * unique
        if unique.shape != (130, 325) or not np.all(counts == 2):
            raise AssertionError("bad antipodal shell reconstruction")
        if not np.all(normalized.sum(axis=0) == 26):
            raise AssertionError("bad normalized shell column sum")
        bad = (normalized < 0).astype(np.int16)
        matrix = np.stack([edge_count, fixed_edge, *bad])
        ranks = {}
        for modulus in (5, 7):
            rank, dependencies = left_dependencies(matrix, modulus)
            if rank != 67 or dependencies.shape != (65, 132):
                raise AssertionError("bad shell rank or nullity")
            if np.any(dependencies @ (matrix % modulus) % modulus):
                raise AssertionError("bad shell left-null witness")
            ranks[str(modulus)] = {
                "rank": rank,
                "left_dependency_dimension": len(dependencies),
            }
        reps[eps] = Y[indices].astype(np.int8)
        bad_rows[eps] = bad
        dimensions[str(eps)] = ranks
        combined_rows.extend(bad)
    combined = np.stack(combined_rows)
    rank, dependencies = left_dependencies(combined, 5)
    if combined.shape != (262, 325) or rank != 113 or dependencies.shape != (149, 262):
        raise AssertionError("bad combined mod-five dimensions")
    if np.any(dependencies @ (combined % 5) % 5):
        raise AssertionError("bad combined left-null witness")
    return reps, {
        "shells": dimensions,
        "combined_mod_5": {
            "equations": 262,
            "edge_variables": 325,
            "rank": rank,
            "left_dependency_dimension": len(dependencies),
        },
        "left_null_audits": True,
    }


def expected_parity_mass(reps: dict, eps: int, c_h: int, boundary: tuple[int, ...]) -> int:
    products = np.prod(reps[eps][:, boundary].astype(np.int16), axis=1)
    return int(np.count_nonzero(-eps * c_h * products == -1))


def audit_batch(path: Path, source_path: Path, source: dict, orbits: list[dict], reps: dict) -> dict:
    payload = json.loads(path.read_text())
    if payload["source_sha256"] != sha256(source_path):
        raise AssertionError("batch source hash mismatch")
    expected = {
        index: (
            tuple(int(value) for value in row["representative_vertices"]),
            int(row["size"]),
        )
        for index, row in enumerate(orbits)
    }
    malformed = []
    valid_excluded = set()
    unknown = set()
    for row in payload["rows"]:
        index = int(row["orbit_index"])
        boundary = tuple(int(value) for value in row["representative_vertices"])
        if index not in expected or expected[index] != (boundary, int(row["orbit_size"])):
            malformed.append(index)
            continue
        for shell in row.get("shell_rows", []):
            eps = int(shell["eps"])
            parity_mass = expected_parity_mass(
                reps, eps, int(source["c_H"]), boundary
            )
            if int(shell["parity_mass"]) != parity_mass:
                malformed.append(index)
            expected_lift = None if parity_mass > 78 else (78 - parity_mass) // 2
            if shell["lift_mass"] != expected_lift:
                malformed.append(index)
        shell_exclusion = any(
            bool(shell.get("mod5_infeasible"))
            and shell.get("solver_status") in {"INFEASIBLE", "PARITY_MASS_INFEASIBLE"}
            for shell in row.get("shell_rows", [])
        )
        combined_exclusion = bool(
            row.get("combined", {}).get("mod5_infeasible")
            and row.get("combined", {}).get("solver_status") == "INFEASIBLE"
        )
        if bool(row["excluded"]) and (shell_exclusion or combined_exclusion):
            valid_excluded.add(index)
        if any(
            shell.get("solver_status") == "UNKNOWN"
            for shell in row.get("shell_rows", [])
        ) or row.get("combined", {}).get("solver_status") == "UNKNOWN":
            unknown.add(index)
    covered = sum(expected[index][1] for index in valid_excluded)
    structurally_valid = bool(
        not malformed
        and len(payload["rows"]) == len(orbits)
        and int(payload["completed"]) == len(orbits)
        and int(payload["excluded"]) == len(valid_excluded)
        and int(payload["unknown"]) == len(unknown)
        and int(payload["covered_boundary_count"]) == covered
    )
    return {
        "path": str(path),
        "sha256": sha256(path),
        "structurally_valid": structurally_valid,
        "orbits": len(orbits),
        "valid_mod5_exclusions": len(valid_excluded),
        "covered_boundaries": covered,
        "unknown_orbits": sorted(unknown),
        "malformed_orbits": sorted(set(malformed)),
    }


def sign_transfer_bijection(minus_rows: list[dict], plus_rows: list[dict]) -> bool:
    q, mul, _add, chi, _frob, _norm, _ia, _ib = field_ctx(5)
    alpha = next(value for value in range(1, q) if chi(value) == -1)
    permutations = stabilizer_permutations(5)

    def canonical(points: tuple[int, ...]) -> tuple[int, ...]:
        return min(
            tuple(sorted(permutation[value] for value in points))
            for permutation in permutations
        )

    plus = {
        canonical(tuple(int(v) for v in row["representative_finite_field"])):
        int(row["size"])
        for row in plus_rows
    }
    images = {}
    for row in minus_rows:
        points = tuple(int(v) for v in row["representative_finite_field"])
        image = canonical(tuple(mul(alpha, value) for value in points))
        images[image] = int(row["size"])
    return images == plus and len(images) == len(minus_rows) == len(plus_rows)


def audit(args: argparse.Namespace) -> dict:
    source_specs = (
        ("minus_no_infinity", args.minus_no_infinity_source, -1, 0),
        ("minus_infinity", args.minus_infinity_source, -1, 1),
        ("plus_infinity", args.plus_infinity_source, 1, 1),
        ("plus_no_infinity", args.plus_no_infinity_source, 1, 0),
    )
    sources = {}
    rows = {}
    fresh_matches = {}
    for name, path, c_h, infinity_value in source_specs:
        sources[name], rows[name], fresh_matches[name] = source_rows(
            path, c_h, infinity_value
        )
    reps, linear_audit = shell_representatives()
    batches = {
        "minus_no_infinity": audit_batch(
            args.minus_no_infinity_result,
            args.minus_no_infinity_source,
            sources["minus_no_infinity"],
            rows["minus_no_infinity"],
            reps,
        ),
        "minus_infinity": audit_batch(
            args.minus_infinity_result,
            args.minus_infinity_source,
            sources["minus_infinity"],
            rows["minus_infinity"],
            reps,
        ),
        "plus_infinity": audit_batch(
            args.plus_infinity_result,
            args.plus_infinity_source,
            sources["plus_infinity"],
            rows["plus_infinity"],
            reps,
        ),
    }
    exception = json.loads(args.exception.read_text())
    exception_valid = bool(
        exception["p"] == 5
        and exception["modulus"] == 7
        and exception["eps"] == 1
        and exception["c_H"] == -1
        and exception["boundary"] == [2, 3, 12, 13]
        and exception["orbit_index_in_c_minus_infinity_zero_source"] == 164
        and exception["orbit_size"] == 24
        and exception["rank"] == 67
        and exception["left_dependency_dimension"] == 65
        and exception["left_null_audit"] is True
        and exception["parity_mass"] == 56
        and exception["lift_mass"] == 11
        and exception["solver_status"] == "INFEASIBLE"
        and exception["mod7_infeasible"] is True
        and batches["minus_no_infinity"]["unknown_orbits"] == [164]
    )
    transfer = sign_transfer_bijection(
        rows["minus_no_infinity"], rows["plus_no_infinity"]
    )
    direct_orbits = sum(batch["orbits"] for batch in batches.values())
    direct_mod5 = sum(batch["valid_mod5_exclusions"] for batch in batches.values())
    direct_boundaries = (
        sum(batch["covered_boundaries"] for batch in batches.values()) + 24
    )
    proved = bool(
        all(fresh_matches.values())
        and all(batch["structurally_valid"] for batch in batches.values())
        and batches["minus_no_infinity"]["valid_mod5_exclusions"] == 488
        and batches["minus_infinity"]["valid_mod5_exclusions"] == 112
        and batches["plus_infinity"]["valid_mod5_exclusions"] == 112
        and batches["minus_infinity"]["unknown_orbits"] == []
        and batches["plus_infinity"]["unknown_orbits"] == []
        and exception_valid
        and transfer
        and direct_orbits == 713
        and direct_mod5 == 712
        and direct_boundaries == 15525
    )
    return {
        "experiment": "p5_size_four_full_shell_audit",
        "status": "independent_structure_coverage_and_sign_transfer_audit",
        "proved": proved,
        "fresh_orbit_reclassification_matches": fresh_matches,
        "source_files": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path, _c_h, _infinity in source_specs
        },
        "linear_system": linear_audit,
        "batches": batches,
        "exception": {
            "path": str(args.exception),
            "sha256": sha256(args.exception),
            "valid": exception_valid,
        },
        "nonsquare_no_infinity_orbit_bijection": transfer,
        "direct_orbits": direct_orbits,
        "direct_mod5_exclusions": direct_mod5,
        "direct_mod7_exclusions": 1 if exception_valid else 0,
        "direct_boundaries_excluded": direct_boundaries,
        "transferred_orbits": 489 if transfer else 0,
        "transferred_boundaries": 10925 if transfer else 0,
        "all_p5_floor_surviving_orbit_sign_cases": 1202,
        "all_p5_floor_surviving_boundary_sign_cases": 26450,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--minus-no-infinity-source", type=Path, required=True)
    parser.add_argument("--minus-no-infinity-result", type=Path, required=True)
    parser.add_argument("--minus-infinity-source", type=Path, required=True)
    parser.add_argument("--minus-infinity-result", type=Path, required=True)
    parser.add_argument("--plus-infinity-source", type=Path, required=True)
    parser.add_argument("--plus-infinity-result", type=Path, required=True)
    parser.add_argument("--plus-no-infinity-source", type=Path, required=True)
    parser.add_argument("--exception", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args)
    atomic_write(args.output, result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
