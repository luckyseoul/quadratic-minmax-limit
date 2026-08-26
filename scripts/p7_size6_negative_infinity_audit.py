#!/usr/bin/env python3
"""Audit the independent p=7 negative infinity-plus-five exclusions."""
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

from residual_size_four_boundary_orbits import stabilizer_permutations  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def atomic_write(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def artifact(path: Path) -> dict:
    return {"file": path.name, "sha256": sha256(path), "bytes": path.stat().st_size}


def stripped(payload: dict, *keys: str) -> dict:
    result = dict(payload)
    for key in keys:
        result.pop(key, None)
    return result


def audit(args: argparse.Namespace) -> dict:
    floor_cuda = load(args.floor_cuda)
    floor_numpy = load(args.floor_numpy)
    for payload, backend in ((floor_cuda, "cuda"), (floor_numpy, "numpy")):
        if (
            payload.get("experiment") != "p7_size6_negative_infinity_floor_gpu"
            or payload.get("status") != "complete_exact_integer_floor_budget_sieve"
            or payload.get("backend") != backend
            or int(payload["p"]) != 7
            or int(payload["c_H"]) != -1
            or int(payload["boundary_size"]) != 6
            or int(payload["infinity_value"]) != 1
            or int(payload["checked_boundaries"]) != 1_906_884
        ):
            raise AssertionError(f"invalid {backend} floor recording")
    floor_ignored = ("backend", "device", "elapsed_seconds")
    if stripped(floor_cuda, *floor_ignored) != stripped(floor_numpy, *floor_ignored):
        raise AssertionError("CUDA and NumPy floor sweeps differ")
    floor_survivors = {
        tuple(int(value) for value in row)
        for row in floor_cuda["survivors_finite_field"]
    }
    if (
        len(floor_survivors) != int(floor_cuda["floor_surviving_boundaries"])
        or len(floor_survivors) != 83_496
    ):
        raise AssertionError("floor survivor list is not a complete set")

    serial = load(args.serial_orbits)
    seeded = load(args.seeded_orbits)
    if stripped(serial, "elapsed_seconds", "upstream_floor_sieve") != stripped(
        seeded, "elapsed_seconds", "upstream_floor_sieve"
    ):
        raise AssertionError("serial and GPU-seeded orbit catalogs differ")
    if (
        int(serial["candidate_boundaries"]) != len(floor_survivors)
        or int(serial["orbit_count"]) != 1_750
        or int(serial["orbit_size_sum"]) != len(floor_survivors)
    ):
        raise AssertionError("orbit catalog summary is inconsistent")
    permutations = stabilizer_permutations(7)
    union: set[tuple[int, ...]] = set()
    for index, row in enumerate(serial["orbits"]):
        representative = tuple(int(value) for value in row["representative_finite_field"])
        orbit = {
            tuple(sorted(permutation[value] for value in representative))
            for permutation in permutations
        }
        if len(orbit) != int(row["size"]):
            raise AssertionError(f"incorrect orbit size at index {index}")
        if union & orbit:
            raise AssertionError("recorded stabilizer orbits overlap")
        union |= orbit
    if union != floor_survivors:
        raise AssertionError("recorded stabilizer orbits do not equal floor survivors")

    serial_sha = sha256(args.serial_orbits)
    affine_nuka = load(args.affine_nuka)
    affine_soulkiller = load(args.affine_soulkiller)
    for payload in (affine_nuka, affine_soulkiller):
        if (
            payload.get("experiment")
            != "p7_size6_negative_infinity_mod7_affine"
            or not payload.get("complete_orbit_enumeration")
            or payload.get("source_sha256") != serial_sha
            or int(payload["checked_boundaries"]) != 1_750
            or int(payload["checked_boundary_weight"]) != 83_496
            or int(payload["elevation_cases"]) != 2_230
            or int(payload["affine_span_surviving_cases"]) != 25
        ):
            raise AssertionError("invalid affine-span recording")
    if stripped(affine_nuka, "elapsed_seconds") != stripped(
        affine_soulkiller, "elapsed_seconds"
    ):
        raise AssertionError("NUKA and Soulkiller affine-span sweeps differ")
    if any(
        tuple(int(value) for value in row["catalog_sizes"]) != (36, 36)
        for row in affine_nuka["surviving_cases"]
    ):
        raise AssertionError("an affine survivor has an unexpected exact catalog size")

    exact_nuka = load(args.exact_nuka)
    exact_soulkiller = load(args.exact_soulkiller)
    affine_hashes = {
        sha256(args.affine_nuka),
        sha256(args.affine_soulkiller),
    }
    for payload in (exact_nuka, exact_soulkiller):
        if (
            payload.get("experiment")
            != "p7_size6_negative_infinity_mod7_exact"
            or payload.get("orbit_source_sha256") != serial_sha
            or payload.get("affine_source_sha256") not in affine_hashes
            or int(payload["source_orbits"]) != 1_750
            or int(payload["affine_rejected_cases"]) != 2_205
            or int(payload["exact_cases"]) != 25
            or int(payload["checked_catalog_pairs"]) != 32_400
            or int(payload["surviving_catalog_pairs"]) != 0
            or payload.get("complete_branch_mod7_infeasible") is not True
        ):
            raise AssertionError("invalid exact catalog-pair recording")
    exact_ignored = (
        "affine_source",
        "affine_source_sha256",
        "elapsed_seconds",
    )
    if stripped(exact_nuka, *exact_ignored) != stripped(
        exact_soulkiller, *exact_ignored
    ):
        raise AssertionError("NUKA and Soulkiller exact pair checks differ")

    paths = (
        args.floor_cuda,
        args.floor_numpy,
        args.serial_orbits,
        args.seeded_orbits,
        args.affine_nuka,
        args.affine_soulkiller,
        args.exact_nuka,
        args.exact_soulkiller,
    )
    return {
        "experiment": "p7_size6_negative_infinity_audit",
        "status": "independent_floor_orbit_affine_and_exact_mod7_audit",
        "proved": True,
        "scope": {
            "p": 7,
            "c_H": -1,
            "boundary_size": 6,
            "infinity_value": 1,
        },
        "all_boundaries": 1_906_884,
        "floor_survivors": 83_496,
        "stabilizer_orbits": 1_750,
        "elevation_cases": 2_230,
        "affine_rejected_cases": 2_205,
        "exact_catalog_cases": 25,
        "checked_exact_catalog_pairs": 32_400,
        "surviving_catalog_pairs": 0,
        "cuda_numpy_survivor_lists_identical": True,
        "serial_and_gpu_seeded_orbit_catalogs_identical": True,
        "nuka_and_soulkiller_affine_results_identical": True,
        "nuka_and_soulkiller_exact_results_identical": True,
        "complete_branch_mod7_infeasible": True,
        "artifacts": [artifact(path) for path in paths],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--floor-cuda", type=Path, required=True)
    parser.add_argument("--floor-numpy", type=Path, required=True)
    parser.add_argument("--serial-orbits", type=Path, required=True)
    parser.add_argument("--seeded-orbits", type=Path, required=True)
    parser.add_argument("--affine-nuka", type=Path, required=True)
    parser.add_argument("--affine-soulkiller", type=Path, required=True)
    parser.add_argument("--exact-nuka", type=Path, required=True)
    parser.add_argument("--exact-soulkiller", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = audit(args)
    atomic_write(args.output, out)
    print(json.dumps({key: value for key, value in out.items() if key != "artifacts"}, indent=2))


if __name__ == "__main__":
    main()
