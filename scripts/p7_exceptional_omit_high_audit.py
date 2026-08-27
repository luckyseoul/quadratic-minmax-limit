#!/usr/bin/env python3
"""Independent audit of all exceptional high-mean omission certificates."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np

from p7_exceptional_omit_high_catalogs import modular_rank
from p7_exceptional_projected_catalogs import sha256
from p7_fixed_boundary_catalog_join import direction_rows, mapped_catalog
from p7_fixed_boundary_modular_cpsat import linear_data


ORBIT_INDICES = (0, 5, 8, 25, 26, 30, 31)
OMITTED_DIRECTIONS = (0, 2, 5, 7)


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def matrix_sha256(matrix: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(matrix).tobytes()).hexdigest()


def pack11(digits: np.ndarray) -> np.ndarray:
    if digits.ndim != 2 or digits.shape[0] > 11:
        raise ValueError("pack11 expects at most eleven rows")
    packed = np.zeros(digits.shape[1], dtype=np.uint64)
    for index, row in enumerate(digits):
        packed |= row.astype(np.uint64) << (4 * index)
    return packed


def packed22(digits: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if digits.shape[0] != 22 or np.any(digits < 0) or np.any(digits >= 7):
        raise AssertionError("invalid mod-seven projection digits")
    return pack11(digits[:11]), pack11(digits[11:])


def high_direction(batch: dict, leaf: dict) -> int | None:
    hits = [
        direction_index
        for direction_index, (row, mean) in enumerate(
            zip(batch["direction_rows"], leaf["scaled_means_direction_order"])
        )
        if int(row["b"]) == 0
        and int(row["phase"]) == 0
        and int(mean) > 16
    ]
    if not hits:
        return None
    if len(hits) != 1 or hits[0] not in OMITTED_DIRECTIONS:
        raise AssertionError("unexpected high-mean pattern")
    return hits[0]


def audit_projection_family(
    projection_dir: Path,
    mean_batches: dict[int, dict],
    dependency: np.ndarray,
    omitted_direction: int,
) -> dict:
    summaries = [
        json.loads(
            (
                projection_dir
                / f"p7_exceptional_omitd{omitted_direction}_p{profile}_all.json"
            ).read_text()
        )
        for profile in range(3)
    ]
    transformed_profiles = []
    selected_basis_sets = []
    catalog_hashes: dict[tuple[int, int], str] = {}
    for profile, summary in enumerate(summaries):
        if (
            summary.get("experiment") != "p7_exceptional_omit_high_catalogs"
            or summary.get("status")
            != "complete_exact_high_direction_eliminating_projection"
            or int(summary.get("p", 0)) != 7
            or int(summary.get("c_H", 0)) != -1
            or int(summary.get("modulus", 0)) != 7
            or int(summary.get("omitted_direction", -1)) != omitted_direction
            or int(summary.get("projection_profile", -1)) != profile
            or int(summary.get("conditioned_dependency_dimension", -1)) != 112
            or int(summary.get("conditioned_dependency_basis_rank", -1)) != 112
            or int(summary.get("omitted_direction_block_rank", -1)) != 23
            or summary.get("conditioned_dependency_block_is_zero") is not True
            or summary.get("large_catalog_treatment")
            != "exactly_eliminated_by_zero_dependency_block"
            or int(summary.get("projection_group_order", 0)) != 7**22
            or int(summary.get("orbit_count", -1)) != 7
        ):
            raise ValueError(
                f"invalid omission summary direction {omitted_direction} profile {profile}"
            )
        coefficients = np.asarray(
            summary["selected_dependency_coefficients_mod7"], dtype=np.int64
        )
        if coefficients.shape != (22, 135):
            raise AssertionError("selected coefficient matrix has wrong shape")
        transformed = coefficients @ dependency % 7
        if modular_rank(coefficients) != 22 or modular_rank(transformed) != 22:
            raise AssertionError("selected dependency profile lost rank")
        high_block = transformed[
            :, 2 + 35 * omitted_direction : 2 + 35 * (omitted_direction + 1)
        ]
        if np.any(high_block):
            raise AssertionError("selected dependency touches omitted direction")
        if matrix_sha256(transformed.astype(np.uint8)) != summary.get(
            "selected_dependency_rows_sha256"
        ):
            raise AssertionError("selected dependency hash mismatch")
        selected_basis = tuple(
            int(value) for value in summary["selected_conditioned_basis_rows"]
        )
        if len(selected_basis) != 22 or len(set(selected_basis)) != 22:
            raise AssertionError("invalid conditioned-basis row selection")
        selected_basis_sets.append(set(selected_basis))
        transformed_profiles.append(transformed.astype(np.uint8))

        orbit_rows = {
            int(row["orbit_index"]): row for row in summary.get("orbits", [])
        }
        if set(orbit_rows) != set(ORBIT_INDICES):
            raise AssertionError("projection summary has wrong orbit coverage")
        for orbit_index in ORBIT_INDICES:
            batch = mean_batches[orbit_index]
            boundary = tuple(int(value) for value in batch["fixed_boundary"])
            leaves = [
                leaf
                for leaf in batch["leaves"]
                if leaf.get("solver_status") != "INFEASIBLE"
                and high_direction(batch, leaf) == omitted_direction
            ]
            orbit_row = orbit_rows[orbit_index]
            if (
                tuple(int(value) for value in orbit_row["fixed_boundary"])
                != boundary
                or set(int(value) for value in orbit_row["high_leaf_indices"])
                != {int(leaf["leaf_index"]) for leaf in leaves}
                or int(orbit_row["high_leaf_count"]) != len(leaves)
            ):
                raise AssertionError("projection high-leaf coverage mismatch")
            catalog_path = projection_dir / (
                f"cminus_exceptional_omitd{omitted_direction}_p{profile}_"
                f"orbit{orbit_index:02d}.npz"
            )
            digest = sha256(catalog_path)
            if digest != orbit_row.get("sha256"):
                raise AssertionError("projection catalog hash mismatch")
            catalog_hashes[(profile, orbit_index)] = digest

            means_by_direction = [
                sorted(
                    {
                        int(leaf["scaled_means_direction_order"][direction_index])
                        for leaf in leaves
                    }
                )
                for direction_index in range(8)
            ]
            rows = direction_rows(-1, boundary)
            expected_keys = {"base_p7_lo", "base_p7_hi"}
            with np.load(catalog_path, allow_pickle=False) as source:
                base_digits = (
                    transformed[:, :2] @ np.asarray([29, 1], dtype=np.int64) % 7
                ).reshape(22, 1)
                expected_lo, expected_hi = packed22(base_digits)
                if not np.array_equal(source["base_p7_lo"], expected_lo) or not np.array_equal(
                    source["base_p7_hi"], expected_hi
                ):
                    raise AssertionError("projection base signature mismatch")
                for direction_index, means in enumerate(means_by_direction):
                    block = transformed[
                        :,
                        2 + 35 * direction_index : 2 + 35 * (direction_index + 1),
                    ]
                    for mean in means:
                        stem = f"d{direction_index}_m{mean}_p7"
                        expected_keys.update((f"{stem}_lo", f"{stem}_hi"))
                        if direction_index == omitted_direction:
                            if mean not in (20, 24, 32) or np.any(block):
                                raise AssertionError("invalid omitted catalog metadata")
                            digits = np.zeros((22, 1), dtype=np.int64)
                        else:
                            values = mapped_catalog(
                                int(rows[direction_index]["b"]),
                                int(rows[direction_index]["phase"]),
                                mean,
                                set(rows[direction_index]["B"]),
                                None,
                            )
                            if len(values) > 1764:
                                raise AssertionError("audit unexpectedly needs a large catalog")
                            digits = block @ ((13 - values.astype(np.int64)).T % 7) % 7
                        expected_lo, expected_hi = packed22(digits)
                        if not np.array_equal(source[f"{stem}_lo"], expected_lo) or not np.array_equal(
                            source[f"{stem}_hi"], expected_hi
                        ):
                            raise AssertionError(
                                f"catalog projection mismatch d={omitted_direction} "
                                f"orbit={orbit_index} profile={profile} key={stem}"
                            )
                if set(source.files) != expected_keys:
                    raise AssertionError("projection catalog key coverage mismatch")

    if any(
        selected_basis_sets[i] & selected_basis_sets[j]
        for i in range(3)
        for j in range(i)
    ):
        raise AssertionError("conditioned projection profiles overlap")
    combined = np.concatenate(transformed_profiles, axis=0)
    if modular_rank(combined) != 66:
        raise AssertionError("combined conditioned projections lost rank")
    expected_block_ranks = [
        0 if direction == omitted_direction else 14 for direction in range(8)
    ]
    actual_block_ranks = [
        modular_rank(combined[:, 2 + 35 * direction : 2 + 35 * (direction + 1)])
        for direction in range(8)
    ]
    if actual_block_ranks != expected_block_ranks:
        raise AssertionError("combined projection block ranks changed")
    family_hash = matrix_sha256(combined.astype(np.uint8))
    if any(
        summary.get("selected_full_dependency_sha256") != family_hash
        or summary.get("selected_full_block_ranks") != expected_block_ranks
        for summary in summaries
    ):
        raise AssertionError("combined projection family metadata mismatch")
    return {
        "omitted_direction": omitted_direction,
        "conditioned_dependency_rows": 112,
        "selected_dependency_rows": 66,
        "selected_block_ranks": expected_block_ranks,
        "selected_family_sha256": family_hash,
        "catalog_hashes": catalog_hashes,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--projection-dir", type=Path, required=True)
    parser.add_argument("--gpu-root", type=Path, required=True)
    parser.add_argument("--mean-dir", type=Path, required=True)
    parser.add_argument("--ordinary-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = time.time()
    mean_batches = {
        orbit_index: json.loads(
            (
                args.mean_dir
                / f"cminus_exceptional_orbit{orbit_index:02d}_means.json"
            ).read_text()
        )
        for orbit_index in ORBIT_INDICES
    }
    _matrix, dependencies, _rows = linear_data((7,))
    dependency = dependencies[7].astype(np.int64)
    if dependency.shape != (135, 282) or modular_rank(dependency) != 135:
        raise AssertionError("common mod-seven dependency basis changed")

    family_rows = [
        audit_projection_family(
            args.projection_dir, mean_batches, dependency, omitted_direction
        )
        for omitted_direction in OMITTED_DIRECTIONS
    ]
    family_by_direction = {row["omitted_direction"]: row for row in family_rows}

    total_high = 0
    orbit_rows = []
    for orbit_index in ORBIT_INDICES:
        batch = mean_batches[orbit_index]
        high = {
            int(leaf["leaf_index"]): high_direction(batch, leaf)
            for leaf in batch["leaves"]
            if leaf.get("solver_status") != "INFEASIBLE"
            and high_direction(batch, leaf) is not None
        }
        output_dir = args.gpu_root / f"orbit{orbit_index:02d}"
        certificates = {
            int(path.stem.removeprefix("leaf")): json.loads(path.read_text())
            for path in sorted(output_dir.glob("leaf*.json"))
        }
        if set(certificates) != set(high):
            raise AssertionError(f"GPU high-leaf coverage mismatch orbit {orbit_index}")
        for leaf_index, omitted_direction in high.items():
            certificate = certificates[leaf_index]
            leaf = batch["leaves"][leaf_index]
            expected_catalog_hashes = [
                family_by_direction[omitted_direction]["catalog_hashes"][
                    (profile, orbit_index)
                ]
                for profile in range(3)
            ]
            if (
                certificate.get("experiment")
                != "p7_exceptional_projected_join_gpu"
                or certificate.get("status")
                != "complete_exact_selected_dependency_gpu_join"
                or certificate.get("projection_mode")
                != "injective_disjoint_mod7_22x3_omitted_high_direction"
                or int(certificate.get("omitted_high_mean_direction", -1))
                != omitted_direction
                or int(certificate.get("leaf_index", -1)) != leaf_index
                or certificate.get("fixed_boundary") != batch["fixed_boundary"]
                or certificate.get("fixed_scaled_means")
                != leaf["scaled_means_direction_order"]
                or int(certificate.get("exact_projected_matches", -1)) != 0
                or certificate.get("projected_modularly_infeasible") is not True
                or certificate.get("finite_mean_allocation_exclusion") is not True
                or certificate.get("projection_catalog_sha256")
                != expected_catalog_hashes
            ):
                raise AssertionError(
                    f"invalid high-mean GPU certificate orbit {orbit_index} leaf {leaf_index}"
                )
        summary_path = args.gpu_root / f"orbit{orbit_index:02d}_summary.json"
        summary = json.loads(summary_path.read_text())
        if (
            summary.get("experiment") != "p7_exceptional_omit_high_gpu_batch"
            or summary.get("status")
            != "complete_exact_high_direction_omission_gpu_batch"
            or int(summary.get("orbit_index", -1)) != orbit_index
            or int(summary.get("selected_high_leaf_count", -1)) != len(high)
            or int(summary.get("eligible_leaf_count", -1)) != len(high)
            or int(summary.get("excluded_leaf_count", -1)) != len(high)
            or int(summary.get("unresolved_leaf_count", -1)) != 0
            or int(summary.get("deferred_leaf_count", -1)) != 0
            or {int(row["leaf_index"]) for row in summary.get("results", [])}
            != set(high)
            or any(
                not row.get("excluded") or int(row.get("matches", -1)) != 0
                for row in summary.get("results", [])
            )
        ):
            raise AssertionError(f"invalid high-mean GPU summary orbit {orbit_index}")
        total_high += len(high)
        orbit_rows.append(
            {
                "orbit_index": orbit_index,
                "fixed_boundary": batch["fixed_boundary"],
                "high_mean_leaves": len(high),
                "gpu_excluded_leaves": len(certificates),
            }
        )
    if total_high != 426:
        raise AssertionError("global high-mean leaf count changed")

    ordinary = json.loads(args.ordinary_audit.read_text())
    if (
        ordinary.get("status") != "passed_independent_all_orbit_coverage_audit"
        or int(ordinary.get("initial_infeasible_leaves", -1)) != 172
        or int(ordinary.get("gpu_excluded_leaves", -1)) != 662
        or int(ordinary.get("high_mean_unknown_leaves_not_claimed", -1)) != 426
    ):
        raise ValueError("ordinary exceptional audit is incompatible")
    if 172 + 662 + total_high != 1260:
        raise AssertionError("exceptional allocation coverage does not close")

    out = {
        "experiment": "p7_exceptional_omit_high_audit",
        "status": "passed_independent_high_mean_and_full_exceptional_audit",
        "p": 7,
        "c_H": -1,
        "exceptional_orbit_count": 7,
        "total_exact_allocations": 1260,
        "initial_cp_infeasible_allocations": 172,
        "ordinary_gpu_excluded_allocations": 662,
        "high_mean_gpu_excluded_allocations": total_high,
        "remaining_allocations": 0,
        "all_seven_exceptional_orbits_excluded_cminus1": True,
        "dependency_family_audits": [
            {key: value for key, value in row.items() if key != "catalog_hashes"}
            for row in family_rows
        ],
        "orbit_rows": orbit_rows,
        "elapsed_seconds": time.time() - started,
    }
    if args.output is not None:
        atomic_json(args.output, out)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
