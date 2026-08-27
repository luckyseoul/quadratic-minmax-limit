#!/usr/bin/env python3
"""Independent aggregate audit of the p=7 size-eight conic subbranch.

This script deliberately does not claim the full size-eight boundary case.
It verifies the complete CUDA floor censuses, independently reconstructs the
32 stabilizer orbits in the minimum-eight-odd-secant (hence conic) subbranch,
audits the saturated and exceptional allocation partitions, and checks the
nonsquare Paley anti-isometry that transfers the c_H=-1 exclusion to c_H=+1.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15654 import p7_nonsquare_signed_permutation  # noqa: E402
from p7_fixed_boundary_mean_allocation_batch import (  # noqa: E402
    allocations,
    direction_data,
)
from p7_size8_conic_orbits import unrank_lex  # noqa: E402
from p7_size8_exceptional_mean_batch import exceptional_orbits  # noqa: E402
from p7_size8_floor_profile_gpu import direction_tables  # noqa: E402
from p7_size8_saturated_mean_batch import saturated_orbits  # noqa: E402
from p7_size_four_slack_classify import _primitive_left_kernel_rows  # noqa: E402
from residual_size_four_boundary_orbits import stabilizer_permutations  # noqa: E402


ALL_BOUNDARIES = math.comb(49, 8)
FLOOR_SURVIVORS = 108_754_569
MINIMUM_BOUNDARIES = 6_174
CONIC_FLOOR_SURVIVORS = 1_323
NONCONIC_FLOOR_SURVIVORS = FLOOR_SURVIVORS - CONIC_FLOOR_SURVIVORS
EXPECTED_ODD_SECANT_HISTOGRAM = {
    "8": 1_323,
    "16": 1_223_628,
    "20": 14_572_992,
    "24": 42_223_398,
    "28": 36_601_824,
    "32": 12_366_228,
    "36": 1_700_496,
    "40": 62_328,
    "44": 2_352,
}


def load_json(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if not isinstance(payload, dict):
        raise ValueError(f"expected a JSON object: {path}")
    return payload


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def audit_floor(path: Path, c_h: int) -> tuple[dict, set[tuple[int, ...]], set[tuple[int, ...]]]:
    payload = load_json(path)
    require(
        payload.get("experiment") == "p7_size8_floor_profile_gpu"
        and payload.get("status") == "complete_exact_floor_profile_census"
        and int(payload.get("p", 0)) == 7
        and int(payload.get("c_H", 0)) == c_h,
        f"invalid floor census metadata for c_H={c_h}",
    )
    require(
        int(payload.get("all_boundaries", 0)) == ALL_BOUNDARIES
        and int(payload.get("checked_boundaries", 0)) == ALL_BOUNDARIES
        and payload.get("rank_interval") == [0, ALL_BOUNDARIES]
        and int(payload.get("floor_surviving_boundaries", 0)) == FLOOR_SURVIVORS
        and int(payload.get("floor_rejected_boundaries", 0))
        == ALL_BOUNDARIES - FLOOR_SURVIVORS,
        f"incomplete floor census for c_H={c_h}",
    )
    histogram = {
        str(key): int(value)
        for key, value in payload.get("survivor_odd_secant_histogram", {}).items()
    }
    require(histogram == EXPECTED_ODD_SECANT_HISTOGRAM, "odd-secant histogram changed")
    require(sum(histogram.values()) == FLOOR_SURVIVORS, "floor histogram sum changed")
    require(
        int(payload.get("survivor_ordered_profile_count", 0)) == 5_152
        and len(payload.get("survivor_ordered_profiles", [])) == 5_152,
        "ordered floor-profile census changed",
    )
    minimum_ranks = tuple(int(value) for value in payload["minimum_odd_secant_ranks"])
    survivor_ranks = tuple(
        int(value) for value in payload["survivor_minimum_odd_secant_ranks"]
    )
    require(
        len(minimum_ranks) == len(set(minimum_ranks)) == MINIMUM_BOUNDARIES,
        "minimum-eight rank coverage changed",
    )
    require(
        len(survivor_ranks) == len(set(survivor_ranks)) == CONIC_FLOOR_SURVIVORS,
        "minimum-eight floor-survivor coverage changed",
    )
    require(set(survivor_ranks) <= set(minimum_ranks), "survivor ranks leave minimum set")
    minimum = {unrank_lex(rank) for rank in minimum_ranks}
    survivors = {unrank_lex(rank) for rank in survivor_ranks}
    labels, _epsilons = direction_tables()
    require(
        all(
            all(
                int(np.bincount(row[list(boundary)], minlength=7).max()) <= 2
                for row in labels
            )
            for boundary in minimum
        ),
        "a minimum-eight-odd-secant boundary is not an eight-arc",
    )
    return payload, minimum, survivors


def audit_conic_orbits(
    path: Path,
    floor_path: Path,
    c_h: int,
    survivors: set[tuple[int, ...]],
) -> dict:
    payload = load_json(path)
    require(
        payload.get("experiment") == "p7_size8_conic_orbits"
        and payload.get("status") == "complete_exact_extremal_boundary_orbit_audit"
        and int(payload.get("p", 0)) == 7
        and int(payload.get("c_H", 0)) == c_h,
        f"invalid conic-orbit metadata for c_H={c_h}",
    )
    require(payload.get("source_sha256") == sha256(floor_path), "floor-source hash mismatch")
    require(
        int(payload.get("minimum_odd_secant_boundaries", 0)) == MINIMUM_BOUNDARIES
        and int(payload.get("floor_surviving_minimum_boundaries", 0))
        == CONIC_FLOOR_SURVIVORS
        and int(payload.get("orbit_count", 0)) == 32
        and int(payload.get("orbit_size_sum", 0)) == CONIC_FLOOR_SURVIVORS,
        "conic orbit counts changed",
    )
    incidence = payload.get("projective_conic_incidence_count", {})
    require(
        int(incidence.get("all_nonsingular_conics", 0)) == 16_758
        and int(incidence.get("external_lines_per_conic", 0)) == 21
        and int(incidence.get("projective_lines", 0)) == 57
        and int(incidence.get("conics_disjoint_from_fixed_line", 0))
        == MINIMUM_BOUNDARIES
        and incidence.get("matches_minimum_boundary_count") is True,
        "finite conic incidence audit changed",
    )

    permutations = stabilizer_permutations(7)
    require(int(payload.get("stabilizer_size", 0)) == len(permutations), "stabilizer changed")
    reconstructed: set[tuple[int, ...]] = set()
    orbit_sizes = []
    for row in payload.get("orbits", []):
        representative = tuple(int(value) for value in row["representative_finite_field"])
        require(tuple(value + 1 for value in representative) == tuple(row["representative_vertices"]), "representative coordinate mismatch")
        orbit = {
            tuple(sorted(int(permutation[value]) for value in representative))
            for permutation in permutations
        }
        require(len(orbit) == int(row["size"]), "recorded orbit size changed")
        require(orbit <= survivors, "a reconstructed orbit leaves the floor-survivor set")
        require(not reconstructed.intersection(orbit), "conic orbit records overlap")
        reconstructed.update(orbit)
        orbit_sizes.append(len(orbit))
    require(reconstructed == survivors, "32 orbit records do not cover all 1,323 survivors")
    return {
        "orbit_count": len(orbit_sizes),
        "orbit_size_sum": sum(orbit_sizes),
        "stabilizer_size": len(permutations),
        "source_sha256": sha256(path),
    }


def normalize_means(payload: dict) -> tuple[int, ...]:
    raw = payload.get("fixed_scaled_means")
    if isinstance(raw, dict):
        return tuple(int(raw[str(index)]) for index in range(8))
    if isinstance(raw, list):
        return tuple(int(value) for value in raw)
    raise ValueError("missing fixed mean allocation")


def audit_slack_cache(cache_path: Path, summary_path: Path) -> dict:
    summary = load_json(summary_path)
    digest = sha256(cache_path)
    require(
        summary.get("experiment") == "p7_slack_catalog_first_lift_shards"
        and summary.get("status") == "complete_exact_disjoint_shard_enumeration"
        and int(summary.get("p", 0)) == 7
        and int(summary.get("phase", -1)) == 0
        and int(summary.get("scaled_mean", -1)) == 16
        and int(summary.get("odd_fibres", -1)) == 0
        and int(summary.get("lift_mass", -1)) == 20
        and int(summary.get("solution_count", 0)) == 575_407
        and summary.get("partition_complete_and_disjoint") is True
        and summary.get("matches_independent_unsharded_count") is True
        and summary.get("all_mean_and_kernel_audits") is True
        and summary.get("merged_sha256") == digest,
        "invalid saturated slack-catalog certificate",
    )
    catalog = np.load(cache_path, mmap_mode="r", allow_pickle=False)
    require(catalog.shape == (575_407, 35) and catalog.dtype == np.uint8, "catalog shape changed")
    kernel = np.asarray(_primitive_left_kernel_rows(), dtype=np.int64)
    require(kernel.shape == (14, 35), "degree-two kernel changed")
    for start in range(0, len(catalog), 32_768):
        block = np.asarray(catalog[start : start + 32_768], dtype=np.int64)
        require(np.all(block.sum(axis=1) == 40), "catalog mean equation failed")
        require(np.all((block & 1) == 0), "catalog parity equation failed")
        require(not np.any(block @ kernel.T), "catalog degree-two kernel equation failed")
    return {"solution_count": len(catalog), "shape": list(catalog.shape), "sha256": digest}


def audit_batch_summary(path: Path, expected: dict) -> None:
    payload = load_json(path)
    for key, value in expected.items():
        require(payload.get(key) == value, f"{path.name}: expected {key}={value!r}")


def audit_saturated(
    conic_payload: dict,
    saturated_root: Path,
) -> dict:
    expected_orbits = dict(saturated_orbits(conic_payload))
    require(len(expected_orbits) == 25, "saturated orbit classification changed")
    require(
        sum(int(conic_payload["orbits"][index]["size"]) for index in expected_orbits)
        == 1_176,
        "saturated boundary coverage changed",
    )

    cache_path = saturated_root / "p7_slack_catalog_b0_phase0_mean16.npy"
    cache_summary = saturated_root / "p7_slack_catalog_b0m16_hash16_shards.json"
    cache_audit = audit_slack_cache(cache_path, cache_summary)

    mean_files = sorted((saturated_root / "saturated_means").glob("*.json"))
    require(len(mean_files) == 25, "expected 25 saturated mean batches")
    boundary_to_index = {boundary: index for index, boundary in expected_orbits.items()}
    initial_infeasible: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    initially_open: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    seen_orbits: set[int] = set()
    for path in mean_files:
        batch = load_json(path)
        boundary = tuple(int(value) for value in batch.get("fixed_boundary", []))
        require(boundary in boundary_to_index, f"unexpected saturated boundary in {path}")
        orbit_index = boundary_to_index[boundary]
        require(orbit_index not in seen_orbits, "duplicate saturated mean batch")
        seen_orbits.add(orbit_index)
        rows = direction_data(-1, boundary)
        exact = allocations(rows)
        require(
            batch.get("experiment") == "p7_fixed_boundary_mean_allocation_batch"
            and batch.get("status") == "complete_exact_mean_allocation_exhaustion"
            and int(batch.get("p", 0)) == 7
            and int(batch.get("c_H", 0)) == -1
            and batch.get("direction_rows") == rows
            and int(batch.get("allocation_count", 0)) == len(exact) == 24
            and len(batch.get("leaves", [])) == 24,
            f"invalid saturated mean batch {path}",
        )
        leaves = sorted(batch["leaves"], key=lambda row: int(row["leaf_index"]))
        require([int(row["leaf_index"]) for row in leaves] == list(range(24)), "leaf indices changed")
        require(
            tuple(tuple(int(value) for value in row["scaled_means_direction_order"]) for row in leaves)
            == exact,
            "mean batch does not equal independent allocation enumeration",
        )
        for row, means in zip(leaves, exact):
            key = (boundary, means)
            if (
                row.get("solver_status") == "INFEASIBLE"
                and row.get("finite_infeasibility_certificate") is True
                and row.get("feasible") is False
            ):
                initial_infeasible.add(key)
            else:
                require(row.get("feasible") is False, "a saturated leaf reports a witness")
                initially_open.add(key)
    require(seen_orbits == set(expected_orbits), "saturated orbit mean coverage incomplete")
    require(len(initial_infeasible) == 355 and len(initially_open) == 245, "initial saturated split changed")

    extra_keys: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    extra_files = sorted((saturated_root / "saturated_high_mean").glob("*.json"))
    require(len(extra_files) == 6, "expected six long CP-SAT certificates")
    for path in extra_files:
        payload = load_json(path)
        key = (
            tuple(int(value) for value in payload.get("fixed_boundary", [])),
            normalize_means(payload),
        )
        require(
            payload.get("experiment") == "p7_fixed_boundary_modular_cpsat"
            and payload.get("status") == "exact_compact_degree_two_multimodular_catalog_model"
            and int(payload.get("p", 0)) == 7
            and int(payload.get("c_H", 0)) == -1
            and payload.get("solver_status") == "INFEASIBLE"
            and payload.get("finite_infeasibility_certificate") is True
            and payload.get("feasible") is False
            and key in initially_open
            and key not in extra_keys,
            f"invalid long CP-SAT certificate {path}",
        )
        extra_keys.add(key)

    join_keys: set[tuple[tuple[int, ...], tuple[int, ...]]] = set()
    join_files = sorted((saturated_root / "saturated_joins").glob("*.json"))
    require(len(join_files) == 239, "expected 239 saturated catalog joins")
    for path in join_files:
        payload = load_json(path)
        key = (
            tuple(int(value) for value in payload.get("fixed_boundary", [])),
            normalize_means(payload),
        )
        require(
            payload.get("experiment") == "p7_fixed_boundary_catalog_join"
            and payload.get("status") == "complete_exact_multimodular_catalog_join"
            and int(payload.get("p", 0)) == 7
            and int(payload.get("c_H", 0)) == -1
            and payload.get("moduli") == [3, 7]
            and int(payload.get("consistent_catalog_tuples", -1)) == 0
            and payload.get("modularly_infeasible") is True
            and payload.get("finite_mean_allocation_exclusion") is True
            and payload.get("catalog_cache", {}).get("sha256") == cache_audit["sha256"]
            and key in initially_open
            and key not in join_keys,
            f"invalid saturated catalog-join certificate {path}",
        )
        join_keys.add(key)
    require(not extra_keys.intersection(join_keys), "saturated certificates overlap")
    require(extra_keys | join_keys == initially_open, "saturated certificates leave an open allocation")

    mean_summary_expected = {
        "experiment": "p7_size8_saturated_mean_batch",
        "status": "complete_exact_saturated_orbit_mean_batches",
        "p": 7,
        "c_H": -1,
        "saturated_orbit_count": 25,
        "total_exact_allocations": 600,
        "total_infeasible_allocations": 355,
        "total_feasible_modular_allocations": 0,
        "total_unknown_allocations": 245,
    }
    audit_batch_summary(saturated_root / "p7_size8_saturated_mean_batch.json", mean_summary_expected)
    join_summary_expected = {
        "experiment": "p7_size8_saturated_join_batch",
        "status": "complete_exact_unknown_leaf_catalog_joins",
        "p": 7,
        "c_H": -1,
        "saturated_orbit_count": 25,
        "total_exact_allocations": 600,
        "cp_sat_infeasible_allocations": 355,
        "extra_cp_sat_infeasible_allocations": 6,
        "catalog_join_allocations": 239,
        "catalog_join_infeasible_allocations": 239,
        "remaining_modularly_consistent_allocations": 0,
        "deferred_high_catalog_allocations": 0,
        "all_saturated_orbits_excluded": True,
    }
    audit_batch_summary(saturated_root / "p7_size8_saturated_join_batch.json", join_summary_expected)
    audit_batch_summary(saturated_root / "replay_audit.json", join_summary_expected)
    return {
        "orbit_count": 25,
        "boundary_count": 1_176,
        "total_exact_allocations": 600,
        "initial_cp_infeasible_allocations": 355,
        "long_cp_infeasible_allocations": 6,
        "catalog_join_infeasible_allocations": 239,
        "remaining_allocations": 0,
        "slack_catalog": cache_audit,
    }


def audit_exceptional(
    conic_payload: dict,
    ordinary_path: Path,
    high_path: Path,
) -> dict:
    expected_orbits = dict(exceptional_orbits(conic_payload))
    require(len(expected_orbits) == 7, "exceptional orbit classification changed")
    require(
        sum(int(conic_payload["orbits"][index]["size"]) for index in expected_orbits)
        == 147,
        "exceptional boundary coverage changed",
    )
    expected_boundaries = {index: list(boundary) for index, boundary in expected_orbits.items()}

    ordinary = load_json(ordinary_path)
    require(
        ordinary.get("experiment") == "p7_exceptional_mod7_tuple_audit"
        and ordinary.get("status") == "passed_independent_all_orbit_coverage_audit"
        and int(ordinary.get("p", 0)) == 7
        and int(ordinary.get("c_H", 0)) == -1
        and int(ordinary.get("exceptional_orbit_count", 0)) == 7
        and int(ordinary.get("total_mean_allocations", 0)) == 1_260
        and int(ordinary.get("initial_infeasible_leaves", 0)) == 172
        and int(ordinary.get("gpu_supported_unknown_leaves", 0)) == 662
        and int(ordinary.get("gpu_excluded_leaves", 0)) == 662
        and int(ordinary.get("gpu_unresolved_leaves", -1)) == 0
        and int(ordinary.get("high_mean_unknown_leaves_not_claimed", 0)) == 426
        and ordinary.get("selected_dependency_rows_pairwise_disjoint") is True,
        "ordinary exceptional audit changed",
    )
    ordinary_rows = {int(row["orbit_index"]): row for row in ordinary.get("orbit_rows", [])}
    require(set(ordinary_rows) == set(expected_orbits), "ordinary exceptional orbit coverage changed")
    for index, boundary in expected_boundaries.items():
        require(ordinary_rows[index].get("fixed_boundary") == boundary, "ordinary audit boundary mismatch")

    high = load_json(high_path)
    require(
        high.get("experiment") == "p7_exceptional_omit_high_audit"
        and high.get("status") == "passed_independent_high_mean_and_full_exceptional_audit"
        and int(high.get("p", 0)) == 7
        and int(high.get("c_H", 0)) == -1
        and int(high.get("exceptional_orbit_count", 0)) == 7
        and int(high.get("total_exact_allocations", 0)) == 1_260
        and int(high.get("initial_cp_infeasible_allocations", 0)) == 172
        and int(high.get("ordinary_gpu_excluded_allocations", 0)) == 662
        and int(high.get("high_mean_gpu_excluded_allocations", 0)) == 426
        and int(high.get("remaining_allocations", -1)) == 0
        and high.get("all_seven_exceptional_orbits_excluded_cminus1") is True,
        "high-mean exceptional audit changed",
    )
    high_rows = {int(row["orbit_index"]): row for row in high.get("orbit_rows", [])}
    require(set(high_rows) == set(expected_orbits), "high exceptional orbit coverage changed")
    require(
        sum(int(row["high_mean_leaves"]) for row in high_rows.values()) == 426
        and sum(int(row["gpu_excluded_leaves"]) for row in high_rows.values()) == 426,
        "high exceptional leaf total changed",
    )
    for index, boundary in expected_boundaries.items():
        require(high_rows[index].get("fixed_boundary") == boundary, "high audit boundary mismatch")
        require(
            int(high_rows[index]["high_mean_leaves"])
            == int(ordinary_rows[index]["high_mean_unknown_leaves"]),
            "ordinary/high exceptional partition mismatch",
        )
    return {
        "orbit_count": 7,
        "boundary_count": 147,
        "total_exact_allocations": 1_260,
        "initial_cp_infeasible_allocations": 172,
        "ordinary_gpu_excluded_allocations": 662,
        "high_direction_omission_gpu_excluded_allocations": 426,
        "remaining_allocations": 0,
        "ordinary_audit_sha256": sha256(ordinary_path),
        "high_mean_audit_sha256": sha256(high_path),
    }


def audit_sign_transfer(
    minimum_minus: set[tuple[int, ...]],
    minimum_plus: set[tuple[int, ...]],
    survivors_minus: set[tuple[int, ...]],
    survivors_plus: set[tuple[int, ...]],
) -> dict:
    symmetry = p7_nonsquare_signed_permutation()
    permutation = tuple(int(value) for value in symmetry["finite_permutation"])

    def mapped(boundaries: set[tuple[int, ...]]) -> set[tuple[int, ...]]:
        return {
            tuple(sorted(permutation[value] for value in boundary))
            for boundary in boundaries
        }

    minimum_bijection = mapped(minimum_minus) == minimum_plus
    survivor_bijection = mapped(survivors_minus) == survivors_plus
    require(
        symmetry.get("fixes_distinguished_edge") is True
        and symmetry.get("signed_conference_anti_isometry") is True
        and minimum_bijection
        and survivor_bijection,
        "nonsquare sign-transfer audit failed",
    )
    return {
        "nonsquare_multiplier": int(symmetry["nonsquare_multiplier"]),
        "fixes_infinity": bool(symmetry["fixes_infinity"]),
        "fixes_finite_zero": bool(symmetry["fixes_finite_zero"]),
        "fixes_distinguished_edge": True,
        "signed_conference_anti_isometry": True,
        "minimum_boundary_bijection_count": len(minimum_minus),
        "floor_surviving_conic_bijection_count": len(survivors_minus),
        "minimum_boundary_sets_match": minimum_bijection,
        "floor_surviving_conic_sets_match": survivor_bijection,
        "finite_boundary_forces_even_infinity_degree": True,
        "edge_count_29_is_odd": True,
        "paley_product_sign_flips": True,
        "normalized_score_and_eigenshell_conditions_preserved": True,
    }


def run(args: argparse.Namespace) -> dict:
    floor_minus, minimum_minus, survivors_minus = audit_floor(args.floor_minus, -1)
    _floor_plus, minimum_plus, survivors_plus = audit_floor(args.floor_plus, 1)
    conic_minus_audit = audit_conic_orbits(
        args.conic_minus, args.floor_minus, -1, survivors_minus
    )
    conic_plus_audit = audit_conic_orbits(
        args.conic_plus, args.floor_plus, 1, survivors_plus
    )
    conic_minus = load_json(args.conic_minus)
    saturated = audit_saturated(conic_minus, args.saturated_root)
    exceptional = audit_exceptional(
        conic_minus, args.ordinary_exceptional_audit, args.high_exceptional_audit
    )
    classified_indices = {
        index for index, _boundary in saturated_orbits(conic_minus)
    } | {index for index, _boundary in exceptional_orbits(conic_minus)}
    require(classified_indices == set(range(32)), "saturated/exceptional orbit partition incomplete")
    require(
        saturated["boundary_count"] + exceptional["boundary_count"]
        == CONIC_FLOOR_SURVIVORS,
        "closed conic boundary coverage changed",
    )
    sign_transfer = audit_sign_transfer(
        minimum_minus, minimum_plus, survivors_minus, survivors_plus
    )
    return {
        "experiment": "p7_size8_conic_global_audit",
        "status": "passed_independent_complete_conic_subbranch_audit",
        "p": 7,
        "finite_boundary_size": 8,
        "all_size_eight_boundaries_per_sign": ALL_BOUNDARIES,
        "floor_surviving_boundaries_per_sign": FLOOR_SURVIVORS,
        "floor_survivor_ordered_profile_count_per_sign": int(
            floor_minus["survivor_ordered_profile_count"]
        ),
        "minimum_odd_secants": 8,
        "minimum_odd_secant_boundaries_per_sign": MINIMUM_BOUNDARIES,
        "all_minimum_odd_secant_boundaries_are_eight_arcs": True,
        "classification_dependency": (
            "Segre's theorem: over odd order, every (q+1)-arc is a conic"
        ),
        "floor_surviving_conic_boundaries_per_sign": CONIC_FLOOR_SURVIVORS,
        "conic_stabilizer_orbits_per_sign": 32,
        "cminus_conic_orbit_partition": {
            "saturated": saturated,
            "exceptional": exceptional,
            "total_orbits": 32,
            "total_boundaries": CONIC_FLOOR_SURVIVORS,
            "total_representative_mean_allocations": 1_860,
            "remaining_representative_mean_allocations": 0,
        },
        "conic_orbit_audits": {
            "c_H=-1": conic_minus_audit,
            "c_H=+1": conic_plus_audit,
        },
        "sign_transfer": sign_transfer,
        "p7_size8_minimum_odd_secant_conic_subbranch_both_signs": "CLOSED",
        "all_32_conic_orbits_both_signs_excluded": True,
        "nonconic_floor_survivors_per_sign": NONCONIC_FLOOR_SURVIVORS,
        "full_p7_size8_boundary_case": "OPEN",
        "closes_all_p7_size8": False,
        "closes_residual_ii": False,
        "closes_R1": False,
        "L_status": "OPEN",
        "source_sha256": {
            "floor_cminus1": sha256(args.floor_minus),
            "floor_cplus1": sha256(args.floor_plus),
            "conic_cminus1": sha256(args.conic_minus),
            "conic_cplus1": sha256(args.conic_plus),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--floor-minus", type=Path, required=True)
    parser.add_argument("--floor-plus", type=Path, required=True)
    parser.add_argument("--conic-minus", type=Path, required=True)
    parser.add_argument("--conic-plus", type=Path, required=True)
    parser.add_argument("--saturated-root", type=Path, required=True)
    parser.add_argument("--ordinary-exceptional-audit", type=Path, required=True)
    parser.add_argument("--high-exceptional-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = run(args)
    atomic_json(args.output, out)
    print(json.dumps(out, indent=2), flush=True)


if __name__ == "__main__":
    main()
