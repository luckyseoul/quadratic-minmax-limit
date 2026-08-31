#!/usr/bin/env python3
"""Exact runner for the proposed support-330 catalog on ``J(13,7)``.

The model reuses Proposition 15.738's 1,638 exact third-difference
identities.  It asks for a Boolean degree-at-most-two function of support
330, anchors one support point, and excludes the 70 anchored members of the
candidate catalog:

* 78 omitted-pair functions ``(1-x_i)(1-x_j)`` in the full catalog;
* 286 all-equal-triple functions in the full catalog.

An ``INFEASIBLE`` result for the unsharded model proves that these 364
supports exhaust the support-330 Boolean quadratics.  Until that exact
status is obtained, this runner makes no classification claim.

Optional shards partition the *unknown Boolean assignment* by a fixed
binary prefix.  Every shard retains all 70 no-goods.  Splitting the no-goods
would not be a valid exhaustion and is deliberately unsupported.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from ortools import __version__ as ORTOOLS_VERSION
from ortools.sat.python import cp_model


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15738 import (  # noqa: E402
    middle_slice_points,
    pair_coordinates,
    selected_third_difference_identities,
    third_difference_rank_certificate,
)
from io_atomic import write_json_atomic  # noqa: E402


EXPERIMENT = "p13_support330_boolean_classifier"
P = 13
N = 13
K = 7
DOMAIN_SIZE = 1716
SUPPORT_SIZE = 330
PAIR_COUNT = 78
TRIPLE_COUNT = 286
CATALOG_SIZE = PAIR_COUNT + TRIPLE_COUNT
ANCHORED_PAIR_COUNT = 15
ANCHORED_TRIPLE_COUNT = 55
ANCHORED_NOGOOD_COUNT = ANCHORED_PAIR_COUNT + ANCHORED_TRIPLE_COUNT
IDENTITY_COUNT = 1638
ANCHOR_INDEX = 0
MAX_SHARD_COUNT = 256
DEFAULT_MAX_TIME_SECONDS = 1800.0
DEFAULT_WORKERS = 32

Support = tuple[int, ...]
SparseRow = tuple[tuple[int, int], ...]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ArithmeticError(message)


def positive_worker_count(raw: str | int) -> int:
    """Parse a strictly positive CP-SAT worker count."""
    value = int(raw)
    if value < 1:
        raise argparse.ArgumentTypeError("workers must be positive")
    return value


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _support_digest(supports: Iterable[Support]) -> str:
    payload = ";".join(",".join(map(str, support)) for support in supports)
    return _sha256_text(payload)


def _proto_hash(model: cp_model.CpModel) -> str:
    return _sha256_text(str(model.Proto()))


@lru_cache(maxsize=1)
def support330_candidate_catalog() -> tuple[
    tuple[Support, ...], tuple[dict[str, object], ...]
]:
    """Build the 364 proposed support-330 Boolean quadratics."""
    points = middle_slice_points()
    _require(len(points) == DOMAIN_SIZE, "J(13,7) point count changed")

    supports: list[Support] = []
    forms: list[dict[str, object]] = []
    for i, j in combinations(range(N), 2):
        support = tuple(
            index
            for index, point in enumerate(points)
            if i not in point and j not in point
        )
        supports.append(support)
        forms.append(
            {
                "family": "omitted_pair",
                "coordinates": [i, j],
                "polynomial": "(1-x_i)(1-x_j)",
                "support_size": len(support),
                "lift_offset_increment": -1,
                "total_b2_target_offset": 3,
            }
        )

    for i, j, k in combinations(range(N), 3):
        support = tuple(
            index
            for index, point in enumerate(points)
            if sum(coordinate in point for coordinate in (i, j, k)) in (0, 3)
        )
        supports.append(support)
        forms.append(
            {
                "family": "all_equal_triple",
                "coordinates": [i, j, k],
                "polynomial": (
                    "1-x_i-x_j-x_k+x_i*x_j+x_i*x_k+x_j*x_k"
                ),
                "support_size": len(support),
                "lift_offset_increment": 1,
                "total_b2_target_offset": 5,
            }
        )

    _require(len(supports) == CATALOG_SIZE, "candidate catalog count changed")
    _require(
        len(set(supports)) == CATALOG_SIZE,
        "candidate catalog supports are not distinct",
    )
    _require(
        all(len(support) == SUPPORT_SIZE for support in supports),
        "candidate support size changed",
    )
    return tuple(supports), tuple(forms)


def _row_vanishes_on_support(row: SparseRow, support: set[int]) -> bool:
    return sum(sign * int(index in support) for index, sign in row) == 0


def _catalog_symmetry_certificate(
    supports: Sequence[Support],
) -> dict[str, object]:
    """Check the catalog under generators of the full coordinate action.

    The twelve adjacent transpositions generate ``S_13``.  We check every
    catalog support under each generator and separately compute the orbit of
    the chosen anchor point.  Thus the support-point anchor is an executable
    symmetry reduction, not an asserted normalization.
    """
    points = middle_slice_points()
    point_to_index = {point: index for index, point in enumerate(points)}
    support_lookup = set(supports)
    generator_maps: list[tuple[int, ...]] = []
    for coordinate in range(N - 1):
        swapped_indices: list[int] = []
        for point in points:
            swapped = tuple(
                sorted(
                    coordinate + 1
                    if value == coordinate
                    else coordinate
                    if value == coordinate + 1
                    else value
                    for value in point
                )
            )
            swapped_indices.append(point_to_index[swapped])
        generator_maps.append(tuple(swapped_indices))

    closure_by_generator = []
    for point_map in generator_maps:
        closed = all(
            tuple(sorted(point_map[index] for index in support)) in support_lookup
            for support in supports
        )
        closure_by_generator.append(closed)

    orbit = {ANCHOR_INDEX}
    frontier = [ANCHOR_INDEX]
    while frontier:
        point_index = frontier.pop()
        for point_map in generator_maps:
            image = point_map[point_index]
            if image not in orbit:
                orbit.add(image)
                frontier.append(image)

    all_closed = all(closure_by_generator)
    anchor_orbit_is_full = len(orbit) == DOMAIN_SIZE
    return {
        "generator_family": "adjacent_coordinate_transpositions",
        "generator_count": len(generator_maps),
        "catalog_images_checked": len(generator_maps) * len(supports),
        "catalog_closed_by_generator": closure_by_generator,
        "catalog_closed_under_generated_S13": all_closed,
        "anchor_orbit_size": len(orbit),
        "anchor_orbit_is_all_J13_7": anchor_orbit_is_full,
        "anchor_reduction_verified": all_closed and anchor_orbit_is_full,
    }


@lru_cache(maxsize=1)
def catalog_arithmetic_certificate() -> dict[str, object]:
    """Verify candidate counts, offsets, anchor filter, and exact identities.

    This certifies the proposed catalog entries and the exact quadratic
    model.  It does not certify that the catalog is exhaustive.
    """
    points = middle_slice_points()
    pairs = pair_coordinates()
    identities, _descriptors, _examined = selected_third_difference_identities()
    rank = third_difference_rank_certificate()
    supports, forms = support330_candidate_catalog()
    support_sets = tuple(map(set, supports))
    symmetry = _catalog_symmetry_certificate(supports)

    family_counts = dict(Counter(str(form["family"]) for form in forms))
    anchored_indices = tuple(
        index for index, support in enumerate(supports) if ANCHOR_INDEX in support
    )
    anchored_family_counts = dict(
        Counter(str(forms[index]["family"]) for index in anchored_indices)
    )

    omitted_target_identity = all(
        4 * (1 - xi) * (1 - xj)
        == 1 - (2 * xi - 1) - (2 * xj - 1)
        + (2 * xi - 1) * (2 * xj - 1)
        for xi, xj in product((0, 1), repeat=2)
    )
    triple_target_identity = all(
        4
        * (
            1
            - xi
            - xj
            - xk
            + xi * xj
            + xi * xk
            + xj * xk
        )
        == 1
        + (2 * xi - 1) * (2 * xj - 1)
        + (2 * xi - 1) * (2 * xk - 1)
        + (2 * xj - 1) * (2 * xk - 1)
        for xi, xj, xk in product((0, 1), repeat=3)
    )
    every_support_satisfies_identities = all(
        _row_vanishes_on_support(row, support)
        for support in support_sets
        for row in identities
    )

    expected_families = {
        "omitted_pair": PAIR_COUNT,
        "all_equal_triple": TRIPLE_COUNT,
    }
    expected_anchored = {
        "omitted_pair": ANCHORED_PAIR_COUNT,
        "all_equal_triple": ANCHORED_TRIPLE_COUNT,
    }
    full_digest = _support_digest(supports)
    anchored_digest = _support_digest(supports[index] for index in anchored_indices)
    verified = bool(
        len(points) == DOMAIN_SIZE
        and len(pairs) == PAIR_COUNT
        and len(identities) == IDENTITY_COUNT
        and rank["proved"]
        and rank["exact_real_rank"] == IDENTITY_COUNT
        and rank["exact_real_nullity"] == PAIR_COUNT
        and family_counts == expected_families
        and anchored_family_counts == expected_anchored
        and len(anchored_indices) == ANCHORED_NOGOOD_COUNT
        and symmetry["generator_count"] == N - 1
        and symmetry["anchor_reduction_verified"]
        and omitted_target_identity
        and triple_target_identity
        and every_support_satisfies_identities
    )
    _require(verified, "support-330 candidate arithmetic changed")
    return {
        "slice": "J(13,7)",
        "point_count": len(points),
        "pair_coordinate_count": len(pairs),
        "third_difference_identity_count": len(identities),
        "third_difference_identity_sha256": rank["selected_identity_sha256"],
        "identity_nullspace_is_exact_degree_at_most_two_space": bool(
            rank["nullspace_equals_degree_at_most_two_evaluation_space"]
        ),
        "support_size": SUPPORT_SIZE,
        "support_density": "330/1716=5/26",
        "full_candidate_count": len(supports),
        "distinct_candidate_count": len(set(supports)),
        "family_counts": family_counts,
        "candidate_catalog_sha256": full_digest,
        "support_point_anchor_index": ANCHOR_INDEX,
        "support_point_anchor": list(points[ANCHOR_INDEX]),
        "support_point_anchor_is_wlog_by_S13_transitivity": symmetry[
            "anchor_reduction_verified"
        ],
        "candidate_catalog_is_S13_invariant": symmetry[
            "catalog_closed_under_generated_S13"
        ],
        "S13_symmetry_certificate": symmetry,
        "anchored_candidate_count": len(anchored_indices),
        "anchored_family_counts": anchored_family_counts,
        "anchored_candidate_catalog_sha256": anchored_digest,
        "omitted_pair_support_formula": "C(11,7)=330",
        "all_equal_triple_support_formula": "C(10,7)+C(10,4)=120+210=330",
        "omitted_pair_anchor_formula": "C(6,2)=15",
        "all_equal_triple_anchor_formula": "C(7,3)+C(6,3)=35+20=55",
        "lift_target_identities_verified": {
            "omitted_pair": omitted_target_identity,
            "all_equal_triple": triple_target_identity,
        },
        "b2_combined_target_offsets": {
            "omitted_pair": 3,
            "all_equal_triple": 5,
        },
        "every_candidate_satisfies_all_1638_identities": (
            every_support_satisfies_identities
        ),
        "candidate_catalog_verified": verified,
        "catalog_exhaustiveness_claimed": False,
    }


def shard_partition(shard_index: int, shard_count: int) -> dict[str, object]:
    """Return one exact binary-prefix partition of the anchored model."""
    if shard_count < 1 or shard_count > MAX_SHARD_COUNT:
        raise ValueError(f"shard-count must lie in [1,{MAX_SHARD_COUNT}]")
    if shard_count & (shard_count - 1):
        raise ValueError("shard-count must be a power of two")
    if not 0 <= shard_index < shard_count:
        raise ValueError("need 0 <= shard-index < shard-count")
    width = shard_count.bit_length() - 1
    point_indices = tuple(range(1, width + 1))
    assignment = tuple((shard_index >> bit) & 1 for bit in range(width))
    return {
        "scheme": "binary_prefix_of_anchored_boolean_assignment",
        "shard_index": shard_index,
        "shard_count": shard_count,
        "prefix_width": width,
        "fixed_point_indices": list(point_indices),
        "fixed_values": list(assignment),
        "all_70_nogoods_retained_in_every_shard": True,
        "nogoods_are_not_sharded": True,
        "partition_is_disjoint_and_exhaustive": True,
    }


def all_shard_assignments(shard_count: int) -> tuple[tuple[int, ...], ...]:
    """Expose the exact prefix partition for focused tests/audits."""
    return tuple(
        tuple(shard_partition(index, shard_count)["fixed_values"])
        for index in range(shard_count)
    )


def build_classifier_model(
    shard_index: int = 0,
    shard_count: int = 1,
) -> tuple[cp_model.CpModel, list[cp_model.IntVar], dict[str, object]]:
    """Build one exact support-330 classifier model without solving it."""
    catalog = catalog_arithmetic_certificate()
    identities, _descriptors, _examined = selected_third_difference_identities()
    supports, forms = support330_candidate_catalog()
    anchored = tuple(
        (support, forms[index])
        for index, support in enumerate(supports)
        if ANCHOR_INDEX in support
    )
    partition = shard_partition(shard_index, shard_count)

    model = cp_model.CpModel()
    values = [model.NewBoolVar(f"f_{index}") for index in range(DOMAIN_SIZE)]
    for ordinal, row in enumerate(identities):
        model.Add(
            sum(sign * values[index] for index, sign in row) == 0
        ).WithName(f"third_difference_{ordinal}")
    model.Add(sum(values) == SUPPORT_SIZE).WithName("support_size_330")
    model.Add(values[ANCHOR_INDEX] == 1).WithName("support_anchor")
    for ordinal, (support, form) in enumerate(anchored):
        family = str(form["family"])
        model.Add(
            sum(values[index] for index in support) <= SUPPORT_SIZE - 1
        ).WithName(f"nogood_{family}_{ordinal}")
    for point_index, fixed in zip(
        partition["fixed_point_indices"], partition["fixed_values"]
    ):
        model.Add(values[int(point_index)] == int(fixed)).WithName(
            f"shard_prefix_f_{point_index}_{fixed}"
        )

    validation = model.Validate()
    _require(not validation, f"invalid support-330 model: {validation}")
    expected_constraints = IDENTITY_COUNT + 2 + ANCHORED_NOGOOD_COUNT + int(
        partition["prefix_width"]
    )
    _require(
        len(model.Proto().variables) == DOMAIN_SIZE
        and len(model.Proto().constraints) == expected_constraints,
        "support-330 model dimensions changed",
    )
    metadata = {
        "boolean_variable_count": len(model.Proto().variables),
        "constraint_count": len(model.Proto().constraints),
        "third_difference_equality_count": IDENTITY_COUNT,
        "support_equation_count": 1,
        "anchor_equation_count": 1,
        "anchored_nogood_count": len(anchored),
        "shard_fixing_count": int(partition["prefix_width"]),
        "model_validation": validation,
        "model_textproto_sha256": _proto_hash(model),
        "candidate_catalog_sha256": catalog["candidate_catalog_sha256"],
        "anchored_candidate_catalog_sha256": catalog[
            "anchored_candidate_catalog_sha256"
        ],
        "third_difference_identity_sha256": catalog[
            "third_difference_identity_sha256"
        ],
        "partition": partition,
        "constraint_scope": [
            "1716 Boolean values",
            "all 1638 exact third-difference identities",
            "exact support size 330",
            "one support-point anchor",
            "all 70 anchored candidate no-goods",
            "optional exact binary-prefix shard fixings",
        ],
    }
    return model, values, metadata


def _witness_audit(
    support: Sequence[int],
    metadata: dict[str, object],
) -> dict[str, object]:
    identities, _descriptors, _examined = selected_third_difference_identities()
    known_supports, _forms = support330_candidate_catalog()
    support_tuple = tuple(sorted(int(index) for index in support))
    support_set = set(support_tuple)
    partition = metadata["partition"]
    all_identities = all(
        _row_vanishes_on_support(row, support_set) for row in identities
    )
    shard_matches = all(
        int(point_index in support_set) == int(fixed)
        for point_index, fixed in zip(
            partition["fixed_point_indices"], partition["fixed_values"]
        )
    )
    outside_catalog = support_tuple not in set(known_supports)
    verified = bool(
        len(support_tuple) == SUPPORT_SIZE
        and len(support_set) == SUPPORT_SIZE
        and ANCHOR_INDEX in support_set
        and all_identities
        and shard_matches
        and outside_catalog
    )
    _require(verified, "solver witness failed direct audit")
    return {
        "support_indices": list(support_tuple),
        "support_sha256": _support_digest((support_tuple,)),
        "support_size": len(support_tuple),
        "anchor_present": ANCHOR_INDEX in support_set,
        "all_1638_identities_hold": all_identities,
        "shard_prefix_matches": shard_matches,
        "outside_364_candidate_catalog": outside_catalog,
        "verified": verified,
    }


def gpu_catalog_cross_check(requested: bool) -> dict[str, object]:
    """Optionally replay catalog sizes/identities on a CUDA device.

    This is diagnostic only.  The CPU integer checks in
    :func:`catalog_arithmetic_certificate` remain the proof premise.
    """
    base = {
        "requested": requested,
        "proof_premise": False,
        "cpu_integer_verification_remains_authoritative": True,
    }
    if not requested:
        return {**base, "status": "NOT_REQUESTED", "backend": None}

    backend = None
    module = None
    diagnostics: list[str] = []
    try:
        import cupy as cp

        if int(cp.cuda.runtime.getDeviceCount()) > 0:
            backend = "cupy"
            module = cp
    except Exception as exc:  # pragma: no cover - environment dependent
        diagnostics.append(f"cupy: {type(exc).__name__}: {exc}")
    if backend is None:
        try:
            import torch

            if torch.cuda.is_available():
                backend = "torch"
                module = torch
        except Exception as exc:  # pragma: no cover - environment dependent
            diagnostics.append(f"torch: {type(exc).__name__}: {exc}")
    if backend is None:
        return {
            **base,
            "status": "SKIPPED_NO_CUDA_BACKEND",
            "backend": None,
            "diagnostics": diagnostics,
        }

    identities, _descriptors, _examined = selected_third_difference_identities()
    supports, _forms = support330_candidate_catalog()
    dense_identities = np.zeros((IDENTITY_COUNT, DOMAIN_SIZE), dtype=np.int32)
    for row_index, row in enumerate(identities):
        for point_index, sign in row:
            dense_identities[row_index, point_index] = sign

    try:
        maximum_residual = 0
        observed_sizes: set[int] = set()
        batch_size = 32
        if backend == "cupy":
            cp = module
            identity_device = cp.asarray(dense_identities)
            for start in range(0, len(supports), batch_size):
                batch = np.zeros(
                    (min(batch_size, len(supports) - start), DOMAIN_SIZE),
                    dtype=np.int32,
                )
                for row_index, support in enumerate(
                    supports[start : start + batch_size]
                ):
                    batch[row_index, list(support)] = 1
                batch_device = cp.asarray(batch)
                residual = identity_device @ batch_device.T
                maximum_residual = max(
                    maximum_residual, int(cp.max(cp.abs(residual)).get())
                )
                observed_sizes.update(
                    int(value) for value in cp.sum(batch_device, axis=1).get()
                )
            arithmetic = "CUDA int32 matrix multiplication"
        else:  # torch CUDA diagnostic; float32 is never a proof premise.
            torch = module
            identity_device = torch.as_tensor(
                dense_identities, device="cuda", dtype=torch.float32
            )
            for start in range(0, len(supports), batch_size):
                batch = np.zeros(
                    (min(batch_size, len(supports) - start), DOMAIN_SIZE),
                    dtype=np.float32,
                )
                for row_index, support in enumerate(
                    supports[start : start + batch_size]
                ):
                    batch[row_index, list(support)] = 1.0
                batch_device = torch.as_tensor(batch, device="cuda")
                residual = identity_device @ batch_device.T
                maximum_residual = max(
                    maximum_residual, int(residual.abs().max().item())
                )
                observed_sizes.update(
                    int(value) for value in batch_device.sum(dim=1).cpu().tolist()
                )
            arithmetic = "CUDA float32 diagnostic matrix multiplication"
        passed = maximum_residual == 0 and observed_sizes == {SUPPORT_SIZE}
        return {
            **base,
            "status": "PASSED_NONPROOF_CROSS_CHECK" if passed else "FAILED",
            "backend": backend,
            "arithmetic": arithmetic,
            "candidate_count_checked": len(supports),
            "identity_count_checked": len(identities),
            "observed_support_sizes": sorted(observed_sizes),
            "maximum_identity_residual": maximum_residual,
            "passed": passed,
            "diagnostics": diagnostics,
        }
    except Exception as exc:  # pragma: no cover - environment dependent
        return {
            **base,
            "status": "AVAILABLE_BUT_CROSS_CHECK_FAILED",
            "backend": backend,
            "diagnostics": diagnostics
            + [f"cross-check: {type(exc).__name__}: {exc}"],
        }


def request_fingerprint(metadata: dict[str, object]) -> str:
    partition = metadata["partition"]
    payload = {
        "experiment": EXPERIMENT,
        "model_textproto_sha256": metadata["model_textproto_sha256"],
        "candidate_catalog_sha256": metadata["candidate_catalog_sha256"],
        "anchored_candidate_catalog_sha256": metadata[
            "anchored_candidate_catalog_sha256"
        ],
        "third_difference_identity_sha256": metadata[
            "third_difference_identity_sha256"
        ],
        "shard_index": partition["shard_index"],
        "shard_count": partition["shard_count"],
        "ortools_version": ORTOOLS_VERSION,
    }
    return _sha256_text(json.dumps(payload, sort_keys=True))


def solve_built_model(
    model: cp_model.CpModel,
    values: Sequence[cp_model.IntVar],
    metadata: dict[str, object],
    *,
    max_time_seconds: float | None,
    workers: int,
    log_search_progress: bool,
    gpu_cross_check_requested: bool,
) -> dict[str, object]:
    """Solve one exact model; one worker additionally makes search deterministic."""
    catalog = catalog_arithmetic_certificate()
    gpu = gpu_catalog_cross_check(gpu_cross_check_requested)
    if workers < 1:
        raise ValueError("workers must be positive")
    solver = cp_model.CpSolver()
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = 0
    solver.parameters.randomize_search = False
    solver.parameters.cp_model_presolve = True
    solver.parameters.symmetry_level = 3
    solver.parameters.log_search_progress = bool(log_search_progress)
    if max_time_seconds is not None:
        if max_time_seconds <= 0:
            raise ValueError("max-time-seconds must be positive")
        solver.parameters.max_time_in_seconds = float(max_time_seconds)

    started = time.perf_counter()
    status = solver.Solve(model)
    elapsed = time.perf_counter() - started
    status_name = solver.StatusName(status)
    infeasible = status == cp_model.INFEASIBLE
    has_witness = status in (cp_model.FEASIBLE, cp_model.OPTIMAL)
    witness = None
    if has_witness:
        support = [
            index for index, value in enumerate(values) if solver.Value(value)
        ]
        witness = _witness_audit(support, metadata)

    partition = metadata["partition"]
    full_unsharded_infeasibility = bool(
        infeasible and int(partition["shard_count"]) == 1
    )
    if full_unsharded_infeasibility:
        result_status = "COMPLETE_EXACT_FULL_MODEL_INFEASIBILITY"
    elif infeasible:
        result_status = "COMPLETE_EXACT_SHARD_INFEASIBILITY"
    elif has_witness:
        result_status = "EXACT_COUNTEREXAMPLE_TO_EXPECTED_CATALOG"
    elif status == cp_model.UNKNOWN:
        result_status = "INCOMPLETE_SOLVER_RUN"
    else:
        result_status = f"UNEXPECTED_SOLVER_STATUS_{status_name}"

    response = solver.ResponseProto()
    return {
        "experiment": EXPERIMENT,
        "scope": "Boolean degree-at-most-two functions of support 330 on J(13,7)",
        "result_status": result_status,
        "request_fingerprint": request_fingerprint(metadata),
        "candidate_catalog": catalog,
        "model": metadata,
        "solver": {
            "name": "OR-Tools CP-SAT",
            "version": ORTOOLS_VERSION,
            "status": status_name,
            "status_code": int(status),
            "num_search_workers": workers,
            "deterministic_one_worker_mode": workers == 1,
            "random_seed": 0,
            "randomize_search": False,
            "cp_model_presolve": True,
            "symmetry_level": 3,
            "max_time_seconds": max_time_seconds,
            "wall_time_seconds": solver.WallTime(),
            "outer_elapsed_seconds": elapsed,
            "deterministic_time": float(response.deterministic_time),
            "branches": solver.NumBranches(),
            "conflicts": solver.NumConflicts(),
            "solution_count_lower_bound": 1 if has_witness else 0,
        },
        "classification": {
            "shard_infeasible": infeasible,
            "full_catalog_exhaustive": full_unsharded_infeasibility,
            "counterexample_found": has_witness,
            "incomplete": status == cp_model.UNKNOWN,
            "theorem_claim_requires_full_unsharded_or_complete_shard_cover": True,
        },
        "witness": witness,
        "gpu_catalog_cross_check": gpu,
        "resume_semantics": {
            "atomic_output": True,
            "terminal_exact_result_may_be_reused": True,
            "unknown_search_state_is_not_checkpointed": True,
            "unknown_result_is_rerun_from_scratch": True,
        },
    }


def build_only_payload(metadata: dict[str, object]) -> dict[str, object]:
    """Return a non-theorem payload for construction-only preflight."""
    return {
        "experiment": EXPERIMENT,
        "scope": "Boolean degree-at-most-two functions of support 330 on J(13,7)",
        "result_status": "MODEL_BUILT_NOT_SOLVED",
        "request_fingerprint": request_fingerprint(metadata),
        "candidate_catalog": catalog_arithmetic_certificate(),
        "model": metadata,
        "solver": None,
        "classification": {
            "shard_infeasible": False,
            "full_catalog_exhaustive": False,
            "counterexample_found": False,
            "incomplete": True,
            "theorem_claim_requires_full_unsharded_or_complete_shard_cover": True,
        },
        "witness": None,
        "gpu_catalog_cross_check": {
            "requested": False,
            "status": "NOT_RUN_IN_BUILD_ONLY_MODE",
            "proof_premise": False,
        },
        "resume_semantics": {
            "atomic_output": True,
            "terminal_exact_result_may_be_reused": True,
            "unknown_search_state_is_not_checkpointed": True,
            "unknown_result_is_rerun_from_scratch": True,
        },
    }


def _terminal_exact_payload(payload: dict[str, object]) -> bool:
    solver = payload.get("solver")
    return isinstance(solver, dict) and solver.get("status") in {
        "INFEASIBLE",
        "FEASIBLE",
        "OPTIMAL",
    }


def run_classifier(
    *,
    output: Path,
    shard_index: int,
    shard_count: int,
    max_time_seconds: float | None,
    workers: int,
    build_only: bool,
    resume: bool,
    gpu_cross_check_requested: bool,
    log_search_progress: bool,
) -> dict[str, object]:
    """Build, optionally resume, solve, and atomically write one shard."""
    if workers < 1:
        raise ValueError("workers must be positive")
    model, values, metadata = build_classifier_model(shard_index, shard_count)
    fingerprint = request_fingerprint(metadata)
    if resume and output.is_file():
        previous = json.loads(output.read_text(encoding="utf-8"))
        if previous.get("experiment") != EXPERIMENT:
            raise ValueError("resume output belongs to another experiment")
        if previous.get("request_fingerprint") != fingerprint:
            raise ValueError("resume output does not match the requested exact model")
        if _terminal_exact_payload(previous):
            return previous

    if build_only:
        payload = build_only_payload(metadata)
    else:
        payload = solve_built_model(
            model,
            values,
            metadata,
            max_time_seconds=max_time_seconds,
            workers=workers,
            log_search_progress=log_search_progress,
            gpu_cross_check_requested=gpu_cross_check_requested,
        )
    write_json_atomic(output, payload)
    return payload


def merge_shard_payloads(payloads: Sequence[dict[str, object]]) -> dict[str, object]:
    """Merge a complete exact binary-prefix cover without rerunning CP-SAT."""
    if not payloads:
        raise ValueError("need at least one shard payload")
    first_model = payloads[0].get("model")
    if not isinstance(first_model, dict):
        raise ValueError("shard payload lacks model metadata")
    first_partition = first_model.get("partition")
    if not isinstance(first_partition, dict):
        raise ValueError("shard payload lacks partition metadata")
    shard_count = int(first_partition["shard_count"])
    if shard_count <= 1:
        raise ValueError("merge mode requires a genuinely sharded run")

    common_keys = (
        "candidate_catalog_sha256",
        "anchored_candidate_catalog_sha256",
        "third_difference_identity_sha256",
        "boolean_variable_count",
        "third_difference_equality_count",
        "anchored_nogood_count",
    )
    expected_common = {key: first_model[key] for key in common_keys}
    by_index: dict[int, dict[str, object]] = {}
    for payload in payloads:
        if payload.get("experiment") != EXPERIMENT:
            raise ValueError("merge input belongs to another experiment")
        model = payload.get("model")
        if not isinstance(model, dict):
            raise ValueError("merge input lacks model metadata")
        if {key: model[key] for key in common_keys} != expected_common:
            raise ValueError("merge inputs use different catalog/identity models")
        partition = model.get("partition")
        if not isinstance(partition, dict):
            raise ValueError("merge input lacks partition metadata")
        if (
            int(partition["shard_count"]) != shard_count
            or partition.get("scheme")
            != "binary_prefix_of_anchored_boolean_assignment"
            or not partition.get("all_70_nogoods_retained_in_every_shard")
            or not partition.get("partition_is_disjoint_and_exhaustive")
        ):
            raise ValueError("merge input has an incompatible shard partition")
        index = int(partition["shard_index"])
        expected_partition = shard_partition(index, shard_count)
        if partition != expected_partition:
            raise ValueError("merge input does not match its exact binary-prefix shard")
        if (
            int(model.get("shard_fixing_count", -1))
            != int(expected_partition["prefix_width"])
            or int(model.get("constraint_count", -1))
            != IDENTITY_COUNT
            + 2
            + ANCHORED_NOGOOD_COUNT
            + int(expected_partition["prefix_width"])
        ):
            raise ValueError("merge input has inconsistent shard model dimensions")
        if index in by_index:
            raise ValueError("duplicate shard index")
        by_index[index] = payload
    if set(by_index) != set(range(shard_count)):
        raise ValueError("shard inputs do not form a complete cover")

    statuses = {}
    for index in range(shard_count):
        solver = by_index[index].get("solver")
        statuses[index] = (
            str(solver.get("status")) if isinstance(solver, dict) else "None"
        )
    all_infeasible = all(status == "INFEASIBLE" for status in statuses.values())
    counterexamples = []
    for index, status in statuses.items():
        if status in {"FEASIBLE", "OPTIMAL"}:
            witness = by_index[index].get("witness")
            if not isinstance(witness, dict) or witness.get("verified") is not True:
                raise ValueError("feasible shard result lacks a verified witness")
            counterexamples.append(witness)
    if all_infeasible:
        result_status = "COMPLETE_EXACT_SHARD_COVER_INFEASIBILITY"
    elif counterexamples:
        result_status = "EXACT_COUNTEREXAMPLE_TO_EXPECTED_CATALOG"
    else:
        result_status = "INCOMPLETE_SHARD_COVER_RESULTS"
    return {
        "experiment": EXPERIMENT,
        "scope": "merged exact binary-prefix shard cover",
        "result_status": result_status,
        "candidate_catalog_sha256": expected_common["candidate_catalog_sha256"],
        "anchored_candidate_catalog_sha256": expected_common[
            "anchored_candidate_catalog_sha256"
        ],
        "third_difference_identity_sha256": expected_common[
            "third_difference_identity_sha256"
        ],
        "partition": {
            "scheme": "binary_prefix_of_anchored_boolean_assignment",
            "shard_count": shard_count,
            "covered_shard_indices": list(range(shard_count)),
            "partition_is_disjoint_and_exhaustive": True,
            "all_70_nogoods_retained_in_every_shard": True,
        },
        "solver_status_by_shard": statuses,
        "classification": {
            "full_catalog_exhaustive": all_infeasible,
            "counterexample_found": bool(counterexamples),
            "incomplete": not all_infeasible and not counterexamples,
        },
        "counterexample_witnesses": counterexamples,
        "source_result_fingerprints": [
            by_index[index]["request_fingerprint"] for index in range(shard_count)
        ],
    }


def _default_output(shard_index: int, shard_count: int) -> Path:
    if shard_count == 1:
        return ROOT / "evidence" / "p13_support330_boolean_classifier.json"
    width = max(2, len(str(shard_count - 1)))
    return (
        ROOT
        / "evidence"
        / (
            "p13_support330_boolean_classifier_"
            f"shard_{shard_index:0{width}d}_of_{shard_count}.json"
        )
    )


def _summary(payload: dict[str, object], output: Path) -> dict[str, object]:
    solver = payload.get("solver")
    model = payload.get("model")
    return {
        "output": str(output),
        "result_status": payload.get("result_status"),
        "solver_status": solver.get("status") if isinstance(solver, dict) else None,
        "wall_time_seconds": (
            solver.get("wall_time_seconds") if isinstance(solver, dict) else None
        ),
        "model_textproto_sha256": (
            model.get("model_textproto_sha256") if isinstance(model, dict) else None
        ),
        "full_catalog_exhaustive": payload.get("classification", {}).get(
            "full_catalog_exhaustive"
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument(
        "--max-time-seconds", type=float, default=DEFAULT_MAX_TIME_SECONDS
    )
    parser.add_argument(
        "--workers",
        type=positive_worker_count,
        default=DEFAULT_WORKERS,
        help=(
            "CP-SAT worker count (default: 32 for speed; use 1 for the "
            "deterministic one-worker exact mode)"
        ),
    )
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--gpu-cross-check", action="store_true")
    parser.add_argument("--log-search-progress", action="store_true")
    parser.add_argument("--merge-shards", type=Path, nargs="+")
    args = parser.parse_args()

    if args.merge_shards:
        if args.shard_index != 0 or args.shard_count != 1:
            parser.error("--merge-shards cannot be combined with shard selectors")
        output = args.output or (
            ROOT / "evidence" / "p13_support330_boolean_classifier_merged.json"
        )
        payloads = [
            json.loads(path.read_text(encoding="utf-8"))
            for path in args.merge_shards
        ]
        result = merge_shard_payloads(payloads)
        write_json_atomic(output, result)
        print(json.dumps(_summary(result, output), indent=2, sort_keys=True))
        return

    output = args.output or _default_output(args.shard_index, args.shard_count)
    result = run_classifier(
        output=output,
        shard_index=args.shard_index,
        shard_count=args.shard_count,
        max_time_seconds=args.max_time_seconds,
        workers=args.workers,
        build_only=args.build_only,
        resume=args.resume,
        gpu_cross_check_requested=args.gpu_cross_check,
        log_search_progress=args.log_search_progress,
    )
    print(json.dumps(_summary(result, output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
