#!/usr/bin/env python3
"""Exact full-quotient direct-semigroup CP-SAT model for z=7 case 0.

For each of the eight directions, the actual catalog row is
``floor + 2*L`` with integer ``L >= 0``.  This model imposes exactly
``K*L=0``, ``sum(L)=5*grade``, and ``L[x]<=grade``.  It derives the full
F3^6 x F7^21 map from the audited quotient context and enforces all 27
modular target equations.  INFEASIBLE is a rigorous full-torsion rejection;
FEASIBLE is necessary only and is accepted only after reconstructing and
auditing all eight complete-catalog rows.  UNKNOWN has no force.
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
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import p7_infinity7_positive_z7_semigroup_case_join as cpu_join  # noqa: E402
import p7_infinity7_positive_z7_semigroup_case_join_gpu as gpu_join  # noqa: E402

EXPERIMENT = "p7_infinity7_positive_z7_semigroup_case_cpsat"
FIXED_CASE_INDEX = 0
DIRECTIONS = 8
AMBIENT = 35
MOD7_COORDINATES = tuple(range(21))
MODULI = (3,) * 6 + (7,) * 21
WIDTH = 27


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def positive_float(value: str) -> float:
    parsed = float(value)
    if not math.isfinite(parsed) or parsed <= 0:
        raise argparse.ArgumentTypeError("value must be a positive finite number")
    return parsed


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return parsed


def nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be nonnegative")
    return parsed


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            digest.update(block)
    return digest.hexdigest()


def array_sha256(array: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def json_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def atomic_write(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    try:
        with temporary.open("w", encoding="utf-8") as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def quotient_block(context: dict[str, Any], direction: int, modulus: int) -> np.ndarray:
    complement = np.asarray(context["quotient_data"][modulus]["complement"], dtype=np.int64)
    block = slice(1 + AMBIENT * direction, 1 + AMBIENT * (direction + 1))
    return np.ascontiguousarray(complement[:, block] % modulus, dtype=np.int64)


def derive_maps(context: dict[str, Any]) -> tuple[tuple[np.ndarray, ...], dict]:
    maps = []
    records = []
    units = np.eye(AMBIENT, dtype=np.int64)
    for direction in range(DIRECTIONS):
        coefficient = np.ascontiguousarray(
            np.vstack(
                (
                    (-2 * quotient_block(context, direction, 3)) % 3,
                    (-2 * quotient_block(context, direction, 7)) % 7,
                )
            ),
            dtype=np.int64,
        )
        direct = cpu_join.semigroup.project_rows(
            -2 * units, direction, context["quotient_data"], MOD7_COORDINATES
        )
        require(
            coefficient.shape == (WIDTH, AMBIENT)
            and np.array_equal(coefficient, direct.T.astype(np.int64)),
            "derived compact map differs from the established full projector",
        )
        maps.append(coefficient)
        records.append(
            {
                "direction": direction,
                "shape": [WIDTH, AMBIENT],
                "sha256_int64": array_sha256(coefficient.astype("<i8", copy=False)),
                "equals_established_projector_on_all_coordinate_units": True,
            }
        )
    return tuple(maps), {
        "quotient": "F3^6 x F7^21",
        "records": records,
        "records_sha256": json_sha256(records),
        "all_27_maps_derived_exactly": True,
    }


def build_case(context: dict[str, Any], maps: tuple[np.ndarray, ...]) -> dict:
    require(len(context["targets"]) == 51, "audited target census changed")
    target_row = context["targets"][FIXED_CASE_INDEX]
    rebuilt = context["rebuilt"]
    orbit, leaf, system, _factory = cpu_join.old_join.validate_parent_survivor(
        target_row, rebuilt
    )
    grades = tuple(cpu_join.leaf_grade(orbit, leaf, d) for d in range(DIRECTIONS))
    require(max(grades) == 3 and sum(grades) == 14, "fixed case grade profile changed")

    anchor_rhs, _raw = cpu_join.affine.anchor_rhs_and_raw_syndromes(
        orbit, leaf, system, rebuilt["anchors"]
    )
    base = cpu_join.project_equation_vector(
        anchor_rhs, context["quotient_data"], MOD7_COORDINATES
    )
    target = np.ascontiguousarray(
        (-base.astype(np.int16)) % np.asarray(MODULI, dtype=np.int16), dtype=np.uint8
    )
    kernel = np.asarray(rebuilt["kernel_rows"], dtype=np.int64)
    floors, offsets, catalogs, records = [], [], [], []
    for direction, grade in enumerate(grades):
        mask = int(orbit["masks"][direction])
        mean = int(leaf["scaled_means"][direction])
        observed_grade, _floor_mean, floor = cpu_join.semigroup.excess_grade(mask, mean)
        require(observed_grade == grade, "grade derivations differ")
        floor = np.ascontiguousarray(floor, dtype=np.int64)
        anchor = rebuilt["anchors"].get(mask, mean)
        catalog = np.ascontiguousarray(
            cpu_join.affine.mapped_catalog(mask, mean), dtype=np.int64
        )
        require(
            len(catalog) == cpu_join.EXPECTED_GRADE_CATALOG_ROWS[grade],
            "complete catalog census changed",
        )
        require(not np.any((catalog - floor[None, :]) % 2), "catalog/floor parity changed")
        catalog_lifts = np.ascontiguousarray((catalog - floor[None, :]) // 2, dtype=np.int64)
        require(
            np.all((0 <= catalog_lifts) & (catalog_lifts <= grade))
            and np.all(catalog_lifts.sum(axis=1) == 5 * grade)
            and not np.any(kernel @ catalog_lifts.T),
            "complete direct catalog does not satisfy compact K/mass/bounds",
        )
        offset_source = np.ascontiguousarray(anchor - floor, dtype=np.int64)
        offset = cpu_join.semigroup.project_rows(
            offset_source[None, :], direction, context["quotient_data"], MOD7_COORDINATES
        )[0]
        direct = cpu_join.semigroup.project_rows(
            anchor[None, :] - catalog,
            direction,
            context["quotient_data"],
            MOD7_COORDINATES,
        )
        compact = np.empty_like(direct)
        compact[:, :6] = (offset[:6] + catalog_lifts @ maps[direction][:6].T) % 3
        compact[:, 6:] = (offset[6:] + catalog_lifts @ maps[direction][6:].T) % 7
        require(np.array_equal(direct, compact), "compact map fails a direct catalog row")
        floors.append(floor)
        offsets.append(np.ascontiguousarray(offset, dtype=np.uint8))
        catalogs.append(catalog)
        records.append(
            {
                "direction": direction,
                "grade": int(grade),
                "mask": mask,
                "scaled_mean": mean,
                "complete_catalog_rows": len(catalog),
                "catalog_sha256_int64": array_sha256(catalog.astype("<i8", copy=False)),
                "offset_sha256_uint8": array_sha256(offset),
                "all_catalog_rows_satisfy_K_mass_bounds_and_compact_map": True,
            }
        )
    key = str(target_row["case_key"])
    return {
        "case_key": key,
        "catalog_pattern": context["current_by_key"][key]["catalog_pattern"],
        "grades": grades,
        "target": target,
        "floors": tuple(floors),
        "offsets": tuple(offsets),
        "catalogs": tuple(catalogs),
        "audit": {
            "fixed_case_index": FIXED_CASE_INDEX,
            "case_key": key,
            "catalog_pattern": context["current_by_key"][key]["catalog_pattern"],
            "directional_grades": list(grades),
            "target_digits": target.tolist(),
            "target_sha256_uint8": array_sha256(target),
            "direction_records": records,
            "direction_records_sha256": json_sha256(records),
        },
    }


def solve(cp_model: Any, context: dict[str, Any], maps: tuple[np.ndarray, ...],
          case: dict, timeout: float, workers: int, seed: int) -> dict:
    kernel = np.asarray(context["rebuilt"]["kernel_rows"], dtype=np.int64)
    model = cp_model.CpModel()
    lifts = [
        [model.NewIntVar(0, int(case["grades"][d]), f"L_{d}_{x}") for x in range(AMBIENT)]
        for d in range(DIRECTIONS)
    ]
    for direction, grade in enumerate(case["grades"]):
        model.Add(sum(lifts[direction]) == 5 * int(grade))
        for row in kernel:
            model.Add(sum(int(row[x]) * lifts[direction][x] for x in range(AMBIENT)) == 0)

    quotient_variables = []
    for coordinate, modulus in enumerate(MODULI):
        wanted = int(case["target"][coordinate])
        constant = sum(int(case["offsets"][d][coordinate]) for d in range(DIRECTIONS))
        terms = []
        maximum = constant
        for direction, grade in enumerate(case["grades"]):
            for x, raw in enumerate(maps[direction][coordinate]):
                coefficient = int(raw)
                if coefficient:
                    terms.append(coefficient * lifts[direction][x])
                    maximum += coefficient * int(grade)
        quotient = model.NewIntVar(0, (maximum - wanted) // modulus, f"q_{coordinate}")
        model.Add(constant + sum(terms) == wanted + modulus * quotient)
        quotient_variables.append(quotient)

    proto = model.Proto()
    scale = {
        "L_variables": DIRECTIONS * AMBIENT,
        "quotient_variables": WIDTH,
        "total_variables": len(proto.variables),
        "kernel_equalities": DIRECTIONS * len(kernel),
        "mass_equalities": DIRECTIONS,
        "modular_equalities": WIDTH,
        "total_constraints": len(proto.constraints),
    }
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = timeout
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = seed
    status = solver.Solve(model)
    status_name = solver.StatusName(status)
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        status_class, rejected, necessary, unknown = "FEASIBLE", False, True, False
    elif status == cp_model.INFEASIBLE:
        status_class, rejected, necessary, unknown = "INFEASIBLE", True, False, False
    elif status == cp_model.UNKNOWN:
        status_class, rejected, necessary, unknown = "UNKNOWN", False, False, True
    else:
        raise AssertionError(f"unsupported CP-SAT status {status_name}")
    result = {
        "solver_status": status_name,
        "status_class": status_class,
        "rigorous_full_torsion_rejection": rejected,
        "necessary_only_survivor": necessary,
        "unknown_has_no_force": unknown,
        "timeout_seconds": timeout,
        "workers": workers,
        "seed": seed,
        "solver_wall_time_seconds": solver.WallTime(),
        "conflicts": solver.NumConflicts(),
        "branches": solver.NumBranches(),
        "model_scale": scale,
        "witness": None,
    }
    if necessary:
        selected = np.asarray([[solver.Value(v) for v in row] for row in lifts], dtype=np.int64)
        catalog_rows = np.asarray(
            [case["floors"][d] + 2 * selected[d] for d in range(DIRECTIONS)], dtype=np.int64
        )
        deltas, projections, membership_indices = [], [], []
        for direction, grade in enumerate(case["grades"]):
            lift = selected[direction]
            require(not np.any(kernel @ lift), "witness violates K")
            require(int(lift.sum()) == 5 * grade, "witness violates mass")
            require(np.all((0 <= lift) & (lift <= grade)), "witness violates bounds")
            hits = np.flatnonzero(
                np.all(case["catalogs"][direction] == catalog_rows[direction][None, :], axis=1)
            )
            require(len(hits) == 1, "SAT lift is absent from complete direct catalog")
            membership_indices.append(int(hits[0]))
            record = case["audit"]["direction_records"][direction]
            anchor = context["rebuilt"]["anchors"].get(record["mask"], record["scaled_mean"])
            delta = np.ascontiguousarray(anchor - catalog_rows[direction], dtype=np.int64)
            direct = cpu_join.semigroup.project_rows(
                delta[None, :], direction, context["quotient_data"], MOD7_COORDINATES
            )[0]
            compact = np.empty(WIDTH, dtype=np.uint8)
            compact[:6] = (case["offsets"][direction][:6] + maps[direction][:6] @ lift) % 3
            compact[6:] = (case["offsets"][direction][6:] + maps[direction][6:] @ lift) % 7
            require(np.array_equal(direct, compact), "witness projection audit failed")
            deltas.append(delta)
            projections.append(direct)
        projected = np.asarray(projections, dtype=np.int16)
        residue = np.empty(WIDTH, dtype=np.uint8)
        residue[:6] = projected[:, :6].sum(axis=0) % 3
        residue[6:] = projected[:, 6:].sum(axis=0) % 7
        require(np.array_equal(residue, case["target"]), "witness misses target")
        witness = {
            "L_rows": selected.tolist(),
            "L_rows_sha256_int64": array_sha256(selected.astype("<i8", copy=False)),
            "catalog_rows": catalog_rows.tolist(),
            "catalog_rows_sha256_int64": array_sha256(catalog_rows.astype("<i8", copy=False)),
            "anchor_relative_deltas": np.asarray(deltas).tolist(),
            "catalog_membership_indices": membership_indices,
            "target_residue": residue.tolist(),
            "quotient_values": [solver.Value(v) for v in quotient_variables],
            "all_K_mass_bounds_catalog_membership_and_27_modular_checks_passed": True,
            "semantics": "necessary full-torsion semigroup membership only",
        }
        result["witness"] = witness
        result["witness_sha256"] = json_sha256(witness)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-input", type=Path, default=cpu_join.DEFAULT_PARENT_INPUT)
    parser.add_argument("--current-join", type=Path, default=cpu_join.DEFAULT_CURRENT_JOIN)
    parser.add_argument("--hilbert-basis", type=Path, default=cpu_join.DEFAULT_HILBERT_BASIS)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--timeout", type=positive_float, default=300.0)
    parser.add_argument("--workers", type=positive_int, default=1)
    parser.add_argument("--seed", type=nonnegative_int, default=1)
    args = parser.parse_args()

    from ortools.sat.python import cp_model

    started = time.time()
    context, context_audit = gpu_join.build_cpu_context(
        args.parent_input, args.current_join, args.hilbert_basis
    )
    maps, map_audit = derive_maps(context)
    case = build_case(context, maps)
    decision = solve(cp_model, context, maps, case, args.timeout, args.workers, args.seed)
    script_path = Path(__file__).resolve()
    result = {
        "experiment": EXPERIMENT,
        "status": "complete_fixed_case_full_quotient_CSP",
        "fixed_case": case["audit"],
        "decision": decision,
        "map_audit": map_audit,
        "context_audit": context_audit,
        "provenance": {
            "script": {"path": str(script_path), "sha256": file_sha256(script_path)},
            "parent_input": {"path": str(args.parent_input.resolve()),
                             "sha256": file_sha256(args.parent_input)},
            "current_join": {"path": str(args.current_join.resolve()),
                             "sha256": file_sha256(args.current_join)},
            "hilbert_basis": {"path": str(args.hilbert_basis.resolve()),
                              "sha256": file_sha256(args.hilbert_basis)},
            "kernel_sha256_int64": array_sha256(
                np.asarray(context["rebuilt"]["kernel_rows"], dtype="<i8")
            ),
            "F3_complement_sha256_int64": array_sha256(
                np.asarray(context["quotient_data"][3]["complement"], dtype="<i8")
            ),
            "F7_complement_sha256_int64": array_sha256(
                np.asarray(context["quotient_data"][7]["complement"], dtype="<i8")
            ),
        },
        "semantics": {
            "INFEASIBLE_is_rigorous_full_torsion_semigroup_rejection": True,
            "FEASIBLE_is_necessary_only_after_witness_audit": True,
            "UNKNOWN_has_no_mathematical_force": True,
            "binary_edge_feasibility_claimed": False,
        },
        "elapsed_seconds": time.time() - started,
        "output_path": str(args.output.resolve()),
    }
    atomic_write(args.output, result)
    print(json.dumps({"status": result["status"], "solver": decision["solver_status"],
                      "output": str(args.output), "elapsed_seconds": result["elapsed_seconds"]},
                     indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
