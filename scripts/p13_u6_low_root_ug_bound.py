#!/usr/bin/env python3
"""Exact separable-energy upper bound for the low p=13,t=4,u=6 U/G cases.

This verifier consumes the exact row tables and finite-field conventions from
the pinned repository U/G table generator, but replaces slow fixed-table CP-SAT
proofs with exhaustive
vectorized enumeration of every coefficient tuple in every fixed model.

For fixed common forms U and G, each nonexact direction has an exact finite
set of permitted row energies.  Therefore the sum of the largest permitted
local energies is an upper bound on every joint energy realization.  If that
upper bound is below Parseval's required nonexact energy, that coefficient
pair is impossible without any convolution or solver search.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import runpy
import subprocess
from pathlib import Path
from typing import Any

import numpy as np


P = 13
ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts" / "p13_u6_joint_ug_tables.py"
SOURCE_SHA256 = "811eb1833d5551f3fadce62ec6a4302296e301a2cd001bedd399df723afaaebe"
HELPER_SHA256 = {
    "src/e1_gmin_m4_prop15740.py":
        "7a2cfbd12a7057971a0cbaaf523d16f65cff989a47d594a79f641062d02439c3",
    "scripts/p13_p5_literal_interpolation.py":
        "31ba186632780d62b3352c214b1112a6beafa48b8a2571882da29e94faf618e5",
}
DEFAULT_OUTPUT = ROOT / "evidence" / "e1_gmin_m4_prop15754_low_root_ug.json"


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def display_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def digest(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def coefficient_grid(degree: int) -> np.ndarray:
    return np.asarray(
        list(itertools.product(range(P), repeat=degree + 1)), dtype=np.int16
    )


def main(output: Path = DEFAULT_OUTPUT) -> None:
    actual_source_sha256 = file_sha256(SOURCE)
    if actual_source_sha256 != SOURCE_SHA256:
        raise ArithmeticError(
            f"input checker changed: {actual_source_sha256} != {SOURCE_SHA256}"
        )
    actual_helper_sha256 = {}
    for relative_path, expected_sha256 in HELPER_SHA256.items():
        actual_sha256 = file_sha256(ROOT / relative_path)
        if actual_sha256 != expected_sha256:
            raise ArithmeticError(
                f"imported helper changed: {relative_path}: "
                f"{actual_sha256} != {expected_sha256}"
            )
        actual_helper_sha256[relative_path] = actual_sha256
    source = runpy.run_path(str(SOURCE), run_name="_qml_u6_joint_ug_input")

    points = tuple(tuple(int(value) for value in point) for point in source["POINTS"])
    signs = tuple(int(value) for value in source["direction_signs"]())
    qr_nonzero = tuple(int(value) for value in source["QR_NONZERO"])
    if len(points) != 14 or signs.count(-1) != 7 or signs.count(1) != 7:
        raise ArithmeticError("projective direction split changed")
    distances = tuple(int(value) for value in source["DISTANCES"])
    safe_w2 = tuple(pow(distance, 2, P) for distance in distances)
    safe_w4 = tuple(pow(distance, 4, P) for distance in distances)
    source_w2 = tuple(int(value) for value in source["W2"])
    source_w4 = tuple(int(value) for value in source["W4"])
    if distances != (1, 2, 3, 4, 5, 6):
        raise ArithmeticError("distance representatives changed")
    if safe_w2 != (1, 4, 9, 3, 12, 10) or source_w2 != safe_w2:
        raise ArithmeticError(f"W2 regression changed: safe={safe_w2}, source={source_w2}")
    if safe_w4 != (1, 3, 3, 9, 1, 9) or source_w4 != safe_w4:
        raise ArithmeticError(f"W4 regression changed: safe={safe_w4}, source={source_w4}")

    tables = source["TABLES"]
    table_maxima: dict[str, np.ndarray] = {}
    table_hashes: dict[str, str] = {}
    for name in ("H1", "H2", "O"):
        maximum = np.full((P, P), -1, dtype=np.int16)
        for u_value, g_value, energy in tables[name]:
            maximum[int(u_value), int(g_value)] = max(
                int(maximum[int(u_value), int(g_value)]), int(energy)
            )
        expected_keys = len({(int(row[0]), int(row[1])) for row in tables[name]})
        if int(np.count_nonzero(maximum >= 0)) != expected_keys:
            raise ArithmeticError(f"{name} maximum-table key loss")
        table_maxima[name] = maximum
        table_hashes[name] = digest(
            [[int(left), int(right), int(maximum[left, right])]
             for left in range(P) for right in range(P)
             if maximum[left, right] >= 0]
        )

    u_basis = np.asarray(source["U_BASIS"], dtype=np.int16)
    u_coefficients = coefficient_grid(2)
    u_evaluations = u_coefficients @ u_basis.T % P
    is_qr_nonzero = np.zeros(P, dtype=np.bool_)
    is_qr_nonzero[list(qr_nonzero)] = True

    def fixed_model_upper(
        *,
        hard_sign: int,
        roots: tuple[int, ...],
        assignment: tuple[int, ...],
    ) -> dict[str, Any]:
        hard = tuple(index for index, sign in enumerate(signs) if sign == hard_sign)
        opposite = tuple(index for index, sign in enumerate(signs) if sign == -hard_sign)
        remaining = tuple(index for index in hard if index not in roots)
        if len(remaining) != len(assignment):
            raise ArithmeticError("assignment length changed")

        eligible_mask = is_qr_nonzero[u_evaluations[:, roots]].all(axis=1)
        eligible_indices = np.flatnonzero(eligible_mask)
        quotient_degree = 4 - len(roots)
        g_coefficients = coefficient_grid(quotient_degree)
        quotient_basis = np.asarray(
            [source["evaluation_basis"](point, quotient_degree) for point in points],
            dtype=np.int16,
        )
        root_values = np.asarray(source["root_product_values"](roots), dtype=np.int16)
        g_evaluations = (g_coefficients @ quotient_basis.T % P) * root_values % P
        if np.any(g_evaluations[:, roots] != 0):
            raise ArithmeticError("root-product parameterization lost a root")

        row_kinds = tuple(
            [(direction, f"H{excess}")
             for direction, excess in zip(remaining, assignment)]
            + [(direction, "O") for direction in opposite]
        )
        coefficient_pairs_checked = 0
        locally_compatible_pairs = 0
        best_upper = -1
        best_u_index: int | None = None
        best_g_index: int | None = None

        for u_index in eligible_indices:
            u_values = u_evaluations[int(u_index)]
            upper = np.zeros(len(g_coefficients), dtype=np.int16)
            compatible = np.ones(len(g_coefficients), dtype=np.bool_)
            for direction, kind in row_kinds:
                local = table_maxima[kind][u_values[direction], g_evaluations[:, direction]]
                present = local >= 0
                compatible &= present
                upper += np.where(present, local, 0)
            coefficient_pairs_checked += len(g_coefficients)
            locally_compatible_pairs += int(np.count_nonzero(compatible))
            if np.any(compatible):
                compatible_indices = np.flatnonzero(compatible)
                local_argmax = int(np.argmax(upper[compatible]))
                g_index = int(compatible_indices[local_argmax])
                value = int(upper[g_index])
                if value > best_upper:
                    best_upper = value
                    best_u_index = int(u_index)
                    best_g_index = g_index

        if best_u_index is None or best_g_index is None:
            witness = None
        else:
            u_values = u_evaluations[best_u_index]
            g_values = g_evaluations[best_g_index]
            local_rows = []
            scalar_sum = 0
            for direction, kind in row_kinds:
                maximum = int(table_maxima[kind][u_values[direction], g_values[direction]])
                if maximum < 0:
                    raise ArithmeticError("maximizing witness has a missing local key")
                scalar_sum += maximum
                local_rows.append({
                    "direction": int(direction),
                    "kind": kind,
                    "U": int(u_values[direction]),
                    "G": int(g_values[direction]),
                    "maximum_energy": maximum,
                })
            if scalar_sum != best_upper:
                raise ArithmeticError("scalar witness sum disagrees with vectorized maximum")
            witness = {
                "U_coefficients": [int(value) for value in u_coefficients[best_u_index]],
                "G_quotient_coefficients": [
                    int(value) for value in g_coefficients[best_g_index]
                ],
                "U_evaluations": [int(value) for value in u_values],
                "G_evaluations": [int(value) for value in g_values],
                "local_maximum_rows": local_rows,
                "local_maximum_energy_sum": scalar_sum,
            }

        return {
            "hard_sign": hard_sign,
            "roots": list(roots),
            "remaining_hard_directions": list(remaining),
            "excess_assignment": list(assignment),
            "eligible_U_coefficient_count": int(len(eligible_indices)),
            "G_quotient_degree": quotient_degree,
            "G_quotient_coefficient_count": int(len(g_coefficients)),
            "coefficient_pairs_checked": coefficient_pairs_checked,
            "locally_compatible_coefficient_pairs": locally_compatible_pairs,
            "maximum_separable_nonexact_energy": best_upper,
            "maximizing_witness": witness,
        }

    cases = (
        {
            "partition": (1, 1, 1, 1, 1),
            "targets": ((0, 303), (1, 329), (2, 355)),
        },
        {
            "partition": (2, 1, 1, 1),
            "targets": ((0, 298), (1, 324), (2, 350)),
        },
    )
    case_results = []
    for case in cases:
        partition = tuple(int(value) for value in case["partition"])
        targets = tuple(
            (int(collision), int(target)) for collision, target in case["targets"]
        )
        if tuple(collision for collision, _target in targets) != (0, 1, 2):
            raise ArithmeticError("low-partition collision target range changed")
        if any(target != targets[0][1] + 26 * collision for collision, target in targets):
            raise ArithmeticError("Parseval collision increment changed")
        exact_count = 7 - len(partition)
        assignments = tuple(sorted(set(itertools.permutations(partition))))
        fixed_models = []
        for hard_sign in (-1, 1):
            hard = tuple(index for index, sign in enumerate(signs) if sign == hard_sign)
            for roots in itertools.combinations(hard, exact_count):
                for assignment in assignments:
                    fixed_models.append(fixed_model_upper(
                        hard_sign=hard_sign,
                        roots=tuple(roots),
                        assignment=tuple(int(value) for value in assignment),
                    ))

        expected_fixed_models = (
            2 * len(tuple(itertools.combinations(range(7), exact_count)))
            * len(assignments)
        )
        if len(fixed_models) != expected_fixed_models:
            raise ArithmeticError("fixed-model coverage count changed")
        compatible_models = [
            row for row in fixed_models
            if int(row["maximum_separable_nonexact_energy"]) >= 0
        ]
        if not compatible_models:
            raise ArithmeticError("no locally compatible fixed model")
        global_maximum = max(
            int(row["maximum_separable_nonexact_energy"])
            for row in compatible_models
        )
        first_maximizer = next(
            row for row in compatible_models
            if int(row["maximum_separable_nonexact_energy"]) == global_maximum
        )
        # Pin the original CP table model to the enumerated maximizing forms.
        # This is a separate scalar/table check that the reported upper bound is
        # attained in the same necessary model, not an artifact of array lookup.
        maximizing_witness = first_maximizer["maximizing_witness"]
        if maximizing_witness is None:
            raise ArithmeticError("global maximizer lost its coefficient witness")
        (
            attainment_model,
            attainment_u_coefficients,
            attainment_g_coefficients,
            _attainment_u_values,
            _attainment_g_values,
            _attainment_remaining,
            _attainment_energy_directions,
        ) = source["build_model"](
            hard_sign=int(first_maximizer["hard_sign"]),
            roots=tuple(int(value) for value in first_maximizer["roots"]),
            assignment=tuple(
                int(value) for value in first_maximizer["excess_assignment"]
            ),
            target=global_maximum,
            extra_collision_budget=0,
        )
        for variable, value in zip(
            attainment_u_coefficients, maximizing_witness["U_coefficients"]
        ):
            attainment_model.Add(variable == int(value))
        for variable, value in zip(
            attainment_g_coefficients,
            maximizing_witness["G_quotient_coefficients"],
        ):
            attainment_model.Add(variable == int(value))
        attainment_solver = source["cp_model"].CpSolver()
        attainment_solver.parameters.num_search_workers = 1
        attainment_solver.parameters.random_seed = 0
        attainment_status = attainment_solver.Solve(attainment_model)
        if attainment_status not in (
            source["cp_model"].OPTIMAL,
            source["cp_model"].FEASIBLE,
        ):
            raise ArithmeticError("original CP table model rejected the maximizer")
        attainment_energies = [
            int(attainment_solver.ResponseProto().solution[index])
            for index, variable in enumerate(attainment_model.Proto().variables)
            if variable.name.startswith("energy_")
        ]
        if sum(attainment_energies) != global_maximum:
            raise ArithmeticError("original CP table model changed the maximum sum")
        attainment_validation = {
            "solver_status": attainment_solver.StatusName(attainment_status),
            "pinned_U_coefficients": [
                int(attainment_solver.Value(value))
                for value in attainment_u_coefficients
            ],
            "pinned_G_quotient_coefficients": [
                int(attainment_solver.Value(value))
                for value in attainment_g_coefficients
            ],
            "row_energies": attainment_energies,
            "energy_sum": sum(attainment_energies),
        }
        target_rows = []
        for collision, target in targets:
            if global_maximum >= int(target):
                raise ArithmeticError(
                    f"separable maximum does not close {partition}, C={collision}"
                )
            target_rows.append({
                "collision": int(collision),
                "target_nonexact_energy": int(target),
                "maximum_separable_nonexact_energy": global_maximum,
                "strict_energy_deficit": int(target) - global_maximum,
                "all_fixed_models_excluded": True,
            })
        case_results.append({
            "partition": list(partition),
            "exact_XNOR_root_count": exact_count,
            "unique_excess_assignment_count": len(assignments),
            "fixed_models_checked": len(fixed_models),
            "expected_fixed_models": expected_fixed_models,
            "coefficient_pairs_checked": sum(
                int(row["coefficient_pairs_checked"]) for row in fixed_models
            ),
            "locally_compatible_coefficient_pairs": sum(
                int(row["locally_compatible_coefficient_pairs"])
                for row in fixed_models
            ),
            "global_maximum_separable_nonexact_energy": global_maximum,
            "first_global_maximizer": first_maximizer,
            "original_CP_table_model_attainment_validation": attainment_validation,
            "targets": target_rows,
            "fixed_model_summary_sha256": digest([
                {
                    key: row[key]
                    for key in (
                        "hard_sign", "roots", "excess_assignment",
                        "eligible_U_coefficient_count",
                        "G_quotient_coefficient_count",
                        "coefficient_pairs_checked",
                        "locally_compatible_coefficient_pairs",
                        "maximum_separable_nonexact_energy",
                    )
                }
                for row in fixed_models
            ]),
            "fixed_models": fixed_models,
        })

    all_targets_excluded = all(
        bool(target["all_fixed_models_excluded"])
        for case in case_results
        for target in case["targets"]
    )
    if not all_targets_excluded:
        raise ArithmeticError("a low-root U/G collision stratum survived")
    result = {
        "title": "p13 t4 u6 low-partition exact joint U/G separable-energy bound",
        "result_status": "exact finite necessary-moment certificate",
        "repository_HEAD": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "source_checker_path": display_path(SOURCE),
        "source_checker_sha256": actual_source_sha256,
        "imported_helper_files": [
            {
                "repository_relative_path": relative_path,
                "sha256": actual_helper_sha256[relative_path],
            }
            for relative_path in HELPER_SHA256
        ],
        "script_path": display_path(Path(__file__)),
        "script_sha256": file_sha256(Path(__file__).resolve()),
        "output_path": display_path(output),
        "numpy_version": np.__version__,
        "field": "F_13",
        "moment_weight_regression_mod_13": {
            "method": "Python builtin pow(distance, degree, 13) before array casting",
            "distances": list(distances),
            "W2": list(safe_w2),
            "W4": list(safe_w4),
            "source_checker_W2_matches": True,
            "source_checker_W4_matches": True,
        },
        "projective_points": [list(point) for point in points],
        "direction_signs": list(signs),
        "nonzero_quadratic_residues_at_exact_XNOR_roots": list(qr_nonzero),
        "table_maximum_sha256": table_hashes,
        "method": (
            "For every hard sign, exact-root subset, excess assignment, eligible "
            "binary-quadratic U coefficient tuple, and root-parameterized binary-"
            "quartic G quotient tuple, sum the exact per-key maxima of all nonexact "
            "row-energy tables. The resulting separable maximum bounds every "
            "joint energy convolution from above."
        ),
        "collision_note": (
            "These two partitions use only H1, H2, and O tables, whose catalogs "
            "are collision-independent here. The listed C=0,1,2 targets advance "
            "by the exact Parseval collision increment 26; exclusion of C>=3 by "
            "the separate raw collision upper rules is outside this artifact."
        ),
        "case_results": case_results,
        "closed_collision_strata": [
            {"partition": case["partition"], "collision": target["collision"]}
            for case in case_results for target in case["targets"]
        ],
        "all_targets_excluded": all_targets_excluded,
        "proved": all_targets_excluded,
        "scope": (
            "Necessary common U/G moment and row-energy compatibility only. "
            "A strict upper-bound deficit excludes the stratum; no graph, cell, "
            "support, or coefficient census is performed."
        ),
    }
    temporary = output.with_suffix(".json.tmp")
    temporary.write_bytes(json.dumps(result, indent=2, sort_keys=True).encode() + b"\n")
    os.replace(temporary, output)
    print(json.dumps({
        "output": str(output),
        "script_sha256": result["script_sha256"],
        "closed_collision_strata": result["closed_collision_strata"],
        "case_summaries": [
            {
                "partition": case["partition"],
                "fixed_models_checked": case["fixed_models_checked"],
                "coefficient_pairs_checked": case["coefficient_pairs_checked"],
                "locally_compatible_coefficient_pairs":
                    case["locally_compatible_coefficient_pairs"],
                "global_maximum": case["global_maximum_separable_nonexact_energy"],
                "targets": case["targets"],
            }
            for case in case_results
        ],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    main(output=arguments.output)
