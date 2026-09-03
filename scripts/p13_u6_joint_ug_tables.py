#!/usr/bin/env python3
"""Exact U=hM2 plus G=hM4-M2^2 tables for p=13,t=4,u=6.

The tables use only the pinned 74 translated-cut catalog and exact finite-field
interpolation.  Each optional CP-SAT replay has a fixed hard sign, exact-XNOR
root set, and excess assignment.  Its tables are the complete six-bin energy
spectra keyed by the common quadratic value U and quartic value G.  No graph,
coefficient-cell, or support census is performed.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
from ortools import __version__ as ortools_version
from ortools.sat.python import cp_model


ROOT = Path(__file__).resolve().parents[1]
CASE_INDEX_TEXT = os.environ.get("QML_U6_CASE_INDEX")
CASE_INDEX = int(CASE_INDEX_TEXT) if CASE_INDEX_TEXT is not None else None
OUTPUT = Path(
    "/tmp/qml_u6_joint_ug_dp.json"
    if CASE_INDEX is None
    else f"/tmp/qml_u6_joint_ug_dp_case{CASE_INDEX}.json"
)
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15740 import translated_cut_vector_catalog  # noqa: E402
from p13_p5_literal_interpolation import (  # noqa: E402
    POINTS,
    direction_signs,
    form_value,
    root_product_values,
)


P = 13
DISTANCE_VALUES = tuple(range(1, 7))
DISTANCES = np.asarray(DISTANCE_VALUES, dtype=np.int16)
QR_NONZERO = tuple(sorted({value * value % P for value in range(1, P)}))
CUTS = np.asarray(translated_cut_vector_catalog()["vectors"], dtype=np.int16)
W2 = np.asarray([pow(value, 2, P) for value in DISTANCE_VALUES], dtype=np.int16)
W4 = np.asarray([pow(value, 4, P) for value in DISTANCE_VALUES], dtype=np.int16)
if tuple(int(value) for value in W2) != (1, 4, 9, 3, 12, 10):
    raise ArithmeticError("the p=13 quadratic moment weights changed")
if tuple(int(value) for value in W4) != (1, 3, 3, 9, 1, 9):
    raise ArithmeticError("the p=13 quartic moment weights changed")


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def digest(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def row_catalog(
    *,
    total: int,
    low: int,
    high: int,
    l1_bound: int,
    cut_upper: int,
    extra_collision_budget: int | None = None,
) -> np.ndarray:
    """Enumerate the exact six-bin necessary row catalog."""
    base = (
        np.indices((high - low + 1,) * 5, dtype=np.int16)
        .reshape(5, -1)
        .T
        + low
    )
    last = total - base.sum(axis=1)
    valid = (last >= low) & (last <= high)
    rows = np.column_stack((base[valid], last[valid]))
    rows = rows[np.abs(rows).sum(axis=1) <= l1_bound]
    rows = rows[np.max(rows @ CUTS.T, axis=1) <= cut_upper]
    if extra_collision_budget is not None:
        overflow = np.maximum(rows - 6, 0) + np.maximum(-7 - rows, 0)
        rows = rows[overflow.sum(axis=1) <= extra_collision_budget]
    return rows


def joint_spectrum(rows: np.ndarray, *, opposite: bool) -> dict[tuple[int, int], tuple[int, ...]]:
    """Return energies keyed by local (N2,G)."""
    n2 = rows @ W2 % P
    n4 = rows @ W4 % P
    g = (-n4 - n2 * n2) % P if opposite else (n4 - n2 * n2) % P
    energy = (rows.astype(np.int32) ** 2).sum(axis=1)
    result: dict[tuple[int, int], tuple[int, ...]] = {}
    for left in range(P):
        for right in range(P):
            values = tuple(
                sorted(
                    {
                        int(value)
                        for value in energy[(n2 == left) & (g == right)]
                    }
                )
            )
            if values:
                result[(left, right)] = values
    return result


def table_for_global_u(
    spectrum: dict[tuple[int, int], tuple[int, ...]], *, opposite: bool
) -> tuple[tuple[int, int, int], ...]:
    rows = []
    for (local_n2, g), energies in sorted(spectrum.items()):
        global_u = (-local_n2) % P if opposite else local_n2
        rows.extend((global_u, g, energy) for energy in energies)
    return tuple(sorted(set(rows)))


H1_ROWS = row_catalog(total=0, low=-3, high=3, l1_bound=56, cut_upper=13)
H2_ROWS = row_catalog(total=-1, low=-5, high=5, l1_bound=55, cut_upper=13)
O_ROWS = row_catalog(total=-9, low=-4, high=1, l1_bound=57, cut_upper=-52)
H3_ROWS = {
    extra: row_catalog(
        total=-2,
        low=-7 - extra,
        high=6 + extra,
        l1_bound=54,
        cut_upper=13,
        extra_collision_budget=extra,
    )
    for extra in range(4)
}

SPECTRA = {
    "H1": joint_spectrum(H1_ROWS, opposite=False),
    "H2": joint_spectrum(H2_ROWS, opposite=False),
    "O": joint_spectrum(O_ROWS, opposite=True),
    **{
        f"H3_D{extra}": joint_spectrum(rows, opposite=False)
        for extra, rows in H3_ROWS.items()
    },
}
TABLES = {
    name: table_for_global_u(spectrum, opposite=(name == "O"))
    for name, spectrum in SPECTRA.items()
}


def evaluation_basis(point: tuple[int, int], degree: int) -> tuple[int, ...]:
    r, s = point
    return tuple(
        pow(r, degree - index, P) * pow(s, index, P) % P
        for index in range(degree + 1)
    )


U_BASIS = tuple(evaluation_basis(point, 2) for point in POINTS)


def add_residue_evaluation(
    model: cp_model.CpModel,
    coefficients: list[cp_model.IntVar],
    basis: tuple[int, ...],
    name: str,
) -> cp_model.IntVar:
    residue = model.NewIntVar(0, P - 1, name)
    quotient = model.NewIntVar(0, 10_000, f"{name}_quotient")
    model.Add(
        sum(value * coefficient for value, coefficient in zip(basis, coefficients))
        == P * quotient + residue
    )
    return residue


def build_model(
    *,
    hard_sign: int,
    roots: tuple[int, ...],
    assignment: tuple[int, ...],
    target: int,
    extra_collision_budget: int,
) -> tuple[
    cp_model.CpModel,
    list[cp_model.IntVar],
    list[cp_model.IntVar],
    list[cp_model.IntVar],
    list[cp_model.IntVar],
    tuple[int, ...],
    tuple[int, ...],
]:
    signs = direction_signs()
    hard = tuple(index for index, sign in enumerate(signs) if sign == hard_sign)
    opposite = tuple(index for index, sign in enumerate(signs) if sign == -hard_sign)
    remaining = tuple(index for index in hard if index not in roots)
    if len(remaining) != len(assignment):
        raise ArithmeticError("root/assignment dimensions disagree")

    root_values = root_product_values(roots)
    quotient_degree = 4 - len(roots)
    model = cp_model.CpModel()
    u_coefficients = [
        model.NewIntVar(0, P - 1, f"U_coefficient_{index}")
        for index in range(3)
    ]
    g_coefficients = [
        model.NewIntVar(0, P - 1, f"G_quotient_coefficient_{index}")
        for index in range(quotient_degree + 1)
    ]

    u_values = [
        add_residue_evaluation(model, u_coefficients, U_BASIS[index], f"U_{index}")
        for index in range(len(POINTS))
    ]
    g_values = []
    for index, point in enumerate(POINTS):
        quotient_basis = evaluation_basis(point, quotient_degree)
        basis = tuple(
            root_values[index] * value % P for value in quotient_basis
        )
        g_values.append(
            add_residue_evaluation(model, g_coefficients, basis, f"G_{index}")
        )

    for root in roots:
        model.AddAllowedAssignments(
            [u_values[root]], [(value,) for value in QR_NONZERO]
        )
        model.Add(g_values[root] == 0)

    energies: list[cp_model.IntVar] = []
    row_metadata: list[tuple[str, int]] = []
    for direction, excess in zip(remaining, assignment):
        table_name = (
            f"H3_D{extra_collision_budget}" if excess == 3 else f"H{excess}"
        )
        table = TABLES[table_name]
        energy = model.NewIntVar(
            min(row[2] for row in table),
            max(row[2] for row in table),
            f"energy_hard_{direction}",
        )
        model.AddAllowedAssignments(
            [u_values[direction], g_values[direction], energy], table
        )
        energies.append(energy)
        row_metadata.append((table_name, direction))
    for direction in opposite:
        table = TABLES["O"]
        energy = model.NewIntVar(
            min(row[2] for row in table),
            max(row[2] for row in table),
            f"energy_opposite_{direction}",
        )
        model.AddAllowedAssignments(
            [u_values[direction], g_values[direction], energy], table
        )
        energies.append(energy)
        row_metadata.append(("O", direction))
    model.Add(sum(energies) == target)
    return (
        model,
        u_coefficients,
        g_coefficients,
        u_values,
        g_values,
        remaining,
        tuple(direction for _name, direction in row_metadata),
    )


CASES = (
    {
        "partition": (1, 1, 1, 1, 1),
        "collision_minimum": 0,
        "targets": ((0, 303), (1, 329)),
    },
    {
        "partition": (2, 1, 1, 1),
        "collision_minimum": 0,
        "targets": ((0, 298), (1, 324)),
    },
    {
        "partition": (2, 2, 1),
        "collision_minimum": 0,
        "targets": ((0, 293), (1, 319), (2, 345)),
    },
    {
        "partition": (3, 1, 1),
        "collision_minimum": 1,
        "targets": ((1, 315), (2, 341), (3, 367), (4, 393)),
    },
)


def solve_target(
    partition: tuple[int, ...],
    collision_minimum: int,
    collision: int,
    target: int,
) -> dict[str, object]:
    exact_count = 7 - len(partition)
    extra_collision_budget = collision - collision_minimum
    signs = direction_signs()
    model_hashes: list[str] = []
    status_counts: dict[str, int] = {}
    attempts = 0
    expected_attempts = 0
    for hard_sign in (-1, 1):
        hard = tuple(
            index for index, sign in enumerate(signs) if sign == hard_sign
        )
        expected_attempts += (
            len(tuple(itertools.combinations(hard, exact_count)))
            * len(set(itertools.permutations(partition)))
        )

    for hard_sign in (-1, 1):
        hard = tuple(
            index for index, sign in enumerate(signs) if sign == hard_sign
        )
        for roots in itertools.combinations(hard, exact_count):
            for assignment in sorted(set(itertools.permutations(partition))):
                attempts += 1
                (
                    model,
                    u_coefficients,
                    g_coefficients,
                    u_values,
                    g_values,
                    remaining,
                    energy_directions,
                ) = build_model(
                    hard_sign=hard_sign,
                    roots=roots,
                    assignment=assignment,
                    target=target,
                    extra_collision_budget=extra_collision_budget,
                )
                model_hashes.append(
                    hashlib.sha256(str(model.Proto()).encode()).hexdigest()
                )
                solver = cp_model.CpSolver()
                solver.parameters.num_search_workers = 1
                solver.parameters.random_seed = 0
                solver.parameters.cp_model_presolve = True
                status = solver.Solve(model)
                status_name = solver.StatusName(status)
                status_counts[status_name] = status_counts.get(status_name, 0) + 1
                if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
                    energy_values = []
                    # Recover named energy variables from the response by rebuilding
                    # their indices from the model proto names.
                    proto = model.Proto()
                    for index, variable in enumerate(proto.variables):
                        if variable.name.startswith("energy_"):
                            energy_values.append(
                                [variable.name, int(solver.ResponseProto().solution[index])]
                            )
                    witness = {
                        "hard_sign": hard_sign,
                        "roots": list(roots),
                        "remaining_hard_directions": list(remaining),
                        "excess_assignment": list(assignment),
                        "U_coefficients": [
                            int(solver.Value(value)) for value in u_coefficients
                        ],
                        "G_quotient_coefficients": [
                            int(solver.Value(value)) for value in g_coefficients
                        ],
                        "U_evaluations": [
                            int(solver.Value(value)) for value in u_values
                        ],
                        "G_evaluations": [
                            int(solver.Value(value)) for value in g_values
                        ],
                        "energy_directions": list(energy_directions),
                        "row_energies": energy_values,
                        "energy_sum": sum(value for _name, value in energy_values),
                    }
                    return {
                        "collision": collision,
                        "target_nonexact_energy": target,
                        "joint_UG_compatible": True,
                        "attempts_before_witness": attempts,
                        "expected_total_fixed_models": expected_attempts,
                        "fixed_model_sha256_aggregate": digest(model_hashes),
                        "solver_status_counts": status_counts,
                        "witness": witness,
                    }
                if status != cp_model.INFEASIBLE:
                    raise ArithmeticError(
                        f"unexpected solver status {status_name} for {partition}, C={collision}"
                    )
    if attempts != expected_attempts:
        raise ArithmeticError("fixed-model coverage count changed")
    return {
        "collision": collision,
        "target_nonexact_energy": target,
        "joint_UG_compatible": False,
        "all_fixed_models_infeasible": True,
        "fixed_models_checked": attempts,
        "expected_total_fixed_models": expected_attempts,
        "fixed_model_sha256_aggregate": digest(model_hashes),
        "solver_status_counts": status_counts,
        "witness": None,
    }


def main() -> None:
    case_rows = []
    selected_cases = CASES if CASE_INDEX is None else (CASES[CASE_INDEX],)
    for case in selected_cases:
        partition = tuple(case["partition"])
        collision_minimum = int(case["collision_minimum"])
        results = []
        for collision, target in case["targets"]:
            print(
                f"start partition={partition} C={collision} target={target}",
                flush=True,
            )
            solved = solve_target(
                partition,
                collision_minimum,
                int(collision),
                int(target),
            )
            results.append(solved)
            print(
                f"done partition={partition} C={collision} "
                f"compatible={solved['joint_UG_compatible']}",
                flush=True,
            )
        case_rows.append(
            {
                "partition": list(partition),
                "exact_XNOR_root_count": 7 - len(partition),
                "collision_minimum": collision_minimum,
                "collision_rows": results,
            }
        )

    row_inputs = {
        "H1": H1_ROWS.tolist(),
        "H2": H2_ROWS.tolist(),
        "O": O_ROWS.tolist(),
        **{f"H3_D{extra}": rows.tolist() for extra, rows in H3_ROWS.items()},
    }
    result = {
        "title": "p13 t4 u6 joint common-quadratic and quartic energy DP",
        "result_status": "exact exploratory necessary-moment reduction",
        "repository_HEAD": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "script_path": str(Path(__file__).resolve()),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "output_path": str(OUTPUT),
        "solver": {
            "name": "OR-Tools CP-SAT",
            "version": ortools_version,
            "workers": 1,
            "random_seed": 0,
            "presolve": True,
        },
        "field": "F_13",
        "projective_points": [list(point) for point in POINTS],
        "direction_signs": list(direction_signs()),
        "hard_sign_classes_checked": [-1, 1],
        "nonzero_quadratic_residues_required_at_exact_XNOR_roots": list(
            QR_NONZERO
        ),
        "common_forms": {
            "U": "U(L)=h*M2(L), homogeneous binary quadratic",
            "G": "G(L)=h*M4(L)-M2(L)^2, homogeneous binary quartic",
            "hard_local_key": "(N2,G)=(U,N4-N2^2)",
            "opposite_local_key": "(N2,G)=(-U,-N4-N2^2)",
        },
        "translated_cut_catalog": {
            "count": len(CUTS),
            "sha256": digest(CUTS.tolist()),
        },
        "row_catalogs": {
            name: {
                "row_count": len(rows),
                "row_sha256": digest(rows),
                "joint_key_count": len(SPECTRA[name]),
                "joint_spectrum_sha256": digest(
                    [
                        [left, right, list(values)]
                        for (left, right), values in sorted(SPECTRA[name].items())
                    ]
                ),
            }
            for name, rows in row_inputs.items()
        },
        "partition_results": case_rows,
        "closed_collision_strata": [
            {
                "partition": row["partition"],
                "collision": result["collision"],
            }
            for row in case_rows
            for result in row["collision_rows"]
            if not result["joint_UG_compatible"]
        ],
        "remaining_collision_strata": [
            {
                "partition": row["partition"],
                "collision": result["collision"],
            }
            for row in case_rows
            for result in row["collision_rows"]
            if result["joint_UG_compatible"]
        ],
        "scope": (
            "Necessary common U/G moment and energy compatibility only; a survivor "
            "is not a graph realization, while an all-INFEASIBLE stratum is excluded "
            "by the stated necessary constraints."
        ),
    }
    temporary = OUTPUT.with_suffix(".json.tmp")
    temporary.write_bytes(json.dumps(result, indent=2, sort_keys=True).encode() + b"\n")
    os.replace(temporary, OUTPUT)
    print(json.dumps({
        "output": str(OUTPUT),
        "closed_collision_strata": result["closed_collision_strata"],
        "remaining_collision_strata": result["remaining_collision_strata"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
