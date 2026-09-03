#!/usr/bin/env python3
"""Exact U/G/J6 necessary-moment checker for the p=13,u=6 four-root strata.

This checks common binary form evaluations and exact row-energy tables.  It is
not a graph, cell, support, or orbit census.
"""
from __future__ import annotations

import argparse
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
DEFAULT_CASE_INDEX = int(os.environ.get("QML_U6_CASE_INDEX", "0"))
CUT_HELPER = ROOT / "src/e1_gmin_m4_prop15740.py"
INTERPOLATION_HELPER = ROOT / "scripts/p13_p5_literal_interpolation.py"
EXPECTED_HELPER_SHA256 = {
    str(CUT_HELPER.relative_to(ROOT)): "7a2cfbd12a7057971a0cbaaf523d16f65cff989a47d594a79f641062d02439c3",
    str(INTERPOLATION_HELPER.relative_to(ROOT)): "31ba186632780d62b3352c214b1112a6beafa48b8a2571882da29e94faf618e5",
}
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from e1_gmin_m4_prop15740 import translated_cut_vector_catalog  # noqa: E402
from p13_p5_literal_interpolation import (  # noqa: E402
    POINTS,
    direction_signs,
    root_product_values,
)


P = 13
QR_NONZERO = tuple(sorted({value * value % P for value in range(1, P)}))
CUTS = np.asarray(translated_cut_vector_catalog()["vectors"], dtype=np.int64)
W2 = np.asarray([pow(a, 2, P) for a in range(1, 7)], dtype=np.int64)
W4 = np.asarray([pow(a, 4, P) for a in range(1, 7)], dtype=np.int64)
W6 = np.asarray([pow(a, 6, P) for a in range(1, 7)], dtype=np.int64)
assert tuple(int(x) for x in W2) == (1, 4, 9, 3, 12, 10)
assert tuple(int(x) for x in W4) == (1, 3, 3, 9, 1, 9)
assert tuple(int(x) for x in W6) == (1, 12, 1, 1, 12, 12)


def canonical_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode()


def digest(value: object) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def display_path(path: Path) -> str:
    resolved = path.resolve()
    try:
        return str(resolved.relative_to(ROOT))
    except ValueError:
        return str(resolved)


ACTUAL_HELPER_SHA256 = {
    name: file_sha256(ROOT / name) for name in EXPECTED_HELPER_SHA256
}
if ACTUAL_HELPER_SHA256 != EXPECTED_HELPER_SHA256:
    raise ArithmeticError(
        f"imported helper changed: {ACTUAL_HELPER_SHA256} != {EXPECTED_HELPER_SHA256}"
    )


def row_catalog(
    *,
    total: int,
    low: int,
    high: int,
    l1_bound: int,
    cut_upper: int,
    extra_collision_budget: int | None = None,
) -> np.ndarray:
    base = (
        np.indices((high - low + 1,) * 5, dtype=np.int16)
        .reshape(5, -1)
        .T
        + low
    )
    last = total - base.sum(axis=1)
    valid = (last >= low) & (last <= high)
    rows = np.column_stack((base[valid], last[valid])).astype(np.int64)
    rows = rows[np.abs(rows).sum(axis=1) <= l1_bound]
    rows = rows[np.max(rows @ CUTS.T, axis=1) <= cut_upper]
    if extra_collision_budget is not None:
        overflow = np.maximum(rows - 6, 0) + np.maximum(-7 - rows, 0)
        rows = rows[overflow.sum(axis=1) <= extra_collision_budget]
    return rows


def joint_table(rows: np.ndarray, *, opposite: bool) -> tuple[tuple[int, ...], ...]:
    """Return global (U,G,J6,energy) tuples from local normalized rows."""
    n2 = rows @ W2 % P
    n4 = rows @ W4 % P
    n6 = rows @ W6 % P
    energy = (rows * rows).sum(axis=1)
    if opposite:
        # Local N_(2r)=(-h)M_(2r).  Thus U=-N2,
        # G=-(N4+N2^2), and J6=-(N6-N2^3).
        u = -n2 % P
        g = -n4 - n2 * n2
        j6 = -n6 + n2 * n2 * n2
    else:
        u = n2
        g = n4 - n2 * n2
        j6 = n6 - n2 * n2 * n2
    return tuple(
        sorted(
            {
                (int(a % P), int(b % P), int(c % P), int(e))
                for a, b, c, e in zip(u, g, j6, energy)
            }
        )
    )


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
ROW_CATALOGS = {
    "H1": H1_ROWS,
    "H2": H2_ROWS,
    "O": O_ROWS,
    **{f"H3_D{extra}": rows for extra, rows in H3_ROWS.items()},
}
TABLES = {
    name: joint_table(rows, opposite=(name == "O"))
    for name, rows in ROW_CATALOGS.items()
}


def evaluation_basis(point: tuple[int, int], degree: int) -> tuple[int, ...]:
    r, s = point
    return tuple(
        pow(r, degree - index, P) * pow(s, index, P) % P
        for index in range(degree + 1)
    )


U_BASIS = tuple(evaluation_basis(point, 2) for point in POINTS)


def add_evaluation(
    model: cp_model.CpModel,
    coefficients: list[cp_model.IntVar],
    basis: tuple[int, ...],
    name: str,
) -> cp_model.IntVar:
    residue = model.NewIntVar(0, P - 1, name)
    # All summands are nonnegative and there are at most three coefficients.
    quotient = model.NewIntVar(0, 100, f"{name}_quotient")
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
) -> tuple[cp_model.CpModel, dict[str, object]]:
    signs = direction_signs()
    hard = tuple(i for i, sign in enumerate(signs) if sign == hard_sign)
    opposite = tuple(i for i, sign in enumerate(signs) if sign == -hard_sign)
    remaining = tuple(i for i in hard if i not in roots)
    if len(roots) != 4 or len(remaining) != len(assignment):
        raise ArithmeticError("this checker is restricted to four exact roots")

    root_values = root_product_values(roots)
    model = cp_model.CpModel()
    u_coefficients = [model.NewIntVar(0, P - 1, f"u_c{i}") for i in range(3)]
    g_coefficient = [model.NewIntVar(0, P - 1, "g_c0")]
    j_coefficients = [model.NewIntVar(0, P - 1, f"j_c{i}") for i in range(3)]
    u_values = [
        add_evaluation(model, u_coefficients, U_BASIS[i], f"u_{i}")
        for i in range(len(POINTS))
    ]
    g_values = [
        add_evaluation(model, g_coefficient, (root_values[i],), f"g_{i}")
        for i in range(len(POINTS))
    ]
    j_values = []
    for i, point in enumerate(POINTS):
        basis = tuple(
            root_values[i] * value % P for value in evaluation_basis(point, 2)
        )
        j_values.append(add_evaluation(model, j_coefficients, basis, f"j_{i}"))

    for root in roots:
        model.AddAllowedAssignments([u_values[root]], [(q,) for q in QR_NONZERO])
        model.Add(g_values[root] == 0)
        model.Add(j_values[root] == 0)

    energies: list[cp_model.IntVar] = []
    energy_names: list[str] = []
    for direction, excess in zip(remaining, assignment):
        table_name = f"H3_D{extra_collision_budget}" if excess == 3 else f"H{excess}"
        table = TABLES[table_name]
        energy = model.NewIntVar(
            min(row[3] for row in table), max(row[3] for row in table),
            f"energy_hard_{direction}",
        )
        model.AddAllowedAssignments(
            [u_values[direction], g_values[direction], j_values[direction], energy],
            table,
        )
        energies.append(energy)
        energy_names.append(energy.Name())
    for direction in opposite:
        table = TABLES["O"]
        energy = model.NewIntVar(
            min(row[3] for row in table), max(row[3] for row in table),
            f"energy_opposite_{direction}",
        )
        model.AddAllowedAssignments(
            [u_values[direction], g_values[direction], j_values[direction], energy],
            table,
        )
        energies.append(energy)
        energy_names.append(energy.Name())
    model.Add(sum(energies) == target)
    return model, {
        "u_coefficients": u_coefficients,
        "g_coefficients": g_coefficient,
        "j_coefficients": j_coefficients,
        "u_values": u_values,
        "g_values": g_values,
        "j_values": j_values,
        "energies": energies,
        "energy_names": energy_names,
        "remaining": remaining,
    }


CASES = (
    {
        "partition": (2, 2, 1),
        "collision_minimum": 0,
        "targets": ((0, 293), (1, 319), (2, 345), (3, 371)),
    },
    {
        "partition": (3, 1, 1),
        "collision_minimum": 1,
        "targets": ((1, 315), (2, 341), (3, 367), (4, 393)),
    },
)


def solve_target(case: dict[str, object], collision: int, target: int) -> dict[str, object]:
    partition = tuple(int(x) for x in case["partition"])
    collision_minimum = int(case["collision_minimum"])
    extra_collision_budget = collision - collision_minimum
    signs = direction_signs()
    attempts = 0
    model_hashes: list[str] = []
    statuses: dict[str, int] = {}
    expected = 0
    for hard_sign in (-1, 1):
        hard = tuple(i for i, sign in enumerate(signs) if sign == hard_sign)
        expected += len(tuple(itertools.combinations(hard, 4))) * len(
            set(itertools.permutations(partition))
        )

    for hard_sign in (-1, 1):
        hard = tuple(i for i, sign in enumerate(signs) if sign == hard_sign)
        for roots in itertools.combinations(hard, 4):
            for assignment in sorted(set(itertools.permutations(partition))):
                attempts += 1
                model, data = build_model(
                    hard_sign=hard_sign,
                    roots=roots,
                    assignment=assignment,
                    target=target,
                    extra_collision_budget=extra_collision_budget,
                )
                model_hashes.append(hashlib.sha256(str(model.Proto()).encode()).hexdigest())
                solver = cp_model.CpSolver()
                solver.parameters.num_search_workers = 1
                solver.parameters.random_seed = 0
                solver.parameters.cp_model_presolve = True
                status = solver.Solve(model)
                status_name = solver.StatusName(status)
                statuses[status_name] = statuses.get(status_name, 0) + 1
                if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
                    def vals(name: str) -> list[int]:
                        return [int(solver.Value(v)) for v in data[name]]

                    return {
                        "collision": collision,
                        "target_nonexact_energy": target,
                        "joint_UGJ6_compatible": True,
                        "attempts_before_witness": attempts,
                        "expected_total_fixed_models": expected,
                        "fixed_model_sha256_aggregate": digest(model_hashes),
                        "solver_status_counts": statuses,
                        "witness": {
                            "hard_sign": hard_sign,
                            "roots": list(roots),
                            "remaining_hard_directions": list(data["remaining"]),
                            "excess_assignment": list(assignment),
                            "U_coefficients": vals("u_coefficients"),
                            "G_quotient_coefficients": vals("g_coefficients"),
                            "J6_quotient_coefficients": vals("j_coefficients"),
                            "U_evaluations": vals("u_values"),
                            "G_evaluations": vals("g_values"),
                            "J6_evaluations": vals("j_values"),
                            "row_energies": list(zip(data["energy_names"], vals("energies"))),
                        },
                    }
                if status != cp_model.INFEASIBLE:
                    raise ArithmeticError(f"unexpected solver status {status_name}")
    if attempts != expected:
        raise ArithmeticError("fixed-model coverage count changed")
    return {
        "collision": collision,
        "target_nonexact_energy": target,
        "joint_UGJ6_compatible": False,
        "all_fixed_models_infeasible": True,
        "fixed_models_checked": attempts,
        "expected_total_fixed_models": expected,
        "fixed_model_sha256_aggregate": digest(model_hashes),
        "solver_status_counts": statuses,
        "witness": None,
    }


COEFFICIENTS_2 = np.asarray(
    list(itertools.product(range(P), repeat=3)), dtype=np.int16
)
QUADRATIC_EVALUATIONS = (
    COEFFICIENTS_2 @ np.asarray(U_BASIS, dtype=np.int16).T % P
)
VALID_KEYS = {name: np.zeros(P**3, dtype=bool) for name in TABLES}
ENERGY_MASKS: dict[str, list[int]] = {
    name: [0] * (P**3) for name in TABLES
}
for table_name, table in TABLES.items():
    for u_value, g_value, j_value, energy in table:
        key = (u_value * P + g_value) * P + j_value
        VALID_KEYS[table_name][key] = True
        ENERGY_MASKS[table_name][key] |= 1 << energy


def mask_values(mask: int) -> list[int]:
    values = []
    while mask:
        bit = mask & -mask
        values.append(bit.bit_length() - 1)
        mask ^= bit
    return values


def convolve_energy_masks(choices: list[tuple[str, int, int]]) -> int:
    reachable = 1
    for _label, _direction, mask in choices:
        next_reachable = 0
        for energy in mask_values(mask):
            next_reachable |= reachable << energy
        reachable = next_reachable
        if reachable == 0:
            break
    return reachable


def recover_energy_sum(
    choices: list[tuple[str, int, int]], target: int
) -> list[tuple[str, int, int]] | None:
    states: dict[int, list[tuple[str, int, int]]] = {0: []}
    for label, direction, mask in choices:
        next_states: dict[int, list[tuple[str, int, int]]] = {}
        for subtotal, path in states.items():
            for energy in mask_values(mask):
                new_total = subtotal + energy
                if new_total <= target and new_total not in next_states:
                    next_states[new_total] = path + [(label, direction, energy)]
        states = next_states
        if not states:
            return None
    return states.get(target)


def coefficient_join(case: dict[str, object]) -> tuple[list[dict[str, object]], dict[str, int]]:
    """Exhaust all U quadratic, scalar G/R4, and quadratic J6/R4 forms."""
    partition = tuple(int(x) for x in case["partition"])
    collision_minimum = int(case["collision_minimum"])
    target_rows = [(int(c), int(t)) for c, t in case["targets"]]
    witnesses: dict[tuple[int, int], dict[str, object] | None] = {
        row: None for row in target_rows
    }
    target_stats: dict[tuple[int, int], dict[str, int]] = {
        row: {
            "hard_locally_compatible_form_assignments": 0,
            "maximum_separable_nonexact_energy": -1,
            "reachable_energy_union_mask": 0,
        }
        for row in target_rows
    }
    signs = direction_signs()
    coverage = {
        "hard_sign_root_sets_checked": 0,
        "U_coefficients_passing_exact_root_QR": 0,
        "UGJ6_coefficient_triples_before_row_pruning": 0,
        "UGJ6_pairs_after_all_seven_opposite_rows": 0,
        "hard_excess_assignment_checks": 0,
    }

    for hard_sign in (-1, 1):
        hard = tuple(i for i, sign in enumerate(signs) if sign == hard_sign)
        opposite = tuple(i for i, sign in enumerate(signs) if sign == -hard_sign)
        for roots in itertools.combinations(hard, 4):
            coverage["hard_sign_root_sets_checked"] += 1
            remaining = tuple(i for i in hard if i not in roots)
            root_values = np.asarray(root_product_values(roots), dtype=np.int16)
            j_evaluations = QUADRATIC_EVALUATIONS * root_values % P
            u_indices = np.flatnonzero(
                np.all(
                    np.isin(
                        QUADRATIC_EVALUATIONS[:, roots],
                        np.asarray(QR_NONZERO, dtype=np.int16),
                    ),
                    axis=1,
                )
            )
            coverage["U_coefficients_passing_exact_root_QR"] += len(u_indices)
            coverage["UGJ6_coefficient_triples_before_row_pruning"] += (
                len(u_indices) * P * len(COEFFICIENTS_2)
            )

            for g_coefficient in range(P):
                g_evaluations = g_coefficient * root_values % P
                first = opposite[0]
                first_keys = (
                    (
                        QUADRATIC_EVALUATIONS[u_indices, first, None] * P
                        + int(g_evaluations[first])
                    )
                    * P
                    + j_evaluations[None, :, first]
                )
                u_positions, j_indices = np.nonzero(VALID_KEYS["O"][first_keys])
                surviving_u = u_indices[u_positions]
                for direction in opposite[1:]:
                    keys = (
                        (
                            QUADRATIC_EVALUATIONS[surviving_u, direction] * P
                            + int(g_evaluations[direction])
                        )
                        * P
                        + j_evaluations[j_indices, direction]
                    )
                    keep = VALID_KEYS["O"][keys]
                    surviving_u = surviving_u[keep]
                    j_indices = j_indices[keep]
                    if len(surviving_u) == 0:
                        break
                coverage["UGJ6_pairs_after_all_seven_opposite_rows"] += len(
                    surviving_u
                )

                for u_index, j_index in zip(surviving_u, j_indices):
                    u_index = int(u_index)
                    j_index = int(j_index)
                    opposite_choices: list[tuple[str, int, int]] = []
                    for direction in opposite:
                        key = (
                            (
                                int(QUADRATIC_EVALUATIONS[u_index, direction]) * P
                                + int(g_evaluations[direction])
                            )
                            * P
                            + int(j_evaluations[j_index, direction])
                        )
                        opposite_choices.append(
                            ("O", direction, ENERGY_MASKS["O"][key])
                        )

                    for assignment in sorted(set(itertools.permutations(partition))):
                        coverage["hard_excess_assignment_checks"] += 1
                        for collision, target in target_rows:
                            if witnesses[(collision, target)] is not None:
                                continue
                            extra = collision - collision_minimum
                            choices = list(opposite_choices)
                            compatible = True
                            for direction, excess in zip(remaining, assignment):
                                table_name = f"H3_D{extra}" if excess == 3 else f"H{excess}"
                                key = (
                                    (
                                        int(QUADRATIC_EVALUATIONS[u_index, direction]) * P
                                        + int(g_evaluations[direction])
                                    )
                                    * P
                                    + int(j_evaluations[j_index, direction])
                                )
                                mask = ENERGY_MASKS[table_name][key]
                                if mask == 0:
                                    compatible = False
                                    break
                                choices.append((table_name, direction, mask))
                            if not compatible:
                                continue
                            stats = target_stats[(collision, target)]
                            stats["hard_locally_compatible_form_assignments"] += 1
                            separable_maximum = sum(
                                mask.bit_length() - 1 for _label, _direction, mask in choices
                            )
                            stats["maximum_separable_nonexact_energy"] = max(
                                stats["maximum_separable_nonexact_energy"],
                                separable_maximum,
                            )
                            reachable = convolve_energy_masks(choices)
                            stats["reachable_energy_union_mask"] |= reachable
                            if ((reachable >> target) & 1) == 0:
                                continue
                            energy_witness = recover_energy_sum(choices, target)
                            if energy_witness is None:
                                raise ArithmeticError("bitset target lost its energy witness")
                            witnesses[(collision, target)] = {
                                "hard_sign": hard_sign,
                                "roots": list(roots),
                                "remaining_hard_directions": list(remaining),
                                "excess_assignment": list(assignment),
                                "U_coefficients": COEFFICIENTS_2[u_index].tolist(),
                                "G_quotient_scalar": g_coefficient,
                                "J6_quotient_coefficients": COEFFICIENTS_2[j_index].tolist(),
                                "U_evaluations": QUADRATIC_EVALUATIONS[u_index].tolist(),
                                "G_evaluations": g_evaluations.tolist(),
                                "J6_evaluations": j_evaluations[j_index].tolist(),
                                "row_energies": [list(row) for row in energy_witness],
                                "energy_sum": sum(row[2] for row in energy_witness),
                            }

            if coverage["hard_sign_root_sets_checked"] % 10 == 0:
                print(
                    "progress root_sets={} survivors={} found={}".format(
                        coverage["hard_sign_root_sets_checked"],
                        coverage["UGJ6_pairs_after_all_seven_opposite_rows"],
                        sum(value is not None for value in witnesses.values()),
                    ),
                    flush=True,
                )
            if all(value is not None for value in witnesses.values()):
                break
        if all(value is not None for value in witnesses.values()):
            break

    exhaustive = coverage["hard_sign_root_sets_checked"] == 70
    results = []
    for collision, target in target_rows:
        witness = witnesses[(collision, target)]
        stats = target_stats[(collision, target)]
        if witness is None and not exhaustive:
            raise ArithmeticError("missing witness without exhaustive coverage")
        reachable_values = mask_values(stats["reachable_energy_union_mask"])
        results.append(
            {
                "collision": collision,
                "target_nonexact_energy": target,
                "joint_UGJ6_compatible": witness is not None,
                "all_form_coefficients_exhausted": witness is None,
                "hard_locally_compatible_form_assignments": stats[
                    "hard_locally_compatible_form_assignments"
                ],
                "maximum_separable_nonexact_energy": stats[
                    "maximum_separable_nonexact_energy"
                ],
                "reachable_energy_value_count": len(reachable_values),
                "reachable_energy_minimum": min(reachable_values, default=None),
                "reachable_energy_maximum": max(reachable_values, default=None),
                "reachable_energy_values_sha256": digest(reachable_values),
                "witness": witness,
            }
        )
    return results, coverage


def main(case_index: int = DEFAULT_CASE_INDEX, output: Path | None = None) -> None:
    if case_index not in range(len(CASES)):
        raise ValueError(f"case index must be in 0..{len(CASES) - 1}")
    case = CASES[case_index]
    if output is None:
        suffix = "221" if case_index == 0 else "311"
        output = ROOT / "evidence" / f"e1_gmin_m4_prop15754_four_root_{suffix}.json"
    print(
        f"catalogs ready; explicit coefficient join; case={case_index} "
        f"partition={case['partition']}",
        flush=True,
    )
    rows, coverage = coefficient_join(case)
    for row in rows:
        print(
            f"done C={row['collision']} compatible={row['joint_UGJ6_compatible']}",
            flush=True,
        )
    all_targets_excluded = all(
        not bool(row["joint_UGJ6_compatible"]) for row in rows
    )
    if not all_targets_excluded:
        raise ArithmeticError("a four-root U/G/J6 collision stratum survived")
    result = {
        "title": "p13 t4 u6 four-root joint U/G/J6 exact necessary-moment checker",
        "result_status": "exact finite necessary-moment certificate",
        "repository_HEAD": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "script_path": display_path(Path(__file__)),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "output_path": display_path(output),
        "algorithm": {
            "name": "explicit finite-field coefficient evaluation join plus exact energy bitsets",
            "U_coefficient_vectors": len(COEFFICIENTS_2),
            "G_quotient_scalars": P,
            "J6_quotient_coefficient_vectors": len(COEFFICIENTS_2),
            "fixed_hard_sign_root_sets_if_exhaustive": 70,
            "coverage": coverage,
        },
        "field": "F_13",
        "projective_points": [list(point) for point in POINTS],
        "direction_signs": list(direction_signs()),
        "hard_sign_classes_checked": [-1, 1],
        "imported_helper_sha256": ACTUAL_HELPER_SHA256,
        "moment_weights": {"W2": W2.tolist(), "W4": W4.tolist(), "W6": W6.tolist()},
        "W6_regression": list((1, 12, 1, 1, 12, 12)),
        "common_forms": {
            "U": "h*M2, binary quadratic",
            "G": "h*M4-M2^2=R4*c, binary quartic",
            "J6": "h*M6-(h*M2)^3=R4*Q2, binary sextic",
            "hard_key": "(U,G,J6)=(N2,N4-N2^2,N6-N2^3)",
            "opposite_key": "(U,G,J6)=(-N2,-N4-N2^2,-N6+N2^3)",
        },
        "translated_cut_catalog": {"count": len(CUTS), "sha256": digest(CUTS.tolist())},
        "row_catalogs": {
            name: {
                "row_count": len(rows_),
                "row_sha256": digest(rows_.tolist()),
                "table_tuple_count": len(TABLES[name]),
                "table_sha256": digest(TABLES[name]),
                "minimum_energy": min(int(row[3]) for row in TABLES[name]),
                "maximum_energy": max(int(row[3]) for row in TABLES[name]),
            }
            for name, rows_ in ROW_CATALOGS.items()
        },
        "partition": list(case["partition"]),
        "exact_XNOR_root_count": 4,
        "collision_minimum": case["collision_minimum"],
        "collision_rows": rows,
        "closed_collision_strata": [
            row["collision"] for row in rows if not row["joint_UGJ6_compatible"]
        ],
        "remaining_collision_strata": [
            row["collision"] for row in rows if row["joint_UGJ6_compatible"]
        ],
        "all_targets_excluded": all_targets_excluded,
        "proved": all_targets_excluded,
        "scope": (
            "Necessary common U/G/J6 form and exact row-energy compatibility only. "
            "A survivor is not a graph realization; exhaustive INFEASIBLE fixed "
            "models exclude the stated stratum."
        ),
    }
    temporary = output.with_suffix(".json.tmp")
    temporary.write_bytes(json.dumps(result, indent=2, sort_keys=True).encode() + b"\n")
    os.replace(temporary, output)
    print(json.dumps({"output": str(output), "closed": result["closed_collision_strata"], "remaining": result["remaining_collision_strata"]}, sort_keys=True), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case-index", type=int, default=DEFAULT_CASE_INDEX, choices=range(len(CASES))
    )
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    main(case_index=arguments.case_index, output=arguments.output)
