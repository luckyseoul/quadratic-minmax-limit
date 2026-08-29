#!/usr/bin/env python3
"""Native-XOR realizability test for the fourteen p=19 survivors.

This is the exact parity-profile problem underlying Proposition 15.689.
It uses both affine Radon equations

    r = A x,             x = A^T r  (mod 2),

where the second identity follows from ``A^T A=I+J`` and ``|x|=16``.
Direction weights and their phase histograms are imposed by guarded exact
cardinality automata.  A pair in a phase-zero b=0 direction is normalized
to field elements 0 and 1.  The optional ``--anchor-second`` mode instead
fixes two phase-zero b=16 directions, translates one selected point to
zero, and uses scalar dilation to fix one missing nonzero fibre.  Sweeping
the nine second-direction slots preserves the complete normalized search.

UNSAT excludes the boundary profile itself.  SAT gives an audited affine
boundary but does not by itself realize a residual edge lift.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import functools
import json
import os
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
)
from e1_gmin_m4_prop15688 import p19_residue_zero_profiles  # noqa: E402


P = 19
SIZE = 16


def atomic_write(path: Path, payload: object) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def survivor_profiles() -> list[dict[str, object]]:
    rows = [
        row
        for row in p19_residue_zero_profiles()["profiles"]
        if int(row["pair_slack"]) >= 16
    ]
    if len(rows) != 14:
        raise ArithmeticError("the p=19 high-slack remainder changed")
    return rows


def solve(
    profile_index: int,
    seconds: float,
    threads: int,
    anchor_second: int | None = None,
    repair_five: bool = False,
    fixed_phase_zero_b16_slots: tuple[int, ...] | None = None,
) -> dict[str, object]:
    from pycryptosat import Solver

    profiles = survivor_profiles()
    if not 0 <= profile_index < len(profiles):
        raise ValueError("profile index must be in 0..13")
    profile = profiles[profile_index]
    directions = projective_directions(P)
    started = time.time()

    next_id = 0

    def new_var() -> int:
        nonlocal next_id
        next_id += 1
        return next_id

    point = [new_var() for _ in range(P * P)]
    parity = [[new_var() for _ in range(P)] for _ in directions]
    selectors: list[dict[int, int]] = []
    phases = []
    records = []
    normalized_index = None
    for index, direction in enumerate(directions):
        eps, labels = field_direction_data(P, direction)
        phase = 0 if eps == 1 else 1  # c_H=-1
        allowed = {
            int(b): int(count)
            for b, count in profile["phase_profiles_b"][str(phase)].items()
        }
        selectors.append({b: new_var() for b in allowed})
        phases.append(phase)
        records.append((direction, eps, labels))
        if labels[0] == labels[1]:
            if normalized_index is not None or eps != 1:
                raise ArithmeticError("normalized direction changed")
            normalized_index = index
    if normalized_index is None:
        raise ArithmeticError("normalized direction missing")

    deleted = [new_var() for _ in point] if repair_five else []
    core = [new_var() for _ in point] if repair_five else []
    core_secant = (
        [[new_var() for _ in range(P)] for _ in directions]
        if repair_five
        else []
    )
    core_occupied = (
        [[new_var() for _ in range(P)] for _ in directions]
        if repair_five
        else []
    )
    semantic_variables = next_id

    solver = Solver(verbose=0, threads=max(1, int(threads)))
    clause_count = 0
    xor_count = 0
    cardinality_count = 0

    def add_clauses(clauses: list[list[int]]) -> None:
        nonlocal clause_count
        if clauses:
            solver.add_clauses(clauses)
            clause_count += len(clauses)

    def add_unit(literal: int) -> None:
        add_clauses([[int(literal)]])

    def add_exact(literals: list[int], bound: int, guard: int | None = None) -> None:
        """Guarded exact cardinality using a deterministic count automaton.

        The clauses are active when ``guard`` is false.  Counts above half
        are encoded through complemented literals to keep the automaton small.
        """
        nonlocal cardinality_count
        literals = [int(literal) for literal in literals]
        bound = int(bound)
        if bound > len(literals) // 2:
            literals = [-literal for literal in literals]
            bound = len(literals) - bound
        if not 0 <= bound <= len(literals):
            add_clauses([[]] if guard is None else [[int(guard)]])
            cardinality_count += 1
            return
        overflow = bound + 1
        states = [
            [new_var() for _ in range(overflow + 1)]
            for _ in range(len(literals) + 1)
        ]

        def gated(clause: list[int]) -> list[int]:
            return clause if guard is None else [int(guard), *clause]

        clauses: list[list[int]] = []
        for row in states:
            clauses.append(gated(list(row)))
            for first in range(len(row)):
                for second in range(first + 1, len(row)):
                    clauses.append(gated([-row[first], -row[second]]))
        clauses.append(gated([states[0][0]]))
        clauses.append(gated([states[-1][bound]]))
        for index, literal in enumerate(literals, start=1):
            previous = states[index - 1]
            current = states[index]
            for count in range(overflow + 1):
                incremented = min(overflow, count + 1)
                clauses.append(gated([-previous[count], literal, current[count]]))
                clauses.append(
                    gated([-previous[count], -literal, current[incremented]])
                )
        add_clauses(clauses)
        cardinality_count += 1

    add_exact(point, SIZE)
    add_unit(point[0])
    if anchor_second is None:
        add_unit(point[1])

    incident_parities: list[list[int]] = [[] for _ in point]
    for index, (_direction, _eps, labels) in enumerate(records):
        fibres = []
        for fibre in range(P):
            fibre_points = [
                point[v] for v, label in enumerate(labels) if label == fibre
            ]
            fibres.append(fibre_points)
            solver.add_xor_clause(
                [*fibre_points, parity[index][fibre]],
                False,
            )
            xor_count += 1
        for v, label in enumerate(labels):
            incident_parities[v].append(parity[index][label])

        direction_selectors = list(selectors[index].values())
        add_clauses([direction_selectors])
        for first in range(len(direction_selectors)):
            for second in range(first + 1, len(direction_selectors)):
                add_clauses(
                    [[-direction_selectors[first], -direction_selectors[second]]]
                )
        for b, selector in selectors[index].items():
            add_exact(parity[index], b, guard=-selector)
            if b == SIZE:
                # With sixteen selected points, parity weight sixteen means
                # sixteen singleton fibres.  Expose that implication in CNF.
                for fibre_points in fibres:
                    for first in range(len(fibre_points)):
                        for second in range(first + 1, len(fibre_points)):
                            add_clauses(
                                [[-selector, -fibre_points[first], -fibre_points[second]]]
                            )

    for v, rows in enumerate(incident_parities):
        if len(rows) != P + 1:
            raise ArithmeticError("affine point degree changed")
        solver.add_xor_clause([point[v], *rows], False)
        xor_count += 1

    # Exact direction-weight histograms within each Paley phase.
    for phase in (0, 1):
        phase_profile = {
            int(b): int(count)
            for b, count in profile["phase_profiles_b"][str(phase)].items()
        }
        phase_indices = [i for i, value in enumerate(phases) if value == phase]
        if len(phase_indices) != 10:
            raise ArithmeticError("quadratic direction split changed")
        for b, target in phase_profile.items():
            add_exact([selectors[i][b] for i in phase_indices], target)

    phase_zero_indices = [i for i, value in enumerate(phases) if value == 0]
    if fixed_phase_zero_b16_slots is not None:
        slots = set(fixed_phase_zero_b16_slots)
        target = int(profile["phase_profiles_b"]["0"].get(16, 0))
        if (
            len(slots) != len(fixed_phase_zero_b16_slots)
            or len(slots) != target
            or not slots <= set(range(len(phase_zero_indices)))
        ):
            raise ValueError(
                "fixed phase-zero b16 slots must be distinct slots matching the profile"
            )
        for slot, index in enumerate(phase_zero_indices):
            add_unit(
                selectors[index][16]
                if slot in slots
                else -selectors[index][16]
            )

    if repair_five:
        # Proposition 15.693 forces every slack-twenty witness to require
        # exactly five repair deletions. If one deleted point had no core
        # secant, restoring it would give a forbidden four-deletion repair.
        # The total slack inequality then forces exactly one core secant
        # through each of the five deleted points.
        add_exact(deleted, 5)
        # Translate a retained core point to zero before applying the
        # direction and missing-fibre normalizations. This is lossless.
        add_unit(core[0])
        for v in range(len(point)):
            add_clauses(
                [
                    [-deleted[v], point[v]],
                    [-core[v], point[v]],
                    [-core[v], -deleted[v]],
                    [-point[v], core[v], deleted[v]],
                ]
            )

        incident_core_secants: list[list[int]] = [[] for _ in point]
        for index, (_direction, _eps, labels) in enumerate(records):
            for fibre in range(P):
                fibre_core = [
                    core[v] for v, label in enumerate(labels) if label == fibre
                ]
                fibre_deleted = [
                    deleted[v]
                    for v, label in enumerate(labels)
                    if label == fibre
                ]
                secant = core_secant[index][fibre]
                occupied = core_occupied[index][fibre]

                # The eleven retained points form an affine arc.
                arc_clauses = []
                for first in range(len(fibre_core)):
                    for second in range(first):
                        for third in range(second):
                            arc_clauses.append(
                                [
                                    -fibre_core[first],
                                    -fibre_core[second],
                                    -fibre_core[third],
                                ]
                            )
                add_clauses(arc_clauses)

                # secant iff this affine line contains two core points.
                add_exact(fibre_core, 2, guard=-secant)
                pair_clauses = []
                for first in range(len(fibre_core)):
                    for second in range(first):
                        pair_clauses.append(
                            [-fibre_core[first], -fibre_core[second], secant]
                        )
                add_clauses(pair_clauses)

                # Slack equality leaves no uncharged bad line. The deleted
                # five-set is therefore an arc. If two deleted points share
                # a line containing a core point, that line must contain two
                # core points and hence be one of the charged core secants.
                add_clauses([[-occupied, *fibre_core]])
                add_clauses([[-literal, occupied] for literal in fibre_core])
                deleted_arc_clauses = []
                deleted_pair_clauses = []
                for first in range(len(fibre_deleted)):
                    for second in range(first):
                        deleted_pair_clauses.append(
                            [
                                -fibre_deleted[first],
                                -fibre_deleted[second],
                                -occupied,
                                secant,
                            ]
                        )
                        for third in range(second):
                            deleted_arc_clauses.append(
                                [
                                    -fibre_deleted[first],
                                    -fibre_deleted[second],
                                    -fibre_deleted[third],
                                ]
                            )
                add_clauses(deleted_pair_clauses)
                add_clauses(deleted_arc_clauses)

            for v, label in enumerate(labels):
                incident_core_secants[v].append(core_secant[index][label])

        for v, rows in enumerate(incident_core_secants):
            if len(rows) != P + 1:
                raise ArithmeticError("affine core point degree changed")
            add_exact(rows, 1, guard=-deleted[v])

    if anchor_second is None:
        add_unit(selectors[normalized_index][0])
        normalization = {
            "mode": "phase-zero-b0-pair",
            "selected_points": [0, 1],
            "direction": list(directions[normalized_index]),
            "phase": 0,
            "b": 0,
            "c_H": -1,
        }
    else:
        if not 1 <= anchor_second < len(phase_zero_indices):
            raise ValueError("anchor_second must be in 1..9")
        if (
            fixed_phase_zero_b16_slots is not None
            and (
                0 not in fixed_phase_zero_b16_slots
                or anchor_second not in fixed_phase_zero_b16_slots
            )
        ):
            raise ValueError("the two anchored directions must be fixed at b=16")
        first_index = phase_zero_indices[0]
        second_index = phase_zero_indices[anchor_second]
        if 16 not in selectors[first_index] or 16 not in selectors[second_index]:
            raise ValueError("profile lacks two phase-zero undetermined directions")
        add_unit(selectors[first_index][16])
        add_unit(selectors[second_index][16])

        # After translating a selected point to zero, F_p^* scales the
        # quotient by the first undetermined direction transitively on its
        # nonzero fibres.  Normalize one of the three missing fibres.
        first_labels = records[first_index][2]
        base_label = first_labels[0]
        missing_label = next(label for label in first_labels if label != base_label)
        fixed_missing_fibre = [
            v for v, label in enumerate(first_labels) if label == missing_label
        ]
        for v in fixed_missing_fibre:
            add_unit(-point[v])
        normalization = {
            "mode": "two-phase-zero-undetermined",
            "selected_point": 0,
            "first_direction": list(directions[first_index]),
            "second_direction": list(directions[second_index]),
            "second_phase_zero_slot": anchor_second,
            "phase": 0,
            "b": 16,
            "fixed_missing_fibre_label": missing_label,
            "c_H": -1,
        }

    build_seconds = time.time() - started
    satisfiable, assignment = solver.solve(time_limit=float(seconds))
    status = (
        "SATISFIABLE"
        if satisfiable is True
        else "UNSATISFIABLE"
        if satisfiable is False
        else "UNKNOWN"
    )
    result: dict[str, object] = {
        "experiment": "p19_second_boundary_profile_cryptominisat",
        "profile_index": profile_index,
        "pair_slack": int(profile["pair_slack"]),
        "phase_profiles_b": profile["phase_profiles_b"],
        "normalization": normalization,
        "fixed_phase_zero_b16_slots": fixed_phase_zero_b16_slots,
        "repair_five_constraints": repair_five,
        "solver": "cryptominisat-native-xor",
        "solver_status": status,
        "feasible_boundary_profile": satisfiable is True,
        "finite_infeasibility_only": satisfiable is False,
        "semantic_variables": semantic_variables,
        "total_variables": next_id,
        "clauses": clause_count,
        "native_xor_constraints": xor_count,
        "cardinality_constraints": cardinality_count,
        "threads": threads,
        "build_seconds": build_seconds,
        "solve_seconds": time.time() - started - build_seconds,
        "elapsed_seconds": time.time() - started,
    }
    if satisfiable is True:
        chosen = [v for v, literal in enumerate(point) if assignment[literal]]
        observed = {0: {}, 1: {}}
        direction_rows = []
        valid = len(chosen) == SIZE and 0 in chosen
        if anchor_second is None:
            valid = valid and 1 in chosen
        for index, (direction, eps, labels) in enumerate(records):
            counts = [0] * P
            for v in chosen:
                counts[labels[v]] += 1
            b = sum(count & 1 for count in counts)
            phase = phases[index]
            observed[phase][b] = observed[phase].get(b, 0) + 1
            direction_rows.append(
                {"direction": list(direction), "eps": eps, "phase": phase, "b": b}
            )
        expected = {
            phase: {
                int(b): int(count)
                for b, count in profile["phase_profiles_b"][str(phase)].items()
            }
            for phase in (0, 1)
        }
        valid = valid and observed == expected
        repair_witness = None
        if repair_five:
            chosen_deleted = [v for v in chosen if assignment[deleted[v]]]
            chosen_core = [v for v in chosen if assignment[core[v]]]
            core_line_counts = []
            deleted_line_counts = []
            for _direction, _eps, labels in records:
                counts = [0] * P
                deleted_counts = [0] * P
                for v in chosen_core:
                    counts[labels[v]] += 1
                for v in chosen_deleted:
                    deleted_counts[labels[v]] += 1
                core_line_counts.append(counts)
                deleted_line_counts.append(deleted_counts)
            deleted_secant_counts = [
                sum(
                    core_line_counts[index][labels[v]] == 2
                    for index, (_direction, _eps, labels) in enumerate(records)
                )
                for v in chosen_deleted
            ]
            core_set = set(chosen_core)
            deleted_set = set(chosen_deleted)
            core_is_arc = max(
                max(counts) for counts in core_line_counts
            ) <= 2
            slack_equality_lines = all(
                deleted_count <= 2
                and not (deleted_count == 2 and core_count == 1)
                for core_counts, deleted_counts in zip(
                    core_line_counts, deleted_line_counts
                )
                for core_count, deleted_count in zip(
                    core_counts, deleted_counts
                )
            )
            valid = valid and (
                len(chosen_deleted) == 5
                and len(chosen_core) == 11
                and core_set | deleted_set == set(chosen)
                and not core_set & deleted_set
                and core_is_arc
                and slack_equality_lines
                and deleted_secant_counts == [1] * 5
            )
            repair_witness = {
                "deleted_points": chosen_deleted,
                "core_points": chosen_core,
                "deleted_core_secant_counts": deleted_secant_counts,
                "core_is_arc": core_is_arc,
                "line_slack_equality_structure": slack_equality_lines,
            }
        if not valid:
            raise AssertionError("CryptoMiniSat p=19 witness failed audit")
        result["boundary"] = chosen
        result["boundary_coordinates"] = [[v % P, v // P] for v in chosen]
        result["direction_rows"] = direction_rows
        result["repair_witness"] = repair_witness
        result["witness_audit_valid"] = True
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--profile", type=int)
    group.add_argument("--all", action="store_true")
    parser.add_argument(
        "--anchor-second",
        type=int,
        choices=range(1, 10),
        help="use the two-undetermined normalization with this phase-zero slot",
    )
    parser.add_argument(
        "--repair-five",
        action="store_true",
        help="encode the forced five-deletion 11-arc repair for slack twenty",
    )
    parser.add_argument(
        "--phase-zero-b16-slots",
        help="comma-separated exact phase-zero slots assigned b=16",
    )
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    fixed_phase_zero_b16_slots = (
        tuple(int(value) for value in args.phase_zero_b16_slots.split(","))
        if args.phase_zero_b16_slots
        else None
    )
    if args.all:
        worker = functools.partial(
            solve,
            seconds=args.seconds,
            threads=args.threads,
            anchor_second=args.anchor_second,
            repair_five=args.repair_five,
            fixed_phase_zero_b16_slots=fixed_phase_zero_b16_slots,
        )
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=max(1, int(args.jobs))
        ) as executor:
            result = list(executor.map(worker, range(14)))
    else:
        result = solve(
            args.profile,
            args.seconds,
            args.threads,
            anchor_second=args.anchor_second,
            repair_five=args.repair_five,
            fixed_phase_zero_b16_slots=fixed_phase_zero_b16_slots,
        )
    rendered = json.dumps(result, indent=2, sort_keys=True)
    print(rendered, flush=True)
    if args.output is not None:
        atomic_write(args.output, result)


if __name__ == "__main__":
    main()
