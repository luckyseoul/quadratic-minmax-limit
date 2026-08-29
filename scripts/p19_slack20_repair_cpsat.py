#!/usr/bin/env python3
"""Exact CP-SAT attack on the four p=19 pair-slack-20 profiles.

Proposition 15.693 forces exactly five repair deletions.  Write the putative
boundary as ``S=A disjoint_union D`` with ``|A|=11`` and ``|D|=5``.  Equality
in the secant/slack bound forces all of the following:

* both ``A`` and ``D`` are affine arcs;
* every point of ``D`` lies on exactly one secant of ``A``;
* a line containing two points of ``D`` contains zero or two points of ``A``;
* the sum of the exact line slacks is twenty.

The model retains the exact phase-labelled odd-fibre profile.  Translation,
phase-preserving field multiplication, and scalar dilation normalize a core
point to zero, one phase-zero undetermined direction to slot zero, and one of
its three missing fibres to label one.  A complete proof requires sweeping
all normalized choices for the remaining undetermined directions; UNKNOWN is
never treated as evidence.
"""
from __future__ import annotations

import argparse
import json
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
SLACK = 20
LIVE_PROFILE_INDICES = (3, 4, 6, 8)


def survivor_profiles() -> list[dict[str, object]]:
    rows = [
        row
        for row in p19_residue_zero_profiles()["profiles"]
        if int(row["pair_slack"]) >= 16
    ]
    if len(rows) != 14:
        raise ArithmeticError("the p=19 high-slack remainder changed")
    return rows


def parse_slots(text: str | None) -> tuple[int, ...] | None:
    if text is None:
        return None
    return tuple(int(value) for value in text.split(","))


def solve(
    profile_index: int,
    anchor_second: int,
    phase_zero_b16_slots: tuple[int, ...] | None,
    first_missing_fibres: tuple[int, ...] | None,
    seconds: float,
    workers: int,
) -> dict[str, object]:
    from ortools.sat.python import cp_model

    if profile_index not in LIVE_PROFILE_INDICES:
        raise ValueError(f"profile must be one of {LIVE_PROFILE_INDICES}")
    profile = survivor_profiles()[profile_index]
    if int(profile["pair_slack"]) != SLACK:
        raise ArithmeticError("selected profile is not a slack-twenty row")

    started = time.time()
    model = cp_model.CpModel()
    selected = [model.new_bool_var(f"x_{v}") for v in range(P * P)]
    core = [model.new_bool_var(f"a_{v}") for v in range(P * P)]
    deleted = [model.new_bool_var(f"d_{v}") for v in range(P * P)]
    model.add(sum(selected) == SIZE)
    model.add(sum(core) == 11)
    model.add(sum(deleted) == 5)
    for v in range(P * P):
        model.add(selected[v] == core[v] + deleted[v])

    # Lossless affine normalization: choose a retained point as the origin.
    model.add(core[0] == 1)

    directions = projective_directions(P)
    records = []
    phases = []
    phase_indices: dict[int, list[int]] = {0: [], 1: []}
    line_rows = []
    incident_secants: list[list[object]] = [[] for _ in selected]
    total_slacks = []
    b_vars = []

    for direction_index, direction in enumerate(directions):
        eps, labels = field_direction_data(P, direction)
        phase = 0 if eps == 1 else 1  # c_H=-1
        records.append((direction, eps, labels))
        phases.append(phase)
        phase_indices[phase].append(direction_index)

        parity_vars = []
        direction_lines = []
        for fibre in range(P):
            vertices = [v for v, label in enumerate(labels) if label == fibre]
            x_occ = model.new_int_var(0, 4, f"xocc_{direction_index}_{fibre}")
            a_occ = model.new_int_var(0, 2, f"aocc_{direction_index}_{fibre}")
            d_occ = model.new_int_var(0, 2, f"docc_{direction_index}_{fibre}")
            model.add(x_occ == sum(selected[v] for v in vertices))
            model.add(a_occ == sum(core[v] for v in vertices))
            model.add(d_occ == sum(deleted[v] for v in vertices))
            # Equality in slack >= 4 sum mu forbids one core plus two
            # deletions; the arc bounds forbid three of either kind.
            model.add_allowed_assignments(
                [a_occ, d_occ, x_occ],
                [
                    (0, 0, 0),
                    (0, 1, 1),
                    (0, 2, 2),
                    (1, 0, 1),
                    (1, 1, 2),
                    (2, 0, 2),
                    (2, 1, 3),
                    (2, 2, 4),
                ],
            )

            secant = model.new_bool_var(f"asec_{direction_index}_{fibre}")
            model.add(a_occ == 2).only_enforce_if(secant)
            model.add(a_occ <= 1).only_enforce_if(secant.Not())
            for v in vertices:
                incident_secants[v].append(secant)

            parity = model.new_bool_var(f"r_{direction_index}_{fibre}")
            line_slack = model.new_int_var_from_domain(
                cp_model.Domain.from_values([0, 4, 8]),
                f"slack_{direction_index}_{fibre}",
            )
            model.add_allowed_assignments(
                [x_occ, parity, line_slack],
                [(0, 0, 0), (1, 1, 0), (2, 0, 0), (3, 1, 4), (4, 0, 8)],
            )
            parity_vars.append(parity)
            total_slacks.append(line_slack)
            direction_lines.append((x_occ, a_occ, d_occ, secant, parity))

        b_var = model.new_int_var(0, SIZE, f"b_{direction_index}")
        model.add(b_var == sum(parity_vars))
        b_vars.append(b_var)
        line_rows.append(direction_lines)

    if any(len(rows) != P + 1 for rows in incident_secants):
        raise ArithmeticError("affine point degree changed")
    for v, rows in enumerate(incident_secants):
        model.add(sum(rows) == 1).only_enforce_if(deleted[v])
    model.add(sum(total_slacks) == SLACK)

    # Exact phase-labelled b histograms.
    indicators: dict[tuple[int, int], object] = {}
    for phase in (0, 1):
        profile_for_phase = {
            int(b): int(count)
            for b, count in profile["phase_profiles_b"][str(phase)].items()
        }
        if len(phase_indices[phase]) != 10:
            raise ArithmeticError("quadratic direction split changed")
        for index in phase_indices[phase]:
            row = []
            for b in profile_for_phase:
                indicator = model.new_bool_var(f"is_{index}_{b}")
                model.add(b_vars[index] == b).only_enforce_if(indicator)
                model.add(b_vars[index] != b).only_enforce_if(indicator.Not())
                indicators[index, b] = indicator
                row.append(indicator)
            model.add_exactly_one(row)
        for b, target in profile_for_phase.items():
            model.add(
                sum(indicators[index, b] for index in phase_indices[phase])
                == target
            )

    phase_zero = phase_indices[0]
    if not 1 <= anchor_second < len(phase_zero):
        raise ValueError("anchor_second must be in 1..9")
    target_b16 = int(profile["phase_profiles_b"]["0"].get(16, 0))
    if phase_zero_b16_slots is not None:
        slots = set(phase_zero_b16_slots)
        if (
            len(slots) != len(phase_zero_b16_slots)
            or len(slots) != target_b16
            or not slots <= set(range(10))
            or 0 not in slots
            or anchor_second not in slots
        ):
            raise ValueError("invalid normalized phase-zero b16 slot set")
        for slot, index in enumerate(phase_zero):
            if slot in slots:
                model.add(b_vars[index] == 16)
            else:
                model.add(b_vars[index] != 16)
    else:
        model.add(b_vars[phase_zero[0]] == 16)
        model.add(b_vars[phase_zero[anchor_second]] == 16)

    # Having fixed the first undetermined direction and translated a core
    # point to zero, scalar dilation is transitive on its nonzero fibres.
    first_labels = records[phase_zero[0]][2]
    base_label = first_labels[0]
    missing_label = next(label for label in first_labels if label != base_label)
    if first_missing_fibres is not None:
        missing_fibres = set(first_missing_fibres)
        if (
            len(missing_fibres) != 3
            or len(missing_fibres) != len(first_missing_fibres)
            or 0 in missing_fibres
            or missing_label not in missing_fibres
            or not missing_fibres <= set(range(P))
        ):
            raise ValueError("first missing fibres must be three distinct nonzero labels including one")
    else:
        missing_fibres = {missing_label}
    for v, label in enumerate(first_labels):
        if label in missing_fibres:
            model.add(selected[v] == 0)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = workers
    solver.parameters.random_seed = 15693000 + 100 * profile_index + anchor_second
    solver.parameters.symmetry_level = 2
    solver.parameters.linearization_level = 0
    status = solver.solve(model)
    feasible = status in (cp_model.OPTIMAL, cp_model.FEASIBLE)
    result: dict[str, object] = {
        "experiment": "p19_slack20_repair_cpsat",
        "profile_index": profile_index,
        "pair_slack": SLACK,
        "phase_profiles_b": profile["phase_profiles_b"],
        "normalization": {
            "core_origin": 0,
            "first_phase_zero_b16_slot": 0,
            "second_phase_zero_b16_slot": anchor_second,
            "fixed_missing_fibre_label": missing_label,
            "fixed_first_missing_fibres": first_missing_fibres,
            "fixed_phase_zero_b16_slots": phase_zero_b16_slots,
        },
        "solver_status": solver.status_name(status),
        "feasible": feasible,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time": solver.wall_time,
        "elapsed_seconds": time.time() - started,
        "rigorous_only_if_infeasible_or_witness_audited": True,
    }

    if feasible:
        chosen = [v for v in range(P * P) if solver.value(selected[v])]
        chosen_core = [v for v in chosen if solver.value(core[v])]
        chosen_deleted = [v for v in chosen if solver.value(deleted[v])]
        observed = {0: {}, 1: {}}
        exact_slack = 0
        deleted_mu = [0] * len(chosen_deleted)
        core_is_arc = True
        deleted_is_arc = True
        equality_lines = True
        direction_rows = []
        for index, (direction, eps, labels) in enumerate(records):
            x_counts = [0] * P
            a_counts = [0] * P
            d_counts = [0] * P
            for v in chosen:
                x_counts[labels[v]] += 1
            for v in chosen_core:
                a_counts[labels[v]] += 1
            for v in chosen_deleted:
                d_counts[labels[v]] += 1
            b = sum(count & 1 for count in x_counts)
            phase = phases[index]
            observed[phase][b] = observed[phase].get(b, 0) + 1
            exact_slack += sum(2 * ((count * (count - 1)) // 2 - count // 2) for count in x_counts)
            core_is_arc = core_is_arc and max(a_counts) <= 2
            deleted_is_arc = deleted_is_arc and max(d_counts) <= 2
            equality_lines = equality_lines and all(
                not (a_count == 1 and d_count == 2)
                for a_count, d_count in zip(a_counts, d_counts)
            )
            for offset, v in enumerate(chosen_deleted):
                deleted_mu[offset] += a_counts[labels[v]] == 2
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
        audited = (
            len(chosen) == 16
            and len(chosen_core) == 11
            and len(chosen_deleted) == 5
            and observed == expected
            and exact_slack == SLACK
            and core_is_arc
            and deleted_is_arc
            and equality_lines
            and deleted_mu == [1] * 5
        )
        if not audited:
            raise AssertionError("CP-SAT p=19 repair witness failed audit")
        result["witness"] = {
            "selected_points": chosen,
            "selected_coordinates": [[v % P, v // P] for v in chosen],
            "core_points": chosen_core,
            "deleted_points": chosen_deleted,
            "deleted_core_secant_counts": deleted_mu,
            "direction_rows": direction_rows,
            "exact_line_slack": exact_slack,
            "audited": True,
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", type=int, choices=LIVE_PROFILE_INDICES, required=True)
    parser.add_argument("--anchor-second", type=int, choices=range(1, 10), required=True)
    parser.add_argument("--phase-zero-b16-slots")
    parser.add_argument("--first-missing-fibres")
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = solve(
        args.profile,
        args.anchor_second,
        parse_slots(args.phase_zero_b16_slots),
        parse_slots(args.first_missing_fibres),
        args.seconds,
        args.workers,
    )
    rendered = json.dumps(result, indent=2, sort_keys=True)
    print(rendered, flush=True)
    if args.output is not None:
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
