#!/usr/bin/env python3
"""Exact native-XOR boundary solver for the p=17 slack-twenty block.

After the Proposition 15.723 replay, the corrected p=17 ledger has 193
pair-slack-twenty profiles.  Nine pairs of those arithmetic rows have
identical phase-labelled odd-fibre histograms, so the finite boundary problem
has 184 distinct signatures.  The former 78/69 counts came from the retracted
blanket floor-plus-two filter.

Proposition 15.707 subsequently excludes the full block without a solver.
This model remains an independent exact-boundary audit path; no solver result
from it is used in that proposition.

For a boundary vector ``x`` and affine line-incidence matrix ``A``, this model
imposes both binary Radon identities

    r = A x,                 x = A^T r  (mod 2).

The second identity is valid because ``wt(x)=16`` is even and, over F_2,
``A^T A=I+J``.  Guarded exact-cardinality automata impose every directional
weight and both Paley-phase histograms.

Every slack-twenty profile has a phase-zero direction of weight zero.  Its
sixteen selected points have even occupancy in every fibre of that direction,
so one fibre contains a selected pair.  A square affine similarity sends the
pair to field elements 0 and 1 and its direction to the canonical direction
``(0,1)`` without changing either phase histogram.  The constraints
``x_0=x_1=1`` and ``b_(0,1)=0`` are therefore lossless.

UNSAT excludes the finite boundary signature.  SAT returns a boundary that is
audited independently from the solver assignment; it does not by itself
provide a residual edge lift.  One invocation solves one signature only so
external process-level sharding remains explicit and crash-safe.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15632 import (  # noqa: E402
    field_direction_data,
    projective_directions,
)
from e1_gmin_m4_prop15700 import (  # noqa: E402
    p17_second_boundary_profile_census,
)


P = 17
SIZE = 16
PHASE_SIZE = 9
SLACK = 20
EXPECTED_PROFILE_COUNT = 193
EXPECTED_SIGNATURE_COUNT = 184
EXPECTED_CENSUS_START = 1364
EXPECTED_CENSUS_STOP = 1557

ProfileKey = tuple[
    tuple[tuple[int, int], ...],
    tuple[tuple[int, int], ...],
]


def atomic_write(path: Path, payload: object) -> None:
    """Atomically render one JSON record, creating its parent if necessary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.{time.time_ns()}.tmp")
    try:
        temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def profile_key(profile: dict[str, Any]) -> ProfileKey:
    """Canonical finite-boundary signature of one arithmetic profile."""
    phases = profile["phase_profiles_b"]
    return tuple(  # type: ignore[return-value]
        tuple(sorted((int(b), int(count)) for b, count in phases[str(phase)].items()))
        for phase in (0, 1)
    )


def _histogram_from_key(key: ProfileKey) -> dict[str, dict[int, int]]:
    return {
        str(phase): {int(b): int(count) for b, count in key[phase]}
        for phase in (0, 1)
    }


def slack20_profiles() -> list[dict[str, Any]]:
    """Return the corrected 193-row block with stable local/census indices."""
    census = p17_second_boundary_profile_census()
    selected = [
        (census_index, row)
        for census_index, row in enumerate(census["profiles"])
        if int(row["pair_slack"]) == SLACK
    ]
    census_indices = [index for index, _row in selected]
    if (
        len(selected) != EXPECTED_PROFILE_COUNT
        or census_indices
        != list(range(EXPECTED_CENSUS_START, EXPECTED_CENSUS_STOP))
    ):
        raise ArithmeticError("the p=17 slack-twenty census changed")

    rows: list[dict[str, Any]] = []
    for profile_index, (census_index, row) in enumerate(selected):
        phases = row["phase_profiles_b"]
        phase_zero_b0 = int(phases["0"].get(0, 0))
        if not 2 <= phase_zero_b0 <= 6:
            raise ArithmeticError("universal phase-zero b=0 normalization changed")
        if any(
            sum(int(count) for count in phases[str(phase)].values()) != PHASE_SIZE
            for phase in (0, 1)
        ):
            raise ArithmeticError("p=17 quadratic direction split changed")
        rows.append(
            {
                "profile_index": profile_index,
                "census_index": census_index,
                "u0": int(row["u0"]),
                "u1": int(row["u1"]),
                "phase_deficits": {
                    str(phase): int(row["phase_deficits"][str(phase)])
                    for phase in (0, 1)
                },
                "phase_profiles_b": {
                    str(phase): {
                        int(b): int(count)
                        for b, count in phases[str(phase)].items()
                    }
                    for phase in (0, 1)
                },
                "pair_slack": int(row["pair_slack"]),
            }
        )
    return rows


def _multinomial(counts: Iterable[int]) -> int:
    values = tuple(int(value) for value in counts)
    result = math.factorial(sum(values))
    for value in values:
        result //= math.factorial(value)
    return result


def _reflection_fixed_colourings(counts: Iterable[int], fixed_slots: int) -> int:
    """Colourings fixed by four swaps and zero/one fixed direction slots."""
    values = tuple(int(value) for value in counts)
    if fixed_slots == 0:
        if any(value & 1 for value in values):
            return 0
        return _multinomial(value // 2 for value in values)
    if fixed_slots != 1:
        raise ValueError("only the p=17 residual reflection is supported")

    total = 0
    for chosen, value in enumerate(values):
        if value % 2 != 1:
            continue
        if any(index != chosen and other % 2 for index, other in enumerate(values)):
            continue
        remaining = list(values)
        remaining[chosen] -= 1
        total += _multinomial(other // 2 for other in remaining)
    return total


def signature_manifest() -> dict[str, Any]:
    """Deduplicate 193 arithmetic profiles into 184 boundary signatures."""
    profiles = slack20_profiles()
    grouped: dict[ProfileKey, list[dict[str, Any]]] = defaultdict(list)
    for profile in profiles:
        grouped[profile_key(profile)].append(profile)

    signatures = []
    for signature_index, (key, members) in enumerate(grouped.items()):
        histogram = _histogram_from_key(key)
        reduced_zero = dict(histogram["0"])
        reduced_zero[0] -= 1  # canonical direction (0,1) is fixed at b=0
        if reduced_zero[0] < 0:
            raise ArithmeticError("normalization consumed a missing b=0 direction")

        raw_assignments = _multinomial(reduced_zero.values()) * _multinomial(
            histogram["1"].values()
        )
        reflection_fixed = _reflection_fixed_colourings(
            reduced_zero.values(), fixed_slots=0
        ) * _reflection_fixed_colourings(
            histogram["1"].values(), fixed_slots=1
        )
        reflection_orbits = (raw_assignments + reflection_fixed) // 2
        signatures.append(
            {
                "signature_index": signature_index,
                "phase_profiles_b": histogram,
                "profile_indices": [int(row["profile_index"]) for row in members],
                "census_indices": [int(row["census_index"]) for row in members],
                "arithmetic_profiles": members,
                "multiplicity": len(members),
                "phase_zero_b0_count": int(histogram["0"].get(0, 0)),
                "direction_assignment_count_after_normalization": raw_assignments,
                "reflection_fixed_assignment_count": reflection_fixed,
                "direction_assignment_orbits_under_residual_reflection": (
                    reflection_orbits
                ),
            }
        )

    multiplicities = Counter(int(row["multiplicity"]) for row in signatures)
    if (
        len(signatures) != EXPECTED_SIGNATURE_COUNT
        or multiplicities != Counter({1: 175, 2: 9})
        or sum(
            int(row["direction_assignment_count_after_normalization"])
            for row in signatures
        )
        != 1_971_382
        or sum(
            int(row["direction_assignment_orbits_under_residual_reflection"])
            for row in signatures
        )
        != 985_730
        or sum(
            int(row["reflection_fixed_assignment_count"])
            for row in signatures
        )
        != 78
    ):
        raise ArithmeticError("the p=17 slack-twenty signature ledger changed")

    return {
        "experiment": "p17_slack20_boundary_signature_manifest",
        "p": P,
        "boundary_size": SIZE,
        "pair_slack": SLACK,
        "profile_count": len(profiles),
        "signature_count": len(signatures),
        "profile_census_indices": [
            int(profile["census_index"]) for profile in profiles
        ],
        "multiplicity_histogram": {
            str(value): count for value, count in sorted(multiplicities.items())
        },
        "normalization": {
            "mode": "phase-zero-b0-pair",
            "selected_points": [0, 1],
            "canonical_direction": [0, 1],
            "phase": 0,
            "b": 0,
            "c_H": -1,
            "lossless_for_every_profile": True,
        },
        "residual_direction_reflection": {
            "phase_zero_fixed": [17],
            "phase_zero_pairs": [[1, 16], [2, 15], [4, 13], [6, 11]],
            "phase_one_fixed": [0],
            "phase_one_pairs": [[3, 14], [5, 12], [7, 10], [8, 9]],
            "raw_assignment_count": 1_971_382,
            "orbit_count": 985_730,
            "fixed_assignment_count": 78,
        },
        "signatures": signatures,
    }


def resolve_signature(
    *, signature_index: int | None = None, profile_index: int | None = None
) -> dict[str, Any]:
    """Resolve exactly one safe CLI identifier to its canonical signature."""
    if (signature_index is None) == (profile_index is None):
        raise ValueError("choose exactly one signature_index or profile_index")
    manifest = signature_manifest()
    signatures = manifest["signatures"]
    if signature_index is not None:
        if not 0 <= signature_index < len(signatures):
            raise ValueError("signature index must lie in 0..68")
        return signatures[signature_index]

    assert profile_index is not None
    if not 0 <= profile_index < EXPECTED_PROFILE_COUNT:
        raise ValueError("profile index must lie in 0..77")
    matches = [
        row for row in signatures if profile_index in row["profile_indices"]
    ]
    if len(matches) != 1:
        raise ArithmeticError("profile-to-signature map is not functional")
    return matches[0]


def radon_geometry() -> dict[str, Any]:
    """Return and audit the canonical p=17 direction/phase geometry."""
    directions = projective_directions(P)
    records = []
    normalized_indices = []
    phase_counts: Counter[int] = Counter()
    for index, direction in enumerate(directions):
        eps, labels = field_direction_data(P, direction)
        phase = 0 if eps == 1 else 1  # normalized c_H=-1 convention
        phase_counts[phase] += 1
        records.append((direction, int(eps), phase, labels))
        if labels[0] == labels[1]:
            normalized_indices.append(index)
    if (
        len(directions) != P + 1
        or phase_counts != Counter({0: PHASE_SIZE, 1: PHASE_SIZE})
        or normalized_indices != [17]
        or directions[17] != (0, 1)
        or records[17][1] != 1
    ):
        raise ArithmeticError("canonical p=17 Radon geometry changed")
    return {
        "directions": directions,
        "records": records,
        "normalized_index": 17,
        "point_variables": P * P,
        "line_parity_variables": P * (P + 1),
        "native_xor_constraints": P * (P + 1) + P * P,
    }


def audit_boundary(
    chosen: Iterable[int], expected: dict[str, dict[int, int]]
) -> dict[str, Any]:
    """Recompute a candidate boundary and all Radon/inverse identities."""
    points = tuple(sorted(int(value) for value in chosen))
    point_set = set(points)
    geometry = radon_geometry()
    observed: dict[str, Counter[int]] = {"0": Counter(), "1": Counter()}
    direction_rows = []
    parity_rows: list[list[int]] = []
    valid_points = (
        len(points) == SIZE
        and len(point_set) == SIZE
        and min(points, default=0) >= 0
        and max(points, default=0) < P * P
        and {0, 1} <= point_set
    )

    for index, (direction, eps, phase, labels) in enumerate(geometry["records"]):
        counts = [0] * P
        for point in points:
            if 0 <= point < P * P:
                counts[labels[point]] += 1
        parity = [count & 1 for count in counts]
        parity_rows.append(parity)
        b = sum(parity)
        observed[str(phase)][b] += 1
        direction_rows.append(
            {
                "direction_index": index,
                "direction": list(direction),
                "eps": eps,
                "phase": phase,
                "b": b,
                "odd_fibres": [i for i, value in enumerate(parity) if value],
            }
        )

    inverse_valid = True
    for point in range(P * P):
        incident_sum = 0
        for index, (_direction, _eps, _phase, labels) in enumerate(
            geometry["records"]
        ):
            incident_sum ^= parity_rows[index][labels[point]]
        if incident_sum != int(point in point_set):
            inverse_valid = False
            break

    expected_counter = {
        str(phase): Counter(
            {int(b): int(count) for b, count in expected[str(phase)].items()}
        )
        for phase in (0, 1)
    }
    normalized_valid = direction_rows[geometry["normalized_index"]]["b"] == 0
    histogram_valid = observed == expected_counter
    valid = valid_points and inverse_valid and normalized_valid and histogram_valid
    return {
        "valid": valid,
        "point_set_valid": valid_points,
        "inverse_radon_valid": inverse_valid,
        "normalization_valid": normalized_valid,
        "phase_histograms_valid": histogram_valid,
        "observed_phase_profiles_b": {
            phase: dict(sorted(counter.items()))
            for phase, counter in observed.items()
        },
        "direction_rows": direction_rows,
    }


def solve_signature(
    signature: dict[str, Any], *, seconds: float, threads: int
) -> dict[str, Any]:
    """Build and solve one exact canonical boundary-signature instance."""
    if seconds <= 0:
        raise ValueError("seconds must be positive")
    if threads <= 0:
        raise ValueError("threads must be positive")
    try:
        from pycryptosat import Solver
    except ModuleNotFoundError as exc:  # pragma: no cover - deployment guard
        raise RuntimeError("pycryptosat is required for exact solving") from exc

    started = time.time()
    geometry = radon_geometry()
    records = geometry["records"]
    expected = {
        str(phase): {
            int(b): int(count)
            for b, count in signature["phase_profiles_b"][str(phase)].items()
        }
        for phase in (0, 1)
    }

    next_id = 0

    def new_var() -> int:
        nonlocal next_id
        next_id += 1
        return next_id

    point = [new_var() for _ in range(P * P)]
    parity = [[new_var() for _ in range(P)] for _ in range(P + 1)]
    selectors: list[dict[int, int]] = []
    for _direction, _eps, phase, _labels in records:
        selectors.append({b: new_var() for b in expected[str(phase)]})
    semantic_variables = next_id

    solver = Solver(verbose=0, threads=int(threads))
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

    def add_exact(
        literals: Iterable[int], bound: int, guard: int | None = None
    ) -> None:
        """Guarded exact count; clauses are active when ``guard`` is false."""
        nonlocal cardinality_count
        terms = [int(literal) for literal in literals]
        target = int(bound)
        if target > len(terms) // 2:
            terms = [-literal for literal in terms]
            target = len(terms) - target
        if not 0 <= target <= len(terms):
            add_clauses([[]] if guard is None else [[int(guard)]])
            cardinality_count += 1
            return
        overflow = target + 1
        states = [
            [new_var() for _ in range(overflow + 1)]
            for _ in range(len(terms) + 1)
        ]

        def gated(clause: list[int]) -> list[int]:
            return clause if guard is None else [int(guard), *clause]

        clauses: list[list[int]] = []
        for state_row in states:
            clauses.append(gated(list(state_row)))
            for first in range(len(state_row)):
                for second in range(first + 1, len(state_row)):
                    clauses.append(gated([-state_row[first], -state_row[second]]))
        clauses.append(gated([states[0][0]]))
        clauses.append(gated([states[-1][target]]))
        for index, literal in enumerate(terms, start=1):
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
    add_unit(point[1])

    incident_parities: list[list[int]] = [[] for _ in point]
    for index, (_direction, _eps, _phase, labels) in enumerate(records):
        fibres = []
        for fibre in range(P):
            fibre_points = [
                point[value] for value, label in enumerate(labels) if label == fibre
            ]
            fibres.append(fibre_points)
            solver.add_xor_clause([*fibre_points, parity[index][fibre]], False)
            xor_count += 1
        for value, label in enumerate(labels):
            incident_parities[value].append(parity[index][label])

        choices = list(selectors[index].values())
        add_clauses([choices])
        for first in range(len(choices)):
            for second in range(first + 1, len(choices)):
                add_clauses([[-choices[first], -choices[second]]])
        for b, selector in selectors[index].items():
            add_exact(parity[index], b, guard=-selector)
            if b == SIZE:
                for fibre_points in fibres:
                    for first in range(len(fibre_points)):
                        for second in range(first + 1, len(fibre_points)):
                            add_clauses(
                                [[-selector, -fibre_points[first], -fibre_points[second]]]
                            )

    for value, rows in enumerate(incident_parities):
        if len(rows) != P + 1:
            raise ArithmeticError("affine point degree changed")
        solver.add_xor_clause([point[value], *rows], False)
        xor_count += 1

    for phase in (0, 1):
        phase_indices = [
            index for index, row in enumerate(records) if int(row[2]) == phase
        ]
        if len(phase_indices) != PHASE_SIZE:
            raise ArithmeticError("p=17 phase size changed")
        for b, target in expected[str(phase)].items():
            add_exact([selectors[index][b] for index in phase_indices], target)

    normalized_index = int(geometry["normalized_index"])
    add_unit(selectors[normalized_index][0])
    if xor_count != 595:
        raise ArithmeticError("p=17 native Radon equation count changed")

    build_seconds = time.time() - started
    satisfiable, assignment = solver.solve(time_limit=float(seconds))
    status = (
        "SATISFIABLE"
        if satisfiable is True
        else "UNSATISFIABLE"
        if satisfiable is False
        else "UNKNOWN"
    )
    result: dict[str, Any] = {
        "experiment": "p17_slack20_boundary_cryptominisat",
        "p": P,
        "boundary_size": SIZE,
        "pair_slack": SLACK,
        "signature_index": int(signature["signature_index"]),
        "profile_indices": list(signature["profile_indices"]),
        "census_indices": list(signature["census_indices"]),
        "phase_profiles_b": expected,
        "normalization": {
            "mode": "phase-zero-b0-pair",
            "selected_points": [0, 1],
            "direction_index": normalized_index,
            "direction": [0, 1],
            "phase": 0,
            "b": 0,
            "c_H": -1,
            "lossless": True,
        },
        "sign_transfer": {
            "nonsquare_dilation_flips_eps_and_c_H_together": True,
            "both_c_H_signs_covered": True,
        },
        "solver": "cryptominisat-native-xor",
        "solver_status": status,
        "feasible_boundary_signature": satisfiable is True,
        "finite_infeasibility_only": satisfiable is False,
        "semantic_variables": semantic_variables,
        "total_variables": next_id,
        "clauses": clause_count,
        "native_xor_constraints": xor_count,
        "cardinality_constraints": cardinality_count,
        "threads": int(threads),
        "build_seconds": build_seconds,
        "solve_seconds": time.time() - started - build_seconds,
        "elapsed_seconds": time.time() - started,
    }

    if satisfiable is True:
        chosen = [value for value, literal in enumerate(point) if assignment[literal]]
        audit = audit_boundary(chosen, expected)
        if not audit["valid"]:
            raise AssertionError("CryptoMiniSat p=17 witness failed independent audit")

        # Also compare every semantic parity/selector literal with the
        # independently reconstructed direction rows.
        semantic_valid = True
        for index, direction_row in enumerate(audit["direction_rows"]):
            odd_fibres = set(direction_row["odd_fibres"])
            for fibre in range(P):
                if bool(assignment[parity[index][fibre]]) != (fibre in odd_fibres):
                    semantic_valid = False
            observed_b = int(direction_row["b"])
            for b, selector in selectors[index].items():
                if bool(assignment[selector]) != (b == observed_b):
                    semantic_valid = False
        if not semantic_valid:
            raise AssertionError("CryptoMiniSat semantic assignment failed audit")

        result["boundary"] = chosen
        result["boundary_coordinates"] = [[value % P, value // P] for value in chosen]
        result["witness_audit"] = audit
        result["semantic_assignment_audit_valid"] = True
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--signature", type=int, help="canonical signature 0..68")
    group.add_argument("--profile", type=int, help="slack-20 profile 0..77")
    group.add_argument(
        "--list-signatures",
        action="store_true",
        help="emit the corrected 193-to-184 manifest without invoking a solver",
    )
    parser.add_argument("--seconds", type=float, default=300.0)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--output", type=Path)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.list_signatures:
        result = signature_manifest()
    else:
        signature = resolve_signature(
            signature_index=args.signature,
            profile_index=args.profile,
        )
        result = solve_signature(
            signature,
            seconds=float(args.seconds),
            threads=int(args.threads),
        )
    if args.output is None:
        print(json.dumps(result, indent=2, sort_keys=True))
    else:
        atomic_write(args.output, result)


if __name__ == "__main__":
    main()
