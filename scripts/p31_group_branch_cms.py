#!/usr/bin/env python3
"""Exact CryptoMiniSat decision model for the p=31 even grouped lemma.

Native XOR clauses express the 480 affine-block parities.  A fixed even
point weight ``s`` is a counterexample exactly when at most ``31-s`` of the
32 direction groups are active.  Any counterexample has at least three silent
directions; PGL(2,31) is triply transitive, so directions 0,1,31 are fixed
silent without loss.  SAT witnesses are reconstructed independently.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from pysat.card import CardEnc, EncType


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_symmetric_halved_mod2 import (  # noqa: E402
    _antipodal_classes,
    _directions,
    _evaluate,
)


P = 31
H = 15
D = 32
N = 480


def block_layout() -> tuple[tuple[tuple[int, ...], ...], ...]:
    points = _antipodal_classes(P)
    squares = tuple(sorted({x * x % P for x in range(1, P)}))
    groups = tuple(
        tuple(
            tuple(
                i + 1
                for i, point in enumerate(points)
                if _evaluate(P, functional, point) ** 2 % P == square
            )
            for square in squares
        )
        for functional in _directions(P)
    )
    if any(len(block) != P for group in groups for block in group):
        raise ArithmeticError("p=31 affine-block layout changed")
    return groups


def verify(point_variables: set[int]) -> dict[str, object]:
    groups = block_layout()
    support = tuple(i - 1 for i in sorted(point_variables) if 1 <= i <= N)
    selected = set(i + 1 for i in support)
    active = tuple(
        direction
        for direction, blocks in enumerate(groups)
        if any(sum(variable in selected for variable in block) & 1 for block in blocks)
    )
    silent = tuple(i for i in range(D) if i not in active)
    return {
        "support_indices": list(support),
        "support_weight": len(support),
        "active_direction_indices": list(active),
        "active_direction_count": len(active),
        "silent_direction_indices": list(silent),
        "silent_direction_count": len(silent),
        "group_branch_weight": len(support) + len(active),
        "counterexample": bool(support) and len(support) + len(active) <= P,
    }


def build_instance(
    point_weight: int,
    required_point: int | None = None,
    fourth_silent: int | None = None,
    fifth_silent: int | None = None,
) -> tuple[list[tuple[int, ...]], list[tuple[int, ...]], int]:
    groups = block_layout()
    parity_start = N + 1
    active_start = parity_start + D * H
    active_variables = [active_start + direction for direction in range(D)]
    top_id = active_start + D - 1
    clauses: list[tuple[int, ...]] = []
    xors: list[tuple[int, ...]] = []

    for direction, blocks in enumerate(groups):
        group_parities: list[int] = []
        for block_number, block in enumerate(blocks):
            parity = parity_start + direction * H + block_number
            group_parities.append(parity)
            # CryptoMiniSat XOR lines have right-hand side true.  Negating y
            # therefore expresses XOR(block points) == y.
            xors.append((*block, -parity))
            clauses.append((-parity, active_variables[direction]))
        clauses.append((-active_variables[direction], *group_parities))

    for direction in (0, 1, 31):
        clauses.append((-active_variables[direction],))
    if fourth_silent is not None:
        if fourth_silent in (0, 1, 31) or not 0 <= fourth_silent < D:
            raise ValueError("fourth silent direction must lie in [0,31] minus {0,1,31}")
        clauses.append((-active_variables[fourth_silent],))
    if fifth_silent is not None:
        if fifth_silent in (0, 1, 31, fourth_silent) or not 0 <= fifth_silent < D:
            raise ValueError("fifth silent direction must be new and lie in [0,31]")
        clauses.append((-active_variables[fifth_silent],))
    if required_point is not None:
        if not 0 <= required_point < N:
            raise ValueError("required point must lie in [0,479]")
        clauses.append((required_point + 1,))

    exact_points = CardEnc.equals(
        lits=list(range(1, N + 1)),
        bound=point_weight,
        top_id=top_id,
        encoding=EncType.seqcounter,
    )
    clauses.extend(tuple(clause) for clause in exact_points.clauses)
    top_id = exact_points.nv
    active_cap = CardEnc.atmost(
        lits=active_variables,
        bound=P - point_weight,
        top_id=top_id,
        encoding=EncType.seqcounter,
    )
    clauses.extend(tuple(clause) for clause in active_cap.clauses)
    top_id = active_cap.nv
    return clauses, xors, top_id


def solve(
    point_weight: int,
    seconds: int,
    threads: int,
    seed: int,
    required_point: int | None = None,
    fourth_silent: int | None = None,
    fifth_silent: int | None = None,
) -> dict[str, object]:
    clauses, xors, variable_count = build_instance(
        point_weight, required_point, fourth_silent, fifth_silent
    )
    with tempfile.NamedTemporaryFile("w", suffix=".cnf") as handle:
        handle.write(f"p cnf {variable_count} {len(clauses) + len(xors)}\n")
        for clause in clauses:
            handle.write(" ".join(map(str, clause)) + " 0\n")
        for xor in xors:
            handle.write("x" + " ".join(map(str, xor)) + " 0\n")
        handle.flush()
        started = time.monotonic()
        completed = subprocess.run(
            [
                "cryptominisat5",
                "--verb", "1",
                "--maxtime", str(seconds),
                "--threads", str(threads),
                "--random", str(seed),
                handle.name,
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        elapsed = time.monotonic() - started

    output = completed.stdout + completed.stderr
    if "s UNSATISFIABLE" in output:
        status = "UNSATISFIABLE"
    elif "s SATISFIABLE" in output:
        status = "SATISFIABLE"
    else:
        status = "UNKNOWN"
    out: dict[str, object] = {
        "p": P,
        "point_weight": point_weight,
        "decision_problem": "wt(f)+active_direction_count<=31",
        "scope": "even point weight; odd weight is symbolically closed",
        "symmetry_break": "directions 0,1,31 silent by PGL(2,31) triple transitivity",
        "solver": "CryptoMiniSat 5 with native XOR clauses",
        "status": status,
        "solver_returncode": completed.returncode,
        "wall_seconds": elapsed,
        "time_limit_seconds": seconds,
        "threads": threads,
        "seed": seed,
        "required_point_index": required_point,
        "fourth_silent_direction_index": fourth_silent,
        "fifth_silent_direction_index": fifth_silent,
        "variables": variable_count,
        "cnf_clauses": len(clauses),
        "xor_clauses": len(xors),
        "proved_no_counterexample_at_weight": status == "UNSATISFIABLE",
    }
    if status == "SATISFIABLE":
        positive: set[int] = set()
        for line in output.splitlines():
            if line.startswith("v "):
                positive.update(int(value) for value in line[2:].split() if int(value) > 0)
        witness = verify(positive)
        if witness["support_weight"] != point_weight or not witness["counterexample"]:
            raise ArithmeticError("CryptoMiniSat witness failed independent verification")
        out["witness"] = witness
    elif status == "UNKNOWN":
        out["solver_log_tail"] = output.splitlines()[-12:]
    # Keep only compact solver counters instead of embedding the full log.
    for line in output.splitlines():
        if line.startswith("c conflicts") or line.startswith("c decisions"):
            key, _, value = line[2:].partition(":")
            out[key.strip().replace(" ", "_")] = value.strip()
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--point-weight", type=int, required=True)
    parser.add_argument("--time-limit", type=int, default=300)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--seed", type=int, default=310037)
    parser.add_argument("--required-point", type=int)
    parser.add_argument("--fourth-silent", type=int)
    parser.add_argument("--fifth-silent", type=int)
    args = parser.parse_args()
    if args.point_weight < 2 or args.point_weight > 30 or args.point_weight % 2:
        parser.error("point weight must be even and lie in [2,30]")
    if args.time_limit <= 0 or args.threads <= 0:
        parser.error("time limit and threads must be positive")
    print(
        json.dumps(
            solve(
                args.point_weight,
                args.time_limit,
                args.threads,
                args.seed,
                args.required_point,
                args.fourth_silent,
                args.fifth_silent,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
