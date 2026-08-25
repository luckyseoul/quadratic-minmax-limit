#!/usr/bin/env python3
"""Find one-active zero-profile vectors at a prescribed exact energy.

This is a finite scout for the first unclassified dual shell.  A profile
``v=(v_s)`` must have integer sum zero and power moments zero modulo ``p``
through degree ``(p-3)/2``.  The CP-SAT model uses one-hot bounded integer
values, so both the energy and all moment congruences are exact.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from ortools.sat.python import cp_model


def find_profile(
    p: int, energy: int, seconds: float, two_double_case: bool = False
) -> dict:
    m = (p - 1) // 2
    bound = 2 if two_double_case else math.isqrt(energy)
    values = tuple(range(-bound, bound + 1))
    model = cp_model.CpModel()
    choose = {
        (s, value): model.new_bool_var(f"x_{s}_{value}")
        for s in range(p)
        for value in values
    }
    for s in range(p):
        model.add_exactly_one(choose[s, value] for value in values)
    model.add(
        sum(value * choose[s, value] for s in range(p) for value in values)
        == 0
    )
    model.add(
        sum(
            value * value * choose[s, value]
            for s in range(p)
            for value in values
        )
        == energy
    )
    if two_double_case:
        # At energy p+3 and mass m, exactly two entries have magnitude two.
        # Affine invariance moves them to 0 and 1; global sign fixes v_0=2.
        model.add(choose[0, 2] == 1)
        model.add(choose[1, -2] + choose[1, 2] == 1)
        model.add(
            sum(choose[s, -2] + choose[s, 2] for s in range(p)) == 2
        )
    else:
        # Translation and global-sign symmetry: move a nonzero entry to zero
        # and orient it positively.
        model.add(sum(choose[0, value] for value in values if value > 0) == 1)
    for degree in range(1, m):
        coefficients = [pow(s, degree, p) for s in range(p)]
        raw_bound = bound * sum(coefficients)
        quotient = model.new_int_var(
            -(raw_bound // p + 1), raw_bound // p + 1, f"q_{degree}"
        )
        model.add(
            sum(
                coefficients[s] * value * choose[s, value]
                for s in range(p)
                for value in values
            )
            == p * quotient
        )
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = seconds
    solver.parameters.num_search_workers = 1
    status = solver.solve(model)
    profile = None
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        profile = [
            next(value for value in values if solver.value(choose[s, value]))
            for s in range(p)
        ]
    return {
        "p": p,
        "energy": energy,
        "two_double_case": two_double_case,
        "status": solver.status_name(status),
        "profile": profile,
        "wall_seconds": solver.wall_time,
        "branches": solver.num_branches,
        "conflicts": solver.num_conflicts,
        "checks": None
        if profile is None
        else {
            "sum": sum(profile),
            "energy": sum(value * value for value in profile),
            "moments_mod_p": [
                sum(value * pow(s, degree, p) for s, value in enumerate(profile))
                % p
                for degree in range(1, m)
            ],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", default="11,13,17")
    parser.add_argument("--seconds", type=float, default=120.0)
    parser.add_argument("--two-double-case", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rows = []
    for p in (int(item) for item in args.primes.split(",")):
        row = find_profile(p, p + 3, args.seconds, args.two_double_case)
        rows.append(row)
        print(json.dumps(row), flush=True)
    result = {
        "experiment": "one_active_profile_at_energy_p_plus_3",
        "finite_scout_only": True,
        "rows": rows,
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
