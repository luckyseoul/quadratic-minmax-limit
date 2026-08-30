#!/usr/bin/env python3
"""Exact quantized type minima for infinity plus p finite boundary points."""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15669 import full_symbolic_floor  # noqa: E402
from e1_gmin_m4_prop15642 import nonbaseline_scaled_cost_floor  # noqa: E402
from e1_gmin_m4_prop15723 import floor_excess_admissible  # noqa: E402


def type_minimum(p: int, phase: int) -> dict[str, object]:
    period = p + 1
    m = period // 2
    lift_floor = nonbaseline_scaled_cost_floor(p)
    best: tuple[int, int, tuple[int, ...], tuple[int, ...]] | None = None
    rows = []
    for u in range(m):
        residue = 2 * u
        target = m - u
        options = []
        for b in range(1, p + 1, 2):
            floor = full_symbolic_floor(p, b, phase)
            for k in range(target + 1):
                excess = residue + period * k - floor
                if floor_excess_admissible(p, b, phase, excess):
                    options.append((k, p - b, b))
        states: dict[int, tuple[int, tuple[int, ...], tuple[int, ...]]] = {
            0: (0, (), ())
        }
        for _ in range(m):
            next_states: dict[
                int, tuple[int, tuple[int, ...], tuple[int, ...]]
            ] = {}
            for quotient, (deficit, profile, quotients) in states.items():
                for k, added_deficit, b in options:
                    new_quotient = quotient + k
                    if new_quotient > target:
                        continue
                    candidate = (
                        deficit + added_deficit,
                        profile + (b,),
                        quotients + (k,),
                    )
                    old = next_states.get(new_quotient)
                    if old is None or candidate[0] < old[0]:
                        next_states[new_quotient] = candidate
            states = next_states
        if target not in states:
            continue
        deficit, profile, quotients = states[target]
        row = {
            "u": u,
            "residue": residue,
            "deficit": deficit,
            "profile": dict(sorted(Counter(profile).items())),
            "quotients": dict(sorted(Counter(quotients).items())),
        }
        rows.append(row)
        candidate = (deficit, u, profile, quotients)
        if best is None or candidate[0] < best[0]:
            best = candidate
    if best is None:
        raise ArithmeticError("no quantized type state")
    return {
        "p": p,
        "phase": phase,
        "nonzero_lift_floor": lift_floor,
        "minimum_deficit": best[0],
        "winning_u": best[1],
        "profile": dict(sorted(Counter(best[2]).items())),
        "quotients": dict(sorted(Counter(best[3]).items())),
        "residue_minima": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", default="17,19,23,29,31,37,41,43")
    args = parser.parse_args()
    for p in (int(value) for value in args.primes.split(",")):
        print(f"p={p} pair_budget={p * (p - 1)}")
        for phase in (0, 1):
            row = type_minimum(p, phase)
            print(
                f"  phase={phase} D={row['minimum_deficit']} "
                f"u={row['winning_u']} profile={row['profile']}"
            )
            print(
                "    residues="
                + str(
                    [
                        (entry["u"], entry["deficit"], entry["profile"])
                        for entry in row["residue_minima"]
                    ]
                )
            )


if __name__ == "__main__":
    main()
