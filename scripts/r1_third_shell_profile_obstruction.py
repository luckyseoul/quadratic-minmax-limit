#!/usr/bin/env python3
"""Exact scout for the only unresolved third-shell profile equality case.

At scaled norm ``2(p+1)``, the difficult zero-common-sum case has one active
profile.  Up to sign, its positive multiset has size ``m=(p-1)/2`` with one
root repeated, its negative multiset has ``m`` distinct roots, and the two
multisets have equal power sums through degree ``m-1``.  This script
enumerates those finite configurations exactly.  Absence is only a finite-p
audit; the intended theorem still needs the uniform polynomial argument.
"""
from __future__ import annotations

import argparse
import itertools
import json
import time
from pathlib import Path


def obstruction_profiles(p: int, stop_after: int = 10) -> dict:
    m = (p - 1) // 2
    field = tuple(range(p))
    checked = 0
    examples = []
    started = time.perf_counter()
    powers = [[pow(x, d, p) for x in field] for d in range(1, m)]

    for repeated in field:
        available = tuple(x for x in field if x != repeated)
        for positive_other in itertools.combinations(available, m - 2):
            positive_set = {repeated, *positive_other}
            remaining = tuple(x for x in field if x not in positive_set)
            target = tuple(
                (2 * powers[d - 1][repeated]
                 + sum(powers[d - 1][x] for x in positive_other))
                % p
                for d in range(1, m)
            )
            for negative in itertools.combinations(remaining, m):
                checked += 1
                signature = tuple(
                    sum(powers[d - 1][x] for x in negative) % p
                    for d in range(1, m)
                )
                if signature != target:
                    continue
                examples.append(
                    {
                        "repeated_positive_root": repeated,
                        "other_positive_roots": list(positive_other),
                        "negative_roots": list(negative),
                    }
                )
                if len(examples) >= stop_after:
                    return {
                        "p": p,
                        "half_mass": m,
                        "checked": checked,
                        "complete": False,
                        "examples": examples,
                        "seconds": time.perf_counter() - started,
                    }
    return {
        "p": p,
        "half_mass": m,
        "checked": checked,
        "complete": True,
        "examples": examples,
        "seconds": time.perf_counter() - started,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", default="11,13,17")
    parser.add_argument("--stop-after", type=int, default=10)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    records = []
    for p in (int(value) for value in args.primes.split(",")):
        record = obstruction_profiles(p, args.stop_after)
        records.append(record)
        print(json.dumps(record), flush=True)
    result = {
        "experiment": "r1_third_shell_profile_obstruction",
        "status": "finite_exact_audit_not_uniform_proof",
        "records": records,
    }
    if args.output is not None:
        args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
