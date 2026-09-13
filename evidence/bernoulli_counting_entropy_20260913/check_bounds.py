"""Exact real-algebra checks for the counting-entropy cutoff arithmetic.

Each query seeks a violation of an inequality used in Section 2 of the
proof note. UNSAT excludes every real violation under its assumptions;
this checks neither Bernstein's inequality nor Shannon's chain rule.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import socket
import time

import z3


def check(name: str, assumptions: list, violation, conclusion: str) -> dict:
    solver = z3.SolverFor("QF_NRA")
    solver.set(timeout=10000)
    solver.add(*assumptions, violation)
    started = time.monotonic()
    outcome = solver.check()
    result = {
        "name": name,
        "outcome": str(outcome),
        "conclusion": conclusion,
        "elapsed_seconds": time.monotonic() - started,
        "query": solver.sexpr(),
    }
    if outcome == z3.sat:
        result["counterexample"] = str(solver.model())
    elif outcome == z3.unknown:
        result["reason_unknown"] = solver.reason_unknown()
    return result


def main() -> int:
    delta, r, u = z3.Reals("delta r u")
    # u=sqrt(n); r is the normalized square-root error. The proof bounds
    # r^2 by delta^2/8 and requires delta*sqrt(n)>=16/3.
    noise = check(
        "normalized_noise_cutoff",
        [delta > 0, u > 0, r >= 0, 8 * r * r <= delta * delta,
         3 * delta * u >= 16],
        3 * r * u + 8 > 3 * delta * u,
        "r + 8/(3u) <= delta",
    )
    n, edges, log_two, entropy = z3.Reals("n edges log_two entropy")
    # Positive denominators are cleared explicitly. log_two and entropy
    # are arbitrary positive reals here, so no floating logarithm enters.
    entropy_loss = check(
        "entropy_loss_cutoff",
        [n >= 2, edges > 0, 2 * edges == n * (n - 1),
         log_two > 0, entropy > 0, n * entropy >= 6 * log_two],
        2 * log_two * (edges + n) > entropy * n * edges,
        "log_two/n + log_two/edges <= entropy/2",
    )
    results = [noise, entropy_loss]
    passed = all(item["outcome"] == "unsat" for item in results)
    print(json.dumps({
        "classification": "symbolic_scalar_corroboration_only",
        "host": socket.gethostname(),
        "solver": z3.get_version_string(),
        "logic": "QF_NRA",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "checks": results,
        "passed": passed,
    }, indent=2, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
