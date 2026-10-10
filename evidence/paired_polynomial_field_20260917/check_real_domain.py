"""Check the full paired-root statement over real variables, not signings.

Separate from the SymPy polynomial identity: ask for an actual violating
real tuple, including the nonnegative square-root branches. Two bounded
queries run on two CPU workers. UNKNOWN is not a passing result.
"""
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import platform


def query(kind):
    import z3

    a, b, q, y, r, u, v, z = z3.Reals("a b q y r u v z")
    plus = a+b*y*y+2*b*r
    minus = a+b*y*y-2*b*r
    solver = z3.SolverFor("QF_NRA")
    solver.set(timeout=60000)
    solver.add(b >= 0, q >= 0, y >= 0, a >= b*q, r*r <= q*y*y)
    if kind == "both_radicands_nonnegative":
        solver.add(z3.Or(plus < 0, minus < 0))
    elif kind == "paired_root_lower":
        solver.add(u >= 0, v >= 0, z >= 0,
                   u*u == plus, v*v == minus, z*z == a, u+v < 2*z)
    else:
        raise ValueError(kind)
    answer = solver.check()
    result = {"query": kind, "result": str(answer)}
    if answer == z3.sat:
        result["counterexample"] = str(solver.model())
    elif answer == z3.unknown:
        result["reason"] = solver.reason_unknown()
    return result


def main():
    import z3

    names = ("both_radicands_nonnegative", "paired_root_lower")
    with ProcessPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(query, names))
    result = {
        "classification": "real-algebra decision checks; analytic imports not formalized",
        "python": platform.python_version(),
        "z3": z3.get_version_string(),
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "workers": 2,
        "timeout_ms_per_query": 60000,
        "queries": results,
        "all_passed": all(row["result"] == "unsat" for row in results),
        "signings_enumerated": 0,
        "convergence_proved": False,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    if not result["all_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
