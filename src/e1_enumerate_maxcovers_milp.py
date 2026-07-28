#!/usr/bin/env python3
"""
Find diverse Max-cover perfect matchings at p=5 via MILP with random objectives.
Verify each: maxR, phi, clique-flip, op-norm.
Writes evidence/e1_milp_maxcovers_verified.json
"""
from __future__ import annotations

import json
import os
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import (  # noqa: E402
    maxR_matching_level,
    paley_conference_prime_power,
    phi_mitm,
)


def main():
    p = 5
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    m = n // 2
    Phi = 0.5 * n * p
    Max = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    Max_full = np.vstack([Max, -Max])
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)

    Chi = np.zeros((len(Max_full), E))
    for yi, y in enumerate(Max_full):
        for k, (i, j) in enumerate(edges):
            Chi[yi, k] = C[i, j] * y[i] * y[j]

    A_eq = np.zeros((n, E))
    for k, (i, j) in enumerate(edges):
        A_eq[i, k] = 1
        A_eq[j, k] = 1

    cons_list = [
        LinearConstraint(A_eq, np.ones(n), np.ones(n)),
        LinearConstraint(Chi, np.ones(len(Max_full)), np.full(len(Max_full), np.inf)),
    ]
    bounds = Bounds(0, 1)

    rng = np.random.default_rng(42)
    found = {}
    for trial in range(80):
        c = rng.normal(size=E)
        res = milp(
            c,
            integrality=np.ones(E),
            bounds=bounds,
            constraints=cons_list,
            options={"time_limit": 20, "disp": False},
        )
        if not res.success:
            print(f"trial {trial}: fail {res.message}", flush=True)
            continue
        x = np.rint(res.x).astype(int)
        M = [edges[k] for k in range(E) if x[k] == 1]
        assert len(M) == m
        key = tuple(sorted(tuple(sorted(e)) for e in M))
        if key in found:
            print(f"trial {trial}: rediscover", flush=True)
            continue
        mr = maxR_matching_level(C, M, p)
        A = C.copy()
        for i, j in M:
            A[i, j] *= -1
            A[j, i] *= -1
        ph = float(phi_mitm(A))
        op = float(np.linalg.norm(A, ord=2))
        rec = {
            "matching": [list(e) for e in M],
            "maxR": int(mr),
            "phi": ph,
            "op": op,
            "undercut": ph < Phi - 0.5,
            "criterion": mr >= 60,
        }
        found[key] = rec
        print(
            f"trial {trial}: NEW cover #{len(found)} maxR={mr} phi={ph} op={op:.4f}",
            flush=True,
        )

    out = {
        "n_unique": len(found),
        "all_no_undercut": all(not r["undercut"] for r in found.values()),
        "all_criterion": all(r["criterion"] for r in found.values()),
        "covers": list(found.values()),
        "status": "OPEN — MILP sample of Max-covers, not exhaustive",
    }
    path = ROOT / "evidence" / "e1_milp_maxcovers_verified.json"
    path.write_text(json.dumps(out, indent=2))
    print("wrote", path, "n=", len(found), flush=True)
    print(
        json.dumps({k: out[k] for k in out if k != "covers"}, indent=2),
        flush=True,
    )


if __name__ == "__main__":
    main()
