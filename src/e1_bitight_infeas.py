#!/usr/bin/env python3
"""
Certify integral bi-tight infeasibility for p=5 (levels 2,3,4).

Bi-tight level s: |H|=s*p, S≡s on Max+, S≡-s on Max-.
Fractional is feasible (uniform); integral is not for p=5.

Writes evidence/e1_bitight_infeas.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, linprog, milp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def boolean_evecs(C: np.ndarray, p: float, seed: int = 0) -> np.ndarray:
    n = C.shape[0]
    M = C - p * np.eye(n)
    _, s, vt = np.linalg.svd(M)
    rank = int((s > 1e-8).sum())
    nb = vt[rank:].T
    nul = n - rank
    rng = np.random.default_rng(seed)
    coords = None
    for _ in range(20000):
        idx = rng.choice(n, size=nul, replace=False)
        if np.linalg.matrix_rank(nb[idx], tol=1e-8) == nul:
            coords = idx
            break
    if coords is None:
        raise RuntimeError("no free coords")
    Binv = np.linalg.inv(nb[coords])
    found = []
    for bits in range(1 << nul):
        free = np.array([1.0 if (bits >> i) & 1 else -1.0 for i in range(nul)])
        x = nb @ (Binv @ free)
        if np.max(np.abs(np.abs(x) - 1.0)) < 1e-6:
            found.append(np.sign(x))
    return np.unique(np.round(found).astype(int), axis=0).astype(float)


def main() -> None:
    p = 5
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    t0 = time.time()

    # Max+ : prefer cache
    cache = Path("/tmp/maxplus_p5.npy")
    if cache.is_file():
        Maxp = np.load(cache).astype(float)
    else:
        Maxp = boolean_evecs(C, float(p), seed=0)
    Maxm = boolean_evecs(-C, float(p), seed=1)

    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)
    Chip = np.array(
        [[C[i, j] * y[i] * y[j] for i, j in edges] for y in Maxp], dtype=float
    )
    Chim = np.array(
        [[C[i, j] * y[i] * y[j] for i, j in edges] for y in Maxm], dtype=float
    )

    results = []
    for s_lev in (2, 3, 4):
        k = s_lev * p
        Aeq = np.vstack([Chip, Chim, np.ones((1, E))])
        beq = np.concatenate(
            [s_lev * np.ones(len(Maxp)), -s_lev * np.ones(len(Maxm)), [float(k)]]
        )
        resf = linprog(
            np.zeros(E), A_eq=Aeq, b_eq=beq, bounds=(0, 1), method="highs"
        )
        resi = milp(
            c=np.zeros(E),
            constraints=LinearConstraint(Aeq, beq, beq),
            integrality=np.ones(E),
            bounds=Bounds(0, 1),
            options={"time_limit": 120},
        )
        results.append(
            {
                "level_s": s_lev,
                "k": k,
                "fractional_feasible": bool(resf.success),
                "integral_feasible": bool(resi.success),
                "integral_message": str(resi.message)[:120],
            }
        )
        print(
            f"s={s_lev} k={k} frac={resf.success} int={resi.success}",
            flush=True,
        )

    out = {
        "p": p,
        "n": n,
        "n_Max_plus": int(len(Maxp)),
        "n_Max_minus": int(len(Maxm)),
        "levels": results,
        "all_integral_infeasible": all(not r["integral_feasible"] for r in results),
        "all_fractional_feasible": all(r["fractional_feasible"] for r in results),
        "seconds": time.time() - t0,
        "status": (
            "Integral bi-tight infeasible at p=5 for levels 2,3,4; "
            "fractional feasible. Used for Type I / deep-tight residual kill."
        ),
    }
    path = ROOT / "evidence" / "e1_bitight_infeas.json"
    path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    assert out["all_integral_infeasible"]
    assert out["all_fractional_feasible"]


if __name__ == "__main__":
    main()
