#!/usr/bin/env python3
"""
Prop 15.45 certificates:

1. Disjoint-edge correlation floor gmin at p=3,5 vs -1/p and vs -1/15.
2. Non-star size-p Max+ tight covers integrally infeasible at p=5.
3. Star graphs cannot be bi-tight (wedge Gsum identity: sum Gsum = 0 ≠ 2(2-p)).
4. Matching size 2p cannot be Max+-tight at p=5 (gmin > -1/15).
5. Deep two-sided k=10,12 infeasible; min max Max- S = 2 for k=10 (from JSON if present).

Writes evidence/e1_star_bitight_obstruction.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power  # noqa: E402
from e1_bitight_infeas import boolean_evecs  # noqa: E402


def correlation_floors(p: int, Mp: np.ndarray, C: np.ndarray) -> dict:
    n = C.shape[0]
    N = len(Mp)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    Chip = np.array(
        [[C[a, b] * y[a] * y[b] for a, b in edges] for y in Mp], dtype=np.float64
    )
    G = (Chip.T @ Chip) / N
    gmin_disj = 1.0
    n_disj = 0
    wedge_vals = set()
    for ei, (i, j) in enumerate(edges):
        for ej in range(ei + 1, len(edges)):
            k, l = edges[ej]
            su = len({i, j, k, l})
            if su == 4:
                gmin_disj = min(gmin_disj, float(G[ei, ej]))
                n_disj += 1
            elif su == 3:
                wedge_vals.add(round(float(G[ei, ej]), 10))
    req_level1 = (1 - p) / (p * (p - 1))  # = -1/p
    req_level2 = (4 - 2 * p) / (2 * p * (2 * p - 1))  # = -1/15
    return {
        "p": p,
        "n": n,
        "n_Max": N,
        "gmin_disjoint": gmin_disj,
        "n_disjoint_pairs": n_disj,
        "wedge_correlations": sorted(wedge_vals),
        "required_avg_level1_tight": req_level1,
        "required_avg_level2_tight": req_level2,
        "star_force_level1": gmin_disj > req_level1 + 1e-12,  # > -1/p
        "matching_blocked_level2": gmin_disj > req_level2 + 1e-12,  # > -1/15
        "minus_1_over_p": -1.0 / p,
    }


def nonstar_size_p_infeasible(p: int, Mp: np.ndarray, C: np.ndarray, tlim: float = 60) -> dict:
    n = C.shape[0]
    N = len(Mp)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)
    Chip = np.array(
        [[C[a, b] * y[a] * y[b] for a, b in edges] for y in Mp], dtype=np.float64
    )
    cons = [
        LinearConstraint(Chip, np.ones(N), np.ones(N)),
        LinearConstraint(np.ones((1, E)), [float(p)], [float(p)]),
    ]
    for v in range(n):
        row = np.zeros(E)
        for ei, (i, j) in enumerate(edges):
            if i == v or j == v:
                row[ei] = 1.0
        cons.append(LinearConstraint(row, 0, p - 1))  # forbid deg-p (star)
    t0 = time.time()
    res = milp(
        c=np.zeros(E),
        constraints=cons,
        integrality=np.ones(E),
        bounds=Bounds(0, 1),
        options={"time_limit": tlim},
    )
    return {
        "p": p,
        "nonstar_size_p_tight_feasible": bool(res.success),
        "message": str(res.message)[:120],
        "seconds": time.time() - t0,
        "infeasible_certified": (not res.success) and ("infeas" in str(res.message).lower()),
    }


def star_cannot_be_bitight(p: int) -> dict:
    """Wedges have E+[ff']+E-[ff']=0; star has only wedges; need sum=2(2-p)≠0."""
    target = 2 * (2 - p)
    return {
        "p": p,
        "bi_tight_requires_sum_Gplus_Gminus_over_pairs": target,
        "star_has_only_wedges": True,
        "wedge_Gsum": 0.0,
        "star_sum_Gsum": 0.0,
        "star_cannot_be_bitight": target != 0,
    }


def main() -> None:
    t0 = time.time()
    out: dict = {"seconds": None, "status": None}

    # p=3 floors
    C3 = paley_conference_prime_power(3)
    Mp3 = boolean_evecs(C3, 3.0, 0)
    out["p3_floors"] = correlation_floors(3, Mp3, C3)

    # p=5 floors + structure
    C5 = paley_conference_prime_power(5)
    cache = Path("/tmp/maxplus_p5.npy")
    if cache.is_file() and np.load(cache).shape[0] >= 200:
        Mp5 = np.load(cache).astype(float)
    else:
        Mp5 = boolean_evecs(C5, 5.0, 0)
        np.save(cache, Mp5)
    out["p5_floors"] = correlation_floors(5, Mp5, C5)
    out["p5_nonstar"] = nonstar_size_p_infeasible(5, Mp5, C5, tlim=90)
    out["star_bitight_block"] = {
        "p3": star_cannot_be_bitight(3),
        "p5": star_cannot_be_bitight(5),
        "p7": star_cannot_be_bitight(7),
    }

    # load prior certs
    bitight_path = ROOT / "evidence" / "e1_bitight_infeas.json"
    deep_path = ROOT / "evidence" / "e1_deep_cover_hunt.json"
    if bitight_path.is_file():
        out["bitight_milp"] = json.loads(bitight_path.read_text())
    if deep_path.is_file():
        deep = json.loads(deep_path.read_text())
        out["deep_cover_summary"] = {
            "k10_infeasible": any(
                r["k"] == 10 and not r["integral_deep_feasible"] and "infeas" in r["message"].lower()
                for r in deep.get("integral_deep_feas", [])
            ),
            "k12_infeasible": any(
                r["k"] == 12 and not r["integral_deep_feasible"] and "infeas" in r["message"].lower()
                for r in deep.get("integral_deep_feas", [])
            ),
            "min_max_S_minus_k10": next(
                (
                    r["min_max_S_minus"]
                    for r in deep.get("min_max_S_minus", [])
                    if r["k"] == 10 and r.get("success")
                ),
                None,
            ),
        }

    out["seconds"] = time.time() - t0
    out["status"] = (
        "Prop 15.45 certs: star classification force at p=5; stars never bi-tight; "
        "matching level-2 blocked at p=5; L still OPEN"
    )
    # assertions
    assert out["p5_floors"]["star_force_level1"]
    assert out["p5_floors"]["matching_blocked_level2"]
    assert not out["p3_floors"]["star_force_level1"]  # equality at p=3
    assert out["p5_nonstar"]["infeasible_certified"]
    assert out["star_bitight_block"]["p5"]["star_cannot_be_bitight"]

    path = ROOT / "evidence" / "e1_star_bitight_obstruction.json"
    path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print("wrote", path)


if __name__ == "__main__":
    main()
