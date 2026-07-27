#!/usr/bin/env python3
"""Census: perfect-matching Max-covers at Paley n=26 and their exact Phi."""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power, phi_mitm  # noqa: E402


def boolean_plus(C: np.ndarray, p: float) -> np.ndarray:
    n = C.shape[0]
    M = C - p * np.eye(n)
    _, s, vt = np.linalg.svd(M)
    rank = int((s > 1e-8).sum())
    nb = vt[rank:].T
    nul = n - rank
    rng = np.random.default_rng(0)
    coords = None
    for _ in range(8000):
        idx = rng.choice(n, size=nul, replace=False)
        if np.linalg.matrix_rank(nb[idx], tol=1e-8) == nul:
            coords = idx
            break
    Binv = np.linalg.inv(nb[coords])
    found = []
    for bits in range(1 << nul):
        free = np.array([1.0 if (bits >> i) & 1 else -1.0 for i in range(nul)])
        x = nb @ (Binv @ free)
        if np.max(np.abs(np.abs(x) - 1.0)) < 1e-6:
            found.append(np.sign(x))
    return np.unique(np.round(found).astype(int), axis=0).astype(float)


def worker(seed: int) -> dict:
    p = 5
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    Maxp = boolean_plus(C, float(p))
    Maxm = boolean_plus(-C, float(p))
    Phi = 0.5 * n * p
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    M = [
        (int(min(perm[2 * i], perm[2 * i + 1])), int(max(perm[2 * i], perm[2 * i + 1])))
        for i in range(n // 2)
    ]

    def evalM(matching):
        sp = np.zeros(len(Maxp))
        sm = np.zeros(len(Maxm))
        for i, j in matching:
            sp += C[i, j] * Maxp[:, i] * Maxp[:, j]
            sm += C[i, j] * Maxm[:, i] * Maxm[:, j]
        return sp.min(), sm.max(), sp, sm

    best = evalM(M)
    bestM = M
    cur = best
    iters = 20000
    for t in range(iters):
        a, b = rng.choice(n // 2, 2, replace=False)
        (u, v), (w, x) = M[a], M[b]
        for pairs in (
            [(min(u, w), max(u, w)), (min(v, x), max(v, x))],
            [(min(u, x), max(u, x)), (min(v, w), max(v, w))],
        ):
            if pairs[0][0] == pairs[0][1] or pairs[1][0] == pairs[1][1]:
                continue
            if len(set([pairs[0][0], pairs[0][1], pairs[1][0], pairs[1][1]])) < 4:
                continue
            newM = M[:]
            newM[a], newM[b] = pairs[0], pairs[1]
            sc = evalM(newM)
            score = (sc[0], -sc[1])
            cscore = (cur[0], -cur[1])
            temp = 2 + 8 * t / iters
            if score > cscore or rng.random() < np.exp((sc[0] - cur[0]) * temp):
                M = newM
                cur = sc
                if (cur[0], -cur[1]) > (best[0], -best[1]):
                    best = cur
                    bestM = M[:]
    mn_sp, mx_sm, sp, sm = best[0], best[1], best[2], best[3]
    if mn_sp < 1 - 1e-9:
        return {"seed": seed, "cover": False, "min_Sp": float(mn_sp), "max_Sm": float(mx_sm)}
    A = C.copy()
    for i, j in bestM:
        A[i, j] *= -1
        A[j, i] *= -1
    max_only = max(max(abs(Phi - 2 * s) for s in sp), max(abs(-Phi - 2 * s) for s in sm))
    ph = float(phi_mitm(A))
    return {
        "seed": seed,
        "cover": True,
        "two_sided": bool(mx_sm <= -1 + 1e-9),
        "min_Sp": float(mn_sp),
        "max_Sm": float(mx_sm),
        "max_on_Max_pm": float(max_only),
        "phi": ph,
        "Phi_C": float(Phi),
        "undercut": bool(ph < Phi - 1e-9),
        "spike_above_Max": bool(ph > max_only + 1e-9),
        "matching": bestM,
    }


def main() -> None:
    workers = max(1, (os.cpu_count() or 2) - 2)
    seeds = list(range(2000, 2050))
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(worker, s) for s in seeds]
        for f in as_completed(futs):
            r = f.result()
            rows.append(r)
            if r.get("cover"):
                print(
                    {k: r[k] for k in r if k != "matching"},
                    flush=True,
                )
    covers = [r for r in rows if r.get("cover")]
    two = [r for r in covers if r.get("two_sided")]
    und = [r for r in covers if r.get("undercut")]
    summary = {
        "n_seeds": len(seeds),
        "n_covers": len(covers),
        "n_two_sided": len(two),
        "n_undercut": len(und),
        "phi_of_covers": sorted(set(r["phi"] for r in covers)),
        "max_on_Max_pm_values": sorted(set(r["max_on_Max_pm"] for r in covers)),
        "all_covers_phi_ge_PhiC": all(r["phi"] + 1e-9 >= r["Phi_C"] for r in covers),
        "all_two_sided_phi_ge_PhiC": all(r["phi"] + 1e-9 >= r["Phi_C"] for r in two),
        "covers": [{k: v for k, v in r.items() if k != "matching"} for r in covers],
        "status": "OPEN — certified numerics, not a general E(1) proof",
    }
    print("SUMMARY", json.dumps({k: summary[k] for k in summary if k != "covers"}, indent=2))
    (ROOT / "evidence" / "e1_n26_matching_covers_census.json").write_text(json.dumps(summary, indent=2))
    if covers:
        (ROOT / "evidence" / "e1_n26_matching_cover_example.json").write_text(
            json.dumps(covers[0], indent=2)
        )
    print("wrote evidence/e1_n26_matching_covers_census.json")


if __name__ == "__main__":
    main()
