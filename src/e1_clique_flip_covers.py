#!/usr/bin/env python3
"""
SA for matching Max-covers at p=5; for each cover try clique-flip spike construction.
Writes evidence/e1_clique_flip_covers.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power, form_Q, phi_mitm  # noqa: E402


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


def clique_flip_spike(C, M, Maxp, p, Phi):
    """Return True if some y in Max+ with S=1 admits r=p clique flip giving QA>=Phi."""
    pairs = list(M)
    m = len(pairs)
    Sigma_needed = (1 + p) // 2  # for s0=1
    r = p
    for y in Maxp:
        s0 = sum(C[i, j] * y[i] * y[j] for i, j in pairs)
        if abs(s0 - 1) > 1e-9:
            continue
        W = np.outer(y, y) * C
        np.fill_diagonal(W, 0)
        edge_chi = [(C[i, j] * y[i] * y[j], i, j) for i, j in pairs]
        for Psel in combinations(range(m), r):
            chis = [edge_chi[b][0] for b in Psel]
            if abs(sum(chis) - Sigma_needed) > 1e-9:
                continue
            for bits in range(1 << r):
                F = []
                for t, b in enumerate(Psel):
                    _ch, i, j = edge_chi[b]
                    F.append(i if ((bits >> t) & 1) == 0 else j)
                ok = True
                for a in range(r):
                    for b2 in range(a + 1, r):
                        if W[F[a], F[b2]] < 0.5:
                            ok = False
                            break
                    if not ok:
                        break
                if not ok:
                    continue
                x = y.copy()
                for i in F:
                    x[i] *= -1
                QC = form_Q(C, x)
                SM = sum(C[i, j] * x[i] * x[j] for i, j in pairs)
                if QC >= Phi - 2 * p - 1e-9 and abs(SM + p) < 1e-9:
                    return True, dict(QC=float(QC), SM=float(SM), QA=float(QC - 2 * SM))
    return False, None


def worker(seed: int) -> dict:
    p = 5
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    Phi = 0.5 * n * p
    Maxp = boolean_plus(C, float(p))
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    M = [
        (int(min(perm[2 * i], perm[2 * i + 1])), int(max(perm[2 * i], perm[2 * i + 1])))
        for i in range(n // 2)
    ]

    def min_S(matching):
        s = np.zeros(len(Maxp))
        for i, j in matching:
            s += C[i, j] * Maxp[:, i] * Maxp[:, j]
        return float(s.min())

    best = min_S(M)
    bestM = M
    cur = best
    for t in range(25000):
        a, b = rng.choice(n // 2, 2, replace=False)
        (u, v), (w, x) = M[a], M[b]
        for pairs in (
            [(min(u, w), max(u, w)), (min(v, x), max(v, x))],
            [(min(u, x), max(u, x)), (min(v, w), max(v, w))],
        ):
            if pairs[0][0] == pairs[0][1] or pairs[1][0] == pairs[1][1]:
                continue
            if len({pairs[0][0], pairs[0][1], pairs[1][0], pairs[1][1]}) < 4:
                continue
            newM = M[:]
            newM[a], newM[b] = pairs[0], pairs[1]
            sc = min_S(newM)
            if sc > cur or rng.random() < np.exp((sc - cur) * (3 + 10 * t / 25000)):
                M = newM
                cur = sc
                if cur > best:
                    best = cur
                    bestM = [tuple(e) for e in M]
    if best < 1 - 1e-9:
        return {"seed": seed, "cover": False, "min_S": best}
    ok, info = clique_flip_spike(C, bestM, Maxp, p, Phi)
    A = C.copy()
    for i, j in bestM:
        A[i, j] *= -1
        A[j, i] *= -1
    ph = float(phi_mitm(A))
    return {
        "seed": seed,
        "cover": True,
        "min_S": best,
        "clique_flip": ok,
        "flip_info": info,
        "phi": ph,
        "Phi_C": float(Phi),
        "undercut": ph < Phi - 1e-9,
        "matching": [list(e) for e in bestM],
    }


def main() -> None:
    workers = max(1, (os.cpu_count() or 2) - 2)
    seeds = list(range(3000, 3100))  # 100 seeds
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = [ex.submit(worker, s) for s in seeds]
        for f in as_completed(futs):
            r = f.result()
            rows.append(r)
            if r.get("cover"):
                print({k: r[k] for k in r if k != "matching"}, flush=True)
    covers = [r for r in rows if r.get("cover")]
    summary = {
        "n_seeds": len(seeds),
        "n_covers": len(covers),
        "n_clique_flip_ok": sum(1 for r in covers if r.get("clique_flip")),
        "n_undercut": sum(1 for r in covers if r.get("undercut")),
        "all_covers_clique_flip": all(r.get("clique_flip") for r in covers) if covers else True,
        "all_covers_phi_ge_PhiC": all(r["phi"] + 1e-9 >= r["Phi_C"] for r in covers) if covers else True,
        "covers": [{k: v for k, v in r.items() if k != "matching"} for r in covers],
        "status": "OPEN — census of clique-flip on SA covers; not a forall-M proof",
    }
    print("SUMMARY", json.dumps({k: summary[k] for k in summary if k != "covers"}, indent=2))
    (ROOT / "evidence" / "e1_clique_flip_covers.json").write_text(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
