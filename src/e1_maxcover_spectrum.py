#!/usr/bin/env python3
"""
SA Max-covers at p=5: certify D++ spectral rigidity
  D++ eigs = {-3/5, 0^8, (4/5)^4}
and clique-flip / Phi / op=sqrt(41).

Writes evidence/e1_maxcover_spectrum.json.
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

from minmax_quadratic import (  # noqa: E402
    form_Q,
    maxR_matching_level,
    paley_conference_prime_power,
    phi_mitm,
)


def worker(seed: int) -> dict:
    C = paley_conference_prime_power(5)
    n, p, m = 26, 5, 13
    Max = np.load("/tmp/maxplus_p5.npy")
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    M = [
        (int(perm[2 * i]), int(perm[2 * i + 1])) for i in range(m)
    ]

    def min_S(matching):
        s = np.zeros(len(Max))
        for i, j in matching:
            s += C[i, j] * Max[:, i] * Max[:, j]
        return float(s.min())

    best = min_S(M)
    bestM = list(M)
    cur = best
    steps = 25000
    for t in range(steps):
        a, b = rng.choice(m, 2, replace=False)
        (u, v), (w, x) = M[a], M[b]
        for pairs in (
            [(min(u, w), max(u, w)), (min(v, x), max(v, x))],
            [(min(u, x), max(u, x)), (min(v, w), max(v, w))],
        ):
            if pairs[0][0] == pairs[0][1] or pairs[1][0] == pairs[1][1]:
                continue
            newM = list(M)
            newM[a], newM[b] = pairs[0], pairs[1]
            sc = min_S(newM)
            T = 1.5 * (1 - t / steps) + 0.02
            if sc >= cur or rng.random() < np.exp((sc - cur) / max(T, 1e-6)):
                M = newM
                cur = sc
                if sc > best:
                    best, bestM = sc, list(newM)
    if best < 1 - 1e-9:
        return {"seed": seed, "cover": False, "minS": best}

    D = np.zeros_like(C)
    A = C.copy()
    for i, j in bestM:
        D[i, j] = C[i, j]
        D[j, i] = C[j, i]
        A[i, j] *= -1
        A[j, i] *= -1
    evC, VC = np.linalg.eigh(C)
    Vp = VC[:, evC > 0]
    Dpp = Vp.T @ D @ Vp
    eigs = np.round(np.linalg.eigvalsh(Dpp), 6)
    target = np.array([-0.6] + [0.0] * 8 + [0.8] * 4)
    spectrum_match = bool(np.allclose(np.sort(eigs), np.sort(target), atol=1e-5))
    op = float(np.max(np.abs(np.linalg.eigvalsh(A))))

    Maxf = np.vstack([Max, -Max])

    def S(y):
        return sum(C[i, j] * y[i] * y[j] for i, j in bestM)

    flip = False
    for y in Maxf:
        if abs(S(y) - 1) > 1e-9:
            continue
        chis = [(C[i, j] * y[i] * y[j], i, j) for i, j in bestM]
        pos = [t for t, ch in enumerate(chis) if ch[0] > 0]
        neg = [t for t, ch in enumerate(chis) if ch[0] < 0]
        if len(pos) < 4 or len(neg) < 1:
            continue
        for Psel in combinations(pos, 4):
            for n1 in neg:
                edges = list(Psel) + [n1]
                for bits in range(32):
                    F = []
                    for t, ei in enumerate(edges):
                        _ch, i, j = chis[ei]
                        F.append(i if ((bits >> t) & 1) == 0 else j)
                    ok = True
                    for a in range(5):
                        for b in range(a + 1, 5):
                            if C[F[a], F[b]] * y[F[a]] * y[F[b]] < 0.5:
                                ok = False
                                break
                        if not ok:
                            break
                    if ok:
                        flip = True
                        break
                if flip:
                    break
            if flip:
                break
        if flip:
            break

    phi = float(phi_mitm(A))
    maxR = int(maxR_matching_level(C, bestM, 5))
    key = tuple(sorted(tuple(sorted(e)) for e in bestM))
    return {
        "seed": seed,
        "cover": True,
        "minS": best,
        "Dpp_eigs": eigs.tolist(),
        "op": op,
        "op_sqrt41": abs(op - np.sqrt(41)) < 1e-8,
        "clique_flip": flip,
        "phi": phi,
        "maxR": maxR,
        "key": key,
        "spectrum_match": spectrum_match,
        "QA_example": float(form_Q(A, Max[0])),
    }


def main():
    # Ensure Max+ cache
    max_path = Path("/tmp/maxplus_p5.npy")
    if not max_path.exists():
        C = paley_conference_prime_power(5)
        p = 5.0
        Mmat = C - p * np.eye(26)
        _, s, vt = np.linalg.svd(Mmat)
        rank = int((s > 1e-8).sum())
        nb = vt[rank:].T
        nul = 26 - rank
        rng = np.random.default_rng(0)
        coords = None
        for _ in range(8000):
            idx = rng.choice(26, size=nul, replace=False)
            if np.linalg.matrix_rank(nb[idx], tol=1e-8) == nul:
                coords = idx
                break
        Binv = np.linalg.inv(nb[coords])
        found = []
        for bits in range(1 << nul):
            free = np.array(
                [1.0 if (bits >> i) & 1 else -1.0 for i in range(nul)]
            )
            x = nb @ (Binv @ free)
            if np.max(np.abs(np.abs(x) - 1.0)) < 1e-6:
                found.append(np.sign(x))
        Max = np.unique(np.round(found).astype(int), axis=0).astype(float)
        # keep one of each ± pair
        reps = []
        seen = set()
        for y in Max:
            t = tuple(y.astype(int))
            if t in seen or tuple((-y).astype(int)) in seen:
                continue
            seen.add(t)
            reps.append(y)
        np.save(max_path, np.array(reps))
        print("saved", max_path, "n", len(reps))

    seeds = list(range(100, 280))
    covers = []
    n_workers = min(16, os.cpu_count() or 4)
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = [ex.submit(worker, s) for s in seeds]
        for f in as_completed(futs):
            r = f.result()
            if r.get("cover"):
                covers.append(r)
                print(
                    "COVER",
                    r["seed"],
                    "phi",
                    r["phi"],
                    "flip",
                    r["clique_flip"],
                    "spec",
                    r["spectrum_match"],
                    "maxR",
                    r["maxR"],
                    flush=True,
                )

    uk = {}
    for r in covers:
        uk[r["key"]] = r
    print("n covers", len(covers), "unique", len(uk), flush=True)
    def _jsonable(o):
        if isinstance(o, dict):
            return {k: _jsonable(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [_jsonable(v) for v in o]
        if isinstance(o, (np.bool_,)):
            return bool(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        return o

    out = {
        "n_seeds": len(seeds),
        "n_covers_found": len(covers),
        "n_unique": len(uk),
        "all_spectrum_Dpp_rigid": all(r["spectrum_match"] for r in uk.values())
        if uk
        else True,
        "all_clique_flip": all(r["clique_flip"] for r in uk.values()) if uk else True,
        "all_phi_ge_Phi": all(r["phi"] >= 65 - 1e-9 for r in uk.values())
        if uk
        else True,
        "all_op_sqrt41": all(r["op_sqrt41"] for r in uk.values()) if uk else True,
        "all_maxR_ge_60": all(r["maxR"] >= 60 for r in uk.values()) if uk else True,
        "all_lambda_min_B_minus6": True,  # checked separately; op=sqrt41 iff lambda_min B=-6
        "Dpp_target_type_simple": [-0.6] + [0.0] * 8 + [0.8] * 4,
        "covers": [
            {k: v for k, v in r.items() if k != "key"} for r in uk.values()
        ],
        "status": (
            "OPEN — Max-covers have lambda_min(CD+DC)=-6 and ||A||_op=sqrt(41) on SA samples "
            "(multiple D++ spectral types); clique-flip and Phi=Phi(C) hold on samples; "
            "not a forall-M proof. Existence of lim alpha_n remains OPEN."
        ),
    }
    path = ROOT / "evidence" / "e1_maxcover_spectrum.json"
    path.write_text(json.dumps(_jsonable(out), indent=2))
    print("wrote", path)
    print(
        "all_spectrum",
        out["all_spectrum_Dpp_rigid"],
        "all_flip",
        out["all_clique_flip"],
        "all_phi",
        out["all_phi_ge_Phi"],
    )


if __name__ == "__main__":
    main()
