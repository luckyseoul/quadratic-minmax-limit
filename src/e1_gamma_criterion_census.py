#!/usr/bin/env python3
"""
Census: Γ-pairing spike criterion maxR >= Phi-p for perfect matchings at p=5.
- Exact maxR on S_M=-p via pairing enum
- SA Max-covers + known example
- Random + min-maxR local search for counterexamples
Writes evidence/e1_gamma_forall_census.json
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power, phi_mitm  # noqa: E402


def build_max_plus(C: np.ndarray, p: int) -> np.ndarray:
    """MITM boolean +p evecs, x[0]=+1 reps then expand."""
    n = C.shape[0]
    half = n // 2
    Cf = C.astype(np.float64)
    C11, C12 = Cf[:half, :half], Cf[:half, half:]
    targets: dict[tuple, list[int]] = {}
    for mask in range(1 << half):
        xL = np.array([1.0 if (mask >> i) & 1 == 0 else -1.0 for i in range(half)])
        t = (p * np.eye(half) - C11) @ xL
        key = tuple(np.rint(t).astype(np.int16))
        targets.setdefault(key, []).append(mask)
    found = set()
    for mask in range(1 << half):
        xR = np.array([1.0 if (mask >> i) & 1 == 0 else -1.0 for i in range(half)])
        key = tuple(np.rint(C12 @ xR).astype(np.int16))
        if key not in targets:
            continue
        for mL in targets[key]:
            xL = np.array([1.0 if (mL >> i) & 1 == 0 else -1.0 for i in range(half)])
            x = np.concatenate([xL, xR])
            if np.allclose(Cf @ x, p * x):
                if x[0] < 0:
                    x = -x
                found.add(tuple(np.rint(x).astype(np.int8)))
    reps = [np.array(t, dtype=np.float64) for t in found]
    return np.vstack(reps + [-r for r in reps])


def maxR_on_level(C: np.ndarray, pairs, p: int) -> int:
    m = len(pairs)
    pairs = [(int(a), int(b)) for a, b in pairs]
    c = np.array([int(C[i, j]) for i, j in pairs], dtype=np.int8)
    U = np.ones((1 << m, m), dtype=np.float64)
    for a in range(m):
        U[:, a] = np.where((np.arange(1 << m) >> a) & 1, -1.0, 1.0)
    best = -10**9
    C8 = C.astype(np.int8)
    for mask in range(1 << m):
        z = U[mask]
        if int(c @ z) != -p:
            continue
        G = np.zeros((m, m), dtype=np.float64)
        for a in range(m):
            i, j = pairs[a]
            za = int(z[a])
            for b in range(a + 1, m):
                k, l = pairs[b]
                zb = int(z[b])
                g = (
                    int(C8[i, k])
                    + int(C8[i, l]) * zb
                    + int(C8[j, k]) * za
                    + int(C8[j, l]) * za * zb
                )
                G[a, b] = G[b, a] = g
        R = 0.5 * np.sum(U * (U @ G), axis=1)
        best = max(best, int(R.max()))
    return int(best)


def min_S(Max: np.ndarray, C: np.ndarray, pairs) -> float:
    s = np.zeros(Max.shape[0])
    for i, j in pairs:
        s += C[i, j] * Max[:, i] * Max[:, j]
    return float(s.min())


def sa_cover(Max, C, n, m, seed: int, steps: int = 30000):
    rng = np.random.default_rng(seed)
    perm = rng.permutation(n)
    pairs = [(int(perm[2 * i]), int(perm[2 * i + 1])) for i in range(m)]
    cur = min_S(Max, C, pairs)
    best, bp = cur, list(pairs)
    for t in range(steps):
        a, b = int(rng.integers(0, m)), int(rng.integers(0, m))
        if a == b:
            continue
        i, j = pairs[a]
        k, l = pairs[b]
        if rng.integers(0, 2) == 0:
            na, nb = (i, k), (j, l)
        else:
            na, nb = (i, l), (j, k)
        if len({na[0], na[1], nb[0], nb[1]}) < 4:
            continue
        old = pairs[a], pairs[b]
        pairs[a], pairs[b] = na, nb
        ns = min_S(Max, C, pairs)
        T = 2.0 * (1 - t / steps) + 1e-6
        if ns >= cur or rng.random() < np.exp((ns - cur) / T):
            cur = ns
            if cur > best:
                best, bp = cur, list(pairs)
                if best >= 1.0 - 1e-9:
                    return bp, best
        else:
            pairs[a], pairs[b] = old
    return bp, best


def main():
    p = 5
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    m = n // 2
    Phi = 0.5 * n * p
    need_R = Phi - p
    print(f"n={n} m={m} Phi={Phi} need_R={need_R}", flush=True)

    cache = Path("/tmp/maxplus_p5.npy")
    if cache.exists():
        reps = np.load(cache).astype(np.float64)
        Max = np.vstack([reps, -reps])
    else:
        Max = build_max_plus(C, p)
    print(f"|Max+|={Max.shape[0]}", flush=True)

    known = [
        [4, 12],
        [5, 14],
        [17, 24],
        [8, 11],
        [1, 2],
        [7, 10],
        [9, 18],
        [19, 25],
        [3, 16],
        [0, 21],
        [20, 22],
        [6, 15],
        [13, 23],
    ]
    covers = [known]
    print(f"known cover minS={min_S(Max, C, known)}", flush=True)

    for seed in range(1000, 1100):
        t0 = time.time()
        bp, best = sa_cover(Max, C, n, m, seed, steps=35000)
        print(f"seed={seed} minS={best:.1f} t={time.time()-t0:.1f}s", flush=True)
        if best >= 1.0 - 1e-9:
            covers.append(bp)
            print(f"  COVER seed={seed}", flush=True)

    cover_rows = []
    for i, pairs in enumerate(covers):
        mr = maxR_on_level(C, pairs, p)
        A = C.copy()
        for a, b in pairs:
            A[a, b] *= -1
            A[b, a] *= -1
        phi = float(phi_mitm(A))
        row = {
            "idx": i,
            "maxR": mr,
            "need_R": need_R,
            "criterion": mr >= need_R - 1e-9,
            "phi": phi,
            "Phi_C": Phi,
            "undercut": phi < Phi - 1e-9,
            "minS": min_S(Max, C, pairs),
            "matching": [[int(a), int(b)] for a, b in pairs],
        }
        cover_rows.append(row)
        print(row, flush=True)

    rng = np.random.default_rng(0)
    random_maxR = []
    worst = 10**9
    for s in range(80):
        perm = rng.permutation(n)
        pairs = [(int(perm[2 * i]), int(perm[2 * i + 1])) for i in range(m)]
        mr = maxR_on_level(C, pairs, p)
        random_maxR.append(mr)
        if mr < worst:
            worst = mr
            print(f"random worst maxR={worst} at s={s}", flush=True)

    # local search minimize maxR from worst random starts
    def sa_min_mr(seed, steps=40):
        rng2 = np.random.default_rng(seed)
        perm = rng2.permutation(n)
        pairs = [(int(perm[2 * i]), int(perm[2 * i + 1])) for i in range(m)]
        cur = maxR_on_level(C, pairs, p)
        best = cur
        for t in range(steps):
            a, b = int(rng2.integers(0, m)), int(rng2.integers(0, m))
            if a == b:
                continue
            i, j = pairs[a]
            k, l = pairs[b]
            if rng2.integers(0, 2) == 0:
                na, nb = (i, k), (j, l)
            else:
                na, nb = (i, l), (j, k)
            if len({na[0], na[1], nb[0], nb[1]}) < 4:
                continue
            old = pairs[a], pairs[b]
            pairs[a], pairs[b] = na, nb
            ns = maxR_on_level(C, pairs, p)
            if ns <= cur or rng2.random() < np.exp((cur - ns) / 3.0):
                cur = ns
                best = min(best, cur)
            else:
                pairs[a], pairs[b] = old
        return best

    sa_worst = min(sa_min_mr(s, 30) for s in range(12))
    print(f"SA min-maxR best (lowest)={sa_worst}", flush=True)

    out = {
        "p": p,
        "n": n,
        "need_R": need_R,
        "n_covers": len(covers),
        "covers": cover_rows,
        "all_covers_criterion": all(r["criterion"] for r in cover_rows),
        "all_covers_no_undercut": all(not r["undercut"] for r in cover_rows),
        "n_random": len(random_maxR),
        "random_maxR_min": int(min(random_maxR)),
        "random_maxR_values": sorted(set(random_maxR)),
        "sa_min_maxR": int(sa_worst),
        "no_counterexample_found": min(random_maxR) >= need_R and sa_worst >= need_R,
        "status": (
            "OPEN — census supports criterion for all tested M at p=5; "
            "not a forall-M proof; matching dichotomy / k_star still open"
        ),
    }
    path = ROOT / "evidence" / "e1_gamma_forall_census.json"
    path.write_text(json.dumps(out, indent=2))
    print("wrote", path, flush=True)
    print(json.dumps({k: out[k] for k in out if k != "covers"}, indent=2), flush=True)


if __name__ == "__main__":
    main()
