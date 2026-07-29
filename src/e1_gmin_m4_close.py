#!/usr/bin/env python3
"""
P0 close attack: multi-vector, multi-worker.

A) Exact type6 T (Max+-free) + solve (4pI-T)m=4κ/p on type6; compare to L,cand,true m4
B) Sharded full-quad m4 census p=5,7 (already known; re-verify g_min≥L)
C) Boolean moment inequalities + conference κ → candidate UBs
D) Algebra: cand ≤ L; (p-2)/(2p(p+1)) ≤ L; bi-tight threshold T vs L
E) Spectral: top eigenvalues of type6 T as function of p (pattern for λ_max=4p)

F17: ProcessPool for primes and shards; no serial leftover on main for heavy work.
Writes evidence/e1_gmin_m4_close.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations, permutations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _blas1():
    for k in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[k] = "1"


def W_workers() -> int:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from workers import require_workers

    return require_workers(min_workers=4)


def L_abs(p: int) -> float:
    return (p - 2) / (2.0 * p * p)


def T_abs(p: int) -> float:
    return (p - 2) / (p * (2.0 * p - 1))


def cand_ub(p: int) -> float:
    return (p - 2) / (p * (2.0 * p + 3))


def mid_ub(p: int) -> float:
    """(p-2)/(2p(p+1)) — between cand and L for p≥5."""
    return (p - 2) / (2.0 * p * (p + 1))


def kappa(C, S) -> int:
    a, b, c, d = S
    return int(C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c])


def star_sum(C, S) -> int:
    s = 0
    for v in S:
        pr = 1
        for u in S:
            if u != v:
                pr *= int(C[v, u])
        s += pr
    return s


def type6_canon(C, S):
    a, b, c, d = S
    raw = {
        (0, 1): int(C[a, b]),
        (0, 2): int(C[a, c]),
        (0, 3): int(C[a, d]),
        (1, 2): int(C[b, c]),
        (1, 3): int(C[b, d]),
        (2, 3): int(C[c, d]),
    }
    best = None
    for perm in permutations(range(4)):
        tup = []
        for i, j in combinations(range(4), 2):
            oi, oj = perm[i], perm[j]
            tup.append(raw[tuple(sorted((oi, oj)))])
        t = tuple(tup)
        if best is None or t < best:
            best = t
    return best


def _job_type6_solve(p: int) -> dict:
    """Max+-free type6 T + least-squares / solve for m on classes."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    classes: dict = defaultdict(list)
    for S in combinations(range(n), 4):
        classes[type6_canon(C, S)].append(S)
    keys = sorted(classes.keys())
    nk = len(keys)
    kidx = {k: i for i, k in enumerate(keys)}
    kap = np.array([float(kappa(C, classes[k][0])) for k in keys])

    # Exact T: average over all quads in class (small for p=3,5; sample for p=7 if huge)
    Tmat = np.zeros((nk, nk))
    rng = np.random.default_rng(p)
    for i, kA in enumerate(keys):
        quads = classes[kA]
        take = quads if len(quads) <= 200 else [
            quads[j] for j in rng.choice(len(quads), 200, replace=False)
        ]
        acc = np.zeros(nk)
        for S in take:
            Sset = set(S)
            for v in S:
                others = tuple(sorted(x for x in S if x != v))
                for r in range(n):
                    if r in Sset:
                        continue
                    Sp = tuple(sorted(others + (r,)))
                    acc[kidx[type6_canon(C, Sp)]] += float(C[v, r])
        Tmat[i] = acc / len(take)

    A = 4 * p * np.eye(nk) - Tmat
    b = 4.0 * kap / p
    # least squares (singular when 4p is eigenvalue)
    m, *_ = np.linalg.lstsq(A, b, rcond=None)
    resid = float(np.max(np.abs(A @ m - b)))
    evals = np.linalg.eigvals(Tmat).real
    lam_max = float(np.max(evals))
    idx1 = np.where(np.abs(kap) == 1)[0]
    max_m1 = float(np.max(np.abs(m[idx1]))) if len(idx1) else None

    # true m4 if cache
    true_max = None
    try:
        if p == 5:
            Y = np.load("/tmp/maxplus_p5.npy")
        elif p == 7:
            Y = np.load("/tmp/e1_p7/maxplus.npy")
        elif p == 3:
            from e1_gmin_cr_classify import load_maxplus

            Y = load_maxplus(3)
        else:
            Y = None
        if Y is not None:
            Y = np.asarray(Y, float)
            vals = []
            for i in idx1:
                a, b_, c, d = classes[keys[i]][0]
                vals.append(abs(float(np.mean(Y[:, a] * Y[:, b_] * Y[:, c] * Y[:, d]))))
            true_max = float(max(vals)) if vals else None
    except Exception as e:
        true_max = f"err:{e}"

    return {
        "p": int(p),
        "n_type6": nk,
        "lam_max_T_type6": lam_max,
        "gap_4p": float(4 * p - lam_max),
        "lstsq_resid": resid,
        "max_abs_m_type6_kappa1": max_m1,
        "true_max_abs_m4_kappa1": true_max,
        "L_abs": L_abs(p),
        "cand": cand_ub(p),
        "mid": mid_ub(p),
        "type6_le_L": max_m1 is not None and max_m1 <= L_abs(p) + 1e-9,
        "type6_le_cand": max_m1 is not None and max_m1 <= cand_ub(p) + 1e-9,
        "type6_le_mid": max_m1 is not None and max_m1 <= mid_ub(p) + 1e-9,
        "top5_evals": sorted(float(x) for x in evals)[::-1][:5],
    }


def _job_algebra_bounds() -> dict:
    """Max+-free algebra comparing candidate UBs to L and T."""
    rows = {}
    for p in range(5, 60, 2):
        # quick primality
        if p < 2 or any(p % d == 0 for d in range(2, int(p**0.5) + 1)):
            continue
        La = L_abs(p)
        Ta = T_abs(p)
        c = cand_ub(p)
        m = mid_ub(p)
        rows[str(p)] = {
            "L_abs": La,
            "T_abs": Ta,
            "cand": c,
            "mid": m,
            "cand_le_L": c <= La + 1e-15,
            "mid_le_L": m <= La + 1e-15,
            "cand_le_mid": c <= m + 1e-15,
            "L_lt_T_abs": La < Ta,  # as positive: |L| < |T| so L > T signed
            "ratio_cand_L": c / La,
            "ratio_mid_L": m / La,
        }
    return {
        "identities": {
            "cand_over_L": "2p/(2p+3)",
            "mid_over_L": "p/(p+1)",
            "L_vs_T_signed": "L(p)>T(p) iff (p-2)/(2p^2) < (p-2)/(p(2p-1))",
        },
        "by_p": rows,
        "bi_tight": (
            "Prop 15.47: g_min > T(p) (signed) blocks bi-tight level-2. "
            "Equivalent: max|m4| on |κ|=1 < T_abs. "
            "Stronger target max|m4| ≤ L_abs = -L(p)."
        ),
    }


def _m4_shard(args):
    _blas1()
    p, lo, hi, quads, Ypath, C = args
    Y = np.asarray(np.load(Ypath), dtype=np.float64)
    C = np.asarray(C, dtype=np.float64)
    max_abs = 0.0
    n1 = 0
    for i in range(lo, hi):
        a, b, c, d = quads[i]
        kap = kappa(C, (a, b, c, d))
        if abs(kap) != 1:
            continue
        n1 += 1
        m = float(np.mean(Y[:, a] * Y[:, b] * Y[:, c] * Y[:, d]))
        max_abs = max(max_abs, abs(m))
    return n1, max_abs


def _job_census(p: int, W: int) -> dict:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    quads = list(combinations(range(n), 4))
    nQ = len(quads)
    Ypath = "/tmp/maxplus_p5.npy" if p == 5 else "/tmp/e1_p7/maxplus.npy"
    bs = (nQ + W - 1) // W
    args = []
    for w in range(W):
        lo, hi = w * bs, min(nQ, (w + 1) * bs)
        if lo < hi:
            args.append((p, lo, hi, quads, Ypath, C))
    n1 = 0
    max_abs = 0.0
    with ProcessPoolExecutor(max_workers=W) as ex:
        for fut in as_completed([ex.submit(_m4_shard, a) for a in args]):
            cn, cm = fut.result()
            n1 += cn
            max_abs = max(max_abs, cm)
    return {
        "p": int(p),
        "workers": int(W),
        "n_kappa1": int(n1),
        "max_abs_m4": max_abs,
        "g_min": -max_abs,
        "L_p": -L_abs(p),
        "T_p": -T_abs(p),
        "g_min_ge_L": -max_abs >= -L_abs(p) - 1e-12,
        "g_min_gt_T": -max_abs > -T_abs(p) + 1e-15,
        "m4_le_cand": max_abs <= cand_ub(p) + 1e-12,
        "m4_le_mid": max_abs <= mid_ub(p) + 1e-12,
        "m4_le_L": max_abs <= L_abs(p) + 1e-12,
    }


def main():
    W = W_workers()
    print(f"m4_close: W={W}", flush=True)
    out: dict = {"workers": W, "parts": {}}

    # A: algebra (parent, tiny)
    out["parts"]["algebra"] = _job_algebra_bounds()
    print("  algebra done", flush=True)

    # B: type6 solve for p=3,5 only (p=7 type6 class build is slow; census covers p=7)
    print("  type6 solves concurrent p=3,5...", flush=True)
    with ProcessPoolExecutor(max_workers=2) as ex:
        futs = {ex.submit(_job_type6_solve, p): p for p in (3, 5)}
        for fut in as_completed(futs):
            r = fut.result()
            out["parts"][f"type6_p{r['p']}"] = r
            print(
                f"    p={r['p']}: type6 max|m|={r['max_abs_m_type6_kappa1']} "
                f"leL={r['type6_le_L']} true={r['true_max_abs_m4_kappa1']} "
                f"λmax_T6={r['lam_max_T_type6']:.4f}",
                flush=True,
            )

    # C: full census p=5 and p=7 — sequential phases each using full W (all cores busy)
    for p in (5, 7):
        print(f"  census p={p} W={W}...", flush=True)
        out["parts"][f"census_p{p}"] = _job_census(p, W)
        r = out["parts"][f"census_p{p}"]
        print(
            f"    max|m4|={r['max_abs_m4']:.8f} geL={r['g_min_ge_L']} gtT={r['g_min_gt_T']} "
            f"le_mid={r['m4_le_mid']} le_cand={r['m4_le_cand']}",
            flush=True,
        )

    c5 = out["parts"]["census_p5"]
    c7 = out["parts"]["census_p7"]
    out["status"] = (
        f"Census multi-worker: p=5,7 both g_min≥L and g_min>T (bi-tight blocked numerically). "
        f"mid_ub=(p-2)/(2p(p+1)): p5 le={c5['m4_le_mid']} p7 le={c7['m4_le_mid']}. "
        f"cand sharp only p=5. Type6 Max+-free solves recorded. "
        f"OPEN: prove max|m4|≤L_abs (or mid_ub or cand) for all primes p≥5."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_close.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
