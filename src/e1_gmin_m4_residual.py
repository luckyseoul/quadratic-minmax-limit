#!/usr/bin/env python3
"""
P0 attack: bound |m4 - κ/p²| on |κ|=1 classes ⇒ g_min ≥ L(p).

Pairing (Prop 15.53): g_min = -max |m4| over |κ|=1 4-sets.
Wick on |κ|=1: m4_Wick = κ/p² (since three pairings give Wick = κ/p²).
If |m4 - κ/p²| ≤ (p-4)/(2p²) for p≥5, then
  |m4| ≤ 1/p² + (p-4)/(2p²) = (p-2)/(2p²) = -L(p),
hence g_min ≥ L(p) > T(p) ⇒ bi-tight empty (Prop 15.47).

Also: free Sym0 Φ-maximizer has ambient zero diag (constraints free).

Multi-worker. Writes evidence/e1_gmin_m4_residual.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _blas1() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


def L_p(p: int) -> float:
    return -(p - 2) / (2.0 * p * p)


def T_p(p: int) -> float:
    return -(p - 2) / (p * (2.0 * p - 1))


def _load_Y(p: int) -> np.ndarray:
    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        return load_maxplus(3).astype(np.float64)
    if p == 5:
        return np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    if p == 7:
        return np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    raise ValueError(p)


def _worker_m4_residuals(p: int) -> dict:
    """Enumerate |κ|=1 classes: m4, Wick, residual, vs L(p) budget."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    # Precompute for speed: for each pair of coords, the vector Y[:,i]*Y[:,j]
    # m4(ijkl) = (1/N) sum_y yi yj yk yl = (1/N) (Yi⊙Yj) · (Yk⊙Yl)

    # Collect all 4-sets with |κ|=1
    # For large n=50, binom(50,4)=230300, feasible
    records = []
    max_abs_m4 = 0.0
    max_abs_resid = 0.0
    worst = None

    # Stream 4-sets
    for pts in combinations(range(n), 4):
        i, j, k, l = pts
        kap = int(
            C[i, j] * C[k, l] + C[i, k] * C[j, l] + C[i, l] * C[j, k]
        )
        if abs(kap) != 1:
            continue
        # m4
        f_ij = Y[:, i] * Y[:, j]
        f_kl = Y[:, k] * Y[:, l]
        m4 = float(np.dot(f_ij, f_kl) / N)
        wick = kap / (p * p)
        resid = m4 - wick
        abs_m4 = abs(m4)
        abs_resid = abs(resid)
        if abs_m4 > max_abs_m4:
            max_abs_m4 = abs_m4
        if abs_resid > max_abs_resid:
            max_abs_resid = abs_resid
            worst = {
                "pts": pts,
                "kappa": kap,
                "m4": m4,
                "wick": wick,
                "resid": resid,
            }

    # Also stratify residuals
    # Sample histogram of (m4, resid) for |κ|=1
    # Re-loop is expensive at p=7; store only aggregates from first pass
    # Second pass for histogram buckets
    resid_list = []
    m4_list = []
    # For p=7, binom(50,4)=230300 is OK single-threaded per worker
    for pts in combinations(range(n), 4):
        i, j, k, l = pts
        kap = int(
            C[i, j] * C[k, l] + C[i, k] * C[j, l] + C[i, l] * C[j, k]
        )
        if abs(kap) != 1:
            continue
        f_ij = Y[:, i] * Y[:, j]
        f_kl = Y[:, k] * Y[:, l]
        m4 = float(np.dot(f_ij, f_kl) / N)
        wick = kap / (p * p)
        resid_list.append(m4 - wick)
        m4_list.append(m4)

    resid_arr = np.array(resid_list)
    m4_arr = np.array(m4_list)
    g_min = -float(np.max(np.abs(m4_arr)))
    Lp = L_p(p)
    Tp = T_p(p)
    eps_budget = (p - 4) / (2.0 * p * p) if p > 4 else float("inf")
    # |m4| ≤ 1/p² + |resid| ≤ 1/p² + max|resid|
    m4_ub_from_resid = 1.0 / (p * p) + float(np.max(np.abs(resid_arr)))
    sufficient = m4_ub_from_resid <= -Lp + 1e-12

    return {
        "p": int(p),
        "N": int(N),
        "n": int(n),
        "n_kappa1_sets": int(len(m4_list)),
        "g_min": g_min,
        "L_p": Lp,
        "T_p": Tp,
        "g_min_ge_L": bool(g_min >= Lp - 1e-12),
        "g_min_gt_T": bool(g_min > Tp + 1e-12),
        "max_abs_m4": float(np.max(np.abs(m4_arr))),
        "max_abs_resid": float(np.max(np.abs(resid_arr))),
        "mean_abs_resid": float(np.mean(np.abs(resid_arr))),
        "eps_budget_(p-4)/(2p^2)": eps_budget,
        "resid_le_eps_budget": bool(
            float(np.max(np.abs(resid_arr))) <= eps_budget + 1e-12
        )
        if p > 4
        else None,
        "m4_ub_from_resid": m4_ub_from_resid,
        "ub_implies_L": sufficient,
        "worst_resid": worst,
        "m4_unique_rounded": sorted({round(float(x), 12) for x in m4_arr})[:30],
        "n_unique_m4": len({round(float(x), 12) for x in m4_arr}),
        "resid_unique_rounded": sorted({round(float(x), 12) for x in resid_arr})[:30],
        "n_unique_resid": len({round(float(x), 12) for x in resid_arr}),
        # exact fractions for small cases
        "target_L_abs": (p - 2) / (2.0 * p * p),
        "wick_abs": 1.0 / (p * p),
    }


def _worker_zero_diag_free(p: int) -> dict:
    """Certify: Φ maximizer on Sym0 has ambient zero diag automatically."""
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    U = (Y @ Vp) / np.sqrt(n)
    rng = np.random.default_rng(p + 1)
    A = rng.standard_normal((d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    A /= np.linalg.norm(A) + 1e-30
    for _ in range(100):
        q = np.einsum("bi,ij,bj->b", U, A, U)
        Ph = (U * q[:, None]).T @ U
        Ph = 0.5 * (Ph + Ph.T)
        Ph -= np.trace(Ph) / d * np.eye(d)
        A = Ph / (np.linalg.norm(Ph) + 1e-30)
    B = Vp @ A @ Vp.T
    diagmax = float(np.max(np.abs(np.diag(B))))
    sAs = np.einsum("bi,ij,bj->b", Y @ Vp, A, Y @ Vp)
    E = float(np.mean(sAs**2))
    return {
        "p": int(p),
        "E_max_Sym0": E,
        "ambient_diag_maxabs": diagmax,
        "zero_diag_automatic": bool(diagmax < 1e-12),
        "qn": E,
        "thr_16": 16.0,
        "le_16": bool(E <= 16.0 + 1e-6),
    }


def _worker_pairing_L_algebra(p: int) -> dict:
    """Algebra: 1/p² + (p-4)/(2p²) = (p-2)/(2p²) = -L(p) for p≥5."""
    wick = 1.0 / (p * p)
    eps = (p - 4) / (2.0 * p * p) if p >= 5 else None
    target = (p - 2) / (2.0 * p * p)
    if eps is not None:
        ok = abs(wick + eps - target) < 1e-15
    else:
        ok = None
    return {
        "p": p,
        "wick": wick,
        "eps_budget": eps,
        "L_abs": target,
        "wick_plus_eps_eq_L_abs": ok,
        "L_gt_T": bool(L_p(p) > T_p(p)),
    }


def main() -> None:
    _blas1()
    W = max(2, (os.cpu_count() or 4) - 2)
    print(f"m4_residual: workers={W}", flush=True)

    rows = []
    free = []
    algs = []
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = {}
        for p in (3, 5, 7):
            futs[ex.submit(_worker_m4_residuals, p)] = ("m4", p)
            futs[ex.submit(_worker_zero_diag_free, p)] = ("zd", p)
        for p in (3, 5, 7, 11, 13, 17, 19):
            futs[ex.submit(_worker_pairing_L_algebra, p)] = ("alg", p)
        for fut in as_completed(futs):
            kind, p = futs[fut]
            r = fut.result()
            if kind == "m4":
                rows.append(r)
                print(
                    f"  m4 p={p}: g_min={r['g_min']:.8f} L={r['L_p']:.8f} "
                    f"ge_L={r['g_min_ge_L']} max|resid|={r['max_abs_resid']:.8f} "
                    f"eps_bud={r['eps_budget_(p-4)/(2p^2)']} "
                    f"resid_ok={r['resid_le_eps_budget']} ub_L={r['ub_implies_L']}",
                    flush=True,
                )
            elif kind == "zd":
                free.append(r)
                print(
                    f"  zd p={p}: E={r['E_max_Sym0']:.6f} diag={r['ambient_diag_maxabs']:.2e} "
                    f"auto0={r['zero_diag_automatic']} le16={r['le_16']}",
                    flush=True,
                )
            else:
                algs.append(r)

    rows.sort(key=lambda r: r["p"])
    free.sort(key=lambda r: r["p"])
    algs.sort(key=lambda r: r["p"])

    # Algebra check all p≥5
    alg_ok = all(r["wick_plus_eps_eq_L_abs"] for r in algs if r["p"] >= 5)
    resid_ok = all(
        r["resid_le_eps_budget"] for r in rows if r["p"] >= 5 and r["resid_le_eps_budget"] is not None
    )
    gmin_ok = all(r["g_min_ge_L"] for r in rows if r["p"] >= 5)
    zd_ok = all(r["zero_diag_automatic"] for r in free)

    out = {
        "m4_residuals": rows,
        "zero_diag_free": free,
        "algebra": algs,
        "sufficient_criterion": (
            "For p≥5: if |m4 - κ/p²| ≤ (p-4)/(2p²) on all |κ|=1 4-sets, "
            "then |m4| ≤ (p-2)/(2p²) = -L(p), hence g_min ≥ L(p) > T(p)."
        ),
        "status": (
            f"Pairing residual criterion: algebra_ok={alg_ok}, "
            f"resid≤eps cert p=5,7={resid_ok}, g_min≥L cert={gmin_ok}, "
            f"zero_diag free={zd_ok}. "
            "OPEN: prove |m4-κ/p²|≤(p-4)/(2p²) for all primes p≥5 (or direct |m4|≤-L(p))."
        ),
        "workers": W,
        "gmin_ge_L_certified_p57": gmin_ok,
        "resid_criterion_holds_p57": resid_ok,
        "algebra_ok": alg_ok,
        "zero_diag_automatic": zd_ok,
    }
    path = ROOT / "evidence" / "e1_gmin_m4_residual.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
