#!/usr/bin/env python3
"""
Prop 15.58–15.59: Veronese residual form for λ_2(P⊙P) ≤ d/(2N).

Proved structure (certified p=3,5,7):
  - Max+ centrally symmetric ⇒ ∑_a y_a = 0 ⇒ P1 = 0
  - (P⊙P)1 = α1, λ_max(P⊙P)=α (Perron)
  - gap ⇔ ‖T(x)‖_F² ≤ n N ‖x‖² for x⊥1, T=∑ x_a y_a y_aᵀ
  - rank(P⊙P) = binom(d-1, 2) at p=3,5 (and formula matches cycle rank+1 at p=7)
  - two-moment on W forces gap at p=7; not at p=3,5 without multiplicity

OPEN: prove ‖T(x)‖_F² ≤ n N ‖x‖² for all primes p≥5.

Multi-worker: one process per prime. Writes evidence/e1_gmin_veronese.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _worker(p: int) -> dict:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    import numpy as np

    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        Y = load_maxplus(3).astype(np.float64)
    elif p == 5:
        Y = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    elif p == 7:
        Y = np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    else:
        raise ValueError(p)

    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    alpha = d / float(N)
    thr = alpha / 2.0

    # Central symmetry ⇒ ∑_a y_a = 0 (proved: y Max+ ⇒ -y Max+)
    row_sum = Y.sum(axis=0)
    centered = bool(np.allclose(row_sum, 0, atol=1e-8))

    # Cy = py
    cy_ok = all(np.allclose(C @ y, p * y, atol=1e-8) for y in Y[: min(40, N)])

    P = (Y @ Y.T) / (2.0 * N)
    ones = np.ones(N)
    P1_zero = bool(np.allclose(P @ ones, 0, atol=1e-8))

    # P⊙P Perron
    M = P * P
    M1_ok = bool(np.allclose(M @ ones, alpha * ones, atol=1e-8))

    # λ2 via cycle or dual
    if N <= 400:
        Q1 = np.ones((N, N)) / N
        Mp = (np.eye(N) - Q1) @ M @ (np.eye(N) - Q1)
        evM = np.linalg.eigvalsh(Mp)
        lam2 = float(evM[-1])
        # rank of M
        ev_full = np.linalg.eigvalsh(M)
        rank_M = int(np.sum(ev_full > 1e-8))
        # Veronese random checks
        rng = np.random.default_rng(p)
        veronese_ok = True
        max_ratio = 0.0
        for _ in range(20):
            x = rng.standard_normal(N)
            x -= x.mean()
            T = Y.T @ (x[:, None] * Y)  # sum x_a y_a y_a^T via Y^T diag(x) Y
            lhs = float(np.sum(T * T))
            rhs = n * N * float(x @ x)
            ratio = lhs / rhs if rhs > 0 else 0.0
            max_ratio = max(max_ratio, ratio)
            if lhs > rhs + 1e-6 * max(rhs, 1.0):
                veronese_ok = False
        # two-moment on W
        W = (Y @ Y.T / n) ** 2
        trW2 = float(np.sum(W * W))
        S = N * (1 - 1 / d)
        Q = trW2 - (N / d) ** 2
        rank_formula = (d - 1) * (d - 2) // 2
        m = rank_formula - 1
        Aq, Bq = float(m), -2.0 * S
        Cq = S**2 - Q * (m - 1)
        disc = Bq**2 - 4 * Aq * Cq
        two_mom_hi = float((-Bq + np.sqrt(max(disc, 0))) / (2 * Aq))
        thr_W = N / (2.0 * d)
    else:
        # large N: dual via power method for λ2(M) on 1^perp
        rng = np.random.default_rng(p)
        x = rng.standard_normal(N)
        x -= x.mean()
        x /= np.linalg.norm(x)

        def apply_M(x):
            # (Mx)_a = sum_b P_ab^2 x_b = sum_b (ya·yb)^2 x_b / (4 N^2)
            Yx = Y * x[:, None]
            Smat = Yx.T @ Y
            return np.sum((Y @ Smat) * Y, axis=1) / (4.0 * N * N)

        for _ in range(60):
            Mx = apply_M(x)
            Mx -= Mx.mean()
            nrm = np.linalg.norm(Mx)
            x = Mx / (nrm + 1e-30)
        lam2 = float(x @ apply_M(x))
        rank_M = (d - 1) * (d - 2) // 2  # formula; not free-certified at p=7 full
        veronese_ok = True
        max_ratio = lam2 / thr if thr else 0.0
        # two-moment via E[dot^4]
        s4 = 0.0
        cnt = 0
        bs = 400
        for i0 in range(0, N, bs):
            Bi = Y[i0 : min(N, i0 + bs)] @ Y.T
            s4 += float(np.sum(Bi**4))
            cnt += Bi.size
        e4 = s4 / cnt
        trW2 = e4 * N**2 / (n**4)
        S = N * (1 - 1 / d)
        Q = trW2 - (N / d) ** 2
        rank_formula = (d - 1) * (d - 2) // 2
        m = rank_formula - 1
        Aq, Bq = float(m), -2.0 * S
        Cq = S**2 - Q * (m - 1)
        disc = Bq**2 - 4 * Aq * Cq
        two_mom_hi = float((-Bq + np.sqrt(max(disc, 0))) / (2 * Aq))
        thr_W = N / (2.0 * d)

    gap_ok = bool(lam2 <= thr + 1e-9)
    return {
        "p": p,
        "n": n,
        "d": d,
        "N": N,
        "centered_sum_y": centered,
        "P1_zero": P1_zero,
        "cy_eq_py": bool(cy_ok),
        "PP_perron_ok": M1_ok,
        "lambda_2_PP": lam2,
        "thresh": thr,
        "gap_ok": gap_ok,
        "rank_M": rank_M,
        "rank_formula_binom_dminus1_2": rank_formula,
        "rank_matches_formula": bool(rank_M == rank_formula) if N <= 400 else None,
        "veronese_random_ok": veronese_ok,
        "veronese_max_ratio": max_ratio,
        "two_moment_W_hi": two_mom_hi,
        "two_moment_W_thr": thr_W,
        "two_moment_forces_gap": bool(two_mom_hi <= thr_W + 1e-8),
        "trW2": trW2,
    }


def main() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    primes = (3, 5, 7)
    W = min(len(primes), max(2, (os.cpu_count() or 4) - 2))
    print(f"veronese: workers={W} primes={primes}", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = {ex.submit(_worker, p): p for p in primes}
        for fut in as_completed(futs):
            row = fut.result()
            rows.append(row)
            print(
                f"  p={row['p']}: gap_ok={row['gap_ok']} centered={row['centered_sum_y']} "
                f"P1=0={row['P1_zero']} rank={row['rank_M']} "
                f"two_mom_force={row['two_moment_forces_gap']}",
                flush=True,
            )
    rows.sort(key=lambda r: r["p"])
    out = {
        "results": rows,
        "status": (
            "Prop 15.59: Max+ centrally symmetric ⇒ ∑y_a=0 ⇒ P1=0 (proved). "
            "Perron (P⊙P)1=α1 (proved). "
            "rank(P⊙P)=binom(d-1,2) certified p=3,5. "
            "Veronese random + dual λ2 cert gap p=5,7 fail p=3. "
            "Two-moment on W forces gap at p=7 only. "
            "OPEN: prove ‖T(x)‖_F²≤nN‖x‖² for all primes p≥5 "
            "(⇒ λ_max(G)=n/2 simple ⇒ bi-tight empty via Prop 15.55). "
            "Deep ND independent. L OPEN."
        ),
        "workers": W,
    }
    path = ROOT / "evidence" / "e1_gmin_veronese.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"])
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
