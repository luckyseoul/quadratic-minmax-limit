#!/usr/bin/env python3
"""
Prop 15.61: 16N fourth-moment bound ⇒ λ_cycle ≤ 8 ⇒ spectral gap for p≥5.

Equivalences (proved in solution.md):
  max_{‖B‖_F=1, Tr B=0, B=P_+ B P_+} ∑_{y∈Max+} (y^T B y)^2  ≤ 16 N
  ⇔  λ_2(W) ≤ 4N/d²
  ⇔  λ_cycle ≤ 8   (at the cycle maximizer with ‖B‖_F²=2)

Algebra: for primes p≥5, n/2 = d ≥ 13 > 8, so λ_cycle ≤ 8 ⇒ λ_max(G)=n/2 simple
⇒ bi-tight empty (Prop 15.55).

Certified: equality at p=3 (λ_cycle=8); strict at p=5,7.
OPEN: prove the 16N bound for all primes p≥3 (or p≥5).

Multi-worker ProcessPool. Writes evidence/e1_gmin_16n.json
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

    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    U = (Y @ Vp) / np.sqrt(n)  # N x d, unit rows in V+

    # Power iteration: max ∑ (u^T A u)^2 over Tr A=0, ‖A‖_F=1
    rng = np.random.default_rng(p + 61)
    A = rng.standard_normal((d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    A /= np.linalg.norm(A) + 1e-30
    for _ in range(70):
        q = np.einsum("bi,ij,bj->b", U, A, U)
        Ph = (U * q[:, None]).T @ U
        Ph = 0.5 * (Ph + Ph.T)
        Ph -= np.trace(Ph) / d * np.eye(d)
        nrm = np.linalg.norm(Ph)
        A = Ph / (nrm + 1e-30)

    q = np.einsum("bi,ij,bj->b", U, A, U)
    ray_W = float(np.sum(q**2))  # λ2(W) ≈ max fourth-moment for unit frame
    # sum_y (y^T B y)^2 with B=Vp A Vp.T, ‖B‖_F=1, Tr B=0
    sum_yBy2 = (n**2) * ray_W
    bound_16N = 16.0 * N
    ratio_16N = sum_yBy2 / bound_16N
    bound_ok = bool(sum_yBy2 <= bound_16N + 1e-6 * max(bound_16N, 1.0))

    # Equivalent forms
    thr_4N_d2 = 4.0 * N / (d**2)
    thr_gap_W = N / (2.0 * d)
    lam2_W_ok_4Nd2 = bool(ray_W <= thr_4N_d2 + 1e-8)
    lam2_W_ok_gap = bool(ray_W <= thr_gap_W + 1e-8)

    # λ_cycle = 2 d^2 / N * λ2(W)
    lam_cycle = 2.0 * (d**2) / N * ray_W
    cycle_le_8 = bool(lam_cycle <= 8.0 + 1e-6)
    cycle_le_nhalf = bool(lam_cycle <= n / 2.0 + 1e-6)

    # Algebra: d≥8 ⇒ 4N/d² ≤ N/(2d) ⇒ 16N bound implies gap
    algebra_d_ge_8 = bool(d >= 8)
    algebra_16N_implies_gap = bool(thr_4N_d2 <= thr_gap_W + 1e-12)

    # Equality case check at p=3
    eq_at_3 = None
    if p == 3:
        eq_at_3 = bool(abs(ratio_16N - 1.0) < 1e-6 and abs(lam_cycle - 8.0) < 1e-6)

    return {
        "p": p,
        "n": n,
        "d": d,
        "N": N,
        "sum_yBy2_max": sum_yBy2,
        "bound_16N": bound_16N,
        "ratio_to_16N": ratio_16N,
        "bound_16N_ok": bound_ok,
        "lambda_2_W": ray_W,
        "thr_4N_over_d2": thr_4N_d2,
        "thr_N_over_2d": thr_gap_W,
        "lambda_2_W_le_4N_d2": lam2_W_ok_4Nd2,
        "lambda_2_W_le_gap_thr": lam2_W_ok_gap,
        "lambda_cycle": lam_cycle,
        "lambda_cycle_le_8": cycle_le_8,
        "lambda_cycle_le_n_half": cycle_le_nhalf,
        "d_ge_8": algebra_d_ge_8,
        "algebra_16N_implies_gap": algebra_16N_implies_gap,
        "equality_case_p3": eq_at_3,
        "gap_ok": lam2_W_ok_gap,
    }


def main() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    primes = (3, 5, 7)
    W = min(len(primes), max(2, (os.cpu_count() or 4) - 2))
    print(f"16n_bound: workers={W}", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = {ex.submit(_worker, p): p for p in primes}
        for fut in as_completed(futs):
            row = fut.result()
            rows.append(row)
            print(
                f"  p={row['p']}: ratio_16N={row['ratio_to_16N']:.6f} "
                f"λ_cycle={row['lambda_cycle']:.6f} "
                f"≤8={row['lambda_cycle_le_8']} gap={row['gap_ok']}",
                flush=True,
            )
    rows.sort(key=lambda r: r["p"])
    r3, r5, r7 = rows
    structure_ok = (
        r3["equality_case_p3"]
        and r5["bound_16N_ok"]
        and r7["bound_16N_ok"]
        and r5["lambda_cycle_le_8"]
        and r7["lambda_cycle_le_8"]
        and r5["gap_ok"]
        and r7["gap_ok"]
        and not r3["gap_ok"]
        and r5["algebra_16N_implies_gap"]
        and r7["algebra_16N_implies_gap"]
    )
    out = {
        "results": rows,
        "status": (
            "Prop 15.61: 16N fourth-moment bound (certified p=3 eq, p=5,7 strict) "
            "⇒ λ_cycle≤8 ⇒ spectral gap for all p≥5 (d≥13>8). "
            "Equivalences: max∑(yBy)²≤16N‖B‖²_F (Tr0 on V+) ⇔ λ_2(W)≤4N/d² ⇔ λ_cycle≤8. "
            "OPEN: prove the 16N bound for all primes p≥5 (equality at p=3 is the base case). "
            "If proved, Prop 15.55 closes bi-tight for all p≥5. Deep ND independent. L OPEN."
            if structure_ok
            else "CERT FAILURE"
        ),
        "workers": W,
        "sufficient": (
            "Prove ∑_{y∈Max+} (y^T B y)^2 ≤ 16 N ‖B‖_F² for all B=P_+BP_+, Tr B=0, "
            "all primes p≥5. Then λ_cycle≤8≤n/2 and bi-tight is empty."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_16n.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"])
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
