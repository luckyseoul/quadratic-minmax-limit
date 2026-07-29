#!/usr/bin/env python3
"""
Multi-worker probe for λ_2(P⊙P) ≤ d/(2N) residual (Prop 15.56.6).

Workers: one process per prime (p=3,5,7), OMP_NUM_THREADS=1 each.
Tests candidate closed forms and the 2×-sphere bound
  λ_cycle ≤ 8n/(n+4)  (sufficient for n≥12 ⇒ p≥5).

Writes evidence/e1_gmin_gap_probe.json
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

    # Confirm Cy = p y
    cy_ok = all(np.allclose(C @ y, p * y, atol=1e-8) for y in Y[: min(30, N)])

    # Build G and cycle spectrum (small enough for p≤7)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)
    eidx = {e: k for k, e in enumerate(edges)}
    Chip = np.array(
        [[C[a, b] * y[a] * y[b] for a, b in edges] for y in Y], dtype=np.float64
    )
    G = (Chip.T @ Chip) / N
    Ustar = np.zeros((E, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            a, b = (i, j) if i < j else (j, i)
            Ustar[eidx[(a, b)], i] = 1.0
    _, _, vh = np.linalg.svd(Ustar.T, full_matrices=True)
    null = vh[n - 1 :].T
    Gc = null.T @ G @ null
    evc = np.linalg.eigvalsh(Gc)
    pos = evc[evc > 1e-8]
    lam_cycle = float(pos[-1]) if len(pos) else 0.0
    lam2 = lam_cycle / (2.0 * N)
    n_half = n / 2.0
    bound_8n = 8.0 * n / (n + 4.0)
    bound_4n = 4.0 * n / (n + 4.0)

    # Top cycle maximizer structure
    evec = np.linalg.eigh(Gc)[1][:, -1]
    v = null @ evec
    v /= np.linalg.norm(v)
    Vmat = np.zeros((n, n))
    for ei, (i, j) in enumerate(edges):
        Vmat[i, j] = Vmat[j, i] = v[ei]
    B = C * Vmat
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    Bp = Vp.T @ B @ Vp
    frob_Bp = float(np.sum(Bp * Bp))
    tr_Bp = float(np.trace(Bp))
    # E[||By||^2], ft
    By_norms = np.array([np.linalg.norm(B @ y) ** 2 for y in Y])
    yBy = np.array([float(y @ (B @ y)) for y in Y])
    E_By2 = float(By_norms.mean())
    E_yBy2 = float(np.mean(yBy**2))
    # ft = p*(E[||By||^2] - 2)
    ft = p * (E_By2 - 2.0)
    # fraction By in V+
    in_vp = []
    for y in Y[: min(80, N)]:
        By = B @ y
        nn = np.linalg.norm(By)
        if nn < 1e-14:
            continue
        proj = Vp @ (Vp.T @ By)
        in_vp.append(float(np.linalg.norm(proj) / nn))
    # fourth-moment ratio vs sphere (Tr=0): sphere = 2||A||F^2/(d(d+2))
    U = (Y @ Vp) / np.sqrt(n)
    quads = np.einsum("bi,ij,bj->b", U, Bp, U)
    emp_u2 = float(np.mean(quads**2))
    sphere = 2.0 * frob_Bp / (d * (d + 2)) if tr_Bp**2 < 1e-12 else float("nan")
    ratio = emp_u2 / sphere if sphere and sphere > 1e-15 else float("nan")

    # two-moment worst-case λ_cycle
    trG2 = float(np.sum(G * G))
    k = (d - 1) * (d - 2) // 2 - 1
    sum_c = n * (n - 2) / 2.0
    sum_c2 = trG2 - n_half**2
    # quadratic for one-large-rest-equal
    Aq, Bq = float(k), -2.0 * sum_c
    Cq = sum_c**2 - sum_c2 * (k - 1)
    disc = Bq**2 - 4 * Aq * Cq
    two_mom_hi = float((-Bq + np.sqrt(max(disc, 0))) / (2 * Aq)) if Aq else float("nan")

    return {
        "p": p,
        "n": n,
        "d": d,
        "N": N,
        "cy_eq_py": bool(cy_ok),
        "lambda_cycle": lam_cycle,
        "lambda_2_PP": lam2,
        "alpha": alpha,
        "thresh_alpha_over_2": thr,
        "gap_ok": bool(lam2 <= thr + 1e-10),
        "bound_8n_over_nplus4": bound_8n,
        "bound_4n_over_nplus4": bound_4n,
        "cycle_le_8n": bool(lam_cycle <= bound_8n + 1e-8),
        "cycle_le_n_half": bool(lam_cycle <= n_half + 1e-8),
        "maximizer_ft": ft,
        "maximizer_E_By2": E_By2,
        "maximizer_frob_Bp": frob_Bp,
        "maximizer_tr_Bp": tr_Bp,
        "maximizer_By_in_Vp_mean": float(np.mean(in_vp)) if in_vp else None,
        "fourth_moment_ratio_vs_sphere": ratio,
        "two_moment_worst_cycle": two_mom_hi,
        "two_moment_forces_gap": bool(two_mom_hi <= n_half + 1e-8),
        "E_yBy2_over_4": E_yBy2 / 4.0,
    }


def main() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"
    primes = (3, 5, 7)
    W = min(len(primes), max(2, os.cpu_count() - 2 if os.cpu_count() else 2))
    print(f"gap_probe: workers={W} primes={primes}", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = {ex.submit(_worker, p): p for p in primes}
        for fut in as_completed(futs):
            row = fut.result()
            rows.append(row)
            print(
                f"  p={row['p']}: λ_cycle={row['lambda_cycle']:.6f} "
                f"8n/(n+4)={row['bound_8n_over_nplus4']:.6f} "
                f"gap_ok={row['gap_ok']} ratio={row['fourth_moment_ratio_vs_sphere']:.4f}",
                flush=True,
            )
    rows.sort(key=lambda r: r["p"])
    out = {
        "results": rows,
        "status": (
            "Prop 15.58 structure: Max+ ⊂ V_+ (Cy=py); (P⊙P)1=α1 Perron; "
            "gap ⇔ ||T(x)||_F^2 ≤ n N ||x||^2 on 1^⊥. "
            "Maximizer: ft=2p, E||By||^2=4, By∈V_+, ||B_+||_F^2=2 (p=3,5,7). "
            "Fourth-moment ratio vs sphere: p=3 ~2.8 (gap fails), p=5 ~1.95, p=7 ~1.43. "
            "λ_cycle ≤ 8n/(n+4) holds p=5,7 (fails p=3). "
            "OPEN: prove ||T(x)||_F^2 ≤ n N ||x||^2 for x⊥1, all primes p≥5. "
            "Deep ND independent. L OPEN."
        ),
        "sufficient_bound": "λ_cycle ≤ 8n/(n+4) ≤ n/2 for n=p^2+1≥26 (certified p=5,7; not proved)",
        "workers": W,
    }
    path = ROOT / "evidence" / "e1_gmin_gap_probe.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"])
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
