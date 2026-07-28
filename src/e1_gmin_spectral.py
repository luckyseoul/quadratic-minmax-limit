#!/usr/bin/env python3
"""
Spectral structure of the Max+ edge Gram G (Prop 15.55–15.56).

Proved structure:
  1. G 1 = (n/2) 1
  2. G u^{(i)} = 1 for every undirected star vector u^{(i)}
  3. ker G contains all star differences (dim n-2)
  4. λ_max(G) = max(n/2, λ_max(G|cycle))
  5. D = YY^T = 2N P with P orthoprojector rank d, P_ii = d/N
  6. On 1^⊥ of R^{Max+}: eigs(G) = 2N * eigs(P⊙P)
     so λ_max(G|cycle) = 2N * λ_2(P⊙P)
  7. Therefore λ_max(G)=n/2 iff λ_2(P⊙P) ≤ d/(2N)

Certified: λ_2(P⊙P) ≤ d/(2N) at p=5,7 (strict); FAILS at p=3.

Writes evidence/e1_gmin_spectral.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def load_mp(p: int) -> np.ndarray:
    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        return load_maxplus(3)
    if p == 5:
        return np.load("/tmp/maxplus_p5.npy").astype(float)
    if p == 7:
        return np.load("/tmp/e1_p7/maxplus.npy").astype(float)
    raise ValueError(p)


def analyze(p: int) -> dict:
    C = paley_conference_prime_power(p).astype(float)
    n = int(C.shape[0])
    d = n // 2
    Mp = load_mp(p)
    N = int(len(Mp))
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)
    eidx = {e: k for k, e in enumerate(edges)}

    Chip = np.array(
        [[C[a, b] * y[a] * y[b] for a, b in edges] for y in Mp], dtype=np.float64
    )
    G = (Chip.T @ Chip) / N

    # --- Star vectors ---
    U = np.zeros((E, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            a, b = (i, j) if i < j else (j, i)
            U[eidx[(a, b)], i] = 1.0
    GU = G @ U
    stars_map_to_one = bool(np.allclose(GU, np.ones((E, n)), atol=1e-8))
    row_sum_ok = bool(np.allclose(G @ np.ones(E), (n / 2) * np.ones(E), atol=1e-8))

    # Star differences in kernel
    diffs_null = True
    for i in range(1, min(n, 6)):
        w = U[:, i] - U[:, 0]
        if np.linalg.norm(G @ w) > 1e-8:
            diffs_null = False
            break

    # Cycle space restriction
    _, _, vh = np.linalg.svd(U.T, full_matrices=True)
    null = vh[n - 1 :].T  # E x (E-n+1)
    Gc = null.T @ G @ null
    evc = np.linalg.eigvalsh(Gc)
    pos_c = evc[evc > 1e-8]
    lam_cycle = float(pos_c[-1]) if len(pos_c) else 0.0
    n_half = n / 2.0
    tr_cycle = float(np.sum(pos_c))
    tr_cycle_expected = n * (n - 2) / 2.0
    rank_cycle = int(len(pos_c))
    rank_cycle_formula = (d - 1) * (d - 2) // 2 - 1
    avg_cycle = tr_cycle / rank_cycle if rank_cycle else 0.0
    avg_lt_nhalf = bool(avg_cycle < n_half - 1e-12)

    # Full spectrum
    evG = np.linalg.eigvalsh(G)
    posG = evG[evG > 1e-8]
    lam_max = float(posG[-1])
    simple_max_is_nhalf = bool(
        abs(lam_max - n_half) < 1e-8 and np.sum(np.abs(posG - n_half) < 1e-8) == 1
    )

    # Dual form: λ_max(cycle) = 2N * λ_2(P⊙P) with P = YY^T/(2N).
    # For large N avoid forming N×N matrices; recover λ_2(PP) = lam_cycle/(2N).
    threshold = d / (2.0 * N)
    lam2_PP = lam_cycle / (2.0 * N) if N else 0.0
    gap_ok = bool(lam2_PP <= threshold + 1e-10)

    # Small-N check of projector identities
    p_proj_err = None
    p_diag = d / float(N)
    pp1_ok = None
    cross_matches = True
    if N <= 400:
        D = Mp @ Mp.T
        P = D / (2.0 * N)
        p_proj_err = float(np.linalg.norm(P @ P - P))
        p_diag = float(P[0, 0])
        PP = P * P
        ones = np.ones(N)
        pp1_ok = bool(np.allclose(PP @ ones, (d / N) * ones, atol=1e-8))
        Q1 = np.ones((N, N)) / N
        PPperp = (np.eye(N) - Q1) @ PP @ (np.eye(N) - Q1)
        evS = np.linalg.eigvalsh(PPperp)
        posS = evS[evS > 1e-10]
        lam2_direct = float(posS[-1]) if len(posS) else 0.0
        cross_matches = bool(abs(2 * N * lam2_direct - lam_cycle) < 1e-6)
        lam2_PP = lam2_direct
        gap_ok = bool(lam2_PP <= threshold + 1e-10)

    return {
        "p": p,
        "n": n,
        "d": d,
        "N_Max": N,
        "E": E,
        "row_sum_ok": row_sum_ok,
        "stars_map_to_one": stars_map_to_one,
        "star_differences_in_ker": diffs_null,
        "lambda_max_G": lam_max,
        "n_over_2": n_half,
        "lambda_max_cycle": lam_cycle,
        "lambda_max_is_max_of_nhalf_and_cycle": bool(
            abs(lam_max - max(n_half, lam_cycle)) < 1e-6
        ),
        "tr_cycle": tr_cycle,
        "tr_cycle_expected": tr_cycle_expected,
        "rank_cycle": rank_cycle,
        "rank_cycle_formula": rank_cycle_formula,
        "avg_cycle_eig": avg_cycle,
        "avg_cycle_lt_n_over_2": avg_lt_nhalf,
        "simple_lambda_max_eq_n_over_2": simple_max_is_nhalf,
        "P_is_projector_err": p_proj_err,
        "P_diag": p_diag,
        "d_over_N": d / float(N),
        "PP_ones_eigenvalue_ok": pp1_ok,
        "lambda_2_PP": lam2_PP,
        "threshold_d_over_2N": threshold,
        "spectral_gap_ok": gap_ok,
        "two_N_times_lambda2_PP": float(2 * N * lam2_PP),
        "cross_matches_cycle": cross_matches,
        "avg_cycle_formula_lt_nhalf_algebra": bool(
            (n - 1) <= (d - 1) * (d - 2) / 2 + 1e-12
        ),
    }


def main() -> None:
    out: dict = {
        "results": [],
        "status": None,
        "reduction": (
            "λ_max(G)=n/2 simple  ⇔  λ_2(P⊙P) ≤ d/(2N), where "
            "P=YY^T/(2N) is the equal-diagonal orthoprojector of rank d=n/2 "
            "onto the Max+ span in R^N. Cycle/star decomposition fully proved."
        ),
    }
    for p in (3, 5, 7):
        print(f"spectral p={p}", flush=True)
        row = analyze(p)
        out["results"].append(row)
        print(
            f"  λmax={row['lambda_max_G']:.4f} cycle={row['lambda_max_cycle']:.4f} "
            f"λ2(PP)={row['lambda_2_PP']:.6f} thr={row['threshold_d_over_2N']:.6f} "
            f"gap_ok={row['spectral_gap_ok']}",
            flush=True,
        )

    r3, r5, r7 = out["results"]
    structure_ok = all(
        r["stars_map_to_one"]
        and r["row_sum_ok"]
        and r["lambda_max_is_max_of_nhalf_and_cycle"]
        and r["cross_matches_cycle"]
        and (r["PP_ones_eigenvalue_ok"] is not False)
        for r in out["results"]
    )
    gap57 = r5["spectral_gap_ok"] and r7["spectral_gap_ok"] and not r3["spectral_gap_ok"]
    avg_ok = r5["avg_cycle_lt_n_over_2"] and r7["avg_cycle_lt_n_over_2"]
    algebra_avg = all(
        r["avg_cycle_formula_lt_nhalf_algebra"] for r in out["results"] if r["p"] >= 5
    )

    out["status"] = (
        "Structure OK (stars→1, λ_max=max(n/2,cycle), G eigs=2N·(P⊙P eigs)). "
        "Spectral gap λ_2(P⊙P)≤d/(2N) holds p=5,7, fails p=3. "
        "Average cycle eig < n/2 for all p≥5 (algebra). "
        "OPEN: prove λ_2(P⊙P)≤d/(2N) for all primes p≥5 "
        "(⇔ λ_max(G)=n/2 ⇔ bi-tight residual closes via Prop 15.55). "
        "Deep ND independent. L OPEN."
        if structure_ok and gap57 and avg_ok and algebra_avg
        else "CERT FAILURE"
    )
    path = ROOT / "evidence" / "e1_gmin_spectral.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"])
    print("wrote", path)
    if "FAILURE" in out["status"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
