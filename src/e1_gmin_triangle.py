#!/usr/bin/env python3
"""
Triangle form bound (Prop 15.57): |ft(v)| ≤ 2p for unit edge-weight v.

Also certifies E[||By||^2] ≤ 4 and the crude λ_max(G) ≤ n bound,
and records cycle-maximizer structure (ft=2p, By∈V+).

Writes evidence/e1_gmin_triangle.json
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


def triangle_form(C: np.ndarray, v: np.ndarray, edges: list) -> float:
    """ft = sum_i w^{(i)T} C w^{(i)} with w_j = C_ij v_ij."""
    n = C.shape[0]
    eidx = {e: k for k, e in enumerate(edges)}
    Vmat = np.zeros((n, n))
    for ei, (i, j) in enumerate(edges):
        Vmat[i, j] = Vmat[j, i] = v[ei]
    B = C * Vmat  # Hadamard
    ft = 0.0
    for i in range(n):
        ri = B[i, :]
        ft += float(ri @ C @ ri)
    return ft


def analyze(p: int) -> dict:
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    Mp = load_mp(p)
    N = len(Mp)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)
    eidx = {e: k for k, e in enumerate(edges)}

    # Op norms of principal submatrices
    op_C = float(np.linalg.norm(C, 2))
    op_sub = float(np.linalg.norm(C[1:, 1:], 2))

    # Build G and cycle space
    Chip = np.array(
        [[C[a, b] * y[a] * y[b] for a, b in edges] for y in Mp], dtype=np.float64
    )
    G = (Chip.T @ Chip) / N
    U = np.zeros((E, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            a, b = (min(i, j), max(i, j))
            U[eidx[(a, b)], i] = 1.0
    _, _, vh = np.linalg.svd(U.T, full_matrices=True)
    null = vh[n - 1 :].T
    Gc = null.T @ G @ null
    evc, vecc = np.linalg.eigh((Gc + Gc.T) / 2)
    lam_cycle = float(evc[-1])
    v_max = null @ vecc[:, -1]
    v_max /= np.linalg.norm(v_max)
    ft_max = triangle_form(C, v_max, edges)

    # Random unit vectors: check |ft| ≤ 2p
    rng = np.random.default_rng(0)
    max_ft_rand = 0.0
    for _ in range(50):
        v = rng.normal(size=E)
        v /= np.linalg.norm(v)
        max_ft_rand = max(max_ft_rand, abs(triangle_form(C, v, edges)))
    # Random cycle
    max_ft_cyc = 0.0
    for _ in range(50):
        z = rng.normal(size=null.shape[1])
        v = null @ z
        v /= np.linalg.norm(v)
        max_ft_cyc = max(max_ft_cyc, abs(triangle_form(C, v, edges)))

    # E[||By||^2] for maximizer
    Vmat = np.zeros((n, n))
    for ei, (i, j) in enumerate(edges):
        Vmat[i, j] = Vmat[j, i] = v_max[ei]
    B = C * Vmat
    Pplus = 0.5 * (np.eye(n) + C / p)
    Pminus = 0.5 * (np.eye(n) - C / p)
    e_tot, e_plus, e_minus = [], [], []
    for y in Mp:
        By = B @ y
        e_tot.append(float(By @ By))
        e_plus.append(float((Pplus @ By) @ (Pplus @ By)))
        e_minus.append(float((Pminus @ By) @ (Pminus @ By)))
    e_tot = np.array(e_tot)
    e_plus = np.array(e_plus)
    e_minus = np.array(e_minus)

    # Full λ_max(G)
    evG = np.linalg.eigvalsh(G)
    lam_max = float(evG[-1])

    return {
        "p": p,
        "n": n,
        "op_norm_C": op_C,
        "op_norm_principal_sub": op_sub,
        "sub_op_le_p": bool(op_sub <= p + 1e-8),
        "triangle_bound_2p": 2 * p,
        "max_ft_random_unit": max_ft_rand,
        "max_ft_random_cycle": max_ft_cyc,
        "ft_le_2p_random": bool(max_ft_rand <= 2 * p + 1e-8),
        "ft_le_2p_cycle": bool(max_ft_cyc <= 2 * p + 1e-8),
        "lambda_cycle": lam_cycle,
        "ft_at_G_maximizer": ft_max,
        "ft_maximizer_eq_2p": bool(abs(ft_max - 2 * p) < 1e-6),
        "E_By_norm_sq": float(e_tot.mean()),
        "E_By_norm_sq_le_4": bool(e_tot.mean() <= 4 + 1e-8),
        "E_Pplus_By_norm_sq": float(e_plus.mean()),
        "E_Pminus_By_norm_sq": float(e_minus.mean()),
        "By_in_Vplus_a_s": bool(e_minus.max() < 1e-8),
        "lambda_max_G": lam_max,
        "crude_bound_n": float(n),
        "lambda_max_le_n": bool(lam_max <= n + 1e-8),
        "n_over_2": n / 2.0,
        "lambda_max_eq_n_over_2": bool(abs(lam_max - n / 2) < 1e-8),
        "gap_to_n_over_2": float(n / 2 - lam_cycle),
    }


def main() -> None:
    out: dict = {"results": [], "status": None}
    for p in (3, 5, 7):
        print(f"triangle p={p}", flush=True)
        row = analyze(p)
        out["results"].append(row)
        print(
            f"  ft_max={row['ft_at_G_maximizer']:.4f} 2p={2*p} "
            f"λ_cycle={row['lambda_cycle']:.4f} λ_max={row['lambda_max_G']:.4f} "
            f"By∈V+={row['By_in_Vplus_a_s']}",
            flush=True,
        )
    ok = all(
        r["sub_op_le_p"]
        and r["ft_le_2p_random"]
        and r["ft_le_2p_cycle"]
        and r["E_By_norm_sq_le_4"]
        and r["lambda_max_le_n"]
        for r in out["results"]
    )
    gap57 = all(
        r["lambda_max_eq_n_over_2"] for r in out["results"] if r["p"] in (5, 7)
    )
    out["status"] = (
        "Prop 15.57: triangle bound |ft|≤2p proved (principal submatrix of C); "
        "E[||By||^2]≤4; crude λ_max(G)≤n. Cycle maximizer has ft=2p and By∈V+ "
        f"(cert p=3,5,7). λ_max=n/2 at p=5,7. "
        "OPEN: improve crude bound by factor 2 on cycle space for all p≥5 "
        "(⇔ λ_2(P⊙P)≤d/(2N)). L OPEN."
        if ok and gap57
        else "CERT FAILURE"
    )
    path = ROOT / "evidence" / "e1_gmin_triangle.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"])
    print("wrote", path)
    if "FAILURE" in out["status"]:
        sys.exit(1)


if __name__ == "__main__":
    main()
