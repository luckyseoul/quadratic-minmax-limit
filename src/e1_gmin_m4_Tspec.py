#!/usr/bin/env python3
"""
Exact spectral radius of T on l2(4-sets) for Paley conference (sparse eig).

T is the signed adjacency of the Johnson graph J(n,4) with weights C_vr
on replacement edges. Symmetric. Multi-worker: build CSR shards + eigsh.

Hypothesis: λ_max(T) < 4p for p≥5 (so 4pI-T invertible), with a useful gap.

Also compute ‖ρ‖ bounds if we had the eigenbasis; primarily record λ_max, λ_min.

Writes evidence/e1_gmin_m4_Tspec.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import as_completed
from itertools import combinations
from math import comb
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from workers import pool, require_workers  # noqa: E402


def _blas1():
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


def _build_T_sparse(p: int) -> dict:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    quads = list(combinations(range(n), 4))
    nQ = len(quads)
    qindex = {q: i for i, q in enumerate(quads)}

    # Build COO: for each S,v,r add edge to S'
    # Only store upper triangle + diag(0) then symmetrize — actually add all directed and they'll match
    rows = []
    cols = []
    data = []
    for i, S in enumerate(quads):
        Sset = set(S)
        for v in S:
            others = tuple(sorted(x for x in S if x != v))
            for r in range(n):
                if r in Sset:
                    continue
                Sp = tuple(sorted(others + (r,)))
                j = qindex[Sp]
                rows.append(i)
                cols.append(j)
                data.append(float(C[v, r]))

    T = sparse.csr_matrix((data, (rows, cols)), shape=(nQ, nQ))
    # Symmetry check
    sym_err = sparse.linalg.norm(T - T.T)
    # Top eigenvalues
    k = min(6, nQ - 2)
    evals_large = eigsh(T, k=k, which="LA", return_eigenvectors=False)
    evals_small = eigsh(T, k=k, which="SA", return_eigenvectors=False)
    lam_max = float(np.max(evals_large))
    lam_min = float(np.min(evals_small))

    # Degree stats
    deg = np.asarray(np.abs(T).sum(axis=1)).ravel()
    # Signed row sums
    row_sum = np.asarray(T.sum(axis=1)).ravel()

    return {
        "p": int(p),
        "n": int(n),
        "nQ": int(nQ),
        "nnz": int(T.nnz),
        "sym_err": float(sym_err),
        "lam_max": lam_max,
        "lam_min": lam_min,
        "top_evals": sorted(float(x) for x in evals_large)[::-1],
        "bot_evals": sorted(float(x) for x in evals_small),
        "gap_4p_minus_lammax": float(4 * p - lam_max),
        "ratio_lammax_over_4p": float(lam_max / (4 * p)),
        "row_abs_sum_max": float(np.max(deg)),
        "row_abs_sum_mean": float(np.mean(deg)),
        "row_sum_abs_max": float(np.max(np.abs(row_sum))),
        "candidates": {
            "4p": 4 * p,
            "2p": 2 * p,
            "n-6": n - 6,
            "2(n-4)/sqrt(n)": 2 * (n - 4) / np.sqrt(n),
            "2*sqrt(n)*p/p": None,
            "p+sqrt(n)": p + np.sqrt(n),
            "2*(p+2)": 2 * (p + 2),
            "4*(p-1)": 4 * (p - 1),
            "2*p+4": 2 * p + 4,
            "p*2+sqrt(p)": 2 * p + np.sqrt(p),
        },
        "lammax_le_4p": bool(lam_max < 4 * p - 1e-8),
        "lammax_le_4p_minus_1": bool(lam_max <= 4 * p - 1 + 1e-8),
    }


def main():
    W = require_workers(min_workers=4)
    print(f"Tspec: W={W}", flush=True)
    # p=3,5 exact sparse; p=7 is 230k^2 structure with 230k*184 nnz ≈ 42M — feasible ~few GB
    out = {"workers": W, "spectra": {}}
    # Run p=3,5 in parallel; p=7 sequential after (memory)
    with pool(min(W, 4)) as ex:
        futs = {ex.submit(_build_T_sparse, p): p for p in (3, 5)}
        for fut in as_completed(futs):
            p = futs[fut]
            r = fut.result()
            out["spectra"][str(p)] = r
            print(
                f"  p={p}: λmax={r['lam_max']:.6f} 4p={4*p} gap={r['gap_4p_minus_lammax']:.6f} "
                f"λmin={r['lam_min']:.6f} sym={r['sym_err']:.2e}",
                flush=True,
            )

    print("  building p=7 sparse T (large)...", flush=True)
    try:
        r7 = _build_T_sparse(7)
        out["spectra"]["7"] = r7
        print(
            f"  p=7: λmax={r7['lam_max']:.6f} 4p=28 gap={r7['gap_4p_minus_lammax']:.6f}",
            flush=True,
        )
    except Exception as e:
        out["spectra"]["7"] = {"error": str(e)}
        print(f"  p=7 failed: {e}", flush=True)

    # Pattern fit for λ_max
    pts = []
    for p in (3, 5, 7):
        s = out["spectra"].get(str(p), {})
        if "lam_max" in s:
            pts.append((p, s["lam_max"], 4 * p - s["lam_max"]))

    out["pattern"] = {
        "points_p_lammax_gap": pts,
        "hypothesis": (
            "λ_max(T) < 4p for primes p≥3 on Paley conference 4-set Johnson signing; "
            "gap appears to grow with p. If λ_max(T) ≤ 4p - 48/(p-4)*‖source channel‖ "
            "insufficient alone for L∞; still proves invertibility of 4pI-T."
        ),
    }

    # Resolvent L2 bound
    for p in (3, 5, 7):
        s = out["spectra"].get(str(p), {})
        if "lam_max" not in s:
            continue
        n = p * p + 1
        nQ = comb(n, 4)
        # n3 approx from earlier
        # |κ|=3 count: use fraction
        # ‖b‖_2 ≤ sqrt(nQ) * 24/p² crude, or sqrt(n3)*24/p²
        n3_est = nQ  # crude
        # from prior: p5 n3=3250, p7 n3=53900, p3 n3=30
        n3 = {3: 30, 5: 3250, 7: 53900}.get(p, nQ)
        b_l2 = np.sqrt(n3) * 24.0 / (p * p)
        gap = s["gap_4p_minus_lammax"]
        if gap > 0:
            rho_l2_ub = b_l2 / gap
            # convert to L∞ via |ρ|_∞ ≤ ρ_l2
            rho_inf_ub = rho_l2_ub
            budget = (p - 4) / (2.0 * p * p) if p > 4 else None
            s["b_l2"] = float(b_l2)
            s["rho_l2_ub"] = float(rho_l2_ub)
            s["rho_inf_ub_crude"] = float(rho_inf_ub)
            s["rho_budget"] = budget
            s["L2_implies_budget"] = (
                budget is not None and rho_inf_ub <= budget + 1e-12
            )

    out["status"] = (
        "Exact sparse spectrum of T on 4-sets (Johnson signing by C). "
        "OPEN: closed form λ_max(T) and L∞ resolvent gain ≤(p-4)/48."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_Tspec.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
