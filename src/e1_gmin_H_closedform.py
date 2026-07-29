#!/usr/bin/env python3
"""
Discover closed form for λ_max(Φ|Z) / residual as function of p,d only.
If the form is ≤ 6+2H(p) with a proved inequality, H is proved.

Also: sample many A in Z, compute E[(sAs)^2], compare to O(d) invariants
Tr(A^2), Tr(A^3), Tr(A^4), sum_i (r_i^T A^2 r_i), etc.

Multi-worker. Writes evidence/e1_gmin_H_closedform.json
"""
from __future__ import annotations

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def _blas1() -> None:
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


def _load_Y(p: int) -> np.ndarray:
    if p == 3:
        from e1_gmin_cr_classify import load_maxplus

        return load_maxplus(3).astype(np.float64)
    if p == 5:
        return np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    if p == 7:
        return np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    raise ValueError(p)


def _A_from_vec(v, d):
    A = np.zeros((d, d))
    k = 0
    for i in range(d):
        for j in range(i, d):
            A[i, j] = A[j, i] = v[k]
            k += 1
    return A


def _worker_sample_invariants(p: int) -> dict:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power

    Y = _load_Y(p)
    C = paley_conference_prime_power(p).astype(np.float64)
    N, n = Y.shape
    d = n // 2
    ew, EV = np.linalg.eigh(C)
    Vp = EV[:, ew > 0]
    S = Y @ Vp
    rows = Vp

    # Nullspace of zero-diag + Tr
    dimS = d * (d + 1) // 2
    M = np.zeros((n + 1, dimS))
    for k in range(dimS):
        e = np.zeros(dimS)
        e[k] = 1.0
        A = _A_from_vec(e, d)
        B = Vp @ A @ Vp.T
        M[:n, k] = np.diag(B)
    kk = 0
    for i in range(d):
        for j in range(i, d):
            if i == j:
                M[n, kk] = 1.0
            kk += 1
    _u, s, vh = np.linalg.svd(M, full_matrices=True)
    rank = int(np.sum(s > 1e-10))
    null = vh[rank:].T
    nullity = null.shape[1]
    rng = np.random.default_rng(p + 99)

    rows_data = []
    for trial in range(80):
        z = rng.standard_normal(nullity)
        c = null @ z
        A = _A_from_vec(c, d)
        nrm = np.linalg.norm(A)
        if nrm < 1e-14:
            continue
        A = A / nrm
        sAs = np.einsum("bi,ij,bj->b", S, A, S)
        E = float(np.mean(sAs**2))
        trA2 = float(np.trace(A @ A))
        trA3 = float(np.trace(A @ A @ A))
        trA4 = float(np.trace(np.linalg.matrix_power(A, 4)))
        # sum_i (r_i^T A^2 r_i)
        A2 = A @ A
        sum_rA2r = float(np.sum([rows[i] @ A2 @ rows[i] for i in range(n)]))
        # sum_i (r_i^T A r_i)^2 = 0 by constraint
        # sum_{i≠j} (r_i^T A r_j)^2 = ‖A‖F^2 = 1
        # Frobenius of A composed with C on ambient
        B = Vp @ A @ Vp.T
        trCB2 = float(np.trace(C @ B @ B))  # = p
        trB4 = float(np.trace(np.linalg.matrix_power(B, 4)))
        sumB4 = float(np.sum(B**4))
        rows_data.append(
            {
                "E": E,
                "trA2": trA2,
                "trA3": trA3,
                "trA4": trA4,
                "sum_rA2r": sum_rA2r,
                "trB4": trB4,
                "sumB4": sumB4,
                "residual": E - 8.0,
            }
        )

    # Also maximizer
    U = S / np.sqrt(n)
    A = rng.standard_normal((d, d))
    A = 0.5 * (A + A.T)
    A -= np.trace(A) / d * np.eye(d)
    A /= np.linalg.norm(A) + 1e-30
    for _ in range(80):
        q = np.einsum("bi,ij,bj->b", U, A, U)
        Ph = (U * q[:, None]).T @ U
        Ph = 0.5 * (Ph + Ph.T)
        Ph -= np.trace(Ph) / d * np.eye(d)
        A = Ph / (np.linalg.norm(Ph) + 1e-30)
    sAs = np.einsum("bi,ij,bj->b", S, A, S)
    Emax = float(np.mean(sAs**2))
    trA4_max = float(np.trace(np.linalg.matrix_power(A, 4)))

    # Linear fit residual ~ a*(trA4 - c) + b*(sum_rA2r - c2) + ...
    X = []
    y = []
    for r in rows_data:
        X.append([1.0, r["trA4"], r["sum_rA2r"], r["trB4"], r["sumB4"], abs(r["trA3"])])
        y.append(r["E"])
    X = np.asarray(X)
    y = np.asarray(y)
    coef, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ coef
    rms = float(np.sqrt(np.mean((pred - y) ** 2)))
    max_err = float(np.max(np.abs(pred - y)))

    H = (p + 2) ** 2 / d
    budget = (p + 1) * (p + 7) / d
    return {
        "p": int(p),
        "d": int(d),
        "N": int(N),
        "nullity": int(nullity),
        "Emax": Emax,
        "H_qn": 6 + 2 * H,
        "residual_max": Emax - 8.0,
        "budget": budget,
        "trA4_at_max": trA4_max,
        "lstsq_coef": [float(c) for c in coef],
        "lstsq_rms": rms,
        "lstsq_max_err": max_err,
        "identity_in_span": bool(max_err < 1e-4),
        "E_min_sample": float(min(r["E"] for r in rows_data)),
        "E_max_sample": float(max(r["E"] for r in rows_data)),
        "residual_min_sample": float(min(r["residual"] for r in rows_data)),
        "residual_max_sample": float(max(r["residual"] for r in rows_data)),
    }


def _worker_formula_screen(p: int) -> dict:
    """Screen closed forms for Emax = λ_max(Φ|Z)."""
    _blas1()
    # Known values
    known = {
        3: 16.0,
        5: 176 / 13,
        7: 10.56234718826407,  # from gen.eig
    }
    if p not in known:
        return {"p": p, "skipped": True}
    E = known[p]
    d = (p * p + 1) // 2
    n = 2 * d
    candidates = {
        "6+2*(p+2)^2/d": 6 + 2 * (p + 2) ** 2 / d,
        "16": 16.0,
        "8": 8.0,
        "4*p": 4 * p,
        "2*(p+1)": 2 * (p + 1),
        "n/d * something": 2.0,
        "8+2*(p+1)*(p+7)/d": 8 + 2 * (p + 1) * (p + 7) / d,  # =6+2H? no 8+budget=8+2(H-1)=6+2H
        "8+(p+1)*(p+7)/d": 8 + (p + 1) * (p + 7) / d,  # =6+2H
        "2*n/(d-1)": 2 * n / (d - 1) if d > 1 else None,
        "4*d/(d-2)": 4 * d / (d - 2) if d > 2 else None,
        "12+4*(d-2)/d": 12 + 4 * (d - 2) / d,
        "6+4*(p+1)/d": 6 + 4 * (p + 1) / d,
        "6+2*(p+1)^2/d": 6 + 2 * (p + 1) ** 2 / d,
        "6+2*(p+3)^2/d": 6 + 2 * (p + 3) ** 2 / d,
        # p=7 measured not H; try other
        "8+4*(d-4)/d": 8 + 4 * (d - 4) / d,
        "8+2*(d-4)/(d-1)": 8 + 2 * (d - 4) / (d - 1) if d > 1 else None,
        "2*(d+3)/(d-1)": 2 * (d + 3) / (d - 1) if d > 1 else None,
        "4*(d+1)/(d-1)": 4 * (d + 1) / (d - 1) if d > 1 else None,
        "16*d/(d+2)": 16 * d / (d + 2),  # 2×sphere full Q/N
        "8*d/(d+2)*2": 16 * d / (d + 2),
    }
    errs = {}
    for name, val in candidates.items():
        if val is None:
            continue
        errs[name] = abs(E - val)
    best = min(errs.items(), key=lambda kv: kv[1])
    f = Fraction(E).limit_denominator(2000)
    return {
        "p": p,
        "E": E,
        "frac": f"{f.numerator}/{f.denominator}",
        "frac_err": abs(E - float(f)),
        "best": {"name": best[0], "err": best[1], "value": candidates[best[0]]},
        "H_qn": 6 + 2 * (p + 2) ** 2 / d,
        "equals_H_qn": bool(abs(E - (6 + 2 * (p + 2) ** 2 / d)) < 1e-6),
        "top5": sorted(errs.items(), key=lambda kv: kv[1])[:8],
    }


def main() -> None:
    _blas1()
    W = max(2, (os.cpu_count() or 4) - 2)
    print(f"H_closedform: workers={W}", flush=True)

    invs, forms = [], []
    with ProcessPoolExecutor(max_workers=W) as ex:
        futs = {}
        for p in (3, 5, 7):
            futs[ex.submit(_worker_sample_invariants, p)] = ("inv", p)
            futs[ex.submit(_worker_formula_screen, p)] = ("form", p)
        for fut in as_completed(futs):
            kind, p = futs[fut]
            r = fut.result()
            if kind == "inv":
                invs.append(r)
                print(
                    f"  inv p={p}: Emax={r['Emax']:.6f} H_qn={r['H_qn']:.6f} "
                    f"resid={r['residual_max']:.6f} budget={r['budget']:.6f} "
                    f"lstsq_err={r['lstsq_max_err']:.3e} in_span={r['identity_in_span']}",
                    flush=True,
                )
            else:
                forms.append(r)
                print(
                    f"  form p={p}: E={r['E']:.6f}={r['frac']} best={r['best']} "
                    f"eq_H={r['equals_H_qn']}",
                    flush=True,
                )

    invs.sort(key=lambda r: r["p"])
    forms.sort(key=lambda r: r["p"])
    out = {
        "invariants": invs,
        "formula_screen": forms,
        "status": (
            "Closed-form screen for λ_max(Φ|Z): H_qn exact at p=3,5; "
            "p=7 strictly below. No low-degree O(d) invariant identity for E. "
            "Uniform H proof still OPEN."
        ),
        "workers": W,
    }
    path = ROOT / "evidence" / "e1_gmin_H_closedform.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
