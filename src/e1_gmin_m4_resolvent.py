#!/usr/bin/env python3
"""
Prop 15.68 residual: exact resolvent gain of T on refined C-classes.

Build class-constant matrix of T on (type6, ext-hist) classes; split |κ|=1 vs |κ|=3;
compute ρ = (4p I - T)^{-1} (Tκ/p²) classwise and report max |ρ| on |κ|=1
vs budget (p-4)/(2p²). Multi-worker class assembly.

Writes evidence/e1_gmin_m4_resolvent.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from concurrent.futures import as_completed
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from workers import pool, require_workers  # noqa: E402


def _blas1():
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[k] = "1"


def _worker_build_T_class(p: int) -> dict:
    _blas1()
    sys.path.insert(0, str(ROOT / "src"))
    from minmax_quadratic import paley_conference_prime_power
    from e1_gmin_moduli import build_classes, kappa as kap_fn

    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    classes, q2k = build_classes(C)
    keys = sorted(classes.keys())
    nk = len(keys)
    key_idx = {k: i for i, k in enumerate(keys)}
    kap_vec = np.array([kap_fn(C, classes[k][0]) for k in keys], dtype=float)

    # Exact T on class-constant functions: sample all quads if small, else many
    # Tmat[i,j] = avg_{S in class i} sum_{v,r: S_vr in class j} C_vr
    samples_per = 80 if p == 5 else 40
    rng = np.random.default_rng(p)
    Tmat = np.zeros((nk, nk))
    for i, kA in enumerate(keys):
        quads = classes[kA]
        take = quads if len(quads) <= samples_per else [
            quads[j] for j in rng.choice(len(quads), samples_per, replace=False)
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
                    kB = q2k[Sp]
                    acc[key_idx[kB]] += float(C[v, r])
        Tmat[i, :] = acc / len(take)

    # Source b = Tκ / p² classwise: Tκ from formula -6*star, constant on type6
    def star(S):
        s = 0
        for v in S:
            prod = 1
            for u in S:
                if u != v:
                    prod *= int(C[v, u])
            s += prod
        return s

    Tkap = np.array([-6.0 * star(classes[k][0]) for k in keys])
    b = Tkap / (p * p)

    # Solve (4p I - T) ρ = b
    A = 4 * p * np.eye(nk) - Tmat
    try:
        rho = np.linalg.solve(A, b)
        solve_ok = True
        resid = float(np.max(np.abs(A @ rho - b)))
    except np.linalg.LinAlgError:
        rho = np.linalg.lstsq(A, b, rcond=None)[0]
        solve_ok = False
        resid = float(np.max(np.abs(A @ rho - b)))

    m4_pred = kap_vec / (p * p) + rho

    # Index |κ|=1 and |κ|=3
    idx1 = np.where(np.abs(kap_vec) == 1)[0]
    idx3 = np.where(np.abs(kap_vec) == 3)[0]

    # Compare to true Max+ m4 on classes
    if p == 5:
        Y = np.load("/tmp/maxplus_p5.npy").astype(np.float64)
    elif p == 7:
        Y = np.load("/tmp/e1_p7/maxplus.npy").astype(np.float64)
    else:
        Y = None

    true_m4 = None
    if Y is not None:
        true_m4 = []
        for k in keys:
            a, b_, c, d = classes[k][0]
            true_m4.append(float(np.mean(Y[:, a] * Y[:, b_] * Y[:, c] * Y[:, d])))
        true_m4 = np.array(true_m4)

    L_abs = (p - 2) / (2.0 * p * p)
    rho_budget = (p - 4) / (2.0 * p * p)
    cand = (p - 2) / (p * (2.0 * p + 3))

    max_abs_rho_1 = float(np.max(np.abs(rho[idx1]))) if len(idx1) else None
    max_abs_m4_pred_1 = float(np.max(np.abs(m4_pred[idx1]))) if len(idx1) else None

    # Same-sign ρ with κ
    same = []
    for i in idx1:
        if np.sign(rho[i]) * np.sign(kap_vec[i]) > 0:
            same.append(abs(rho[i]))
    max_same_rho = float(max(same)) if same else 0.0

    # Effective gain: max |ρ| on |κ|=1 / (24/p²)
    source_inf = 24.0 / (p * p)
    gain = max_abs_rho_1 / source_inf if max_abs_rho_1 is not None else None
    gain_budget = (p - 4) / 48.0

    out = {
        "p": int(p),
        "n_classes": nk,
        "n_kappa1_classes": int(len(idx1)),
        "n_kappa3_classes": int(len(idx3)),
        "solve_ok": solve_ok,
        "solve_resid": resid,
        "lam_max_T": float(np.max(np.linalg.eigvals(Tmat).real)),
        "lam_min_T": float(np.min(np.linalg.eigvals(Tmat).real)),
        "gap_4p_minus_lammax": float(4 * p - np.max(np.linalg.eigvals(Tmat).real)),
        "max_abs_rho_kappa1": max_abs_rho_1,
        "max_abs_rho_same_sign": max_same_rho,
        "rho_budget_(p-4)/(2p^2)": rho_budget,
        "same_sign_rho_le_budget": max_same_rho <= rho_budget + 1e-9,
        "max_abs_m4_pred_kappa1": max_abs_m4_pred_1,
        "L_abs": L_abs,
        "pred_m4_le_L": max_abs_m4_pred_1 is not None
        and max_abs_m4_pred_1 <= L_abs + 1e-9,
        "candidate_UB": cand,
        "pred_m4_le_candidate": max_abs_m4_pred_1 is not None
        and max_abs_m4_pred_1 <= cand + 1e-9,
        "effective_gain": gain,
        "gain_budget_(p-4)/48": gain_budget,
        "gain_le_budget": gain is not None and gain <= gain_budget + 1e-9,
        "source_inf_24/p2": source_inf,
    }
    if true_m4 is not None:
        out["max_abs_true_m4_kappa1"] = float(np.max(np.abs(true_m4[idx1])))
        out["pred_vs_true_m4_max_err"] = float(np.max(np.abs(m4_pred - true_m4)))
        out["pred_matches_true"] = bool(np.max(np.abs(m4_pred - true_m4)) < 1e-3)
    return out


def main():
    W = require_workers(min_workers=4)
    print(f"m4_resolvent: W={W}", flush=True)
    out = {"workers": W, "certs": {}}
    # p=5,7 in parallel (each heavy on class T)
    with pool(W) as ex:
        futs = {ex.submit(_worker_build_T_class, p): p for p in (5, 7)}
        for fut in as_completed(futs):
            r = fut.result()
            out["certs"][str(r["p"])] = r
            print(
                f"  p={r['p']}: classes={r['n_classes']} "
                f"max|ρ|_1={r['max_abs_rho_kappa1']:.6f} bud={r['rho_budget_(p-4)/(2p^2)']:.6f} "
                f"same_le={r['same_sign_rho_le_budget']} "
                f"gain={r['effective_gain']} gain_bud={r['gain_budget_(p-4)/48']:.6f} "
                f"pred_le_L={r['pred_m4_le_L']} "
                f"true_err={r.get('pred_vs_true_m4_max_err')}",
                flush=True,
            )

    out["status"] = (
        "Class-level resolvent of (4pI-T) on refined C-classes. "
        "If pred matches true m4 and same-sign ρ ≤ budget, structure is sound. "
        "OPEN: closed-form gain ≤(p-4)/48 for all primes p≥5."
    )
    path = ROOT / "evidence" / "e1_gmin_m4_resolvent.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"], flush=True)
    print(f"wrote {path}", flush=True)


if __name__ == "__main__":
    main()
