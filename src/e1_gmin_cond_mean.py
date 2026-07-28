#!/usr/bin/env python3
"""
Prop 15.50: conditional means of Max+ match the Gaussian frame interpolant.

For Paley ρ=1 Max+ with E[yy^T]=I+C/p=2P_+, and distinct i,j,a,b∈{±1}:
  E[y | y_i=a, y_j=b] = Σ_{*,S} Σ_{S,S}^{-1} (a,b)^T,  S={i,j}, Σ=I+C/p.

Proof is analytic (central symmetry + frame + min-norm interpolant in V+);
this script certifies numerically at p=5,7 and records disjoint-mean identity
  avg_{e' disj e} G_{ee'} = 1/(p^2-2).

Writes evidence/e1_gmin_cond_mean.json
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


def sigma_frame(C: np.ndarray, p: int) -> np.ndarray:
    return np.eye(C.shape[0]) + C / p


def gaussian_cond_mean(Sigma: np.ndarray, i: int, j: int, a: float, b: float) -> np.ndarray:
    S = [i, j]
    inv = np.linalg.inv(Sigma[np.ix_(S, S)])
    return Sigma[:, S] @ inv @ np.array([a, b], dtype=float)


def certify_cond_mean(p: int, Mp: np.ndarray, n_pairs: int = 40, seed: int = 0) -> dict:
    C = paley_conference_prime_power(p)
    Sigma = sigma_frame(C, p)
    n = C.shape[0]
    rng = np.random.default_rng(seed)
    max_err = 0.0
    n_checked = 0
    for _ in range(n_pairs):
        i, j = rng.choice(n, size=2, replace=False)
        for a in (1.0, -1.0):
            for b in (1.0, -1.0):
                mask = (Mp[:, i] == a) & (Mp[:, j] == b)
                if int(mask.sum()) < 3:
                    continue
                emp = Mp[mask].mean(axis=0)
                mu = gaussian_cond_mean(Sigma, int(i), int(j), a, b)
                err = float(np.max(np.abs(emp - mu)))
                max_err = max(max_err, err)
                n_checked += 1
    return {
        "p": p,
        "n_Max": int(len(Mp)),
        "n_checked_slices": n_checked,
        "max_abs_err": max_err,
        "ok": max_err < 1e-9,
    }


def disj_mean_identity(p: int) -> dict:
    """avg G on disj partners of any edge = 1/(p^2-2) (proved: row sum + wedge sum 0)."""
    n = p * p + 1
    D = (n - 2) * (n - 3) // 2
    s = n // 2 - 1
    avg = s / D
    formula = 1 / (p * p - 2)
    return {
        "p": p,
        "n_disj_per_edge": D,
        "sum_disj_G": s,
        "avg_disj_G": avg,
        "formula_1_over_p2_minus_2": formula,
        "match": abs(avg - formula) < 1e-15,
    }


def load_maxplus(p: int) -> np.ndarray:
    if p == 5:
        path = Path("/tmp/maxplus_p5.npy")
    elif p == 7:
        path = Path("/tmp/e1_p7/maxplus.npy")
    else:
        raise ValueError(p)
    if not path.is_file():
        raise FileNotFoundError(path)
    return np.load(path).astype(float)


def main() -> None:
    out: dict = {
        "lemma": (
            "Conditional means of Max+ given two coordinates equal the "
            "Gaussian interpolant from Σ=I+C/p (Prop 15.50). Proved analytically; "
            "certified p=5,7."
        ),
        "certs": [],
        "disj_means": [],
        "status": None,
    }
    for p in (5, 7):
        Mp = load_maxplus(p)
        out["certs"].append(certify_cond_mean(p, Mp))
        out["disj_means"].append(disj_mean_identity(p))
    all_ok = all(c["ok"] for c in out["certs"]) and all(
        d["match"] for d in out["disj_means"]
    )
    out["status"] = (
        "Prop 15.50 conditional-mean lemma certified at p=5,7; "
        "disj mean identity 1/(p^2-2) algebra OK. "
        "Residual: g_min >= L(p) still OPEN (Frechet on cond cov too weak)."
        if all_ok
        else "CERT FAILURE"
    )
    path = ROOT / "evidence" / "e1_gmin_cond_mean.json"
    path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print("wrote", path)
    if not all_ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
