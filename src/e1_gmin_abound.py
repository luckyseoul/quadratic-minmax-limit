#!/usr/bin/env python3
"""
Slice a-bound analysis for Prop 15.51 residual.

On f_e=1: a(e') = E[f_{e'}|f_e=1] = (1+p G)/(p+1).
Need min a >= 1/(2p-1) <=> g_min >= T.

Compares:
  - true Max+ a_min
  - Wick/Gaussian fourth-moment a
  - mean-only product C_kl mu_k mu_l (residual R=0)

Writes evidence/e1_gmin_abound.json
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

from e1_gmin_moduli import L_p, T_p, load_maxplus  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def analyze(p: int) -> dict:
    C = paley_conference_prime_power(p).astype(float)
    n = C.shape[0]
    Sigma = np.eye(n) + C / p
    i, j = 0, 1
    S = [i, j]
    inv = np.linalg.inv(Sigma[np.ix_(S, S)])
    cij = float(C[i, j])
    states = [np.array([1.0, cij]), np.array([-1.0, -cij])]
    a_T = 1.0 / (2 * p - 1)

    Mp = load_maxplus(p)
    fij = C[i, j] * Mp[:, i] * Mp[:, j]
    mask = np.abs(fij - 1.0) < 1e-9

    pstates = []
    for ys in states:
        m = (np.abs(Mp[:, i] - ys[0]) < 1e-9) & (np.abs(Mp[:, j] - ys[1]) < 1e-9)
        pstates.append(int(m.sum()))
    w = np.array(pstates, dtype=float)
    w /= w.sum()
    mus = [Sigma[:, S] @ inv @ ys for ys in states]

    min_a_true = 1.0
    min_a_wick = 1.0
    min_a_mean = 1.0
    n_disj = 0
    for k in range(n):
        for l in range(k + 1, n):
            if len({i, j, k, l}) < 4:
                continue
            n_disj += 1
            fkl = C[k, l] * Mp[:, k] * Mp[:, l]
            a_true = float(np.mean(fkl[mask]))
            min_a_true = min(min_a_true, a_true)

            m4w = (
                C[i, j] * C[k, l] + C[i, k] * C[j, l] + C[i, l] * C[j, k]
            ) / (p * p)
            Gw = float(C[i, j] * C[k, l] * m4w)
            aw = (1.0 + p * Gw) / (p + 1)
            min_a_wick = min(min_a_wick, aw)

            am = 0.0
            for wt, mu in zip(w, mus):
                am += wt * float(C[k, l] * mu[k] * mu[l])
            min_a_mean = min(min_a_mean, am)

    G_true = (min_a_true * (p + 1) - 1.0) / p
    return {
        "p": p,
        "n": n,
        "n_disj": n_disj,
        "a_T": a_T,
        "min_a_true": min_a_true,
        "min_a_wick": min_a_wick,
        "min_a_mean_only": min_a_mean,
        "true_ge_aT": bool(min_a_true >= a_T - 1e-15),
        "wick_ge_aT": bool(min_a_wick >= a_T - 1e-15),
        "mean_ge_aT": bool(min_a_mean >= a_T - 1e-15),
        "G_true_from_a": G_true,
        "T": T_p(p),
        "L": L_p(p),
        "G_gt_T": bool(G_true > T_p(p)),
        "G_ge_L": bool(G_true >= L_p(p) - 1e-15),
        "state_weights": w.tolist(),
        "note": (
            "Wick a_min uses Gaussian 4th moments on the frame (often > true). "
            "mean-only ignores residual R; Fréchet on R is F15 (too weak)."
        ),
    }


def main() -> None:
    out: dict = {"results": [], "status": None}
    for p in (5, 7):
        print(f"abound p={p}", flush=True)
        row = analyze(p)
        out["results"].append(row)
        print(
            f"  aT={row['a_T']:.6f} true={row['min_a_true']:.6f} "
            f"wick={row['min_a_wick']:.6f} mean={row['min_a_mean_only']:.6f} "
            f"G={row['G_true_from_a']:.6f}",
            flush=True,
        )
    ok = all(r["true_ge_aT"] and r["G_gt_T"] for r in out["results"])
    out["status"] = (
        "a_min certs OK at p=5,7 (true > a_T). Wick/mean-only do not furnish a "
        "uniform LB for all p (mean-only undershoots residual). OPEN: prove "
        "a_min>=1/(2p-1) for all primes p>=5."
        if ok
        else "CERT FAILURE"
    )
    path = ROOT / "evidence" / "e1_gmin_abound.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"])
    print("wrote", path)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
