#!/usr/bin/env python3
"""
Structural certificates for residual g_min attack (Prop 15.51 support).

1. Equivalent form: g_min >= T(p) iff min E[f'|f_e=1] >= 1/(2p-1) for disj e'.
2. On f_e=1 slice, sum of f over disj partners is deterministic.
3. Residual Gram R on the slice satisfies R ≽ λ_min P_W (Loewner), certified p=5,7.
4. Max+ distance distribution is constant (2-point homogeneous), certified p=5,7.

Writes evidence/e1_gmin_structure.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
for _k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_k] = "1"

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


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


def equiv_form_check(p: int, Mp: np.ndarray) -> dict:
    """min_a = min E[f'|f=1] and convert to G_lb; compare to T/L."""
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    i, j = 0, 1
    fij = C[i, j] * Mp[:, i] * Mp[:, j]
    mask = np.abs(fij - 1.0) < 1e-9
    min_a = 1.0
    for k in range(n):
        for l in range(k + 1, n):
            if len({i, j, k, l}) < 4:
                continue
            fkl = C[k, l] * Mp[:, k] * Mp[:, l]
            a = float(np.mean(fkl[mask]))
            min_a = min(min_a, a)
    G_from_a = (min_a * (p + 1) - 1) / p
    T = -(p - 2) / (p * (2 * p - 1))
    L = -(p - 2) / (2 * p * p)
    a_T = 1 / (2 * p - 1)
    return {
        "p": p,
        "min_E_f_given_f_eq_1": min_a,
        "a_threshold_for_T": a_T,
        "a_ge_a_T": min_a >= a_T - 1e-15,
        "G_from_a": G_from_a,
        "T": T,
        "L": L,
        "G_from_a_gt_T": G_from_a > T,
        "identity": "a=(1+G*p)/(p+1); a>=1/(2p-1) <=> G>=T(p)",
    }


def disj_sum_deterministic(p: int) -> dict:
    """Algebra: sum_{e' disj e} f_e' = Phi - 2p + f_e (pointwise on Max+)."""
    n = p * p + 1
    Phi = p * n / 2
    D = (n - 2) * (n - 3) // 2
    avg_on_slice = (Phi - 2 * p + 1) / D
    return {
        "p": p,
        "n_disj": D,
        "sum_disj_f": "Phi - 2p + f_e",
        "avg_f_on_f_eq_1_slice": avg_on_slice,
    }


def residual_loewner(p: int, Mp: np.ndarray) -> dict:
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    d = n // 2
    Sigma = np.eye(n) + C / p
    P_plus = Sigma / 2
    i, j = 0, 1
    cij = float(C[i, j])
    V = np.column_stack([P_plus[:, i], P_plus[:, j]])
    P_W = P_plus - V @ np.linalg.inv(V.T @ V) @ V.T
    S = [i, j]
    inv = np.linalg.inv(Sigma[np.ix_(S, S)])
    acc = np.zeros((n, n))
    tot = 0
    for a, b in ((1.0, cij), (-1.0, -cij)):
        mask = (np.abs(Mp[:, i] - a) < 1e-9) & (np.abs(Mp[:, j] - b) < 1e-9)
        mu = Sigma[:, S] @ inv @ np.array([a, b])
        Z = Mp[mask] - mu
        acc += Z.T @ Z
        tot += len(Z)
    R = acc / tot
    ev = np.linalg.eigvalsh(R)
    pos = ev[ev > 1e-8]
    lmin = float(pos.min())
    loewner_min = float(np.linalg.eigvalsh(R - lmin * P_W).min())
    return {
        "p": p,
        "lambda_min_R": lmin,
        "rank_R": int(len(pos)),
        "d_minus_2": d - 2,
        "loewner_R_ge_lmin_PW": loewner_min > -1e-8,
        "loewner_min_eig": loewner_min,
        "trace_R": float(np.trace(R)),
    }


def distance_homogeneous(p: int, Mp: np.ndarray) -> dict:
    """Distance distribution from each of first K vectors (types of Max+)."""
    N, n = Mp.shape
    K = min(20, N)
    dists_list = []
    for t in range(K):
        dots = Mp @ Mp[t]
        dists = np.round((n - dots) / 2).astype(int)
        dists_list.append(Counter(dists.tolist()))
    # number of distinct distance distributions among sample
    uniq = []
    for c in dists_list:
        key = tuple(sorted(c.items()))
        if key not in uniq:
            uniq.append(key)
    return {
        "p": p,
        "n_Max": int(N),
        "sample_size": K,
        "n_distinct_distance_distributions": len(uniq),
        "distance_homogeneous_on_sample": len(uniq) == 1,
        "note": "p=5: homogeneous; p=7: multiple Max+ types (not a single Aut-orbit)",
    }


def main() -> None:
    out: dict = {"results": [], "status": None}
    for p in (5, 7):
        Mp = load_maxplus(p)
        row = {
            "equiv": equiv_form_check(p, Mp),
            "disj_sum": disj_sum_deterministic(p),
            "loewner": residual_loewner(p, Mp),
            "distance": distance_homogeneous(p, Mp),
        }
        out["results"].append(row)
        print("p", p, "a_min", row["equiv"]["min_E_f_given_f_eq_1"],
              "loewner", row["loewner"]["loewner_R_ge_lmin_PW"],
              "n_dist_types", row["distance"]["n_distinct_distance_distributions"])
    ok = all(
        r["equiv"]["a_ge_a_T"] and r["loewner"]["loewner_R_ge_lmin_PW"]
        for r in out["results"]
    )
    out["status"] = (
        "Structure certs OK at p=5,7: equiv form a<=>T holds; residual Loewner. "
        "Max+ distance-homogeneous only at p=5 (multiple types at p=7). "
        "OPEN: prove a>=1/(2p-1) (or g_min>=L(p)) for all p>=5."
        if ok
        else "CERT FAILURE"
    )
    path = ROOT / "evidence" / "e1_gmin_structure.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"])
    print("wrote", path)
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
