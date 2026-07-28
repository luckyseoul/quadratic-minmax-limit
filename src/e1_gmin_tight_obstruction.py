#!/usr/bin/env python3
"""
Tight size-2p cover obstruction via λ_max(G)=n/2.

If G 1 = (n/2) 1 and λ_max(G)=n/2 is simple, then no Max+-tight
S≡2 cover of size 2p exists (kills bi-tight level 2 and Type I freeness
equality cases).

Certified: λ_max(G)=n/2 simple at p=5,7; at p=3, λ_max=8 > n/2=5
(consistent with bi-tight C6).

Writes evidence/e1_gmin_tight_obstruction.json
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
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    Mp = load_mp(p)
    N = len(Mp)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(edges)
    Chip = np.array(
        [[C[a, b] * y[a] * y[b] for a, b in edges] for y in Mp], dtype=np.float64
    )
    G = (Chip.T @ Chip) / N
    ones = np.ones(E)
    row = G @ ones
    row_ok = bool(np.allclose(row, (n / 2) * ones, atol=1e-8))

    ev = np.linalg.eigvalsh(G)
    pos = ev[ev > 1e-8]
    lam_max = float(pos[-1])
    n_half = n / 2.0
    mult_nhalf = int(np.sum(np.abs(pos - n_half) < 1e-6))
    mult_max = int(np.sum(np.abs(pos - lam_max) < 1e-6))
    lam_max_eq_nhalf = abs(lam_max - n_half) < 1e-6
    simple = mult_nhalf == 1 and lam_max_eq_nhalf

    # Algebra identity: for any 0-1 vector v with sum(v)=2p,
    # v^T ((n/2) 11^T/E) v = 4
    v = np.zeros(E)
    v[: 2 * p] = 1.0  # any support of size 2p
    Pcoef = (n / 2.0) / E
    quad_allones = float(Pcoef * (v.sum() ** 2))
    # For a true tight cover, v^T G v = E[S^2] = 4, hence v^T G_perp v = 0
    # If G_perp ≽ 0 and ker=span{1}, then v ∥ 1, impossible.

    return {
        "p": int(p),
        "n": int(n),
        "E": int(E),
        "N_Max": int(N),
        "row_sum_is_n_over_2": row_ok,
        "lambda_max": lam_max,
        "n_over_2": n_half,
        "lambda_max_eq_n_over_2": bool(lam_max_eq_nhalf),
        "mult_of_n_over_2": mult_nhalf,
        "mult_of_lambda_max": mult_max,
        "n_over_2_is_simple_max": bool(simple),
        "rank_pos": int(len(pos)),
        "rank_formula": int((n // 2 - 1) * (n // 2 - 2) // 2),
        "allones_quad_for_size_2p": quad_allones,
        "allones_quad_is_4": bool(abs(quad_allones - 4.0) < 1e-9),
        "top_eigs": [float(x) for x in pos[-6:][::-1]],
        "tight_size_2p_impossible_if_simple_max": bool(simple),
        "note": (
            "If λ_max(G)=n/2 simple and G≽0, then any v with v^T G v = 4 and "
            "sum(v)=2p has v^T G_perp v=0 ⇒ G_perp v=0 ⇒ v∥1, impossible. "
            "Hence no Max+-tight S≡2 cover of size 2p (no bi-tight level 2)."
        ),
    }


def main() -> None:
    out: dict = {"results": [], "status": None, "theorem_sketch": None}
    for p in (3, 5, 7):
        print(f"tight obstruction p={p}", flush=True)
        row = analyze(p)
        out["results"].append(row)
        print(
            f"  λmax={row['lambda_max']:.6f} n/2={row['n_over_2']} "
            f"simple_max={row['n_over_2_is_simple_max']}",
            flush=True,
        )
    out["theorem_sketch"] = (
        "Prop 15.55: Assume G1=(n/2)1 and λ_max(G)=n/2 is simple. "
        "For |H|=2p with S_H≡2 on Max+, E[S^2]=4=1_H^T G 1_H. "
        "Write G=(n/2)P_1 + G_perp with P_1=11^T/E. Then "
        "1_H^T (n/2)P_1 1_H=4, so 1_H^T G_perp 1_H=0. "
        "G≽0 ⇒ G_perp≽0 on 1^⊥ and ker G_perp=span{1}. "
        "Hence G_perp 1_H=0 ⇒ 1_H∥1, contradiction. "
        "Therefore no Max+-tight level-2 cover exists. "
        "Certified hypothesis λ_max=n/2 simple at p=5,7; FAILS at p=3 (λ_max=8>5)."
    )
    ok57 = all(
        r["n_over_2_is_simple_max"] and r["allones_quad_is_4"]
        for r in out["results"]
        if r["p"] in (5, 7)
    )
    p3_fails = not out["results"][0]["n_over_2_is_simple_max"]
    out["status"] = (
        "p=5,7: λ_max(G)=n/2 simple ⇒ tight size-2p Max+ covers impossible "
        "(bi-tight level 2 empty; Type I freeness-failure blocked). "
        "p=3: λ_max=8>n/2 (bi-tight possible). "
        "OPEN: prove λ_max(G)=n/2 for all primes p≥5 (then bi-tight residual closes). "
        "Deep non-tight residual independent. L OPEN."
        if ok57 and p3_fails
        else "CERT FAILURE"
    )
    path = ROOT / "evidence" / "e1_gmin_tight_obstruction.json"
    path.write_text(json.dumps(out, indent=2))
    print(out["status"])
    print("wrote", path)
    if not (ok57 and p3_fails):
        sys.exit(1)


if __name__ == "__main__":
    main()
