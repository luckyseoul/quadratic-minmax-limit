#!/usr/bin/env python3
"""Historical spectrum census for the retracted Prop 15.55 obstruction.

The numerical spectra remain valid. The former discrete conclusion is false:
ker(G-(n/2)P1) also contains ker G, including star differences. Proposition
15.720 is the valid replacement for required bi-tight levels 2 and 3.
Its bi-tight level-4 corollary does not exclude one-sided tightness.
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
    # The old next step was invalid: ker(G_perp) also contains ker(G).

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
        "tight_size_2p_impossible_if_simple_max": False,
        "spectral_hypothesis_holds": bool(simple),
        "retracted": True,
        "note": (
            "The spectrum can be simple at n/2 while G_perp has ker(G) in "
            "its kernel. Tightness only puts the centered indicator in ker(G); "
            "it does not make the indicator constant."
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
        "Retraction: ker(G_perp)=span{1}+(ker G intersect 1-perp), not "
        "span{1}. The p=5,7 spectral hypothesis is certified but has no "
        "tight-cover conclusion. Use Prop 15.720 for required bi-tight levels."
    )
    ok57 = all(
        r["n_over_2_is_simple_max"] and r["allones_quad_is_4"]
        for r in out["results"]
        if r["p"] in (5, 7)
    )
    p3_fails = not out["results"][0]["n_over_2_is_simple_max"]
    out["status"] = (
        "p=5,7: λ_max(G)=n/2 simple is certified, but the former "
        "tight-cover implication is retracted. "
        "p=3: λ_max=8>n/2 (bi-tight possible). "
        "Required bi-tight levels are handled by Prop 15.720. L OPEN."
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
