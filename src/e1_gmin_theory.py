#!/usr/bin/env python3
"""
Theory probes for g_min (disjoint edge correlation floor) and spike lemma.

1. Closed-form candidates for g_min at p=3,5,7
2. Linear identity for 4-point moments from Cy=py
3. Spike: multi-bit / local search from S=2 Max+ with low sigma
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from minmax_quadratic import paley_conference_prime_power  # noqa: E402


def g_min_table() -> dict:
    """Known certified g_min values and comparisons."""
    rows = []
    # (p, g_min, N_Max, n, Phi)
    data = [
        (3, -1 / 3, 12, 10, 15),
        (5, -3 / 65, 260, 26, 65),
        (7, -436 / 11452, 11452, 50, 175),
    ]
    for p, g, N, n, Phi in data:
        rows.append(
            {
                "p": p,
                "g_min": g,
                "N": N,
                "n": n,
                "Phi": Phi,
                "g_times_N": g * N,
                "minus_1_over_p": -1 / p,
                "strict_gt_minus_1_over_p": g > -1 / p + 1e-15,
                "minus_1_over_15": -1 / 15,
                "matching_level2_blocked": g > -1 / 15 + 1e-15,
            }
        )
    return {"rows": rows}


def moment_identity_check(p: int = 5) -> dict:
    """Verify p * m_ijkl = wick_num/p + sum_{r notin} C_ir m_rjkl."""
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    if p == 5:
        Mp = np.load("/tmp/maxplus_p5.npy").astype(float)
    else:
        raise ValueError("need cached Max+")
    N = len(Mp)

    def m4(a, b, c, d):
        return float(np.mean(Mp[:, a] * Mp[:, b] * Mp[:, c] * Mp[:, d]))

    errs = []
    rng = np.random.default_rng(0)
    for _ in range(200):
        i, j, k, l = rng.choice(n, 4, replace=False)
        left = p * m4(i, j, k, l)
        wick = (
            C[i, j] * C[k, l] + C[i, k] * C[j, l] + C[i, l] * C[j, k]
        ) / p
        residual = 0.0
        for r in range(n):
            if r in (i, j, k, l):
                continue
            residual += C[i, r] * m4(r, j, k, l)
        # when r in {j,k,l} the identity uses reduced moments:
        reduced = (
            C[i, j] * C[k, l] / p
            + C[i, k] * C[j, l] / p
            + C[i, l] * C[j, k] / p
        )
        # full RHS = reduced + residual, but reduced already is wick/p? wick/p = reduced
        # left = reduced*p? Let's use full sum over all r including j,k,l:
        full = 0.0
        for r in range(n):
            if r == i:
                continue
            # E[y_r y_j y_k y_l] with possible repeats
            if r == j:
                full += C[i, r] * float(np.mean(Mp[:, k] * Mp[:, l]))  # y_j^2=1
            elif r == k:
                full += C[i, r] * float(np.mean(Mp[:, j] * Mp[:, l]))
            elif r == l:
                full += C[i, r] * float(np.mean(Mp[:, j] * Mp[:, k]))
            else:
                full += C[i, r] * m4(r, j, k, l)
        errs.append(abs(left - full))
    return {
        "max_abs_err": float(max(errs)),
        "mean_abs_err": float(np.mean(errs)),
        "identity_holds": max(errs) < 1e-8,
    }


def gmin_from_psd_bound(p: int) -> dict:
    """
    Attempt analytic lower bound on G for disjoint pairs.

    Centered features u_e = f_e - 1/p live in a space of rank r = Phi = p n / 2
    (certified rank at p=3,5). Number of edges M = binom(n,2).

    For ANY set of unit vectors, min inner product can be -1.
    Restrict to a matching of size m = floor(n/2): m vectors in R^r.
    If they were equiangular simplex, <u,v> = -1/(m-1) for centered? 

    Actual: Gram of raw f_e has diag 1, off >= g_min for disjoint.
    PSD of matching Gram: 1 + (m-1) g_avg >= 0 if constant off-diag,
    so g_avg >= -1/(m-1).

    This constrains averages on large matchings, not the global min of a single pair.
    """
    n = p * p + 1
    Phi = p * n / 2
    m = n // 2
    return {
        "p": p,
        "rank_Phi": Phi,
        "matching_size": m,
        "simplex_bound_on_avg": -1 / (m - 1),
        "minus_1_over_p": -1 / p,
        "simplex_beats_minus_1_over_p": -1 / (m - 1) > -1 / p,
    }


def main() -> None:
    out = {
        "g_min_table": g_min_table(),
        "moment_identity_p5": moment_identity_check(5),
        "psd_bounds": [gmin_from_psd_bound(p) for p in (3, 5, 7, 11)],
    }
    # Key observation: for p>=5, m=n/2 = (p^2+1)/2, -1/(m-1) = -2/(p^2-1)
    # Compare -2/(p^2-1) vs -1/p: 2/(p^2-1) < 1/p iff 2p < p^2-1 iff p^2 - 2p -1 > 0 iff p>=3.
    # So simplex avg bound -1/(m-1) > -1/p for p>=3.
    # At p=3: -1/(5-1)= -0.25 > -1/3, but actual g_min=-1/3 < -0.25!
    # So some matching pairs MUST have G < simplex constant-equiangular average —
    # meaning G is non-constant on matchings; the bound is on the average of a PSD Gram,
    # and min can be more negative than -1/(m-1) if other entries compensate.
    path = ROOT / "evidence" / "e1_gmin_theory.json"
    path.write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))
    print("wrote", path)


if __name__ == "__main__":
    main()
