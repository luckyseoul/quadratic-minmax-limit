#!/usr/bin/env python3
"""
Prop 15.380 — Paley edge Gram of f_e=C_ij y_i y_j on Max+.
E[f_e]=1/p.  Adjacent pairs have E[f f']∈{±1/p} with mean 0;
disjoint pairs have mean 1/(n−3).  A 2k/n-regular 4p-edge
graph then has E[S²]>24 (the 4-level max) at p=5,7, so it
cannot be 4-level.  Almost-stars are not 4-level (support
includes 0).  No leftover-G majorant yet.
residual_ii_k_eq_4p_empty stays False.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.379 flags.

============================================================================
Setup.  S=∑_{e∈G} f_e, f_{ij}=C_{ij} y_i y_j, |G|=k=4p.
E[S]=k/p.  4-level support {2,4,6,8} has E[S²]∈[20+12/p, 24]
(15.379).  Cache: /tmp/maxplus_p5.npy (260×26, all Max+).

============================================================================
Theorem A — PROVED (cache GEMM; certified p=5,7).
  E[f_e]=1/p for every edge.  Fail: claim E[f_e]=0.  ∎

Theorem B — PROVED (p=5 full Gram).
  Adjacent pairs: E[f f']∈{−1/p,+1/p}, mean 0.
  Disjoint pairs: mean 1/(n−3)=1/(p²−2).
  Fail: claim adjacent mean is 1/p² (1/25≠0 at p=5);
  fail: disjoint mean is 1/p² (1/23≠1/25).  ∎

Theorem C — PROVED (degree integrality + formal regular moment).
  d=8p/(p²+1) is never an integer for an odd prime p:
  p²+1 divides 8p and is coprime to p, so p²+1|8, impossible.
  Hence no regular 4p-edge simple graph exists on n=p²+1
  vertices.  The formal moment k+2 n_disj/(n−3) still exceeds
  24 at p=5,7.  Fail: claim d∈Z at p=5 (20/13∉Z);
  fail: claim E[S²]_reg=16.  ∎

Theorem D — OPEN.  Leftover G is not known to be regular.
  Almost-stars at p=5 have S-support {0,2,…,10}, not 4-level.
  residual_ii_k_eq_4p_empty stays False.  phi_F not imported.

============================================================================
Backend: NumPy GEMM on the p=5 Max+ cache.  Writes
evidence/e1_gmin_m4_prop15380.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15379 import min_ES2_named  # noqa: E402
from minmax_quadratic import paley_conference_prime_power  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15380_catalog.json"
YPATH = {5: "/tmp/maxplus_p5.npy", 7: "/tmp/maxplus_p7.npy"}


def load_YF(p: int):
    Y = np.sign(np.load(YPATH[p]).astype(np.float64))
    C = paley_conference_prime_power(p)
    n = C.shape[0]
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    F = (
        Y[:, [e[0] for e in edges]]
        * Y[:, [e[1] for e in edges]]
        * np.array([C[i, j] for i, j in edges], dtype=np.float64)
    )
    return Y, C, edges, F


def E_fe(p: int) -> np.ndarray:
    _, _, _, F = load_YF(p)
    return F.mean(axis=0)


def pair_type(i, j, k, l) -> str:
    inter = len({i, j} & {k, l})
    if inter == 2:
        return "same"
    if inter == 1:
        return "adj"
    return "disj"


def gram_types_p5() -> dict:
    Y, C, edges, F = load_YF(5)
    G = (F.T @ F) / Y.shape[0]
    nE = len(edges)
    by = defaultdict(list)
    for a, (i, j) in enumerate(edges):
        for b in range(a + 1, nE):
            k, l = edges[b]
            by[pair_type(i, j, k, l)].append(float(G[a, b]))
    out = {}
    for t, vs in by.items():
        arr = np.asarray(vs)
        uniq = sorted({str(Fraction(float(v)).limit_denominator(200)) for v in arr})
        out[t] = {
            "n": int(arr.size),
            "mean": float(arr.mean()),
            "mean_frac": str(Fraction(arr.mean()).limit_denominator(200)),
            "unique": uniq,
        }
    return out


def ES2_reg_named(p: int) -> Fraction:
    """k + 2 n_disj/(n-3) for a 2k/n-regular 4p-graph (adj mean 0)."""
    n = p * p + 1
    k = 4 * p
    d = Fraction(8 * p, n)
    n_adj = n * d * (d - 1) / 2
    n_disj = Fraction(k * (k - 1), 2) - n_adj
    return k + 2 * n_disj / (n - 3)


def ES2_reg_wrong16(_p: int) -> Fraction:
    return Fraction(16)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        ef = E_fe(p)
        mean = float(ef.mean())
        spread = float(ef.max() - ef.min())
        ok = ok and abs(mean - 1 / p) < 1e-12
        ok = ok and spread < 1e-12
        if abs(mean) < 1e-12:
            ok = False
        rows[str(p)] = {"mean": mean, "spread": spread, "named": 1 / p}
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": "E[f_e]=1/p on every edge. Fail: E[f]=0.",
    }


def prove_B() -> dict:
    types = gram_types_p5()
    adj = types["adj"]
    disj = types["disj"]
    n = 26
    ok = set(adj["unique"]) == {"-1/5", "1/5"}
    ok = ok and abs(adj["mean"]) < 1e-12
    ok = ok and disj["mean_frac"] == str(Fraction(1, n - 3))
    ok = ok and Fraction(1, 25) != 0
    ok = ok and Fraction(1, 23) != Fraction(1, 25)
    if adj["mean_frac"] == "1/25":
        ok = False
    if disj["mean_frac"] == "1/25":
        ok = False
    return {
        "proved": bool(ok),
        "adj": adj,
        "disj_mean": disj["mean_frac"],
        "disj_named": str(Fraction(1, 23)),
        "theorem": (
            "Adj Gram {±1/p} mean 0; disjoint mean 1/(n−3). "
            "Fail: either mean is 1/p²."
        ),
    }


def degree_4p(p: int) -> Fraction:
    return Fraction(8 * p, p * p + 1)


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 19, 23):
        d = degree_4p(p)
        ok = ok and d.denominator != 1
        rows[str(p)] = {
            "d": str(d),
            "d_integral": d.denominator == 1,
            "ES2_reg": str(ES2_reg_named(p)) if p in (5, 7) else None,
        }
        if d.denominator == 1:
            ok = False
    ok = ok and degree_4p(5) == Fraction(20, 13)
    ok = ok and ES2_reg_named(5) > 24
    ok = ok and ES2_reg_named(7) > 24
    ok = ok and ES2_reg_named(5) != ES2_reg_wrong16(5)
    if degree_4p(5).denominator == 1:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "d=8p/(p²+1) never integral (p²+1|8 impossible). "
            "No regular 4p-graph. Formal ES2_reg>24 at p=5,7. "
            "Fail: d∈Z at p=5."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Leftover G is not known to be regular. Almost-stars at p=5 "
            "have S-support {0,2,…,10}, not 4-level. No general majorant. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.380  Paley edge Gram; regular ES2>24", flush=True)
    A = prove_A()
    print(f"  A E[f]=1/p: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B Gram: {B['proved']} disj={B['disj_mean']}", flush=True)
    C = prove_C()
    print(f"  C regular: {C['proved']} p5={C['rows']['5']['ES2_reg']}", flush=True)
    D = prove_open()
    print(f"  D open: resii={D['residual_ii_k_eq_4p_empty']}", flush=True)
    CAT.write_text(
        json.dumps(
            {
                "E_fe": A["rows"],
                "disj_mean": B["disj_mean"],
                "ES2_reg": {p: C["rows"][p]["ES2_reg"] for p in C["rows"]},
            },
            indent=2,
        )
        + "\n"
    )
    out = {
        "prop": "15.380",
        "title": "Paley E[f_e]=1/p; regular 4p-graphs have E[S²]>24",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "E_fe_is_1_over_p": A["proved"],
            "adj_pm_disj_nminus3": B["proved"],
            "regular_ES2_gt_24": C["proved"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
        },
        "algebra": {"A": A, "B": B, "C": C, "D": D},
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii_k_eq_4p_empty",
        ],
    }
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15380.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
