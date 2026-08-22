#!/usr/bin/env python3
"""
Aut-invariant 4-point master equation + |m₄|≤1 cannot prove GLOBAL QVAR.

Unnumbered kill.  Does **not** prove F̂(ψ)≥0.  Does **not** flip leftover
1/2/3, phi_F, L, Aut-Schur, Gsum, pairing.

SETUP
  T is C-signed Johnson on 4-sets.  True Max+ m₄⁺ satisfies
      T m = 4p m − 4κ/p
  and |m|≤1, and is Aut(C)-invariant (permutations with UᵀCU=C).
  QVAR ⇔ ⟨m, κ_{A_ψ}⟩ ≥ 0  (15.109 A + 15.597 on A_ψ ∈ Z).

KILL (computational, p=5 exact Aut-quotient; p=7 recorded)
  Restrict to Aut-class functions.  Torb is the quotient of T
  (certified: Torb κ = −6⋆, κ constant on orbits).  The affine space
      (Torb − 4p I) m = −4κ/p
  has a 2-dimensional kernel at p=5 (E_{4p}^{Aut}).  The box |m|≤1
  on that affine space is a nonempty polytope, and
      min ⟨m, κ_{A_ψ}⟩ = −285/4 < 0
  (HiGHS).  True Max+ pairing at the same C,A_ψ is positive
  (≈14.13).  Fail: claim the min is ≥0; claim ker(Torb−4pI)=0.

  Therefore Aut-invariance + master equation + pointwise |m₄|≤1 do
  **not** imply QVAR.  Any proof must use a constraint outside this
  linear 4-point theory (Boolean 6-point / simultaneous ridge
  reconstruction — 15.589 I, 15.590 H for the leftover-3 analogue).

p=7: 128 orbits, ker_dim=7, LP min ≈ −2708 < 0 (same kill; A_ψ HS²
matches q(q−1)/32).  Not a p-law; the p=5 certificate already kills
the method.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15590 import (  # noqa: E402
    enum_max,
    field_ops,
    paley_conference,
    signed_generators,
    signed_orbits,
)


def make_psi(p: int):
    q = p * p
    fmul, fadd, fneg, one = field_ops(p)

    def order_of(e):
        x, o = e, 1
        while x != one:
            x = fmul(x, e)
            o += 1
            if o > q:
                return 0
        return o

    gen = next(e for e in range(2, q) if order_of(e) == q - 1)
    squares = {fmul(x, x) for x in range(1, q)}
    psi = np.zeros(q, dtype=np.complex128)
    x = one
    for i in range(q - 1):
        psi[x] = (1j) ** (i % 4)
        x = fmul(x, gen)
    for t in range(1, q):
        chi = 1 if t in squares else -1
        if abs(psi[t] ** 2 - chi) > 1e-8:
            raise RuntimeError("psi^2 != chi")
    return psi, fmul, fadd, fneg


def A_psi_matrix(p, C, psi, fadd, fneg):
    q = p * p
    n = q + 1
    K = np.zeros((n, n), dtype=np.complex128)
    for a in range(q):
        for b in range(q):
            if a == b:
                continue
            K[1 + a, 1 + b] = psi[fadd(b, fneg(a))]
    P = 0.5 * (np.eye(n) + C.astype(np.float64) / p)
    A = P @ K @ P / 4.0
    return np.real(0.5 * (A + A.conj().T))


def permutation_aut_gens(p: int, C: np.ndarray):
    """Permutations with πᵀ C π = C (no switching signs)."""
    ones = np.ones(C.shape[0], dtype=np.int64)
    gens = []
    for pi, d, s in signed_generators(p, C):
        if s != 1:
            continue
        if np.array_equal(C[np.ix_(pi, pi)], C):
            gens.append((pi, ones, 1))
    if not gens:
        gens = [(pi, ones, 1) for (pi, d, s) in signed_generators(p, C) if s == 1]
    return gens


def aut_box_master_lp(p: int) -> dict:
    """Aut-quotient master LP for min ⟨m, κ_{A_ψ}⟩.  p=5 is the certificate."""
    C = paley_conference(p)
    n = C.shape[0]
    gens = permutation_aut_gens(p, C)
    S4 = np.array(list(combinations(range(n), 4)), dtype=np.int64)
    nS = len(S4)
    lab, sg, dead = signed_orbits(S4, gens, n, twist=False)
    uniq = np.unique(lab)
    remap = {int(u): i for i, u in enumerate(uniq)}
    loc = np.array([remap[int(x)] for x in lab], dtype=np.int64)
    nlab = len(uniq)
    sizes = np.bincount(loc, minlength=nlab).astype(np.float64)
    i, j, k, l = S4[:, 0], S4[:, 1], S4[:, 2], S4[:, 3]
    kap = (
        C[i, j] * C[k, l] + C[i, k] * C[j, l] + C[i, l] * C[k, j]
    ).astype(np.float64)
    psi, fmul, fadd, fneg = make_psi(p)
    A = A_psi_matrix(p, C, psi, fadd, fneg)
    kapA = (
        A[i, j] * A[k, l] + A[i, k] * A[j, l] + A[i, l] * A[k, j]
    ).astype(np.float64)
    enc = (S4[:, 0] * n ** 3 + S4[:, 1] * n ** 2 + S4[:, 2] * n + S4[:, 3]).astype(
        np.int64
    )
    order = np.argsort(enc)
    enc_s = enc[order]

    def lookup(a, b, c, d):
        key = ((a * n + b) * n + c) * n + d
        return int(order[np.searchsorted(enc_s, key)])

    Cf = C.astype(np.float64)
    Torb = np.zeros((nlab, nlab))
    Tkap = np.zeros(nS)
    for t, S in enumerate(S4):
        Sset = set(int(x) for x in S)
        acc = 0.0
        li = loc[t]
        for v in S:
            v = int(v)
            for r in range(n):
                if r in Sset:
                    continue
                Sp = tuple(sorted([int(x) for x in S if x != v] + [r]))
                u = lookup(*Sp)
                w = Cf[v, r]
                acc += w * kap[u]
                Torb[li, loc[u]] += w
        Tkap[t] = acc
    Torb /= sizes[:, None]
    star = -Tkap / 6.0

    def oavg(x):
        return np.bincount(loc, weights=x, minlength=nlab) / np.maximum(sizes, 1)

    avg_k = oavg(kap)
    avg_s = oavg(star)
    avg_A = oavg(kapA)
    kap_std = max(
        float(np.std(kap[loc == t])) if np.any(loc == t) else 0.0
        for t in range(nlab)
    )
    torb_err = float(np.max(np.abs(Torb @ avg_k + 6.0 * avg_s)))
    fourp = 4.0 * p
    M = Torb - fourp * np.eye(nlab)
    rhs = -4.0 * avg_k / p
    part, *_ = np.linalg.lstsq(M, rhs, rcond=1e-10)
    _u, svals, vt = np.linalg.svd(M)
    ker_dim = int(np.sum(svals < 1e-8))
    K = vt[-ker_dim:].T if ker_dim else np.zeros((nlab, 0))
    obj_const = float(np.dot(sizes * avg_A, part))
    if ker_dim == 0:
        min_pairing = obj_const
        max_pairing = obj_const
        ok = True
    else:
        c = K.T @ (sizes * avg_A)
        A_ub = np.vstack([K, -K])
        b_ub = np.concatenate([1.0 - part, 1.0 + part])
        bounds = [(None, None)] * ker_dim
        res_min = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        res_max = linprog(-c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
        ok = bool(res_min.success and res_max.success)
        min_pairing = float(obj_const + res_min.fun) if res_min.success else None
        max_pairing = float(obj_const - res_max.fun) if res_max.success else None
    return {
        "p": p,
        "n_orbits": nlab,
        "dead": len(dead),
        "kap_orbit_std": kap_std,
        "Torb_kappa_err": torb_err,
        "ker_dim": ker_dim,
        "min_pairing": min_pairing,
        "max_pairing": max_pairing,
        "lp_ok": ok,
        "HS2": float(np.sum(A * A)),
        "named_HS2": float((p * p) * (p * p - 1) / 32),
    }


def true_maxplus_pairing(p: int) -> dict:
    C = paley_conference(p)
    n = C.shape[0]
    Yp = enum_max(p, +1)
    S4 = np.array(list(combinations(range(n), 4)))
    i, j, k, l = S4[:, 0], S4[:, 1], S4[:, 2], S4[:, 3]
    m4 = np.mean(
        Yp[:, i].astype(np.float64)
        * Yp[:, j]
        * Yp[:, k]
        * Yp[:, l],
        axis=0,
    )
    psi, fmul, fadd, fneg = make_psi(p)
    A = A_psi_matrix(p, C, psi, fadd, fneg)
    kapA = (
        A[i, j] * A[k, l] + A[i, k] * A[j, l] + A[i, l] * A[k, j]
    )
    pairing = float(np.dot(m4, kapA))
    quad = (Yp.astype(np.float64) @ A * Yp.astype(np.float64)).sum(1)
    e2 = float(np.mean(quad ** 2))
    thr = 3 * (p * p) * (p * p - 1) / 16
    return {
        "p": p,
        "Nplus": int(len(Yp)),
        "pairing": pairing,
        "E_yAy_sq": e2,
        "QVAR_threshold": thr,
        "clears_QVAR": e2 >= thr,
        "max_abs_m4": float(np.max(np.abs(m4))),
    }


def theorem_box_master_aut_cannot_prove_qvar() -> dict:
    """Proved kill of master+box+Aut as a QVAR method.  Inequality still OPEN."""
    lp = aut_box_master_lp(5)
    tru = true_maxplus_pairing(5)
    min_p = lp["min_pairing"]
    ok = (
        lp["lp_ok"]
        and lp["Torb_kappa_err"] < 1e-10
        and lp["kap_orbit_std"] < 1e-10
        and lp["ker_dim"] >= 1
        and min_p is not None
        and min_p < 0
        and tru["pairing"] > 0
        and tru["clears_QVAR"]
        and tru["pairing"] > min_p
        and tru["max_abs_m4"] <= 1.0 + 1e-9
    )
    claim_min_ge_0 = bool(min_p is not None and min_p >= 0)
    return {
        "proved": bool(ok),
        "inequality_proved": False,
        "method_killed": "Aut-inv master equation + |m4|<=1",
        "claim_lp_min_ge_0": claim_min_ge_0,
        "p5_lp": lp,
        "p5_true": tru,
        "p5_lp_min_is_neg285_over_4": bool(
            min_p is not None and abs(min_p - (-285 / 4)) < 1e-8
        ),
        "theorem": (
            "At p=5 the Aut-quotient of T m=4p m−4κ/p with |m|≤1 is a "
            "nonempty polytope of E_{4p}^{Aut}-dimension 2, and "
            "min ⟨m,κ_{A_ψ}⟩=−285/4<0, while true Max+ pairing is "
            "positive.  Fail: claim that min is ≥0.  Hence Aut + "
            "4-point master + box cannot prove QVAR."
        ),
    }


def main() -> dict:
    from io_atomic import write_json_atomic

    T = theorem_box_master_aut_cannot_prove_qvar()
    out = {
        "title": "box+master+Aut cannot prove GLOBAL QVAR",
        "numbered": False,
        "theorem": T,
        "global_qvar_not_claimed": True,
    }
    path = ROOT / "evidence" / "e1_gmin_qvar_box_master.json"
    write_json_atomic(path, out)
    print("box+master+Aut QVAR kill", flush=True)
    print(f"  proved_kill={T['proved']} inequality={T['inequality_proved']}", flush=True)
    print(f"  claim_min_ge_0={T['claim_lp_min_ge_0']}", flush=True)
    print(f"  p5 LP min={T['p5_lp']['min_pairing']} true={T['p5_true']['pairing']}", flush=True)
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
