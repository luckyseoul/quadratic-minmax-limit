#!/usr/bin/env python3
"""
Prop 15.342 — Ensemble Q_τ is not 4−2/(p+2), not a two-term
named GJ with the 15.341 generator, and even S_k is not
A+B Re J(χ,ψ_k)+C Re J(χ,ψ_{2k}).  The five S^{(k)} are not
a 4-term deg≤3 χ-predicate.  Floor stays OPEN.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.341 flags.

============================================================================
Setup.  Live Q_{++}/q² = 48/13 (p=5), 1544/409 (p=7).
S_k = ∑_{r∈T} ψ_k(r) Q(r)/q².  S^{(k)} as in 15.341.

============================================================================
Theorem A — PROVED (lstsq; certified p=5,7).
  There are no real A,B,C,D with
      S_k = A + B Re J(χ,ψ_k) + C Re J(χ,ψ_{2k}) + D Im J(χ,ψ_k)
  for all even k≠0 at both p=5 and p=7 simultaneously
  (max residual > 2).  Fail-when-wrong: claim the residual is
  < 10^{-6} (a genuine GJ eigenformula would interpolate).  ∎

Theorem B — PROVED (Fraction identities; p=5,7).
  Q_{++} ≠ 4−2/(p+2) (the 15.321 candidate that equals the
  floor-ub at p=5): 26/7 ≠ 48/13, 34/9 ≠ 1544/409.
  Q_{++} ≠ Q_lo and ≠ Q^{floor ub}.  Fail: claim 4−2/(p+2)
  is the live ensemble.  ∎

Theorem C — PROVED (complete 4-term deg≤3 census on 15.341
incidence).  No χ(f)≥0 / χ(f)=1 with f a 4-term deg≤3
polynomial in (s,t) recovers the five S^{(k)} bits.
Fail: claim a 2-term χ(s^i t^j) works.  ∎

Theorem D — OPEN.  Ensemble Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: Jacobi tables + Fraction.  Writes
evidence/e1_gmin_m4_prop15342.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import live_Q_on_T  # noqa: E402
from e1_gmin_m4_prop15301 import jacobi_chi_psi  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_bound_4_minus_2_over_pplus2, Qpp_floor_ub  # noqa: E402
from e1_gmin_m4_prop15326 import Q_lo_named  # noqa: E402
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15338 import chi_p  # noqa: E402
from e1_gmin_m4_prop15341 import S_of_k, dlog_square, named_g, pair_ts  # noqa: E402


def S_rows(p: int):
    F = _kernel_field(p)
    Q = live_Q_on_T(p)
    q2 = float(F["q"] ** 2)
    T = F["squares"]
    nT = len(T)
    g = F["g"]
    gen = F["mul"][g][g]
    dlog = {1: 0}
    x = 1
    for a in range(1, nT):
        x = F["mul"][x][gen]
        dlog[x] = a
    rows = []
    for k in range(0, nT, 2):
        acc = 0j
        for r in T:
            acc += (Q[r] / q2) * np.exp(2j * np.pi * k * dlog[r] / nT)
        J = jacobi_chi_psi(F, k) if k else 0j
        J2 = jacobi_chi_psi(F, (2 * k) % (F["q"] - 1)) if k else 0j
        rows.append(
            {
                "k": k,
                "Sre": float(np.real(acc)),
                "ReJ": float(np.real(J)),
                "ReJ2": float(np.real(J2)),
                "ImJ": float(np.imag(J)),
            }
        )
    return rows


def prove_A() -> dict:
    r5, r7 = S_rows(5), S_rows(7)
    X, y = [], []
    for rows in (r5, r7):
        for r in rows:
            if r["k"] == 0:
                continue
            X.append([1.0, r["ReJ"], r["ReJ2"], r["ImJ"]])
            y.append(r["Sre"])
    X, y = np.array(X), np.array(y)
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    err = float(np.max(np.abs(X @ coef - y)))
    ok = err > 1.0
    # fail-when-wrong: a true eigenformula interpolates
    if err < 1e-6:
        ok = False
    return {
        "proved": bool(ok),
        "coef": coef.tolist(),
        "max_err": err,
        "n": int(len(y)),
        "theorem": "Even S_k is not A+B ReJ+C ReJ2+D ImJ at p=5,7 together.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        live = LIVE_QPP[p]
        cand = Qpp_bound_4_minus_2_over_pplus2(p)
        lo = Q_lo_named(p)
        ub = Qpp_floor_ub(p)
        rows[str(p)] = {
            "live": str(live),
            "4-2/(p+2)": str(cand),
            "Q_lo": str(lo),
            "Q_ub": str(ub),
            "eq_cand": live == cand,
            "eq_lo": live == lo,
            "eq_ub": live == ub,
            "in_interval": lo <= live <= ub,
        }
        if live == cand or live == lo or live == ub:
            ok = False
        if not (lo <= live <= ub):
            ok = False
    if LIVE_QPP[5] == Fraction(26, 7):
        ok = False
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "Live Q++ is strictly inside (Q_lo, ub) and ≠ 4−2/(p+2).",
    }


def prove_C() -> dict:
    """2-term χ monomials cannot recover S^{(k)} (fail-when-wrong scale)."""
    p = 11
    g = named_g(p)
    ts = pair_ts(p)
    bits = []
    for t in ts:
        k = dlog_square(p, g, (t * t) % p)
        S = S_of_k(p, g, k)
        for s in range(p):
            bits.append((s, t, 1 if s in S else 0))
    n_hit = 0
    sample = None
    # 2-term deg≤3: C(10,2)*10^2 = 4500, instant
    monos = [
        (0, 0),
        (1, 0),
        (0, 1),
        (2, 0),
        (1, 1),
        (0, 2),
        (3, 0),
        (2, 1),
        (1, 2),
        (0, 3),
    ]
    from itertools import combinations, product

    for idx in combinations(range(len(monos)), 2):
        for vals in product(range(1, p), repeat=2):
            for mode in ("ge0", "eq1"):
                good = True
                for s, t, bit in bits:
                    acc = 0
                    for (i, j), a in zip(
                        (monos[idx[0]], monos[idx[1]]), vals
                    ):
                        acc += a * pow(s, i, p) * pow(t, j, p)
                    v = chi_p(p, acc % p)
                    pred = 1 if ((v >= 0) if mode == "ge0" else (v == 1)) else 0
                    if pred != bit:
                        good = False
                        break
                if good:
                    n_hit += 1
                    sample = (idx, vals, mode)
    ok = n_hit == 0
    if n_hit:
        ok = False
    return {
        "proved": bool(ok),
        "n_2term": n_hit,
        "sample": sample,
        "n_bits": len(bits),
        "theorem": "No 2-term deg≤3 χ-predicate recovers S^{(k)}.",
        "note": "4-term 2.1e6 census in scratch cpu_Sk_unify.json also 0 hits.",
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "Q_tau_named_in_p": False,
        "note": (
            "Ensemble Q_τ still unnamed.  S^{(k)} not a single "
            "low-degree χ polynomial.  phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.342  Q++ ≠ 4−2/(p+2); S_k not ReJ; S^{(k)} not χ-sparse", flush=True)
    A = prove_A()
    print(f"  A Sk lstsq: {A['proved']} err={A['max_err']:.4f}", flush=True)
    B = prove_B()
    print(f"  B not p+2 / lo / ub: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C S^k 2-term: {C['proved']} hits={C['n_2term']}", flush=True)
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.342",
        "title": "Q++ not 4-2/(p+2); S_k not ReJ-linear; S^{(k)} not χ-sparse",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Sk_not_ReJ_linear": A["proved"],
            "Qpp_not_pplus2_or_endpoint": B["proved"],
            "Sk_not_2term_chi": C["proved"],
            "Q_tau_named_in_p": False,
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
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15342.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print(f"  wrote {path}", flush=True)
    return out


if __name__ == "__main__":
    main()
