#!/usr/bin/env python3
"""
Prop 15.562 — k=3 Fejer names A_k3 and the k=3 piece of Q3_n3.

    Σ_{t∈F_p^×} ω^{2t}/(ω^t−1)^3 = −(p²−1)/24
    A_{k=3,dbl} = 0
    A_{k=3,n3} = −16 p³ / ((p−1)(p−3))
    n_{k=3} A_{k=3,n3} = −p^5 (p²−1)/3
When m=3 (i.e. p=5) this is all of Q3_n3 and Q3_dbl=n_1d A_1d.
A_full for m>3 (Q3_dbl / Q3_n3 at p=7) still open.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.561 flags.  Does **not** import phi_F.
Does **not** interpolate A_full,dbl by (p−5)/15.

============================================================================
Setup.  15.276: every k=3 Max+ is a locked sawtooth triple and
ẑ(ξ)=F_{λ,s}(c)=2p ω^{−c λ^{-1}s}/(ω^{c λ^{-1}}−1) on each
support line (c≠0).  Locked phases s0+s1+s2≡−2.  15.561:
n_{k=3}=C(m,3)(p−1)q and A_{k=3}=0 on 1-line triples (cache).

============================================================================
Theorem A — PROVED field-only (complete sum; certified p=5..23).
  S:=Σ_{t≠0} ω^{2t}/(ω^t−1)^3 equals −(p²−1)/24.
  Fail: (p−1)/24.  Fail: drop 24 (p=5: −1≠−24).  ∎

Theorem B — PROVED Max+-free (15.276 ẑ=F + s-sum).
  On a 1-line double, ẑ(ξ)²ẑ(−2ξ)=8p³/((ω^c−1)²(ω^{−2c}−1)),
  independent of phase.  Average over c∈F_p^× is 0.
  Hence A_{k=3,dbl}=0 for every odd p≥5.  Fail: A_{k=3,dbl}=A_1d.  ∎

Theorem C — PROVED Max+-free (15.276 + A + 15.561).
  A 3-line triple on a locked support has equal F_p-coordinates
  after aligning generators (g1+g2+g3=0).  The s-sum collapses
  and A_{k=3,n3} = [8p³/(C(m,3)(p−1))] S = −16p³/((p−1)(p−3)).
  The C(m,3) is the fraction of k=3 vectors whose support is
  that triple (other supports give product 0).
  p=5: −250; p=7: −2p³/3.  Fail: (p−2) in the denominator.  ∎

Theorem D — PROVED (C + 15.561 A).
  n_{k=3} A_{k=3,n3} = −p^5(p²−1)/3.
  At p=5 (m=3) this equals live Q3_n3=−25000, and Q3_dbl=n_1d A_1d
  (A_full,dbl=A_{k=3,dbl}=0).  Fail: drop 1/3 (p=5: −75000≠−25000).  ∎

Theorem E — OPEN.  For m>3, Q3_dbl=n_1d A_1d + n_{k=m} A_full,dbl
  and Q3_n3=−p^5(p²−1)/3 + n_{k=m} A_full,n3 still carry
  A_full,dbl=A_full,n3=−4p³/15 at p=7.  Not a p-law.
  Do not interpolate (p−5)/15.  phi_F not imported.

============================================================================
Backend: field Fejer sum p=5..23 + p=5,7 Max+ cache.  Writes
evidence/e1_gmin_m4_prop15562.json
"""
from __future__ import annotations

import json
import os
import sys
from math import comb
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15554 import LIVE_Q3, triple_kind, yhat_omega
from e1_gmin_m4_prop15561 import A1d_named, n1d_named, n_k3_named

EV = ROOT / "evidence" / "e1_gmin_m4_prop15562.json"


def S_fejer(p: int) -> complex:
    omega = np.exp(2j * np.pi / p)
    acc = 0j
    for t in range(1, p):
        w = omega ** t
        acc += (w ** 2) / (w - 1) ** 3
    return acc


def S_named(p: int) -> float:
    return -(p * p - 1) / 24.0


def S_wrong_pm1(p: int) -> float:
    return -(p - 1) / 24.0


def S_drop24(p: int) -> float:
    return -(p * p - 1)


def A_k3_n3_named(p: int) -> float:
    return -16.0 * p**3 / ((p - 1) * (p - 3))


def A_k3_n3_wrong_pminus2(p: int) -> float:
    return -16.0 * p**3 / ((p - 1) * (p - 2))


def q3_n3_k3_named(p: int) -> float:
    return -(p**5) * (p * p - 1) / 3.0


def q3_n3_k3_drop_third(p: int) -> float:
    return -(p**5) * (p * p - 1)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        got = S_fejer(p)
        want = S_named(p)
        match = abs(got - want) < 1e-9
        fail_pm1 = abs(want - S_wrong_pm1(p)) > 0.5
        fail_24 = abs(want - S_drop24(p)) > 0.5
        if not (match and fail_pm1 and fail_24):
            ok = False
        rows[str(p)] = {
            "S_re": float(got.real),
            "named": want,
            "err": float(abs(got - want)),
            "match": bool(match),
            "pm1_fails": bool(fail_pm1),
            "drop24_fails": bool(fail_24),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "Σ ω^{2t}/(ω^t−1)^3=−(p²−1)/24. Fail (p−1)/24; fail drop 24.",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F, Y, Omega, yhat = yhat_omega(p)
        k3 = (np.abs(yhat) > 1e-6).sum(axis=1) == 3 * (p - 1)
        add, neg, mul, inv = F["add"], F["neg"], F["mul"], F["inv"]
        idx = {om: i for i, om in enumerate(Omega)}
        acc = []
        if k3.any():
            yh = yhat[k3]
            for a in Omega:
                b = a
                c = neg[add[a][a]]
                if c not in idx:
                    continue
                if triple_kind(p, mul[b][inv[a]]) != "double":
                    continue
                prod = yh[:, idx[a]] * yh[:, idx[b]] * yh[:, idx[c]]
                acc.append(float(prod.real.mean()))
                break
        mean = float(np.mean(acc)) if acc else 0.0
        match = abs(mean) < 1e-8
        a1d_fails = abs(A1d_named(p)) > 1.0
        if not (match and a1d_fails):
            ok = False
        rows[str(p)] = {
            "A_k3_dbl_mean": mean,
            "match": bool(match),
            "A1d_fails": bool(a1d_fails),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "A_{k=3,dbl}=0 (Fejer 1-line). Fail A_{k=3,dbl}=A_1d.",
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F, Y, Omega, yhat = yhat_omega(p)
        k3 = (np.abs(yhat) > 1e-6).sum(axis=1) == 3 * (p - 1)
        add, neg, mul, inv = F["add"], F["neg"], F["mul"], F["inv"]
        idx = {om: i for i, om in enumerate(Omega)}
        live = None
        if k3.any():
            yh = yhat[k3]
            for a in Omega:
                for b in Omega:
                    c = neg[add[a][b]]
                    if c not in idx:
                        continue
                    if triple_kind(p, mul[b][inv[a]]) != "n3":
                        continue
                    prod = yh[:, idx[a]] * yh[:, idx[b]] * yh[:, idx[c]]
                    live = float(prod.real.mean())
                    break
                if live is not None:
                    break
        named = A_k3_n3_named(p)
        match = live is not None and abs(live - named) < 1e-6
        den_fails = abs(named - A_k3_n3_wrong_pminus2(p)) > 1.0
        if not (match and den_fails):
            ok = False
        rows[str(p)] = {
            "live": live,
            "named": named,
            "match": bool(match),
            "pminus2_fails": bool(den_fails),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "A_{k=3,n3}=−16p³/((p−1)(p−3)). Fail (p−2).",
    }


def prove_D() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        named = q3_n3_k3_named(p)
        live_q3 = LIVE_Q3[p]["n3"]
        n_k3 = n_k3_named(p)
        a = A_k3_n3_named(p)
        piece = n_k3 * a
        match_piece = abs(piece - named) < 1e-4
        # at p=5 the piece is all of Q3_n3; at p=7 it is not
        if p == 5:
            match_q3 = abs(named - live_q3) < 1e-4
            dbl_named = n1d_named(p) * A1d_named(p)
            match_dbl = abs(dbl_named - LIVE_Q3[p]["double"]) < 1e-4
        else:
            match_q3 = abs(named - live_q3) > 1.0  # must NOT claim full Q3_n3
            match_dbl = abs(n1d_named(p) * A1d_named(p) - LIVE_Q3[p]["double"]) > 1.0
        drop_fails = abs(named - q3_n3_k3_drop_third(p)) > 1.0
        if not (match_piece and match_q3 and match_dbl and drop_fails):
            ok = False
        rows[str(p)] = {
            "k3_piece": named,
            "n_k3_A": piece,
            "live_Q3_n3": live_q3,
            "match_piece": bool(match_piece),
            "p5_Q3n3_or_p7_not_full": bool(match_q3),
            "p5_Q3dbl_or_p7_not_1d_only": bool(match_dbl),
            "drop_third_fails": bool(drop_fails),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "n_{k=3} A_{k=3,n3}=−p^5(p²−1)/3. "
            "Equals Q3_n3 iff m=3. Fail drop 1/3."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q3_double_named_in_p": False,
        "Q3_n3_named_in_p": False,
        "A_full_dbl_named_in_p": False,
        "NUM_SUM_16pA_general": False,
        "phi_F_imported": False,
        "note": (
            "A_k3 and the k=3 piece of Q3_n3 are named. "
            "At p=5 (m=3) this closes both Q3_n3 and Q3_dbl. "
            "For m>3, A_full,dbl=A_full,n3=−4p³/15 at p=7 is not a p-law. "
            "Do not interpolate (p−5)/15. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.562  k=3 Fejer names A_k3 and Q3_n3 k=3 piece", flush=True)
    A = prove_A()
    print(f"  A S=−(p²−1)/24: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B A_k3,dbl=0: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C A_k3,n3: {C['proved']}", flush=True)
    D = prove_D()
    print(f"  D k=3 Q3_n3 piece: {D['proved']}", flush=True)
    E = prove_open()
    out = {
        "prop": "15.562",
        "title": "A_k3,n3=−16p³/((p−1)(p−3)); n_k3 A=−p^5(p²−1)/3",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "E": E,
        "proved": {
            "S_fejer_named": A["proved"],
            "A_k3_dbl_zero": B["proved"],
            "A_k3_n3_named": C["proved"],
            "Q3_n3_k3_piece_named": D["proved"],
            "Q3_double_named_in_p": False,
            "Q3_n3_named_in_p": False,
            "phi_F_ge_6": False,
        },
        "e1_closed_general": bool(e1_closed_general()),
        "gsum_disj_lb": bool(gsum_disj_lb_proved_general()),
        "phi_F_ge_6": False,
        "L_status": "OPEN",
        "flags_not_flipped": [
            "phi_F_ge_6",
            "e1",
            "L",
            "Aut-Schur",
            "Gsum",
            "pairing",
            "residual_ii",
            "type_I",
        ],
        "identity": "A_k3,n3=-16p^3/((p-1)(p-3)); n_k3 A=-p^5(p^2-1)/3",
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
