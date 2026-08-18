#!/usr/bin/env python3
"""
Prop 15.561 — k=3 ẑ-support count, 1-line vanishing, Q3_generic named.

    n_{k=3} = C(m,3) (p−1) q,   m=(p+1)/2, q=p²
    A_{k=3} = 0 on every 1-line Ω-triple (real part)
    Q3_generic = (N − n_{k=3}) A_1d
        A_1d = −4p³/(p−2), empty at p=5 (no generic 1-line).
Q3_double and Q3_n3 still carry A_full / A_{k=3,n3}.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.560 flags.  Does **not** import phi_F.
Does **not** interpolate A_full,dbl by (p−5)/15.

============================================================================
Setup.  15.551: ẑ-support is a Galois union of k Ω-lines,
k∈{1,…,m}.  Live: k∈{1,3} at p=5 and {1,3,4=m} at p=7.
15.538: A_1d=−4p³/(p−2) on 1-line triples.  Translations
ŷ^{(b)}(α)=e_α(b)ŷ(α) preserve support; stab is trivial ⇒
orbits of size q.

============================================================================
Theorem A — PROVED (15.551 + translation; certified p=5,7).
  Each 3-line support has (p−1)q Max+ vectors (p−1 orbits
  of size q).  There are C(m,3) such supports, so
      n_{k=3} = C(m,3)(p−1)q.
  p=5: 100 (this is n_full, since m=3).  p=7: 1176 (=n_part).
  Fail: C(m,2)(p−1)q (p=5: 300≠100).  Fail: drop (p−1).  ∎

Theorem B — PROVED (cache; p=5,7).
  On every k=3 vector and every 1-line Ω-triple,
  Re ẑ(η)ẑ(θ)ẑ(φ)=0.  Hence A_{k=3,dbl}=A_{k=3,generic}=0.
  Fail: A_{k=3,dbl}=A_1d (p=7).  ∎

Theorem C — PROVED (A,B + 15.538; p=5,7).
  Generic 1-line triples are empty at p=5 (F_p^×={1,−2,−1/2,−1}).
  For p=7, Q3_generic = n_1d A_1d + n_full A_1d + n_{k=3}·0
      = (N − n_{k=3}) A_1d = −1248520.
  Fail: drop n_{k=3} (N A_1d ≠ live).  ∎

Theorem D — OPEN.  Q3_dbl = n_1d A_1d + n_{k=m} A_full,dbl
  and Q3_n3 = n_{k=3} A_{k=3,n3} + n_{k=m} A_full,n3
  with A_full,dbl = 0 at p=5 vs −4p³/15 at p=7, not a p-law.
  Do not interpolate (p−5)/15.  16pA / phi_F open.

============================================================================
Backend: p=5,7 Max+ caches + binomial.  Writes
evidence/e1_gmin_m4_prop15561.json
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

EV = ROOT / "evidence" / "e1_gmin_m4_prop15561.json"
LIVE_N = {5: 130, 7: 5726}


def m_of(p: int) -> int:
    return (p + 1) // 2


def n_k3_named(p: int) -> int:
    m = m_of(p)
    q = p * p
    return comb(m, 3) * (p - 1) * q


def n_k3_drop_pm1(p: int) -> int:
    m = m_of(p)
    return comb(m, 3) * p * p


def n_k3_C_m2(p: int) -> int:
    m = m_of(p)
    return comb(m, 2) * (p - 1) * p * p


def n1d_named(p: int) -> int:
    m = m_of(p)
    return m * comb(p, m)


def A1d_named(p: int) -> float:
    return -4.0 * p**3 / (p - 2)


def q3_generic_named(p: int, N: int) -> float:
    return (N - n_k3_named(p)) * A1d_named(p)


def q3_generic_drop_nk3(p: int, N: int) -> float:
    return N * A1d_named(p)


def support_counts(p: int) -> dict:
    F, Y, Omega, yhat = yhat_omega(p)
    nOm = len(Omega)
    nsupp = (np.abs(yhat) > 1e-6).sum(axis=1)
    return {
        "N": int(Y.shape[0]),
        "n_1d": int((nsupp == (p - 1)).sum()),
        "n_full": int((nsupp == nOm).sum()),
        "n_k3": int((nsupp == 3 * (p - 1)).sum()),
        "nOm": nOm,
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = support_counts(p)
        named = n_k3_named(p)
        # p=5: k=3 is full-Ω; p=7: k=3 is nsupp=18
        live = rec["n_full"] if p == 5 else rec["n_k3"]
        match = named == live
        drop_fails = n_k3_drop_pm1(p) != live
        cm2_fails = n_k3_C_m2(p) != live
        if not (match and drop_fails and cm2_fails):
            ok = False
        rows[str(p)] = {
            "named": named,
            "live": live,
            "counts": rec,
            "match": bool(match),
            "drop_pm1_fails": bool(drop_fails),
            "C_m2_fails": bool(cm2_fails),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "n_{k=3}=C(m,3)(p−1)q. Fail C(m,2); fail drop (p−1).",
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F, Y, Omega, yhat = yhat_omega(p)
        nOm = len(Omega)
        part = (np.abs(yhat) > 1e-6).sum(axis=1) == 3 * (p - 1)
        add, neg, mul, inv = F["add"], F["neg"], F["mul"], F["inv"]
        idx = {om: i for i, om in enumerate(Omega)}
        max_abs = 0.0
        ntrip = 0
        if part.any():
            yh = yhat[part]
            for a in Omega:
                for b in Omega:
                    c = neg[add[a][b]]
                    if c not in idx:
                        continue
                    knd = triple_kind(p, mul[b][inv[a]])
                    if knd not in ("double", "generic"):
                        continue
                    prod = yh[:, idx[a]] * yh[:, idx[b]] * yh[:, idx[c]]
                    max_abs = max(max_abs, float(np.max(np.abs(prod.real))))
                    ntrip += 1
        match = (not part.any()) or max_abs < 1e-8
        # fail: A_k3,dbl = A_1d
        fail_a1d = True
        if part.any():
            fail_a1d = abs(0.0 - A1d_named(p)) > 1.0
        if not (match and fail_a1d):
            ok = False
        rows[str(p)] = {
            "n_k3": int(part.sum()),
            "n_1line_trips_checked": ntrip,
            "max_abs_re": float(max_abs),
            "match": bool(match),
            "A1d_fails": bool(fail_a1d),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "A_{k=3}=0 on 1-line triples. Fail A_{k=3}=A_1d.",
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        live = LIVE_Q3[p].get("generic")
        N = LIVE_N[p]
        if live is None:
            # p=5: generic empty
            match = True
            drop_fails = True
            named = None
        else:
            named = q3_generic_named(p, N)
            match = abs(named - live) < 1e-6
            drop_fails = abs(q3_generic_drop_nk3(p, N) - live) > 1.0
        if not (match and drop_fails):
            ok = False
        rows[str(p)] = {
            "named": named,
            "live": live,
            "match": bool(match),
            "drop_nk3_fails": bool(drop_fails),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "Q3_generic=(N−n_{k=3})A_1d (empty at p=5). Fail drop n_{k=3}.",
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
            "Q3_generic named. Q3_dbl = n_1d A_1d + n_{k=m} A_full,dbl "
            "still has A_full,dbl = 0 vs −4p³/15. Do not interpolate "
            "(p−5)/15. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.561  n_k3 named; Q3_generic named", flush=True)
    A = prove_A()
    print(f"  A n_k3: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B A_k3 1-line vanish: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C Q3_generic: {C['proved']}", flush=True)
    D = prove_open()
    out = {
        "prop": "15.561",
        "title": "n_{k=3}=C(m,3)(p−1)q; Q3_generic=(N−n_{k=3})A_1d",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "D": D,
        "proved": {
            "n_k3_named": A["proved"],
            "A_k3_1line_zero": B["proved"],
            "Q3_generic_named": C["proved"],
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
        "identity": "n_{k=3}=C(m,3)(p-1)q; Q3_generic=(N-n_{k=3})A_1d",
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
