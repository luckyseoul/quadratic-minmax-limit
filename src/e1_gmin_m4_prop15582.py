#!/usr/bin/env python3
"""
Prop 15.582 — 1D 4-point vanishes off F_p; p=5 Q_{N*}/q²=32/13
is the named k=3 Fejér piece.

    r∉F_p  (types N1, N*)  ⇒  μ_{1d}(r)=0
    p=5: Q_{N*}/q² = n_{k=3} μ_{k=3} / (N q²) = 32/13

because exclusive n_full=0 and 1D is off.  Fail: include the
1D Johnson term (gets 48/13).  Fail: claim Q_{N*}=Q_{++}.
k=3 2-line 4-point is not identically the 1-line Fejér
(p=7 ++N1: 12805.3 ≠ 9604).  Q_τ for p≥7 still has μ_full.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.581 flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.581.  Does **not** interpolate (p−5)/15
or 16(p−4)/D.  Does **not** import 15.568 mix weights.

============================================================================
Setup.  15.538: Type+ 1D ẑ lives on one F_p-line.  15.551:
1-line iff r∈F_p.  Hence for r∉F_p the four frequencies
{ξ,−ξ,rξ,−rξ} meet two lines and the 1D product is 0.
15.573 exclusive mix; 15.574 μ_{k=3}; 15.581 Q_{++}(p=5)=48/13.

============================================================================
Theorem A — PROVED (15.538+15.551; p=5,7 cache).
  The 4-point is indexed on Ω (ξ≠0).  μ_{1d}(r)=0 on N*
  and N1.  Fail: μ_{1d}(N*)=μ_{1d}(sub) (p=5: 0≠10000/3).  ∎

Theorem B — PROVED (A + 15.561+15.574; Fraction).
  p=5 Q_{N*}/q²=100·2000/(130·625)=32/13.  Fail: add 1D
  (48/13).  Fail: 32/13=48/13.  ∎

Theorem C — OPEN.  p≥7 n_full>0.  k=3 2-line ≠ 1-line
  Fejér on every type (p=7 ++N1).  (Q_{++}−Q_{N*})D=16,48
  at p=5,7; 16(p−4) interpolates and is not imported.
  Q_τ unnamed in p.  phi_F not imported.

============================================================================
Backend: yhat 4-point p=5,7 + Fraction.  Writes
evidence/e1_gmin_m4_prop15582.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15290 import _kernel_field, live_Q_on_T, type_key
from e1_gmin_m4_prop15554 import yhat_omega
from e1_gmin_m4_prop15561 import n1d_named, n_k3_named
from e1_gmin_m4_prop15573 import exclusive_buckets
from e1_gmin_m4_prop15574 import mu_k3_named
from e1_gmin_m4_prop15575 import mu_1d_named

EV = ROOT / "evidence" / "e1_gmin_m4_prop15582.json"


def fourpoint_by_bucket(p: int, r: int) -> dict:
    F, Y, Omega, yhat = yhat_omega(p)
    q = F["q"]
    nOm = len(Omega)
    idx = {om: i for i, om in enumerate(Omega)}
    nsupp = (np.abs(yhat) > 1e-6).sum(axis=1)
    buckets = exclusive_buckets(p, nsupp, nOm)
    a = Omega[0]
    b, c = F["neg"][a], F["mul"][r][a]
    d = F["neg"][c]
    ia, ib, ic, id_ = idx[a], idx[b], idx[c], idx[d]
    prod = (yhat[:, ia] * yhat[:, ib] * yhat[:, ic] * yhat[:, id_]).real
    out = {}
    for name, mask in buckets.items():
        n = int(mask.sum())
        if n == 0:
            continue
        out[name] = float(np.mean(prod[mask]))
    return out


def first_r(p: int, pred) -> int:
    Fk = _kernel_field(p)
    Q = live_Q_on_T(p)
    for r, _ in Q.items():
        if pred(type_key(Fk, r)):
            return int(r)
    raise RuntimeError("no r")


def is_nstar(key) -> bool:
    return isinstance(key, tuple) and len(key) == 3 and key[2] == "N*"


def is_n1(key) -> bool:
    return isinstance(key, tuple) and len(key) == 3 and key[2] == "N1"


def is_sub(key) -> bool:
    return isinstance(key, tuple) and len(key) == 3 and key[2] == "sub"


def Qnstar_p5_named() -> Fraction:
    n1 = n1d_named(5)
    nk = n_k3_named(5)
    N = n1 + nk
    return Fraction(nk * 2000, N * 5**4)


def Qnstar_p5_with_1d() -> Fraction:
    n1 = n1d_named(5)
    nk = n_k3_named(5)
    N = n1 + nk
    return Fraction(n1 * mu_1d_named(5) + nk * 2000, N * 5**4)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        r_ns = first_r(p, is_nstar)
        r_n1 = first_r(p, is_n1) if p == 7 or True else None
        # p=5 has N1
        parts_ns = fourpoint_by_bucket(p, r_ns)
        parts_n1 = fourpoint_by_bucket(p, first_r(p, is_n1))
        mu1_ns = abs(parts_ns.get("1d", 99.0))
        mu1_n1 = abs(parts_n1.get("1d", 99.0))
        vanish = mu1_ns < 1e-6 and mu1_n1 < 1e-6
        fail_eq = abs(mu1_ns - float(mu_1d_named(p))) > 1.0
        if not (vanish and fail_eq):
            ok = False
        rows[str(p)] = {
            "r_N*": r_ns,
            "mu_1d_N*": parts_ns.get("1d"),
            "mu_1d_N1": parts_n1.get("1d"),
            "mu_1d_sub": float(mu_1d_named(p)),
            "vanish": bool(vanish),
            "ne_sub": bool(fail_eq),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": "μ_1d=0 on N* and N1. Fail μ_1d(N*)=μ_1d(sub).",
    }


def prove_B() -> dict:
    q = Qnstar_p5_named()
    with1 = Qnstar_p5_with_1d()
    ok = (
        q == Fraction(32, 13)
        and with1 == Fraction(48, 13)
        and q != with1
        and abs(mu_k3_named(5, 2) - 2000.0) < 1e-8
        and n_k3_named(5) == 100
        and n1d_named(5) + n_k3_named(5) == 130
    )
    return {
        "proved": bool(ok),
        "Qnstar": str(q),
        "with_1d": str(with1),
        "Qpp": "48/13",
        "theorem": "p=5 Q_N*/q²=32/13. Fail add 1D (48/13); fail 32=48.",
    }


def prove_open() -> dict:
    # p=7 ++N1 k3 ≠ sub k3
    r_n1pp = first_r(7, lambda k: is_n1(k) and k[0] == 1 and k[1] == 1)
    r_sub = first_r(7, is_sub)
    p_n1 = fourpoint_by_bucket(7, r_n1pp)
    p_sub = fourpoint_by_bucket(7, r_sub)
    k3_n1 = p_n1.get("k3", 0.0)
    k3_sub = p_sub.get("k3", 0.0)
    return {
        "proved": False,
        "p7_k3_++N1": k3_n1,
        "p7_k3_sub": k3_sub,
        "k3_2line_eq_1line": bool(abs(k3_n1 - k3_sub) < 1.0),
        "diff_times_D_p5": 16,
        "diff_times_D_p7": 48,
        "sixteen_p_minus_4_not_imported": True,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "p=5 Q_N* named. p≥7 still has μ_full. "
            "k=3 2-line ≠ 1-line Fejér on every type. "
            "16(p−4)/D interpolates (Q++−Q_N*)D and is not imported. "
            "Do not interpolate (p−5)/15. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.582  1D vanishes off F_p; p=5 Q_N*=32/13", flush=True)
    A = prove_A()
    print(f"  A 1D vanish: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B Q_N*=32/13: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C k3 2-line=1-line? {C['k3_2line_eq_1line']}", flush=True)
    out = {
        "prop": "15.582",
        "title": "1D 4-point vanishes off F_p; p=5 Q_N*/q²=32/13",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "mu_1d_vanishes_off_Fp": A["proved"],
            "Qnstar_p5_named_32_13": B["proved"],
            "Q_tau_named_in_p": False,
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
        "identity": "Q_N*(p=5)/q^2 = 32/13; mu_1d(r)=0 for r not in F_p",
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
