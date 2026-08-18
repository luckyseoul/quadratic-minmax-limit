#!/usr/bin/env python3
"""
Prop 15.573 — Q_sub is the exclusive 1D / k=3 / full-Ω mix
of the 1-line 4-point of ẑ.

    Q(r) = E[ẑ(ξ)ẑ(−ξ)ẑ(rξ)ẑ(−rξ)]   (r∈F_p^×, type sub)
    = (n_1d μ_1d + n_{k=3} μ_{k=3} + n_{k>3} μ_full) / N
with exclusive support buckets (at p=5, k=3=full, do not
double-count).  Reconstructs live Q(++,sub)/q² at p=5,7.
Fail: count k=3 twice at p=5 (6.15≠48/13).

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.572 flags.  Does **not** import phi_F.
Does **not** overwrite 15.550–15.572.  Does **not** interpolate (p−5)/15.

============================================================================
Setup.  15.290: Q constant on Paley×norm.  15.298: Q is the
4-linear of ŷ on {ξ,−ξ,rξ,−rξ}.  For r∈F_p those four
frequencies lie on one Ω-line.  15.561: live supports k=1,3
(p=5) and k=1,3,m (p=7).  ẑ(−ξ)=conj ẑ(ξ), so the 4-point
equals u(ξ)u(rξ).

============================================================================
Theorem A — PROVED (15.298 + exclusive support; p=5,7 cache).
  Exclusive mix reconstructs live Q(++,sub).  Fail double-count
  of k=3=full at p=5.  ∎

Theorem B — OPEN.  μ_1d, μ_{k=3}, μ_full are not yet a
  Gauss/Jacobi integer in p.  Q_N* / Q_{--N1} unnamed.
  Do not interpolate (p−5)/15.  phi_F not imported.

============================================================================
Backend: p=5,7 yhat 4-point.  Writes
evidence/e1_gmin_m4_prop15573.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15290 import _kernel_field, live_Q_on_T, type_key
from e1_gmin_m4_prop15554 import yhat_omega

EV = ROOT / "evidence" / "e1_gmin_m4_prop15573.json"


def exclusive_buckets(p: int, nsupp: np.ndarray, n_omega: int) -> dict[str, np.ndarray]:
    m = (p + 1) // 2
    oned = nsupp == (p - 1)
    k3 = nsupp == 3 * (p - 1)
    full = nsupp == n_omega
    if m == 3:
        return {"1d": oned, "k3": k3}
    return {"1d": oned, "k3": k3, "full": full}


def sub_r(p: int) -> int:
    Fk = _kernel_field(p)
    Q = live_Q_on_T(p)
    for r in Q:
        key = type_key(Fk, r)
        if len(key) == 3 and key[2] == "sub":
            return int(r)
    raise RuntimeError("no sub r")


def q_sub_mix(p: int, *, double_count: bool = False) -> dict:
    F, Y, Omega, yhat = yhat_omega(p)
    q = F["q"]
    N = yhat.shape[0]
    nOm = len(Omega)
    idx = {om: i for i, om in enumerate(Omega)}
    nsupp = (np.abs(yhat) > 1e-6).sum(axis=1)
    buckets = exclusive_buckets(p, nsupp, nOm)
    if double_count and p == 5:
        buckets["full_dup"] = nsupp == nOm
    r = sub_r(p)
    a = Omega[0]
    b = F["neg"][a]
    c = F["mul"][r][a]
    d = F["neg"][c]
    ia, ib, ic, id_ = idx[a], idx[b], idx[c], idx[d]
    prod = yhat[:, ia] * yhat[:, ib] * yhat[:, ic] * yhat[:, id_]
    q2 = float(q * q)
    mix = 0.0
    parts = {}
    for name, mask in buckets.items():
        n = int(mask.sum())
        if n == 0:
            continue
        mean = float(np.mean(prod[mask].real))
        parts[name] = {"n": n, "mu": mean, "contrib": n * mean / (N * q2)}
        mix += n * mean / N
    live = float(live_Q_on_T(p)[r] / q2)
    return {
        "r": r,
        "mix": mix / q2,
        "live": live,
        "err": abs(mix / q2 - live),
        "parts": parts,
        "N": int(N),
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        hit = q_sub_mix(p)
        bad = q_sub_mix(p, double_count=True) if p == 5 else None
        match = hit["err"] < 1e-9
        fail_dup = True
        if bad is not None:
            fail_dup = abs(bad["mix"] - hit["live"]) > 1.0
        if not (match and fail_dup):
            ok = False
        rows[str(p)] = {
            "mix": hit["mix"],
            "live": hit["live"],
            "err": hit["err"],
            "parts": hit["parts"],
            "dup_mix": None if bad is None else bad["mix"],
            "match": bool(match),
            "dup_fails": bool(fail_dup),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "Exclusive 1D/k=3/full mix names Q_sub at p=5,7. "
            "Fail double-count k=3 at p=5."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Q_sub_mu_named_in_p": False,
        "Q_Nstar_named_in_p": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "Q_sub is the exclusive support mix of a 1-line 4-point. "
            "The three μ pieces are not a p-law. Do not interpolate "
            "(p−5)/15. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.573  Q_sub is the exclusive 1D/k=3/full 4-point mix", flush=True)
    A = prove_A()
    print(f"  A mix: {A['proved']}", flush=True)
    B = prove_open()
    out = {
        "prop": "15.573",
        "title": "Q_sub exclusive 1D/k=3/full mix of the 1-line 4-point",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "proved": {
            "Q_sub_exclusive_mix": A["proved"],
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
        "identity": "Q_sub = (n_1d mu_1d + n_k3 mu_k3 + n_full mu_full)/N",
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
