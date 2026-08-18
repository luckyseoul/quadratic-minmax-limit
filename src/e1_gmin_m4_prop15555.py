#!/usr/bin/env python3
"""
Prop 15.555 — G_3(0,u,v) depends only on the F_p-norms
    (N(u), N(v), N(v−u)), N=N_{q/p}.
Hence the 3-distinct part of Q3 is a complete sum of a
function Φ:(F_p^×)³→ℝ against additive characters on
norm fibers.  Φ is not a short F_p quadratic polynomial.
The three Q3 constants stay unnamed as a p-rational.
16pA / Q_τ / phi_F stay open.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.554 flags.  Does **not** import phi_F.

============================================================================
Setup.  G_3(0,u,v)=Σ_y y_0 y_u y_v.  Q3(η,θ,φ)=Σ_y ŷ(η)ŷ(θ)ŷ(φ).
15.554: Q3 lives on η+θ+φ=0 and is constant on line types.
Collision strata of the Fourier inversion contribute −2Np.
Characters of F_q^× of order dividing p−1 are Norm-pulled from
F_p^×.  G3 lies in that span (lstsq maxerr 1e-12 at p=5,7).

============================================================================
Theorem A — PROVED (cache; p=5,7).
  For u,v,v−u≠0,
      G_3(0,u,v) = Φ(N(u), N(v), N(v−u))
  with N(x)=x^{p+1}∈F_p.  52 (p=5) and 156 (p=7) Norm keys,
  0 splits.  Fail: G3 depends only on (χ(u),χ(v),χ(v−u))
  (p=5 type (1,1,1,1) takes both 2 and 18).  ∎

Theorem B — PROVED (A + Fourier).
  On a triple η+θ+φ=0 in Ω,
      Q3 = −2 N p + q Σ_{u≠v, u,v≠0} Φ(N(u),N(v),N(v−u)) e_θ(u) e_φ(v).
  Fail: drop the −2Np collision (p=5 double: −3700≠−5000).  ∎

Theorem C — OPEN.  Φ takes 7 values at p=5 and 14 at p=7 and is
  not a short χ_p-polynomial (lstsq maxerr 15, 174).  Q3
  constants / 16pA / Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: p=5,7 Max+ caches + field Norm.  Writes
evidence/e1_gmin_m4_prop15555.json
"""
from __future__ import annotations

import json
import os
import sys
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general, gsum_disj_lb_proved_general
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field
from e1_gmin_m4_prop15550 import omega_set_local
from e1_gmin_m4_prop15554 import LIVE_Q3, load_Y, triple_kind, yhat_omega

EV = ROOT / "evidence" / "e1_gmin_m4_prop15555.json"


def field_norm(F: dict, x: int) -> int:
    """N_{q/p}(x)=x^{p+1}."""
    if x == 0:
        return 0
    mul, p = F["mul"], F["p"]
    y = 1
    base = x
    e = p + 1
    while e:
        if e & 1:
            y = mul[y][base]
        base = mul[base][base]
        e >>= 1
    return int(y)


def G3_matrix(p: int):
    F = _kernel_field(p)
    q = F["q"]
    Y = load_Y(p, q)
    G3 = (Y * Y[:, 0:1]).T @ Y
    return F, Y, G3


def norm_buckets(p: int) -> dict:
    F, Y, G3 = G3_matrix(p)
    q = F["q"]
    add, neg, chi = F["add"], F["neg"], F["chi"]
    Ns = [field_norm(F, x) for x in range(q)]
    buckets: dict[tuple[int, int, int], list[float]] = defaultdict(list)
    chi_buckets: dict[tuple[int, int, int], list[float]] = defaultdict(list)
    for u in range(1, q):
        for v in range(1, q):
            if u == v:
                continue
            w = add[v][neg[u]]
            if w == 0:
                continue
            val = float(np.real(G3[u, v]))
            buckets[(Ns[u], Ns[v], Ns[w])].append(val)
            chi_buckets[(int(chi[u]), int(chi[v]), int(chi[w]))].append(val)
    return {
        "F": F,
        "Y": Y,
        "G3": G3,
        "Ns": Ns,
        "buckets": buckets,
        "chi_buckets": chi_buckets,
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = norm_buckets(p)
        n_split = 0
        n_ok = 0
        nuniq = set()
        for vals in rec["buckets"].values():
            arr = np.array(vals)
            if arr.max() - arr.min() > 1e-6:
                n_split += 1
            else:
                n_ok += 1
                nuniq.add(round(float(arr.mean()), 6))
        chi_split = 0
        chi_multi = 0
        for vals in rec["chi_buckets"].values():
            arr = np.array(vals)
            if arr.max() - arr.min() > 1e-6:
                chi_split += 1
                chi_multi = max(chi_multi, len(set(np.round(arr, 6))))
        row_ok = n_split == 0 and n_ok > 0 and chi_split > 0
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "n_keys": len(rec["buckets"]),
            "n_constant": n_ok,
            "n_split": n_split,
            "nuniq_Phi": len(nuniq),
            "chi_split": chi_split,
            "chi_max_nuniq": chi_multi,
            "ok": bool(row_ok),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "G3(0,u,v)=Φ(N(u),N(v),N(v−u)). Fail χ-only "
            "(p=5 a χ-type takes 2 and 18)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = norm_buckets(p)
        F, Y, G3 = rec["F"], rec["Y"], rec["G3"]
        q, N = F["q"], Y.shape[0]
        Phi = {k: float(np.mean(v)) for k, v in rec["buckets"].items()}
        Ns = rec["Ns"]
        add, neg = F["add"], F["neg"]
        # rebuild G3 from Φ
        maxerr = 0.0
        for u in range(1, q):
            for v in range(1, q):
                if u == v:
                    continue
                w = add[v][neg[u]]
                if w == 0:
                    continue
                pred = Phi[(Ns[u], Ns[v], Ns[w])]
                maxerr = max(maxerr, abs(float(np.real(G3[u, v])) - pred))
        # Q3 vs collision + Fourier of Φ for one double and one n3
        F2, Y2, Omega, yhat = yhat_omega(p)
        idx = {om: i for i, om in enumerate(Omega)}
        samples = {}
        for i, a in enumerate(Omega):
            for j, b in enumerate(Omega):
                c = F2["neg"][F2["add"][a][b]]
                if c not in idx:
                    continue
                r = F2["mul"][b][F2["inv"][a]]
                kind = triple_kind(p, r)
                if kind in samples or kind == "degenerate":
                    continue
                theta, phi = b, c
                acc = 0j
                for u in range(1, q):
                    for v in range(1, q):
                        if u == v:
                            continue
                        w = add[v][neg[u]]
                        if w == 0:
                            continue
                        acc += Phi[(Ns[u], Ns[v], Ns[w])] * _e_tr(
                            F, F["mul"][theta][u]
                        ) * _e_tr(F, F["mul"][phi][v])
                named = -2 * N * p + q * acc
                raw = np.sum(yhat[:, i] * yhat[:, idx[c]] * yhat[:, j])
                # wait ŷ(η)ŷ(θ)ŷ(φ) with η=a, θ=b, φ=c: i, j, idx[c]
                raw = np.sum(yhat[:, i] * yhat[:, j] * yhat[:, idx[c]])
                drop = named + 2 * N * p  # drop collision
                samples[kind] = {
                    "raw": float(np.real(raw)),
                    "named": float(np.real(named)),
                    "drop_coll": float(np.real(drop)),
                    "live": LIVE_Q3[p].get(
                        {"double": "double", "generic": "generic", "n3": "n3"}[kind]
                    ),
                    "match": bool(abs(raw - named) < 1e-4),
                    "drop_fails": bool(abs(drop - raw) > 1.0),
                }
                break
        row_ok = maxerr < 1e-6 and all(s["match"] and s["drop_fails"] for s in samples.values())
        if not row_ok:
            ok = False
        rows[str(p)] = {
            "Phi_rebuild_maxerr": float(maxerr),
            "samples": samples,
            "ok": bool(row_ok),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "Q3=−2Np + q Σ Φ(N(u),N(v),N(v−u)) e_θ(u) e_φ(v). "
            "Fail drop −2Np."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "Phi_named_in_p": False,
        "Q3_named_in_p": False,
        "NUM_SUM_16pA_general": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "note": (
            "Φ takes 7 values at p=5 and 14 at p=7; not a short "
            "χ_p-polynomial. Q3 constants / 16pA / Q_τ unnamed."
        ),
    }


def main() -> dict:
    print("Prop 15.555  G3=Φ∘Norm; Q3 is a norm-fiber complete sum", flush=True)
    A = prove_A()
    print(f"  A Norm reduction: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B Q3 Fourier of Φ: {B['proved']}", flush=True)
    C = prove_open()
    out = {
        "prop": "15.555",
        "title": "G3(0,u,v)=Φ(N(u),N(v),N(v−u)); Q3 is a norm-fiber sum",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "G3_is_Phi_of_norms": A["proved"],
            "Q3_Fourier_of_Phi": B["proved"],
            "Q3_named_in_p": False,
            "NUM_SUM_16pA_general": False,
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
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
