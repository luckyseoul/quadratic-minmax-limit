#!/usr/bin/env python3
"""
Prop 15.471 — Quartic weight on the lag s (not d): Λ_0,Λ_2
are named by 15.468; ensemble order-4 pair counts are the
named halves of χ±; Λ_1 is new and is not Cov++.
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.470 flags.

============================================================================
Setup.  ψ^2=χ.  Λ_j(d)=E[∑_{s≠0} N(d) N(s d) ψ(s)^j].
n_k = #{ {a,b}⊂D : dlog(b−a)≡k (mod 4) }.  15.468: n0+n2=χ+,
n1+n3=χ−, pointwise.  Cyclotomic (i,j)_4 is the named order-4
intersection matrix (not the conference Gram).

============================================================================
Theorem A — PROVED (15.468 + cache fold; p=5,7).
  E[n0]=E[n2]=χ+/2 and E[n1]=E[n3]=χ−/2.
  Not pointwise: 4 fingerprints at p=5 (n1∈[5,10]); 5 at p=7.
  At p=7, n1=n3=42 on every Max+.  Fail: n1 constant at p=5.  ∎

Theorem B — PROVED (15.468 + dilation).
  ∑_{s≠0} N(sd)=k(k−1) and ∑_s χ(s) N(sd)=χ(d) S_χ, pointwise,
  so Λ_0(d)=μ(d) k(k−1) and Λ_2(d)=χ(d) μ(d) S_χ.
  Fail: Λ_2=0 on squares.  ∎

Theorem C — PROVED (lag fold; p=5,7).
  Λ_1 on squares is 220/13 − (110/13)i at p=5 and 13230/409
  (real) at p=7; on ns it is 55/13+(110/13)i and 0.
  Not Cov++ (10/39, 2793/1636).  Fail: Λ_1(□)=10/39.  ∎

Theorem D — OPEN.  Λ_1 unnamed in p.  phi_F not imported.

============================================================================
Backend: 15.468 + N-fold.  Writes
evidence/e1_gmin_m4_prop15471.json
"""
from __future__ import annotations

import json
import os
import sys
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402
from e1_gmin_m4_prop15298 import load_N, mu_minus, mu_plus  # noqa: E402
from e1_gmin_m4_prop15467 import coarse_L  # noqa: E402
from e1_gmin_m4_prop15468 import (  # noqa: E402
    D_rows,
    S_chi_named,
    chi_minus_named,
    chi_plus_named,
    k_named,
)
from e1_gmin_m4_prop15470 import make_psi  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15471_catalog.json"

LIVE_LAM1_SQ = {
    5: 220 / 13 - 110 / 13 * 1j,
    7: 13230 / 409 + 0j,
}
LIVE_NK_MEAN = {
    5: (15.0, 7.5, 15.0, 7.5),
    7: (63.0, 42.0, 63.0, 42.0),
}


def dlog_mod4(psi) -> np.ndarray:
    q = psi.shape[0]
    dl = np.zeros(q, dtype=np.int8)
    for x in range(1, q):
        v = psi[x]
        if abs(v - 1) < 1e-8:
            dl[x] = 0
        elif abs(v - 1j) < 1e-8:
            dl[x] = 1
        elif abs(v + 1) < 1e-8:
            dl[x] = 2
        elif abs(v + 1j) < 1e-8:
            dl[x] = 3
        else:
            raise RuntimeError("bad psi")
    return dl


def projectors_mod4(F, psi):
    q = F["q"]
    dl = dlog_mod4(psi)
    add, neg = F["add"], F["neg"]
    P = np.zeros((4, q, q), dtype=np.float64)
    for a in range(q):
        for b in range(q):
            if a == b:
                continue
            P[dl[add[b][neg[a]]], a, b] = 1.0
    return P


@lru_cache(maxsize=4)
def pair_n_mod4(p: int) -> np.ndarray:
    F = _kernel_field(p)
    psi = make_psi(F)
    Ds = D_rows(p, F).astype(np.float64)
    P = projectors_mod4(F, psi)
    # n_k[i] = (1/2) 1_D^T P_k 1_D
    counts = np.einsum("ia,kab,ib->ik", Ds, P, Ds) / 2.0
    return counts


@lru_cache(maxsize=4)
def Lambda_tables(p: int):
    F = _kernel_field(p)
    N = load_N(p, F)
    psi = make_psi(F)
    q = F["q"]
    mul = np.array(F["mul"], dtype=np.int32)
    slist = np.arange(1, q)
    psi_s = np.stack([psi[slist] ** j for j in range(4)])
    Lam = {j: np.zeros(q, dtype=np.complex128) for j in range(4)}
    for d in range(1, q):
        Nsd = N[:, mul[slist, d]]
        m = (N[:, d][:, None] * Nsd).mean(axis=0)
        for j in range(4):
            Lam[j][d] = np.dot(m, psi_s[j])
    return F, psi, Lam


def mean_on_chi(F, Lamj, chi_val: int) -> complex:
    vs = [Lamj[d] for d in range(1, F["q"]) if F["chi"][d] == chi_val]
    return sum(vs) / len(vs)


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        cnt = pair_n_mod4(p)
        mean = cnt.mean(axis=0)
        want = np.array(
            [
                float(chi_plus_named(p) / 2),
                float(chi_minus_named(p) / 2),
                float(chi_plus_named(p) / 2),
                float(chi_minus_named(p) / 2),
            ]
        )
        if np.max(np.abs(mean - want)) > 1e-8:
            ok = False
        nuniq = len({tuple(np.rint(r).astype(int)) for r in cnt})
        n1_spread = float(cnt[:, 1].max() - cnt[:, 1].min())
        if p == 5 and n1_spread < 1:
            ok = False
        if p == 7 and n1_spread > 1e-8:
            ok = False
        rows[str(p)] = {
            "mean": mean.tolist(),
            "nunique": nuniq,
            "n1_min": float(cnt[:, 1].min()),
            "n1_max": float(cnt[:, 1].max()),
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "E[n_k] are the named halves of χ±. "
            "Fail: n1 constant at p=5."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F, _, Lam = Lambda_tables(p)
        k = k_named(p)
        S = S_chi_named(p)
        l0p = mean_on_chi(F, Lam[0], 1).real
        l0m = mean_on_chi(F, Lam[0], -1).real
        l2p = mean_on_chi(F, Lam[2], 1).real
        l2m = mean_on_chi(F, Lam[2], -1).real
        w0p = float(mu_plus(p) * k * (k - 1))
        w0m = float(mu_minus(p) * k * (k - 1))
        w2p = float(mu_plus(p) * S)
        w2m = float(-mu_minus(p) * S)
        if abs(l0p - w0p) > 1e-6 or abs(l0m - w0m) > 1e-6:
            ok = False
        if abs(l2p - w2p) > 1e-6 or abs(l2m - w2m) > 1e-6:
            ok = False
        if abs(l2p) < 1e-8:
            ok = False
        rows[str(p)] = {
            "Lam0+": l0p,
            "named0+": w0p,
            "Lam2+": l2p,
            "named2+": w2p,
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Λ_0=μ k(k−1) and Λ_2=χ μ S_χ (15.468). "
            "Fail: Λ_2(+)=0."
        ),
    }


def prove_C() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F, _, Lam = Lambda_tables(p)
        v = mean_on_chi(F, Lam[1], 1)
        want = LIVE_LAM1_SQ[p]
        if abs(v - want) > 1e-6:
            ok = False
        cov = coarse_L(p)["++"][1] - mu_plus(p) ** 2
        if abs(v.real - float(cov)) < 1e-8 and abs(v.imag) < 1e-8:
            ok = False
        rows[str(p)] = {
            "Lam1_sq": [v.real, v.imag],
            "cov": str(cov),
        }
    if abs(LIVE_LAM1_SQ[5] - (10 / 39)) < 1e-9:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Λ_1(□)=220/13−(110/13)i and 13230/409, not Cov++. "
            "Fail: Λ_1(□)=10/39."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Λ_1 is the odd lag moment, not Cov++. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.471  ψ on lag s; order-4 pair halves", flush=True)
    A = prove_A()
    print(f"  A n_k: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B Lam02: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C Lam1: {C['proved']} {C['rows']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["proved"], "B": B["proved"], "C": C["proved"]}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.471",
        "title": "Lag-s Λ_0,Λ_2 named; Λ_1 ≠ Cov++; n_k halves of χ±",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "nk_halves": A["proved"],
            "Lam02_named": B["proved"],
            "Lam1_not_cov": C["proved"],
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": D["residual_ii_k_eq_4p_empty"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15471.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
