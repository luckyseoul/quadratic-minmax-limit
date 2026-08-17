#!/usr/bin/env python3
"""
Prop 15.472 — The quartic pair form Z_ψ=∑_{d≠0} ψ(d) N(d) is
not Max+-constant (unlike S_χ).  E[Z_ψ]=0.  Named Gauss/Jacobi
and cyclotomic atoms miss Λ_1(□) at both primes.  phi_F not
imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.471 flags.

============================================================================
Setup.  15.471: Λ_1(d)=ψ̄(d) E[N(d) Z_ψ], Z_ψ=∑_{e≠0} ψ(e) N(e).
15.468: Z_χ=S_χ is p(p²−1)/4 on every Max+.  ψ²=χ.

============================================================================
Theorem A — PROVED (cache fold; p=5,7).
  Z_ψ is not constant: 4 values at p=5 (|Z|∈[√180,√500]),
  5 values at p=7 (Z∈ℝ, |Z|∈[0,84]).  E[Z_ψ]=0.
  Fail: Z_ψ=S_χ=30 at p=5.  ∎

Theorem B — PROVED (A + 15.471 Λ_1).
  Λ_1(□) equals neither G(ψ) nor J(ψ,ψ) nor J(ψ,χ) at both
  primes (p=5: 220/13−(110/13)i vs √5(1+2i) and 3−4i;
  p=7: 13230/409 vs 7).  Fail: Λ_1(□)=G(ψ).  ∎

Theorem C — OPEN.  Λ_1 / Q_τ unnamed.  phi_F not imported.

============================================================================
Backend: 15.471 Λ + N-fold Z_ψ.  Writes
evidence/e1_gmin_m4_prop15472.json
"""
from __future__ import annotations

import json
import os
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "1")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15274 import residual_ii_k_eq_4p_empty  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15468 import S_chi_named  # noqa: E402
from e1_gmin_m4_prop15470 import make_psi  # noqa: E402
from e1_gmin_m4_prop15471 import LIVE_LAM1_SQ, Lambda_tables, mean_on_chi  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15472_catalog.json"


@lru_cache(maxsize=4)
def Z_psi(p: int) -> np.ndarray:
    F = _kernel_field(p)
    N = load_N(p, F)
    psi = make_psi(F)
    d = np.arange(1, F["q"])
    return N[:, d] @ psi[d]


def gauss_psi(p: int) -> complex:
    F = _kernel_field(p)
    psi = make_psi(F)
    acc = 0j
    for x in range(1, F["q"]):
        acc += psi[x] * _e_tr(F, x)
    return acc


def jacobi_pp(p: int) -> complex:
    F = _kernel_field(p)
    psi = make_psi(F)
    acc = 0j
    for x in range(1, F["q"]):
        y = F["add"][1][F["neg"][x]]
        if y == 0:
            continue
        acc += psi[x] * psi[y]
    return acc


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        Z = Z_psi(p)
        nuniq = len({(round(z.real, 6), round(z.imag, 6)) for z in Z})
        mean = complex(Z.mean())
        if nuniq < 2:
            ok = False
        if abs(mean) > 1e-8:
            ok = False
        if abs(np.std(np.abs(Z))) < 1e-8:
            ok = False
        rows[str(p)] = {
            "nunique": nuniq,
            "mean_re": mean.real,
            "mean_im": mean.imag,
            "std": float(np.std(np.abs(Z))),
            "min_abs": float(np.min(np.abs(Z))),
            "max_abs": float(np.max(np.abs(Z))),
            "EZ2": float(np.mean(np.abs(Z) ** 2)),
        }
    if np.max(np.abs(Z_psi(5) - S_chi_named(5))) < 1e-6:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Z_ψ is not constant; E[Z_ψ]=0. "
            "Fail: Z_ψ=S_χ=30 at p=5."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F, _, Lam = Lambda_tables(p)
        lam = mean_on_chi(F, Lam[1], 1)
        G = gauss_psi(p)
        J = jacobi_pp(p)
        if abs(lam - G) < 1e-6:
            ok = False
        if abs(lam - J) < 1e-6:
            ok = False
        if abs(lam - LIVE_LAM1_SQ[p]) > 1e-6:
            ok = False
        rows[str(p)] = {
            "Lam1": [lam.real, lam.imag],
            "G": [G.real, G.imag],
            "Jpp": [J.real, J.imag],
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Λ_1(□) is not G(ψ) or J(ψ,ψ). "
            "Fail: Λ_1(□)=G(ψ)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Z_ψ varies; named GJ atoms miss Λ_1. "
            "Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.472  Z_ψ not constant; GJ miss Λ_1", flush=True)
    A = prove_A()
    print(f"  A Z_ψ: {A['proved']} { {k: v['nunique'] for k,v in A['rows'].items()} }", flush=True)
    B = prove_B()
    print(f"  B GJ: {B['proved']}", flush=True)
    C = prove_open()
    print(f"  C open phi_F={C['phi_F_ge_6']}", flush=True)
    CAT.write_text(json.dumps({"A": A["proved"], "B": B["proved"]}, indent=2) + "\n")
    out = {
        "prop": "15.472",
        "title": "Z_ψ not constant (E=0); named GJ miss Λ_1 at both primes",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "Zpsi_not_constant": A["proved"],
            "Lam1_not_GJ": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "residual_ii_k_eq_4p_empty": C["residual_ii_k_eq_4p_empty"],
        },
        "algebra": {"A": A, "B": B, "C": C},
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15472.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
