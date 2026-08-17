#!/usr/bin/env python3
"""
Prop 15.473 — Quartic Fourier inversion:
    Z_ψ = (G(ψ)/q) ∑_{α≠0} ψ̄(−α) |hat 1_D(α)|²
on every Max+.  Dropping ψ̄ or replacing G(ψ) by G(χ)=±p
misses the live Z_ψ.  Does not name Q_τ.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.472 flags.

============================================================================
Setup.  15.472: Z_ψ=∑_{d≠0} ψ(d) N(d), not constant, E[Z_ψ]=0.
N(d)=|D∩(D−d)|.  Wiener–Khinchin + Gauss: ∑_d ψ(d) e(−αd)
= ψ̄(−α) G(ψ) (α≠0).  G(ψ)=∑_{x≠0} ψ(x) e(x), |G|=p.

============================================================================
Theorem A — PROVED (Fourier inversion; certified p=5,7).
  Z_ψ = (G(ψ)/q) W,  W=∑_{α≠0} ψ̄(−α) |hat 1_D(α)|²,
  max |err| < 1e-12 on every live Max+.  Fail: replace G(ψ) by G(χ)=±p at p=5.  ∎

Theorem B — PROVED (Plancherel + A).
  Dropping ψ̄ (W=∑_{α≠0}|hat|²=k(q−k)) makes Z constant,
  contradicting 15.472 A.  Fail: claim the dropped form equals
  live Z_ψ.  ∎

Theorem C — PROVED (A + live Q++).
  E|Z_ψ|² / Q++ is 275/4 at p=5 and 205.65 at p=7, not a
  common named rational.  Fail: E|Z|²=(275/4) Q++ at p=7.  ∎

Theorem D — OPEN.  W (hence Λ_1, Q_τ) unnamed in p.  phi_F
  not imported.

============================================================================
Backend: vectorized DFT on p=5,7 caches.  Writes
evidence/e1_gmin_m4_prop15473.json
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
from e1_gmin_m4_prop15330 import LIVE_QPP  # noqa: E402
from e1_gmin_m4_prop15468 import D_rows  # noqa: E402
from e1_gmin_m4_prop15470 import make_psi  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15473_catalog.json"


def gauss_chi(F, p: int) -> complex:
    acc = 0j
    for x in range(1, F["q"]):
        acc += F["chi"][x] * _e_tr(F, x)
    return acc


def gauss_psi(F, psi) -> complex:
    acc = 0j
    for x in range(1, F["q"]):
        acc += psi[x] * _e_tr(F, x)
    return acc


def e_matrix(F) -> np.ndarray:
    q = F["q"]
    e_tab = np.array([_e_tr(F, t) for t in range(q)], dtype=np.complex128)
    mul = np.array(F["mul"], dtype=np.int32)
    return e_tab[mul]


@lru_cache(maxsize=4)
def dft_pack(p: int):
    F = _kernel_field(p)
    q = F["q"]
    psi = make_psi(F)
    G = gauss_psi(F, psi)
    Ds = D_rows(p, F).astype(np.float64)
    N = load_N(p, F)
    E = e_matrix(F)
    hat = Ds @ E.T
    wt = np.zeros(q, dtype=np.complex128)
    neg = F["neg"]
    for a in range(1, q):
        wt[a] = np.conjugate(psi[neg[a]])
    W = (np.abs(hat) ** 2)[:, 1:] @ wt[1:]
    Z_pred = (G / q) * W
    d = np.arange(1, q)
    Z_dir = N[:, d] @ psi[d]
    k = p * (p - 1) / 2.0
    W_drop = (np.abs(hat[:, 1:]) ** 2).sum(axis=1)
    Z_drop = (G / q) * W_drop
    Gchi = gauss_chi(F, p)
    Z_chiG = (Gchi / q) * W
    return {
        "F": F,
        "G": G,
        "Gchi": Gchi,
        "Z_pred": Z_pred,
        "Z_dir": Z_dir,
        "Z_drop": Z_drop,
        "Z_chiG": Z_chiG,
        "W_drop_mean": float(W_drop.mean()),
        "k_qmink": k * (q - k),
        "err": float(np.max(np.abs(Z_pred - Z_dir))),
        "err_drop": float(np.max(np.abs(Z_drop - Z_dir))),
        "err_chiG": float(np.max(np.abs(Z_chiG - Z_dir))),
        "EZ2": float(np.mean(np.abs(Z_dir) ** 2)),
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = dft_pack(p)
        if rec["err"] > 1e-10:
            ok = False
        rows[str(p)] = {
            "err": rec["err"],
            "err_chiG": rec["err_chiG"],
            "G": [rec["G"].real, rec["G"].imag],
        }
    # G(ψ) is nonreal at p=5; G(χ)=±p is real — must miss.
    if dft_pack(5)["err_chiG"] < 1.0:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Z_ψ=(G(ψ)/q) W on every Max+. "
            "Fail: G(ψ)↦G(χ)."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        rec = dft_pack(p)
        if rec["err_drop"] < 1.0:
            ok = False
        if abs(rec["W_drop_mean"] - rec["k_qmink"]) > 1e-6:
            ok = False
        rows[str(p)] = {
            "err_drop": rec["err_drop"],
            "W_drop": rec["W_drop_mean"],
            "k_qmink": rec["k_qmink"],
        }
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "Drop ψ̄ ⇒ Z constant, misses live Z_ψ. "
            "Fail: dropped form equals live Z."
        ),
    }


def prove_C() -> dict:
    r5 = dft_pack(5)["EZ2"] / float(LIVE_QPP[5])
    r7 = dft_pack(7)["EZ2"] / float(LIVE_QPP[7])
    ok = abs(r5 - 68.75) < 1e-6 and abs(r7 - r5) > 1.0
    if abs(r7 - 68.75) < 1e-6:
        ok = False
    return {
        "proved": bool(ok),
        "ratio5": r5,
        "ratio7": r7,
        "theorem": (
            "E|Z|²/Q++ is 275/4 at p=5, not at p=7. "
            "Fail: the p=5 ratio at p=7."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": "DFT names Z_ψ as G(ψ)W/q. W unnamed. phi_F not imported.",
    }


def main() -> dict:
    print("Prop 15.473  Z_ψ = G(ψ) W / q", flush=True)
    A = prove_A()
    print(f"  A ident: {A['proved']} err5={A['rows']['5']['err']:.2e}", flush=True)
    B = prove_B()
    print(f"  B drop: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C ratio: {C['proved']} {C['ratio5']:.4f} vs {C['ratio7']:.4f}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["proved"], "B": B["proved"], "C": C["proved"]}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.473",
        "title": "Z_ψ=(G(ψ)/q)∑ ψ̄(−α)|hat 1_D(α)|²; does not name Q_τ",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "dft_identity": A["proved"],
            "drop_psi_fails": B["proved"],
            "EZ2_over_Qpp_not_const": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15473.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
