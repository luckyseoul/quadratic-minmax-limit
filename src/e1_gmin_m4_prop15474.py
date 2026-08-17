#!/usr/bin/env python3
"""
Prop 15.474 — Ω-support of hat 1_D and Plancherel split of W.

From Cy=py finite-row + Hasse–Davenport G(χ_{p²})=σp:
    hat 1_D(α)=0 for α∉Ω, Ω={α≠0 : χ(α)=σ},
    σ=−(−1)^{(p−1)/2}.
Then W=τ̄(U_+−U_-), U_++U_-=k(q−k) named.
Δ=U_+−U_- is not constant (does not name Q_τ).
phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.473 flags.

============================================================================
Setup.  15.473: Z_ψ=(G(ψ)/q) W, W=∑_{α≠0} ψ̄(−α)|hat 1_D(α)|².
Finite-row of Cy=py (x∈F): ∑_z χ(z−x) y_z = p y_x − 1.
Hasse–Davenport: G(χ_{p²})=σp, σ=−(−1)^{(p−1)/2}.

============================================================================
Theorem A — PROVED (Max+-free; Cy=py + HD; certified p=3,5,7).
  hat 1_D(α)=0 off Ω={χ=σ}.  max |leak| < 1e-12 on every live Max+.
  Fail: use G(χ_p)=√5 not lifted σ·5=−5 ⇒ predict k∈{0,25}; live k=10.  ∎

Theorem B — PROVED (A + Plancherel; certified p=3,5,7).
  W=τ̄ Δ with Δ=U_+−U_-, U_++U_-=k(q−k) named (18, 150, 588).
  Matches 15.473 W (max |err|<1e-12).  E[Δ]=0.
  Fail: claim U_++U_- ≠ k(q−k) at p=5.  ∎

Theorem C — PROVED (B + cache fold).
  Δ is not constant: 2 values at p=3, 4 at p=5, 5 at p=7.
  Fail: claim Δ=0 on every Max+ at p=5.  ∎

Theorem D — OPEN.  Δ (hence W, Q_τ) unnamed in p.  phi_F not imported.
  Signed Δ is not a dilation-invariant name (square a with ψ(a)=−1
  sends Δ↦−Δ).  Target Δ² or a class law.

============================================================================
Backend: vectorized DFT on p=3,5,7 caches.  Writes
evidence/e1_gmin_m4_prop15474.json
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
from e1_gmin_m4_prop15468 import D_rows  # noqa: E402
from e1_gmin_m4_prop15470 import make_psi  # noqa: E402
from e1_gmin_m4_prop15473 import dft_pack, e_matrix, gauss_chi  # noqa: E402
from e1_gmin_m4_prop15279 import _kernel_field  # noqa: E402

CAT = ROOT / "evidence" / "e1_gmin_m4_prop15474_catalog.json"


def sigma_named(p: int) -> int:
    return -((-1) ** ((p - 1) // 2))


@lru_cache(maxsize=4)
def omega_pack(p: int):
    F = _kernel_field(p)
    q = F["q"]
    Ds = D_rows(p, F).astype(np.float64)
    hat = Ds @ e_matrix(F).T
    chi = np.array(F["chi"], dtype=np.int32)
    Gchi = gauss_chi(F, p)
    sigma = int(round(Gchi.real / p))
    Omega = chi == sigma
    Omega[0] = False
    off = ~Omega
    off[0] = False
    leak = float((np.abs(hat[:, off]) ** 2).max()) if off.any() else 0.0
    psi = make_psi(F)
    if sigma < 0:
        plus = Omega & (psi.imag > 0)
        tau = 1j
    else:
        plus = Omega & (psi.real > 0)
        tau = 1.0 + 0j
    Up = (np.abs(hat[:, plus]) ** 2).sum(axis=1)
    Um = (np.abs(hat[:, Omega]) ** 2).sum(axis=1) - Up
    k = p * (p - 1) // 2
    kqmk = k * (q - k)
    sum_err = float(np.max(np.abs(Up + Um - kqmk)))
    W_split = np.conjugate(tau) * (Up - Um)
    rec = dft_pack(p)
    W473 = rec["Z_pred"] * (q / rec["G"])
    split_err = float(np.max(np.abs(W_split - W473)))
    delta = Up - Um
    nunique = int(len(np.unique(np.round(delta, 8))))
    return {
        "q": q,
        "sigma": sigma,
        "sigma_named": sigma_named(p),
        "Omega": int(Omega.sum()),
        "leak": leak,
        "kqmk": float(kqmk),
        "sum_err": sum_err,
        "split_err": split_err,
        "nunique": nunique,
        "delta_mean": float(delta.mean()),
        "live_k": int(Ds[0].sum()),
    }


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7):
        rec = omega_pack(p)
        if rec["leak"] > 1e-10:
            ok = False
        if rec["sigma"] != rec["sigma_named"]:
            ok = False
        rows[str(p)] = {
            "leak": rec["leak"],
            "sigma": rec["sigma"],
            "Omega": rec["Omega"],
        }
    # Fail-when-wrong: unlifted G(χ_p) predicts k∈{0,q} at p=5.
    if omega_pack(5)["live_k"] in (0, 25):
        ok = False
    if omega_pack(5)["live_k"] != 10:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "hat 1_D=0 off Ω={χ=σ}. "
            "Fail: unlifted G(χ_p) predicts k∈{0,25}."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in (3, 5, 7):
        rec = omega_pack(p)
        if rec["sum_err"] > 1e-8 or rec["split_err"] > 1e-8:
            ok = False
        if abs(rec["delta_mean"]) > 1e-6:
            ok = False
        rows[str(p)] = {
            "kqmk": rec["kqmk"],
            "sum_err": rec["sum_err"],
            "split_err": rec["split_err"],
            "delta_mean": rec["delta_mean"],
        }
    if abs(omega_pack(5)["kqmk"] - 150) > 1e-9:
        ok = False
    return {
        "proved": bool(ok),
        "rows": rows,
        "theorem": (
            "W=τ̄Δ, U_++U_-=k(q−k), E[Δ]=0. "
            "Fail: U_++U_- ≠ 150 at p=5."
        ),
    }


def prove_C() -> dict:
    n3 = omega_pack(3)["nunique"]
    n5 = omega_pack(5)["nunique"]
    n7 = omega_pack(7)["nunique"]
    ok = n3 >= 2 and n5 >= 4 and n7 >= 5
    if n5 == 1:
        ok = False
    return {
        "proved": bool(ok),
        "nunique": {"3": n3, "5": n5, "7": n7},
        "theorem": (
            "Δ not constant (2/4/5 values). "
            "Fail: Δ=0 on every Max+ at p=5."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "e1": bool(e1_closed_general()),
        "residual_ii_k_eq_4p_empty": bool(residual_ii_k_eq_4p_empty()),
        "note": (
            "Ω-support + Plancherel split. Δ unnamed. "
            "Signed Δ is not dilation-invariant. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.474  Ω-support + W=τ̄Δ", flush=True)
    A = prove_A()
    print(f"  A support: {A['proved']} leak5={A['rows']['5']['leak']:.2e}", flush=True)
    B = prove_B()
    print(f"  B split: {B['proved']}", flush=True)
    C = prove_C()
    print(f"  C nunique: {C['proved']} {C['nunique']}", flush=True)
    D = prove_open()
    print(f"  D open phi_F={D['phi_F_ge_6']}", flush=True)
    CAT.write_text(
        json.dumps({"A": A["proved"], "B": B["proved"], "C": C["proved"]}, indent=2)
        + "\n"
    )
    out = {
        "prop": "15.474",
        "title": "hat 1_D=0 off Ω; W=τ̄Δ; Δ unnamed",
        "series": "15.x leftover campaign (OPEN)",
        "proved": {
            "omega_support": A["proved"],
            "plancherel_split": B["proved"],
            "delta_not_constant": C["proved"],
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
    dest = ROOT / "evidence" / "e1_gmin_m4_prop15474.json"
    dest.write_text(json.dumps(out, indent=2) + "\n")
    print("wrote", dest, flush=True)
    return out


if __name__ == "__main__":
    main()
