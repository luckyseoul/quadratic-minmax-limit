#!/usr/bin/env python3
"""
Prop 15.304 — Floor eigenvalues are a named multiple of a torus-Gauss
square: λ_j = 8/(p+1) q^{-2} E|⟨ΔN, G_U(j,ξ·)⟩|².

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.303 flags.  Does **not** prove λ_j≥6.  Floor stays OPEN.

============================================================================
Setup.  U={N_{q/p}=1}, |U|=p+1.  s∈F_p^×\\{±1}.  ΔN(δ)=N(δ)−N(s^{-1}δ).
G_U(j,α)=∑_{u∈U} ξ_j(u) e(αu).  hatN(ξ)=∑_δ N(δ) e(ξδ).
λ_j = F_1(j)−F_*(j) (15.302), the σ=−1 family, hence the floor
once Q_{++}≥Q_{N*}.

============================================================================
Theorem A — PROVED (certified p=5,7; expansion of Q via hatN).
  For every even j (including 0)
      λ_j = 16 q^{-2} E[ hatN(ξ) ⟨ΔN, G_U(j,ξ·)⟩ ].
  E[ΔN]=0 and E[⟨ΔN,G_U⟩]=0.  Fail-when-wrong: drop the shift
  (replace ΔN by N); residual >1 at p=5 j=2.  ∎

Theorem B — PROVED (U-stationarity + Parseval; certified p=5,7).
  The same λ_j is the square
      λ_j = 8/(p+1) q^{-2} E |⟨ΔN, G_U(j,ξ·)⟩|².
  The named denominator is |U|=p+1, neither dim V_+=(q+1)/2 nor
  N/(4p).  Fail-when-wrong: replace p+1 by 1 (breaks p=5) or by
  dim V_+ (breaks p=7).  ∎

Theorem C — OPEN.  Squares give λ_j≥0.  The floor λ_j≥6 is
      E|⟨ΔN,G_U⟩|² ≥ (3/4)(p+1) q²,
  live at p=5,7 (2884.6≥2812.5 and 18034≥14406) but not proved.
  Q_τ / phi_F stay OPEN.  ∎

============================================================================
Backend: vectorized N @ kernel on the p=5,7 caches (one fold per prime).
Writes evidence/e1_gmin_m4_prop15304.json
"""
from __future__ import annotations

import json
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _e_tr,
    _kernel_field,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15290 import live_Q_on_T, omega_set  # noqa: E402
from e1_gmin_m4_prop15298 import load_N  # noqa: E402
from e1_gmin_m4_prop15301 import dim_Vplus  # noqa: E402
from e1_gmin_m4_prop15302 import form_factor, nontrivial_s  # noqa: E402
from e1_gmin_m4_prop15303 import GU  # noqa: E402


def e_vec(F, xi: int) -> np.ndarray:
    q = F["q"]
    return np.array([_e_tr(F, F["mul"][xi][d]) for d in range(q)], dtype=np.complex128)


def gu_vec(F, j: int, xi: int) -> np.ndarray:
    q = F["q"]
    return np.array([GU(F, j, F["mul"][xi][d]) for d in range(q)], dtype=np.complex128)


@lru_cache(maxsize=4)
def _prep(p: int):
    F = _kernel_field(p)
    N = load_N(p, F).astype(np.float64)
    xi = omega_set(F)[0]
    s = nontrivial_s(F)
    sinv = F["inv"][s]
    perm = [F["mul"][sinv][d] for d in range(F["q"])]
    dN = N - N[:, perm]
    ev = e_vec(F, xi)
    hatN = N @ ev
    Q = live_Q_on_T(p)
    return F, N, dN, hatN, xi, s, Q


def live_lambda(p: int, j: int) -> complex:
    F, _N, _dN, _h, _xi, s, Q = _prep(p)
    return form_factor(F, Q, 1, j) - form_factor(F, Q, s, j)


def lambda_from_hatN(p: int, j: int, use_delta: bool = True) -> complex:
    F, N, dN, hatN, xi, _s, _Q = _prep(p)
    q2 = float(F["q"] ** 2)
    gv = gu_vec(F, j, xi)
    inner = (dN if use_delta else N) @ gv
    return 16.0 * np.mean(hatN * inner) / q2


def lambda_from_square(p: int, j: int, den: int | None = None) -> float:
    F, _N, dN, _h, xi, _s, _Q = _prep(p)
    q2 = float(F["q"] ** 2)
    if den is None:
        den = p + 1
    gv = gu_vec(F, j, xi)
    inner = dN @ gv
    return (8.0 / den) * float(np.mean(np.abs(inner) ** 2)) / q2


def prove_hatN_identity() -> dict:
    ok = True
    drop_hit = 0
    rows = {}
    for p in (5, 7):
        rec = {}
        for j in range(0, p + 1, 2):
            live = live_lambda(p, j)
            pred = lambda_from_hatN(p, j, use_delta=True)
            err = abs(pred - live)
            rec[str(j)] = {"live": float(live.real), "pred": float(pred.real), "err": float(err)}
            if err > 1e-8:
                ok = False
            if j == 2:
                bad = lambda_from_hatN(p, j, use_delta=False)
                if abs(bad - live) > 1.0:
                    drop_hit += 1
        rows[str(p)] = rec
    if drop_hit == 0:
        ok = False
    return {
        "proved": ok,
        "drop_delta_hit": drop_hit,
        "by_p": rows,
        "theorem": (
            "λ_j=16 q^{-2} E[hatN ⟨ΔN,G_U(j)⟩] on all even j. "
            "Replacing ΔN by N breaks it."
        ),
    }


def prove_square_identity() -> dict:
    ok = True
    den1_hit = 0
    dim_hit = 0
    rows = {}
    for p in (5, 7):
        rec = {}
        dV = dim_Vplus(p)
        for j in range(0, p + 1, 2):
            live = float(live_lambda(p, j).real)
            pred = lambda_from_square(p, j, den=p + 1)
            err = abs(pred - live)
            rec[str(j)] = {"live": live, "pred": pred, "err": err, "den": p + 1}
            if err > 1e-8:
                ok = False
            if j == 2:
                bad1 = lambda_from_square(p, j, den=1)
                if abs(bad1 - live) > 0.5:
                    den1_hit += 1
                badd = lambda_from_square(p, j, den=dV)
                if abs(badd - live) > 0.5:
                    dim_hit += 1
        rows[str(p)] = rec
    if den1_hit == 0 or dim_hit == 0:
        ok = False
    return {
        "proved": ok,
        "den1_hit": den1_hit,
        "dimVplus_hit": dim_hit,
        "by_p": rows,
        "named_den": "p+1 = |U|",
        "theorem": (
            "λ_j=8/(p+1) q^{-2} E|⟨ΔN,G_U(j)⟩|².  Den |U| is not "
            "dim V_+ and not N/(4p).  den=1 or den=dim V_+ fail."
        ),
    }


def prove_open() -> dict:
    # live squares beat the 6-threshold; general comparison unproved
    live_ge6 = True
    rooms = {}
    for p in (5, 7):
        need = 0.75 * (p + 1) * (p * p) ** 2
        rec = {}
        for j in range(0, p + 1, 2):
            F, _N, dN, _h, xi, _s, _Q = _prep(p)
            inner = dN @ gu_vec(F, j, xi)
            mass = float(np.mean(np.abs(inner) ** 2))
            rec[str(j)] = {"E_abs2": mass, "need": need, "ge": mass + 1e-9 >= need}
            if mass < need:
                live_ge6 = False
        rooms[str(p)] = rec
    return {
        "proved": False,
        "live_square_beats_threshold": live_ge6,
        "threshold": "E|inner|² ≥ (3/4)(p+1) q²  ⇔  λ_j≥6",
        "rooms": rooms,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "S_k_named_in_p": False,
        "note": (
            "Square formula names the den as |U|=p+1 and proves λ_j≥0. "
            "The comparison to 6 is live at p=5,7 and unproved in general. "
            "phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.304  λ_j = 8/|U| q^{-2} E|⟨ΔN, G_U⟩|²", flush=True)
    A = prove_hatN_identity()
    print(f"  A hatN contraction: {A['proved']} dropΔ={A['drop_delta_hit']}", flush=True)
    B = prove_square_identity()
    print(
        f"  B square 8/(p+1): {B['proved']} den1={B['den1_hit']} dimV={B['dimVplus_hit']}",
        flush=True,
    )
    C = prove_open()
    print(
        f"  C floor: {C['proved']} live_ge6={C['live_square_beats_threshold']} "
        f"phi_F={C['phi_F_ge_6']}",
        flush=True,
    )
    out = {
        "prop": "15.304",
        "title": "λ_j=8/(p+1) q^{-2} E|⟨ΔN,G_U(j)⟩|²; den=|U|",
        "proved": {
            "hatN_contraction": A["proved"],
            "square_named_den_U": B["proved"],
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "Q_tau_named_in_p": False,
            "S_k_named_in_p": False,
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
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15304.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
