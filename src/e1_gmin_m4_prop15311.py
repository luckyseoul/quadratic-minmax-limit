#!/usr/bin/env python3
"""
Prop 15.311 — Q on the AP⊙NC orbit is the Ω-average of
|ẑ_hs−2Γ|², Γ the signed norm-circle Gauss.  At p=5 this is
the unique NL type and names full-ensemble Q_{++}.  Extra NL
orbits at p≥7 stay unnamed.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.310 flags.  Floor stays OPEN.

============================================================================
Setup.  y_* = hs ⊙ z_C from 15.138 (lex-first Hoffman norm circle C
of C_hs).  Γ(ξ)=∑_{u∈C} z_hs(u) e_ξ(u).  u_*=|ẑ_hs−2Γ|² on Ω.

============================================================================
Theorem A — PROVED (Fourier; Max+-free).
  ẑ_{y_*}(ξ) = ẑ_hs(ξ) − 2 Γ(ξ) on Ω.  Fail-when-wrong: drop the 2
  (use ẑ_hs−Γ).  ∎

Theorem B — PROVED (certified p=5,7).
  Q^{ystar}(r) := |Ω|⁻¹ ∑_{ξ∈Ω} u_*(ξ) u_*(rξ)
  equals the live G-orbit average of Q on AP⊙NC.  Class average on
  ++ is 64/15 at p=5 and 8/3 at p=7.  Fail-when-wrong: replace u_*
  by the AP-Fejer Q_AP (15.310); or use a single ξ0 at p=7.  ∎

Theorem C — PROVED at p=5 only (named weights).
  H_+ = AP-hs ∪ AP⊙NC, sizes |G|/(4p) and |G|/(p+1).
      Q_{++} = (n_AP Q_{++}^{AP} + n_* Q_{++}^{ystar}) / (n_AP+n_*)
  with Q_{++}^{AP} the 15.310 Fejer class average.  Recovers 48/13.
  Fail-when-wrong: drop n_* or put Q_AP in the * slot.  ∎

Theorem D — OPEN.  p≥7 has extra NL orbits (QR0⊙NC, gen2, …).
  Full-ensemble Q_τ unnamed for all p≥5.  phi_F not imported.  ∎

============================================================================
Backend: CPU construct_ystar + Ω Fourier ( |Ω|=(q−1)/2 ).
Writes evidence/e1_gmin_m4_prop15311.json
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15138 import construct_ystar  # noqa: E402
from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field, rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15290 import in_subfield, omega_set, zhat_omega  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15308 import G_order  # noqa: E402
from e1_gmin_m4_prop15310 import Q_AP_Fejer  # noqa: E402


def circle_points(ystar: dict) -> list[int]:
    return [i - 1 for i in ystar["S"] if i >= 1]


def recover_zhs_zc(ystar: dict, q: int) -> tuple[np.ndarray, np.ndarray]:
    z_y = np.asarray(ystar["y"][1 : 1 + q], dtype=np.float64)
    zc = np.ones(q, dtype=np.float64)
    for u in circle_points(ystar):
        zc[u] = -1.0
    z_hs = z_y / zc
    return z_hs, zc


def Gamma_on_Omega(z_hs: np.ndarray, circle: list[int], F, Omega) -> np.ndarray:
    Gam = np.zeros(len(Omega), dtype=np.complex128)
    for i, xi in enumerate(Omega):
        acc = 0j
        for u in circle:
            acc += z_hs[u] * _e_tr(F, F["mul"][xi][u])
        Gam[i] = acc
    return Gam


def prove_fourier() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        ys = construct_ystar(p)
        F = _kernel_field(p)
        q = F["q"]
        Omega = omega_set(F)
        z_hs, zc = recover_zhs_zc(ys, q)
        z_y = np.asarray(ys["y"][1 : 1 + q], dtype=np.float64)
        zh_hs = zhat_omega(z_hs, F, Omega)
        zh_y = zhat_omega(z_y, F, Omega)
        Gam = Gamma_on_Omega(z_hs, circle_points(ys), F, Omega)
        err2 = float(np.max(np.abs(zh_y - (zh_hs - 2 * Gam))))
        err1 = float(np.max(np.abs(zh_y - (zh_hs - Gam))))
        if err2 > 1e-9:
            ok = False
        if err1 < 1e-6:
            ok = False  # fail-when-wrong: ẑ_hs−Γ would also match
        rows[str(p)] = {"err_minus2": err2, "err_minus1": err1, "n_Omega": len(Omega)}
    return {
        "proved": ok,
        "by_p": rows,
        "theorem": "ẑ_{y_*}=ẑ_hs−2Γ on Ω.  Dropping the 2 fails.",
    }


def Q_ystar_Omega(p: int) -> dict[str, float]:
    """Class-average of Ω-mean u_*(ξ)u_*(rξ)/q² on coarse types."""
    ys = construct_ystar(p)
    F = _kernel_field(p)
    q = F["q"]
    q2 = float(q * q)
    Omega = omega_set(F)
    z_hs, _ = recover_zhs_zc(ys, q)
    Gam = Gamma_on_Omega(z_hs, circle_points(ys), F, Omega)
    zh_hs = zhat_omega(z_hs, F, Omega)
    u = np.abs(zh_hs - 2 * Gam) ** 2
    om = {xi: i for i, xi in enumerate(Omega)}
    acc: dict[str, list[float]] = defaultdict(list)
    for r in F["squares"]:
        sm = 0.0
        for xi in Omega:
            sm += u[om[xi]] * u[om[F["mul"][r][xi]]]
        acc[merge_cls(F, r)].append(sm / len(Omega) / q2)
    return {c: float(np.mean(vs)) for c, vs in acc.items()}


def Q_ystar_xi0_only(p: int) -> dict[str, float]:
    """Fail-when-wrong: one ξ0 instead of Ω-average."""
    ys = construct_ystar(p)
    F = _kernel_field(p)
    q = F["q"]
    q2 = float(q * q)
    Omega = omega_set(F)
    z_hs, _ = recover_zhs_zc(ys, q)
    z_y = np.asarray(ys["y"][1 : 1 + q], dtype=np.float64)
    u = np.abs(zhat_omega(z_y, F, Omega)) ** 2
    om = {xi: i for i, xi in enumerate(Omega)}
    xi0 = Omega[0]
    acc: dict[str, list[float]] = defaultdict(list)
    for r in F["squares"]:
        eta = F["mul"][r][xi0]
        acc[merge_cls(F, r)].append(float(u[om[xi0]] * u[om[eta]] / q2))
    return {c: float(np.mean(vs)) for c, vs in acc.items()}


def Q_AP_class(p: int) -> dict[str, float]:
    F = _kernel_field(p)
    acc: dict[str, list[float]] = defaultdict(list)
    for r in F["squares"]:
        pred = Q_AP_Fejer(p, r % p) if in_subfield(F, r) else 0.0
        acc[merge_cls(F, r)].append(pred)
    return {c: float(np.mean(vs)) for c, vs in acc.items()}


def prove_ystar_Q() -> dict:
    ok = True
    rows = {}
    expect = {5: Fraction(64, 15), 7: Fraction(8, 3)}
    for p in (5, 7):
        Qs = Q_ystar_Omega(p)
        Q0 = Q_ystar_xi0_only(p)
        QA = Q_AP_class(p)
        qpp = Fraction(Qs["++"]).limit_denominator(300)
        if qpp != expect[p]:
            ok = False
        if abs(Q0["++"] - float(expect[p])) < 1e-6 and p == 7:
            ok = False  # ξ0-only must miss at p=7
        if abs(QA["++"] - float(expect[p])) < 1e-6:
            ok = False  # Fejer is not the ystar value
        rows[str(p)] = {
            "Qpp_ystar": str(qpp),
            "Qpp_xi0": Q0["++"],
            "Qpp_Fejer": QA["++"],
            "Q_Nstar": Qs.get("N*"),
        }
    return {
        "proved": ok,
        "by_p": rows,
        "theorem": (
            "Ω-average of |ẑ_hs−2Γ|² gives Q_{++}^{ystar}=64/15, 8/3. "
            "Fejer and ξ0-only fail at p=7."
        ),
    }


def prove_p5_mixture() -> dict:
    p = 5
    n_AP = G_order(p) // (4 * p)
    n_star = G_order(p) // (p + 1)
    QA = Q_AP_class(p)["++"]
    Qs = Q_ystar_Omega(p)["++"]
    mix = (n_AP * QA + n_star * Qs) / (n_AP + n_star)
    mix_f = Fraction(mix).limit_denominator(50)
    drop = QA  # drop n_*
    wrong_slot = (n_AP * QA + n_star * QA) / (n_AP + n_star)
    ok = mix_f == Fraction(48, 13)
    if abs(drop - 48 / 13) < 1e-9:
        ok = False
    if abs(wrong_slot - 48 / 13) < 1e-9:
        ok = False
    return {
        "proved": ok,
        "n_AP": n_AP,
        "n_star": n_star,
        "Q_AP": QA,
        "Q_ystar": Qs,
        "mix": str(mix_f),
        "drop_star": drop,
        "fejer_in_star_slot": wrong_slot,
        "theorem": (
            "p=5 Q_{++}=(n_AP Q_AP + n_* Q_ystar)/(n_AP+n_*)=48/13. "
            "Dropping n_* or using Fejer in the * slot fails."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "p5_Qpp_named": True,
        "note": (
            "AP⊙NC Q_τ named by ẑ_hs−2Γ. p=5 full Q_{++} named. "
            "p≥7 extra NL orbits unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.311  Q_ystar = Ω-avg |ẑ_hs−2Γ|²; p=5 mixture 48/13", flush=True)
    A = prove_fourier()
    print(f"  A ẑ_y=ẑ_hs−2Γ: {A['proved']}", flush=True)
    B = prove_ystar_Q()
    print(f"  B Q_ystar ++: {B['proved']} { {k: v['Qpp_ystar'] for k,v in B['by_p'].items()} }", flush=True)
    C = prove_p5_mixture()
    print(f"  C p=5 mix: {C['proved']} {C['mix']}", flush=True)
    D = prove_open()
    print(f"  D floor: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.311",
        "title": "AP⊙NC Q is Ω-avg |ẑ_hs−2Γ|²; p=5 Q++ named; p≥7 NL extra",
        "proved": {
            "fourier_minus_2Gamma": A["proved"],
            "Q_ystar_Omega": B["proved"],
            "p5_mixture_48_13": C["proved"],
            "Q_tau_named_in_p": False,
            "lambda_ge_6": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
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
        ],
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15311.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
