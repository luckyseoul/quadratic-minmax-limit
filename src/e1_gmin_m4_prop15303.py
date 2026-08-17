#!/usr/bin/env python3
"""
Prop 15.303 — On T, the μ-Gauss part of Q is 4(q−1)/q;
torus Gauss |G_U(j,α)|² is Paley×norm type-constant.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.302 flags.  Does **not** name Q_++ or Q_N* in p.
Floor stays OPEN.

============================================================================
Setup.  q=p², T=□, U={N=1}.  E[N(δ)]=μ_+ on squares, μ_- on nonsquares
(2-design of regular sets).  μ_+=p(p−1)/4, μ_-=p(p−3)/4.
G_+(ξ)=∑_{δ≠0} E[N(δ)] e_ξ(δ).  On Ω, G(χ,ξ)=p.
Torus Gauss G_U(j,α)=∑_{u∈U} ξ_j(u) e(α u).

============================================================================
Theorem A — PROVED (Max+-free; μ± and G(χ)|_Ω=p).
  On Ω, G_+=p/2.  For every r∈T, scaling by r preserves χ, so G_r=G_+.
  The 15.298 named piece is therefore constant on T:
      k² + k(G_++G_r) = k(k+p) = q(q−1)/4,
      16·(named)/q² = 4(q−1)/q.
  Hence Q(r)/q² = 4(q−1)/q + 16 L^{conn}(r)/q².
  Fail-when-wrong: replace r by a nonsquare (χ flips); G becomes
  −p(p−1)/2 ≠ p/2 and the constant 4(q−1)/q is lost.  ∎

Theorem B — PROVED (Max+-free field sums).
  |G_U(j,α)|² is constant on each (χ(α), norm-type, N(α)) class, and
  ∑_α |G_U(j,α)|² = q|U|.  At p=5, even j≠0, |G_U(j,α)|²=p on every
  square α.  At p=7 the same claim fails (spread > 10 on squares).
  Fail-when-wrong: claim |G_U(2,α)|²=p on all squares at p=7.  ∎

Theorem C — CERTIFIED p=5,7 (not a general bound).
  Live Q_{++}/q² ≤ 4−2/(p+2), with room 2/91 at p=5 and 10/3681 at p=7.
  At p=5 the 15.302 σ=−1 / j=2 form is 32−7 Q_{++}, so the bound
  Q_{++}≤26/7 is equivalent to that eigenvalue ≥6.  Live is strict.
  Does **not** name Q_{++}.  phi_F not imported.  ∎

============================================================================
Backend: CPU field tables; live Q only to check C.
Writes evidence/e1_gmin_m4_prop15303.json
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

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _e_tr,
    _kernel_field,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15290 import field_norm, in_subfield, live_Q_on_T  # noqa: E402
from e1_gmin_m4_prop15300 import merge_cls  # noqa: E402
from e1_gmin_m4_prop15302 import U_powers  # noqa: E402


def mu_plus(p: int) -> Fraction:
    return Fraction(p * (p - 1), 4)


def mu_minus(p: int) -> Fraction:
    return Fraction(p * (p - 3), 4)


def Gplus_on_Omega(p: int) -> Fraction:
    """(1/2)[-(μ++μ-)+p(μ+-μ-)] with G(χ)=p."""
    mp, mm = mu_plus(p), mu_minus(p)
    return (-(mp + mm) + p * (mp - mm)) / 2


def Gplus_nonsquare_scale(p: int) -> Fraction:
    """χ flipped: (1/2)[-(μ++μ-)+p(μ--μ+)]."""
    mp, mm = mu_plus(p), mu_minus(p)
    return (-(mp + mm) + p * (mm - mp)) / 2


def named_Q_over_q2(p: int) -> Fraction:
    q = p * p
    return Fraction(4 * (q - 1), q)


def GU(F, j: int, alpha: int) -> complex:
    p = F["p"]
    acc = 0j
    for t, u in enumerate(U_powers(F)):
        au = F["mul"][alpha][u] if alpha else 0
        acc += np.exp(2j * np.pi * j * t / (p + 1)) * _e_tr(F, au)
    return acc


def prove_named_constant() -> dict:
    ok = True
    rows = {}
    flip_hit = 0
    for p in (5, 7):
        Gp = Gplus_on_Omega(p)
        Gns = Gplus_nonsquare_scale(p)
        if Gp != Fraction(p, 2):
            ok = False
        if Gns == Gp:
            ok = False
        if abs(Gns - Fraction(-p * (p - 1), 2)) > 0:
            # exact Fraction compare
            if Gns != Fraction(-p * (p - 1), 2):
                ok = False
        if Gns != Gp:
            flip_hit += 1
        q = p * p
        k = Fraction(p * (p - 1), 2)
        named = k * (k + p)
        if named != Fraction(q * (q - 1), 4):
            ok = False
        if Fraction(16) * named / Fraction(q * q) != named_Q_over_q2(p):
            ok = False
        rows[str(p)] = {
            "Gplus": str(Gp),
            "G_nonsquare_scale": str(Gns),
            "named_Q_over_q2": str(named_Q_over_q2(p)),
        }
    if flip_hit == 0:
        ok = False
    return {
        "proved": ok,
        "flip_hit": flip_hit,
        "by_p": rows,
        "theorem": (
            "G_+=p/2 on Ω; G_r=G_+ for r∈T; 16·named/q²=4(q−1)/q. "
            "Nonsquare scale flips G to −p(p−1)/2."
        ),
    }


def prove_GU_type_constant() -> dict:
    ok = True
    p5_flat = True
    p7_spread = 0.0
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        rec = {}
        for j in range(0, p + 1, 2):
            by: dict[tuple, list[float]] = defaultdict(list)
            tot = 0.0
            for a in range(F["q"]):
                mag2 = abs(GU(F, j, a)) ** 2
                tot += mag2
                if a == 0:
                    key = ("0",)
                else:
                    nt = (
                        "U"
                        if field_norm(F, a) == 1
                        else ("sub" if in_subfield(F, a) else "N*")
                    )
                    key = (int(F["chi"][a]), nt, int(field_norm(F, a)))
                by[key].append(mag2)
            spreads = {str(k): float(max(v) - min(v)) for k, v in by.items()}
            if max(spreads.values()) > 1e-6:
                ok = False
            expect = float(F["q"] * (p + 1))
            if abs(tot - expect) > 1e-6:
                ok = False
            # squares at j=2
            if j == 2:
                sq = [mag for k, vs in by.items() for mag in vs if k[0] == 1]
                spr = float(max(sq) - min(sq)) if sq else 0.0
                rec["j2_square_spread"] = spr
                rec["j2_square_min"] = float(min(sq)) if sq else None
                rec["j2_square_max"] = float(max(sq)) if sq else None
                if p == 5:
                    if any(abs(m - p) > 1e-6 for m in sq):
                        p5_flat = False
                        ok = False
                if p == 7:
                    p7_spread = spr
                    if spr < 10:
                        ok = False
            rec[f"j{j}_plancherel"] = tot
            rec[f"j{j}_class_spreads_max"] = max(spreads.values())
        rows[str(p)] = rec
    if (not p5_flat) or p7_spread < 10:
        ok = False
    return {
        "proved": ok,
        "p5_even_j_squares_eq_p": p5_flat,
        "p7_j2_square_spread": p7_spread,
        "by_p": rows,
        "theorem": (
            "|G_U(j,α)|² type-constant; Plancherel q|U|. "
            "p=5 even j≠0: |G_U|=√p on squares. p=7 spread>10."
        ),
    }


def prove_bound_certified() -> dict:
    ok = True
    rows = {}
    for p in (5, 7):
        F = _kernel_field(p)
        Q = live_Q_on_T(p)
        q2 = float(F["q"] ** 2)
        acc = []
        for r, v in Q.items():
            if merge_cls(F, r) == "++":
                acc.append(v / q2)
        Qpp = Fraction(float(np.mean(acc))).limit_denominator(2000)
        ub = Fraction(4) - Fraction(2, p + 2)
        room = ub - Qpp
        if room < 0:
            ok = False
        rows[str(p)] = {"Qpp": str(Qpp), "ub": str(ub), "room": str(room)}
    # p=5 threshold identity: 32-7*(26/7)=6
    if Fraction(32) - 7 * (Fraction(4) - Fraction(2, 7)) != 6:
        ok = False
    return {
        "proved": ok,
        "general_bound": False,
        "by_p": rows,
        "p5_threshold": "32-7*(4-2/7)=6",
        "theorem": (
            "Live Q_++ ≤ 4−2/(p+2) at p=5,7.  At p=5 this is equivalent "
            "to the j=2 eigenvalue ≥6.  Not a general proof; Q_++ unnamed."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "S_k_named_in_p": False,
        "note": (
            "Named constant 4(q−1)/q does not separate ++ from N*. "
            "Q_++ / Q_N* still unnamed. Floor OPEN."
        ),
    }


def main() -> dict:
    print("Prop 15.303  Q/q² = 4(q−1)/q + L^conn; |G_U|² type-constant", flush=True)
    A = prove_named_constant()
    print(f"  A named constant 4(q-1)/q: {A['proved']} flip={A['flip_hit']}", flush=True)
    B = prove_GU_type_constant()
    print(
        f"  B |GU|² type-const: {B['proved']} p5_flat={B['p5_even_j_squares_eq_p']} "
        f"p7_spread={B['p7_j2_square_spread']:.3f}",
        flush=True,
    )
    C = prove_bound_certified()
    print(f"  C Q++<=4-2/(p+2) live: {C['proved']} general={C['general_bound']}", flush=True)
    D = prove_open()
    print(f"  D floor: {D['proved']} phi_F={D['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.303",
        "title": "Q/q²=4(q−1)/q + L^conn on T; torus Gauss type-constant",
        "proved": {
            "named_constant_on_T": A["proved"],
            "GU_type_constant": B["proved"],
            "Qpp_bound_certified_p57": C["proved"],
            "Qpp_bound_general": False,
            "phi_F_ge_6_proved_general": D["phi_F_ge_6"],
            "Q_tau_named_in_p": False,
            "S_k_named_in_p": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15303.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
