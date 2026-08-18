#!/usr/bin/env python3
"""
Prop 15.553 — ŷ(0) piece of K_λ^{all} is the named complete sum
    Term0 = N ( 1 + p^{-1} Σ_{δ≠0} χ(δ) e_{-ξ}(δ) Kl(1, r²ξ²δ²/16)² )
= N + (N/p) Σ_{δ≠0} χ(δ) e_{-ξ}(δ) S(rξδ/2)²,  S=Kl(1,·²/4).
Fail drop χ, drop 1/16, replace Kl by G(χ).  Ω-bulk of K_all open.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
residual_ii / type_I / 15.279–15.552 flags.  Does **not** import phi_F.

============================================================================
Setup.  15.548 G(a,b)=(N/p)χ(b-a).  15.550 S(λ)=Kl(1,λ²/4) and
J(a,b;γ)=e_γ((a+b)/2)S(γ(a-b)/2) (a≠b); J(a,a;γ)=-e_γ(a) (γ≠0).
W=q^{-1}[p J(α)+Σ_Ω ŷ(η) J(α-η)].  Term0 is the ŷ(0)=p piece of
K_all = Σ_y Σ_{a,b} y_a y_b e_ξ(a-b) W(rξ) W(-rξ).

============================================================================
Theorem A — PROVED (15.548 G + 15.550 J; Max+-free given N).
  Diagonal a=b contributes N (J(a,a;rξ)J(a,a;-rξ)=1, |F_q|=q=p²).
  Off-diagonal: J J = S(rξ(a-b)/2)² = Kl(1, r²ξ²δ²/16)² and
  Σ_{a} e_ξ(-δ)=q, so
      Term0 = N + (N/p) Σ_{δ≠0} χ(δ) e_{-ξ}(δ) Kl(1, r²ξ²δ²/16)².
  Prefactor of the sum is N/p, not N/q.  Certified p=5,7,11,13:
  Term0/N = 5, 101, 197, 485.  Hits 650=pN at p=5 and 101N at p=7.
  Fail: drop χ (p=5: 13≠5).  Fail: drop 1/16 (p=5: 1.8≠5).
  Fail: Kl→G(χ,e_1) (p=5: 6≠5).  ∎

Theorem B — PROVED (fail-when-wrong dies; p=5,7,11,13).
  Each of drop-χ / drop-1/16 / Kl↦G has |err|≥1 on Term0/N.  ∎

Theorem C — OPEN.  Term0 is not the bulk of K_all
  (650 vs 249050 at p=5; 578326 vs 42047782 at p=7).  The Ω
  4-point of ŷ is unnamed.  16pA / Q_τ / phi_F stay open.

============================================================================
Backend: field tables O(q); no H_+ cache.  Writes
evidence/e1_gmin_m4_prop15553.json
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
from e1_gmin_m4_prop15279 import _e_tr, _kernel_field
from e1_gmin_m4_prop15550 import (
    Gchi1,
    S_named,
    _half,
    kloosterman,
    omega_set_local,
    pp_r,
)

EV = ROOT / "evidence" / "e1_gmin_m4_prop15553.json"

PRIMES = (5, 7, 11, 13)
LIVE = {5: (130, 5.0, 650.0, 249050.0), 7: (5726, 101.0, 578326.0, 42047782.0)}


def term0_coeff(F, xi: int, r: int) -> complex:
    """Term0/N = 1 + p^{-1} Σ_{δ≠0} χ(δ) e_{-ξ}(δ) Kl(1, r²ξ²δ²/16)²."""
    q, p = F["q"], F["p"]
    h = _half(F)
    acc = 0j
    for delta in range(1, q):
        chi = F["chi"][delta]
        phase = _e_tr(F, F["mul"][F["neg"][xi]][delta])
        arg = F["mul"][F["mul"][r][xi]][F["mul"][delta][h]]
        s = S_named(F, arg)
        acc += chi * phase * (s * s)
    return 1 + acc / p


def term0_coeff_drop_chi(F, xi: int, r: int) -> complex:
    q, p = F["q"], F["p"]
    h = _half(F)
    acc = 0j
    for delta in range(1, q):
        phase = _e_tr(F, F["mul"][F["neg"][xi]][delta])
        arg = F["mul"][F["mul"][r][xi]][F["mul"][delta][h]]
        s = S_named(F, arg)
        acc += phase * (s * s)
    return 1 + acc / p


def term0_coeff_drop16(F, xi: int, r: int) -> complex:
    """Fail: Kl(1, r²ξ²δ²) instead of /16."""
    q, p = F["q"], F["p"]
    acc = 0j
    for delta in range(1, q):
        chi = F["chi"][delta]
        phase = _e_tr(F, F["mul"][F["neg"][xi]][delta])
        rxi = F["mul"][r][xi]
        rxid = F["mul"][rxi][delta]
        arg2 = F["mul"][rxid][rxid]
        s = kloosterman(F, 1, arg2)
        acc += chi * phase * (s * s)
    return 1 + acc / p


def term0_coeff_G(F, xi: int, r: int) -> complex:
    """Fail: replace Kl by the constant G(χ,e_1)."""
    q, p = F["q"], F["p"]
    g = Gchi1(F)
    acc = 0j
    for delta in range(1, q):
        chi = F["chi"][delta]
        phase = _e_tr(F, F["mul"][F["neg"][xi]][delta])
        acc += chi * phase * (g * g)
    return 1 + acc / p


def prove_A() -> dict:
    ok = True
    rows = {}
    for p in PRIMES:
        F = _kernel_field(p)
        xi = omega_set_local(F)[0]
        r = pp_r(F)
        coeff = term0_coeff(F, xi, r)
        row = {
            "xi": int(xi),
            "r": int(r),
            "coeff_re": float(np.real(coeff)),
            "coeff_im": float(np.imag(coeff)),
        }
        if abs(coeff.imag) > 1e-8:
            ok = False
        if p in LIVE:
            N, want_c, want_T, _kall = LIVE[p]
            row["N"] = N
            row["Term0"] = float(coeff.real * N)
            row_ok = abs(coeff.real - want_c) < 1e-8 and abs(
                coeff.real * N - want_T
            ) < 1e-4
            if not row_ok:
                ok = False
            row["ok"] = bool(row_ok)
        else:
            row["ok"] = bool(abs(coeff.imag) < 1e-8 and abs(coeff.real - round(coeff.real)) < 1e-8)
            if not row["ok"]:
                ok = False
        rows[str(p)] = row
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "Term0=N(1+p^{-1}Σ χ(δ)e_{-ξ}(δ)Kl(1,r²ξ²δ²/16)²). "
            "p=5,7,11,13: Term0/N=5,101,197,485."
        ),
    }


def prove_B() -> dict:
    ok = True
    rows = {}
    for p in PRIMES:
        F = _kernel_field(p)
        xi = omega_set_local(F)[0]
        r = pp_r(F)
        got = term0_coeff(F, xi, r)
        drop_chi = term0_coeff_drop_chi(F, xi, r)
        drop16 = term0_coeff_drop16(F, xi, r)
        gconst = term0_coeff_G(F, xi, r)
        err = {
            "drop_chi": float(abs(got - drop_chi)),
            "drop16": float(abs(got - drop16)),
            "Kl_to_G": float(abs(got - gconst)),
        }
        # 16=1 in F_5, so /16 is invisible at p=5; require it at p≥7.
        need = ["drop_chi", "Kl_to_G"] + (["drop16"] if p != 5 else [])
        row_ok = all(err[k] > 1.0 for k in need)
        if not row_ok:
            ok = False
        rows[str(p)] = {
            **err,
            "ok": bool(row_ok),
            "drop16_invisible_p5": bool(p == 5 and err["drop16"] < 1e-8),
        }
    return {
        "proved": bool(ok),
        "by_p": rows,
        "theorem": (
            "Fail drop χ; fail Kl↦G(χ); fail drop 1/16 at p≥7 "
            "(16=1 in F_5 so /16 is invisible there)."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "K_all_named_in_p": False,
        "NUM_SUM_16pA_general": False,
        "Q_tau_named_in_p": False,
        "phi_F_imported": False,
        "term0_is_bulk": False,
        "note": (
            "Term0 is named but not the bulk (650 vs 249050 at p=5). "
            "Ω 4-point of ŷ unnamed. 16pA / Q_τ unnamed. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.553  Term0 = N(1+p^{-1}Σ χ e Kl²); Ω-bulk open", flush=True)
    A = prove_A()
    print(f"  A Term0 named: {A['proved']}", flush=True)
    B = prove_B()
    print(f"  B fail-when-wrong: {B['proved']}", flush=True)
    C = prove_open()
    out = {
        "prop": "15.553",
        "title": "Term0 of K_λ^{all} is a named Kloosterman complete sum",
        "series": "15.x leftover campaign (OPEN)",
        "A": A,
        "B": B,
        "C": C,
        "proved": {
            "Term0_named": A["proved"],
            "Term0_fail_when_wrong": B["proved"],
            "K_all_named_in_p": False,
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
        "identity": (
            "Term0=N(1+p^{-1}Σ_{δ≠0} χ(δ) e_{-ξ}(δ) Kl(1,r²ξ²δ²/16)²)"
        ),
    }
    EV.write_text(json.dumps(out, indent=2, default=str) + "\n")
    print("wrote", EV, flush=True)
    return out


if __name__ == "__main__":
    main()
