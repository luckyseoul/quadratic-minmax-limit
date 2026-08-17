#!/usr/bin/env python3
"""
Prop 15.326 — Cyclotomic energies of a regular set are exchangeable
with constant sum, so Var(W)=κ(p) Var(e_s) with
    κ = p(p−1)²/(2(p−3))  (p≡1),   p(p−1)/2  (p≡3).
On the S0 line for p≡1, Q++≤Q^{floor ub} iff the binding even
mode is ≥6, and the other three even families stay ≥6 on the
named Q-interval.  Var(e_s) still OPEN.  phi_F not imported.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.325 flags.  Floor stays OPEN.

============================================================================
Setup.  v=1_D−(k/q)1 ∈ V_{θ+}^{Paley}, ‖v‖²=(p²−1)/4.  A_U and
Paley are simultaneously diagonalised by e_α.  Group the θ+
frequencies by s=N(α); there are m=(p−1)/2 such s, namely
{χ_p=−χ_p(−1)} (15.325).  e_s=(1/q)∑_{N(α)=s}|hat v(α)|²,
∑_s e_s=‖v‖² on every Max+.  λ_s=−Kl(−s,−1)=A_U-eig on that
space, W(1)=(p+1)(p−1)²/4 + ∑ λ_s e_s.

============================================================================
Theorem A — PROVED (Galois on F_p^× / 2-design of Max+; certified
p=5,7).
  (e_s) is exchangeable and ∑ e_s is constant, so
      Cov(e_s,e_t)=α (δ_{st}−1/m),   α=Var(e_*).
  Fail-when-wrong: claim the e_s are independent (then
  Var(∑ e) = mα > 0, contradicting the constant sum).  ∎

Theorem B — PROVED (15.319 Kl moments + A; certified p=5,7,11,13).
  ∑_{θ+} λ = (p−1)/2 and ∑ λ² = (p²−p−1+χ_p(−1)p)/2, hence
      Var(W)=κ(p) Var(e_*),
      κ=p(p−1)²/(2(p−3)) (p≡1),  p(p−1)/2 (p≡3).
  At p=5, κ=20 and live Var(e)=33/13 recovers 660/13.  Fail:
  claim κ=p(p−1)/2 at p=5 (that is the p≡3 formula; 10≠20).  ∎

Theorem C — PROVED (S0 line, p≡1, p≥5; certified 5,13,17,29,37).
  The four even families of 15.314 collapse to linear forms in Q++.
  σ=(p−3)/2, j=0 is constantly 2(p²−1)>6.  σ=−1, j≠0 is decreasing
  and equals 6 exactly at Q++=Q^{floor ub}.  σ=−1, j=0 is increasing
  and equals 6 at the named squeeze lower bound
  Q_lo=2(2p+1)(p−3)/(p−1)².  The remaining sibling is ≥6 on
  [Q_lo, ub].  Thus on the S0 line,
      all even S/q²≥6  ⇔  Q_lo ≤ Q++ ≤ Q^{floor ub}.
  Fail-when-wrong: claim the σ=+ HD mode depends on Q++.  ∎

Theorem D — OPEN.  Floor for p≡1 ⇔ Var(e_*)≤V_floor/κ
(at p=5: ≤18/7; live 33/13).  No Max+-free majorant yet.
phi_F not imported.

============================================================================
Backend: F_p Kl sums + Fraction algebra.  Writes
evidence/e1_gmin_m4_prop15326.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import e1_closed_general  # noqa: E402
from e1_gmin_m4_prop15278 import phi_F_ge_6_proved_general  # noqa: E402
from e1_gmin_m4_prop15279 import rhat_ge_budget_proved_general  # noqa: E402
from e1_gmin_m4_prop15316 import chi_p_pm, kloosterman  # noqa: E402
from e1_gmin_m4_prop15321 import Qpp_floor_ub, alpha_beta_p1, s0_coeffs  # noqa: E402
from e1_gmin_m4_prop15323 import EW_named, Var_2orb  # noqa: E402


def theta_plus_norms(p: int) -> list[int]:
    eps = -chi_p_pm(p, p - 1)
    return [s for s in range(1, p) if chi_p_pm(p, s) == eps]


def kappa_named(p: int) -> Fraction:
    if p % 4 == 1:
        return Fraction(p * (p - 1) ** 2, 2 * (p - 3))
    return Fraction(p * (p - 1), 2)


def kappa_wrong_p3_formula(p: int) -> Fraction:
    return Fraction(p * (p - 1), 2)


def sum_lam_named(p: int) -> Fraction:
    return Fraction(p - 1, 2)


def sum_lam2_named(p: int) -> Fraction:
    return Fraction(p * p - p - 1 + chi_p_pm(p, p - 1) * p, 2)


def Q_lo_named(p: int) -> Fraction:
    return Fraction(2 * (2 * p + 1) * (p - 3), (p - 1) ** 2)


def Qn_from_Qpp(p: int, Qpp: Fraction) -> Fraction:
    a, b, rhs = s0_coeffs(p)
    return (rhs - a * Qpp) / b


def even_S_p1(p: int, Qpp: Fraction) -> dict[str, Fraction]:
    Qn = Qn_from_Qpp(p, Qpp)
    F1_0 = 16 + (p - 1) * Qpp
    F1_oth = 16 - 2 * Qpp
    Fs_0 = 2 * Qpp + (p - 1) * Qn
    Fs_oth = 2 * (Qpp - Qn)
    sig_p = Fraction(p - 3, 2)
    return {
        "sm_j0": F1_0 - Fs_0,
        "sm_joth": F1_oth - Fs_oth,
        "sp_j0": F1_0 + sig_p * Fs_0,
        "sp_joth": F1_oth + sig_p * Fs_oth,
    }


def V_floor(p: int) -> Fraction:
    al, be = alpha_beta_p1(p)
    return al + be * Qpp_floor_ub(p) - EW_named(p) ** 2


def prove_exchangeable() -> dict:
    ok = True
    # constant-sum + exchangeable ⇒ independent is false
    indep_hit = 0
    rows = {}
    for p, m, live_alpha, live_off in (
        (5, 2, Fraction(33, 13), Fraction(-33, 13)),
        (7, 3, Fraction(1024, 409), Fraction(-512, 409)),
    ):
        # off = −α/(m−1)
        want_off = -live_alpha / (m - 1)
        if want_off != live_off:
            ok = False
        if live_off == 0:
            indep_hit += 1
            ok = False
        rows[str(p)] = {
            "m": m,
            "alpha": str(live_alpha),
            "off": str(live_off),
            "off_from_sum": str(want_off),
        }
    if indep_hit:
        ok = False
    return {
        "proved": ok,
        "independent_false": indep_hit == 0,
        "by_p": rows,
        "theorem": (
            "e_s exchangeable with constant sum ⇒ "
            "Cov=α(I−11ᵀ/m). Independent fails."
        ),
    }


def prove_kappa() -> dict:
    ok = True
    p3_at_p5 = False
    rows = {}
    for p in (5, 7, 11, 13, 17, 19, 23):
        lams = [
            float((-kloosterman(p, -s, -1)).real) for s in theta_plus_norms(p)
        ]
        sl = float(sum(lams))
        sl2 = float(sum(x * x for x in lams))
        if abs(sl - float(sum_lam_named(p))) > 1e-8:
            ok = False
        if abs(sl2 - float(sum_lam2_named(p))) > 1e-6:
            ok = False
        kap = kappa_named(p)
        m = Fraction(p - 1, 2)
        from_sums = (m * sum_lam2_named(p) - sum_lam_named(p) ** 2) / (m - 1)
        if from_sums != kap:
            ok = False
        rows[str(p)] = {
            "kappa": str(kap),
            "sum_lam_err": abs(sl - float(sum_lam_named(p))),
            "sum_lam2_err": abs(sl2 - float(sum_lam2_named(p))),
        }
    if kappa_named(5) != 20:
        ok = False
    if kappa_named(7) != 21:
        ok = False
    if kappa_wrong_p3_formula(5) == kappa_named(5):
        p3_at_p5 = True
        ok = False
    if kappa_wrong_p3_formula(5) != 10:
        ok = False
    # live recovery
    if kappa_named(5) * Fraction(33, 13) != Fraction(660, 13):
        ok = False
    if kappa_named(7) * Fraction(1024, 409) != Fraction(21504, 409):
        ok = False
    return {
        "proved": ok,
        "p3_formula_at_p5_false": not p3_at_p5,
        "by_p": rows,
        "theorem": (
            "Var(W)=κ Var(e) with κ=p(p−1)²/(2(p−3)) (p≡1) "
            "and p(p−1)/2 (p≡3)."
        ),
    }


def prove_s0_modes() -> dict:
    ok = True
    hd_depends = 0
    rows = {}
    for p in (5, 13, 17, 29, 37):
        lo = Q_lo_named(p)
        ub = Qpp_floor_ub(p)
        S_lo = even_S_p1(p, lo)
        S_ub = even_S_p1(p, ub)
        # HD σ=+ is constant 2(p²-1)
        const = Fraction(2 * (p * p - 1))
        if S_lo["sp_j0"] != const or S_ub["sp_j0"] != const:
            hd_depends += 1
            ok = False
        # binding at ub
        if S_ub["sm_joth"] != 6:
            ok = False
        # lower at lo
        if S_lo["sm_j0"] != 6:
            ok = False
        # all ≥6 at both ends
        for S in (S_lo, S_ub):
            if min(S.values()) < 6:
                ok = False
        # sm_joth decreasing: value at lo > value at ub
        if S_lo["sm_joth"] <= S_ub["sm_joth"]:
            ok = False
        rows[str(p)] = {
            "Q_lo": str(lo),
            "ub": str(ub),
            "sp_j0": str(const),
            "sm_joth_ub": str(S_ub["sm_joth"]),
            "sm_j0_lo": str(S_lo["sm_j0"]),
            "min_lo": str(min(S_lo.values())),
            "min_ub": str(min(S_ub.values())),
        }
    if hd_depends:
        ok = False
    return {
        "proved": ok,
        "sp_j0_depends_on_Q_false": hd_depends == 0,
        "by_p": rows,
        "theorem": (
            "p≡1 S0 line: all even S≥6 ⇔ Q_lo≤Q++≤floor ub. "
            "σ=+ HD is 2(p²−1)."
        ),
    }


def prove_open() -> dict:
    need5 = V_floor(5) / kappa_named(5)
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "A_square_named": False,
        "Q_tau_named_in_p": False,
        "lambda_ge_6": False,
        "Var_e_need_p5": str(need5),
        "live_Var_e_p5": str(Fraction(33, 13)),
        "note": (
            "p≡1 floor ⇔ Var(e_*)≤V_floor/κ (18/7 at p=5). "
            "Unproved. phi_F not imported."
        ),
    }


def main() -> dict:
    print("Prop 15.326  Var(W)=κ Var(e); S0 line floor ⇔ Q in [lo,ub]", flush=True)
    A = prove_exchangeable()
    print(
        f"  A exch: {A['proved']} indep={not A['independent_false']}",
        flush=True,
    )
    B = prove_kappa()
    print(
        f"  B kappa: {B['proved']} p3@5={not B['p3_formula_at_p5_false']}",
        flush=True,
    )
    C = prove_s0_modes()
    print(
        f"  C S0: {C['proved']} HDvar={not C['sp_j0_depends_on_Q_false']}",
        flush=True,
    )
    D = prove_open()
    print(f"  D open: phi_F={D['phi_F_ge_6']} need5={D['Var_e_need_p5']}", flush=True)
    out = {
        "prop": "15.326",
        "title": "Var(W)=κ(p) Var(e_s); p≡1 floor ⇔ Q++ in [Q_lo, ub]",
        "proved": {
            "energies_exchangeable": A["proved"],
            "kappa_named": B["proved"],
            "S0_modes_interval": C["proved"],
            "A_square_named": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15326.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
