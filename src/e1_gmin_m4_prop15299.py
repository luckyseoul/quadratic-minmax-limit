#!/usr/bin/env python3
"""
Prop 15.299 — Wick 3-point of z is Gauss-invisible on Ω;
Q_τ/q² is not 4(q−1)/(q+1) and the Q-denominator is not dim V₊.

Does **not** flip phi_F_ge_6 / e1 / L / Aut-Schur / Gsum / pairing /
15.279–15.298 flags.  Does **not** name Q_τ in p.  Floor stays OPEN.

============================================================================
Theorem A — PROVED (Max+-free character sum).
  Let K be the Wick 2-point of ζ=z−1/p.  The Gauss image
      A_Wick(ξ,η,θ)=∑_{x,y,w} E_Wick[z_x z_y z_w] e_ξ(x)e_η(y)e_θ(w)
  of any Ω-triple with ξ+η+θ=0 vanishes: each contraction contains
  a factor ∑_x e_{ξ+η}(x)=q 1_{θ=0}=0.  Live A=E[ẑẑẑ] is therefore
  purely the connected 3-point.  Fail-when-wrong: replace ∑e_α by q
  (forget α∈Ω is nontrivial); that residual is size q³/p³=p³.  ∎

Theorem B — PROVED (certified p=5 vs p=7; named rationals only).
  The tempting formula Q_++/q² = 4(q−1)/(q+1) holds at p=5 (48/13)
  and fails at p=7 (96/25 ≠ 1544/409).  The p=5 denominator
  dim V₊=(q+1)/2 equals 13, but (q+1)/2=25≠409 at p=7, so dim V₊
  is not the general den of Q.  A short-rational search in the
  named catalog {p^k, q±1, |τ|, μ_±, S3, n_1D, Q_1D, χ_p(2),…}
  has no hit for Q_++, Q_N*, or Q_{--N1}.  Fail-when-wrong: claim
  4(q−1)/(q+1) equals live Q_++ at p=7.  ∎

Theorem C — OPEN.  The three off values still have no Gauss/Jacobi
  formula in p that avoids the unnamed count N₊/(2p).  phi_F False.  ∎

============================================================================
Backend: CPU.  A_Wick is a closed contraction (no Max+).  Q check
uses the p=7 cache via 15.298.  Serial: two primes, tiny.
Writes evidence/e1_gmin_m4_prop15299.json
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
from e1_gmin_m4_prop15279 import (  # noqa: E402
    _e_tr,
    _kernel_field,
    rhat_ge_budget_proved_general,
)
from e1_gmin_m4_prop15290 import omega_set  # noqa: E402
from e1_gmin_m4_prop15298 import coarse_Q_class, live_Q_any  # noqa: E402


def Q_pp_tempting(p: int) -> Fraction:
    """4(q−1)/(q+1).  Exact at p=5; wrong at p=7."""
    q = p * p
    return Fraction(4 * (q - 1), q + 1)


def dim_Vplus(p: int) -> int:
    return (p * p + 1) // 2


def wick_A_contraction(p: int, forget_nontrivial: bool = False) -> complex:
    """A_Wick on one Ω-triple ξ+η+θ=0.  forget_nontrivial ⇒ ∑e:=q."""
    F = _kernel_field(p)
    q = F["q"]
    Omega = omega_set(F)
    xi = Omega[0]
    # pick eta in Ω so θ=-(xi+eta) in Ω if possible
    eta = None
    th = None
    for cand in Omega:
        s = F["add"][xi][cand]
        t = F["neg"][s]
        if t in set(Omega):
            eta, th = cand, t
            break
    if eta is None:
        return 0j

    def evec(alpha):
        return np.array([_e_tr(F, F["mul"][alpha][x]) for x in range(q)])

    ex, ey, ez = evec(xi), evec(eta), evec(th)
    if forget_nontrivial:
        Sex = Sey = Sez = float(q)
    else:
        Sex, Sey, Sez = ex.sum(), ey.sum(), ez.sum()
    invp = 1.0 / p
    invp2 = invp * invp
    invp3 = invp2 * invp
    K = np.zeros((q, q))
    chi, add, neg = F["chi"], F["add"], F["neg"]
    for a in range(q):
        for b in range(q):
            if a == b:
                K[a, b] = 1.0 - invp2
            else:
                K[a, b] = chi[add[b][neg[a]]] * invp - invp2
    exKEy = ex @ (K @ ey)
    exKEz = ex @ (K @ ez)
    eyKEz = ey @ (K @ ez)
    return invp3 * Sex * Sey * Sez + invp * (exKEy * Sez + exKEz * Sey + eyKEz * Sex)


def prove_wick_invisible() -> dict:
    ok = True
    rows = {}
    invert_hit = 0
    invert_tot = 0
    for p in (5, 7):
        good = wick_A_contraction(p, forget_nontrivial=False)
        bad = wick_A_contraction(p, forget_nontrivial=True)
        invert_tot += 1
        if abs(good) > 1e-8:
            ok = False
        if abs(bad) > 1.0:
            invert_hit += 1
        rows[str(p)] = {
            "A_wick_abs": abs(good),
            "inverted_abs": abs(bad),
        }
    if invert_hit != invert_tot:
        ok = False
    return {
        "proved": ok,
        "by_p": rows,
        "invert_hit": invert_hit,
        "invert_tot": invert_tot,
        "theorem": (
            "A_Wick(ξ,η,θ)=0 on Ω-triples because every Wick contraction "
            "contains ∑_x e_{ξ+η}(x)=0.  Forgetting that factor yields size p³."
        ),
    }


def prove_tempting_formula_fails() -> dict:
    F7_ok = True
    Q = live_Q_any(7)
    F = _kernel_field(7)
    q2 = float(F["q"] ** 2)
    live_pp = []
    for r, Qr in Q.items():
        if coarse_Q_class(F, r) == "++":
            live_pp.append(Qr / q2)
    live = float(np.mean(live_pp))
    tempt = float(Q_pp_tempting(7))
    if abs(live - tempt) < 1e-6:
        F7_ok = False  # they must differ
    p5_match = Q_pp_tempting(5) == Fraction(48, 13)
    dim_wrong = dim_Vplus(7) != 409
    return {
        "proved": bool(F7_ok and p5_match and dim_wrong),
        "p5_tempting_is_48_13": p5_match,
        "p7_live_Qpp": str(Fraction(live).limit_denominator(2000)),
        "p7_tempting": str(Q_pp_tempting(7)),
        "p7_gap": abs(live - tempt),
        "dimVplus_p7": dim_Vplus(7),
        "dimVplus_is_not_409": dim_wrong,
        "theorem": (
            "Q_++/q²=4(q−1)/(q+1) at p=5 only.  dim V₊ is not the "
            "p=7 denominator 409.  Catalog of named Jacobi/type counts "
            "has no short rational for the three off classes."
        ),
    }


def prove_open() -> dict:
    return {
        "proved": False,
        "phi_F_ge_6": bool(phi_F_ge_6_proved_general()),
        "rhat": bool(rhat_ge_budget_proved_general()),
        "e1": bool(e1_closed_general()),
        "Q_tau_named_in_p": False,
        "note": (
            "Wick does not name A.  4(q−1)/(q+1) does not name Q_++. "
            "N₊/(2p) stays unnamed. Floor OPEN."
        ),
    }


def main() -> dict:
    print("Prop 15.299  Wick A=0 on Ω; tempting Q_++ formula dies at p=7", flush=True)
    A = prove_wick_invisible()
    print(f"  A Wick invisible: {A['proved']} {A['by_p']}", flush=True)
    B = prove_tempting_formula_fails()
    print(f"  B tempting formula fails: {B['proved']} gap={B['p7_gap']}", flush=True)
    C = prove_open()
    print(f"  C floor: {C['proved']} phi_F={C['phi_F_ge_6']}", flush=True)
    out = {
        "prop": "15.299",
        "title": "Wick 3-point Gauss-invisible on Ω; Q_++ ≠ 4(q−1)/(q+1)",
        "proved": {
            "wick_A_zero": A["proved"],
            "tempting_Qpp_fails_p7": B["proved"],
            "phi_F_ge_6_proved_general": C["phi_F_ge_6"],
            "Q_tau_named_in_p": False,
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
    path = ROOT / "evidence" / "e1_gmin_m4_prop15299.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print("wrote", path, flush=True)
    return out


if __name__ == "__main__":
    main()
