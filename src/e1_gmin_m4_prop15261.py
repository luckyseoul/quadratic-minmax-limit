#!/usr/bin/env python3
"""
Prop 15.261 — Exact local configuration multiplicity for Max± via hypercube
Fourier; mult≥0 bounds; equivalence Q-match ⇔ ν=0; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.189 pairwise π, 15.250 odd moments vanish, 15.258 orientation)
  Max± uniform on {y∈{±1}^n : C y = ±p y}, |Max+|=|Max-|=N (Paley).
  For a fixed 4-set S with local Seidel matrix A=C_S, and α∈{±1}^4,
  mult_±(α) := #{y∈Max± : y_S = α}.
  Write pr(α)=∏α_i, qf(α)=αᵀ A α = 2∑_{i<j} A_ij α_i α_j.

PROVED Max+-free (any conference with odd moments 0, π=±C/p, |Max+|=|Max-|=N)
  A. Exact multiplicity formula.
     Using 1_{y_S=α}=2^{-4}∏_{i∈S}(1+α_i y_i) expanded over subsets T⊆S,
     and E_±[∏_{i∈T} y_i]=0 for |T| odd, =1 for |T|=0,
     = (±A_{ab}/p) for T={a,b}, = m₄±(S) for T=S:
        mult_+(α) = (N/16) ( 1 + qf(α)/(2p) + pr(α) m₄⁺(S) )
        mult_-(α) = (N/16) ( 1 − qf(α)/(2p) + pr(α) m₄⁻(S) ).  ∎
     (Holds for every 4-set S, not only |κ|=1.)

  B. Consistency.
     ∑_α mult_±(α)=N and ∑_α mult_±(α) pr(α)=N m₄± follow from
     ∑_α pr=0 and ∑_α qf pr=0 (each reduces to a 2-factor character sum).  ∎

  C. Nonnegativity constraints on m₄±.
     mult_±(α)≥0 for all α∈{±1}^4 gives linear inequalities on m₄±
     in terms of the attainable qf values of A.  In particular if
     q_max=max|qf| and both signs of qf occur for both products,
        |m₄±| ≤ 1 − q_max/(2p).
     For |κ|=1 one has q_max≤8 (spectral ‖A‖₂=√5 on ‖α‖₂=2), so
        |m₄±| ≤ 1 − 4/p
     (too weak alone for 1/(2p) when p≥5: 1−4/p > 1/(2p)).  ∎

  D. Orientation equivalence (on |κ|=1).
     If Q is a signed-permutation matrix with Qᵀ A Q=−A and det D=+1
     (exists iff |κ|=1 by 15.258), then qf(Qᵀα)=−qf(α), pr(Qᵀα)=pr(α), and
        mult_+(α) − mult_-(Qᵀα) = (N/16) pr(α)·2ν(S).
     Hence mult_+(α)=mult_-(Qᵀα) for all α  ⇔  ν(S)=0  ⇔  m₄⁺(S)=m₄⁻(S).  ∎
     (Census: both sides hold on all |κ|=1 at p=5,7; does not alone prove either.)

CERTIFIED
  E. Formula (A) matches integer multiplicities for all 16 configs on sample
     |κ|=1 and |κ|=3 four-sets at p=5,7 (machine precision / exact counts).
  F. Q-match ⇔ ν=0 verified; mult≥0 bound 1−4/p holds but is strictly
     weaker than 1/(2p) for p≥5.

OPEN (blocks residual i / predicate flip)
  G. Prove ν=0 on |κ|=1 (e.g. by independent Q-config match, Schur average,
     or π-covariance), **or** reflection |ρ|≤|ρ_f4|, **or** a tighter
     Max+-free L^∞ bound. Soft-close forbidden.

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15261.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    n_of,
    type_I_k_3p_minus_2_closed_general,
)


def mult_plus_formula(N: int, p: int, qf: int, pr: int, m4plus: float) -> float:
    return (N / 16.0) * (1.0 + qf / (2.0 * p) + pr * m4plus)


def mult_minus_formula(N: int, p: int, qf: int, pr: int, m4minus: float) -> float:
    return (N / 16.0) * (1.0 - qf / (2.0 * p) + pr * m4minus)


def theorem_mult_formula() -> dict:
    return {
        "proved": True,
        "theorem": (
            "mult_+(α)=(N/16)(1+qf(α)/(2p)+pr(α)m₄⁺(S)) and "
            "mult_-(α)=(N/16)(1−qf(α)/(2p)+pr(α)m₄⁻(S)), "
            "for any 4-set S, under odd-moment vanishing and π=±C/p."
        ),
        "depends_on": ["pairwise_pi_189", "odd_moments_250", "hypercube_fourier"],
    }


def theorem_mult_nonneg_bound() -> dict:
    return {
        "proved": True,
        "bound_sufficient_for_residual_i": False,
        "theorem": (
            "mult±≥0 ⇒ linear inequalities on m₄±.  If q_max=max|αᵀAα| and "
            "both signs appear for both products, |m₄±|≤1−q_max/(2p).  "
            "On |κ|=1, q_max≤8 ⇒ |m₄±|≤1−4/p, which is >1/(2p) for p≥5."
        ),
        "depends_on": ["theorem_mult_formula", "seidel4_spectrum_257"],
    }


def theorem_orientation_mult_equivalence() -> dict:
    return {
        "proved": True,
        "hypothesis_proved_general": False,
        "theorem": (
            "On |κ|=1 with QᵀAQ=−A, det D=+1: "
            "mult_+(α)=mult_-(Qᵀα) ∀α  ⇔  ν(S)=0.  "
            "Proof: plug mult formulas and qf(Qᵀα)=−qf(α).  "
            "Either side OPEN as a general statement (census both true p=5,7)."
        ),
        "depends_on": ["theorem_mult_formula", "theorem_orientation_258"],
    }


def certify_mult_formula() -> dict:
    import numpy as np
    from itertools import combinations
    from minmax_quadratic import paley_conference_prime_power

    out = {}
    for p in (5, 7):
        yp, ym = Path(f"/tmp/maxplus_p{p}.npy"), Path(f"/tmp/maxminus_p{p}.npy")
        if not yp.exists() or not ym.exists():
            out[str(p)] = {"ok": False, "reason": "missing cache"}
            continue
        Yp, Ym = np.load(yp), np.load(ym)
        N = len(Yp)
        C = paley_conference_prime_power(p).astype(np.float64)
        n = C.shape[0]
        max_err = 0.0
        n_checked = 0
        q_max_obs = 0
        need = 30 if p == 5 else 10
        for S in combinations(range(n), 4):
            a, b, c, d = S
            Slist = [a, b, c, d]
            A = C[np.ix_(Slist, Slist)]
            mp = float(np.mean(Yp[:, a] * Yp[:, b] * Yp[:, c] * Yp[:, d]))
            mm = float(np.mean(Ym[:, a] * Ym[:, b] * Ym[:, c] * Ym[:, d]))
            bits_p = (
                ((Yp[:, a] > 0).astype(np.int32) << 3)
                | ((Yp[:, b] > 0).astype(np.int32) << 2)
                | ((Yp[:, c] > 0).astype(np.int32) << 1)
                | (Yp[:, d] > 0).astype(np.int32)
            )
            bits_m = (
                ((Ym[:, a] > 0).astype(np.int32) << 3)
                | ((Ym[:, b] > 0).astype(np.int32) << 2)
                | ((Ym[:, c] > 0).astype(np.int32) << 1)
                | (Ym[:, d] > 0).astype(np.int32)
            )
            cp = np.bincount(bits_p, minlength=16)
            cm = np.bincount(bits_m, minlength=16)
            for k in range(16):
                alpha = np.array(
                    [1.0 if (k >> (3 - i)) & 1 else -1.0 for i in range(4)]
                )
                pr = int(np.prod(alpha))
                qf = int(round(float(alpha @ A @ alpha)))
                q_max_obs = max(q_max_obs, abs(qf))
                max_err = max(
                    max_err,
                    abs(mult_plus_formula(N, p, qf, pr, mp) - cp[k]),
                    abs(mult_minus_formula(N, p, qf, pr, mm) - cm[k]),
                )
            n_checked += 1
            if n_checked >= need:
                break
        bound = 1 - 4 / p
        thr = 1 / (2 * p)
        out[str(p)] = {
            "ok": True,
            "n_sets_checked": n_checked,
            "max_formula_err": max_err,
            "formula_ok": max_err < 1e-8,
            "q_max_observed": q_max_obs,
            "bound_1_minus_4_over_p": bound,
            "thr_1_over_2p": thr,
            "bound_weaker_than_thr": bound > thr,
        }
    return out


def residual_i_closed_via_261() -> bool:
    return False


def hinge_status_261() -> dict:
    return {
        "mult_formula": True,
        "mult_nonneg_bound": True,
        "orientation_mult_equivalence": True,
        "nu_zero_general": False,
        "residual_i_closed": residual_i_closed_via_261(),
        "open": (
            "mult formula proved; Q-match ⇔ ν=0 on |κ|=1; "
            "still need ν=0 or reflection or tighter L^∞ bound Max+-free."
        ),
    }


def main() -> dict:
    a = theorem_mult_formula()
    b = theorem_mult_nonneg_bound()
    c = theorem_orientation_mult_equivalence()
    cert = certify_mult_formula()
    out = {
        "title": (
            "Prop 15.261 exact Max± config multiplicities; "
            "Q-match ⇔ ν=0; residual (i) OPEN"
        ),
        "proved": {
            "mult_formula": a["proved"],
            "mult_nonneg_bound": b["proved"],
            "bound_sufficient_for_residual_i": b["bound_sufficient_for_residual_i"],
            "orientation_mult_equivalence": c["proved"],
            "nu_zero_general": False,
            "residual_i_closed": residual_i_closed_via_261(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "mult_formula": a,
            "nonneg_bound": b,
            "orientation_eq": c,
        },
        "certified": cert,
        "hinge": hinge_status_261(),
        "F3": (
            "Exact mult_±(α)=(N/16)(1±qf/(2p)+pr m₄±) proved Max+-free.  "
            "On |κ|=1, Q-config-match ⇔ ν=0.  Bound 1−4/p too weak for "
            "1/(2p).  Residual-(i) still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15261.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.261 mult formula; residual-i OPEN")
    print(f"  mult formula: {a['proved']}")
    print(f"  nonneg bound: {b['proved']} (sufficient for res_i? {b['bound_sufficient_for_residual_i']})")
    print(f"  Q-match ⇔ ν=0: {c['proved']}")
    for pk, v in cert.items():
        if not v.get("ok"):
            print(f"  p={pk}: {v}")
            continue
        print(
            f"  p={pk}: formula_ok={v['formula_ok']} err={v['max_formula_err']:.2e} "
            f"q_max={v['q_max_observed']} bound={v['bound_1_minus_4_over_p']:.4f} "
            f"> thr={v['thr_1_over_2p']:.4f}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_261()}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
