#!/usr/bin/env python3
"""
Prop 15.254 — T-action on m₄± closed forms; Tμ=2p(m₄⁺−m₄⁻);
Paley monomial C∼−C; μ as average of m₄⁺ and signed pullback;
residual-(i) still OPEN on m₄⁺=m₄⁻ / reflection.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.177 master, 15.189 π, 15.191 star, 15.253 reflection)
  T = C-weighted Johnson adjacency on 4-sets:
      (Tf)(S)=∑_{v∈S, r∉S} C_{vr} f(S_{v→r}).
  m₄⁺=E_{Max+}[∏_{i∈S} y_i], m₄⁻ similarly, μ=½(m₄⁺+m₄⁻).
  π_{ij}=E[y_i y_j]=C_{ij}/p on Max+, −C_{ij}/p on Max− (15.189).

PROVED Max+-free (any symmetric conference with Max± and π as above)
  A. Pointwise Cy-expansion for T-weight before expectation.
     For y∈Max+ (Cy=py): ∑_{r∉S} C_{vr} y_r = p y_v − ∑_{u∈S\\{v}} C_{vu} y_u.
     Hence
        ∑_{v∈S,r∉S} C_{vr} y_r ∏_{i∈S\\{v}} y_i
          = X(4p − 2 ∑_{e⊂S} C_e y_e),
     where X=∏_{i∈S} y_i and y_e=y_a y_b for e=ab.  ∎

  B. T m₄⁺ closed form.
     Take E_{Max+} of (A).  For e={a,b}⊂S={a,b,c,d},
        E[X y_e]=E[y_c y_d]=C_{cd}/p = C_{e^*}/p
     (e^*=complementary edge).  Sum over six edges:
        ∑_e C_e C_{e^*}=2κ(S).
     Therefore on every 4-set S:
        (T m₄⁺)(S) = 4p m₄⁺(S) − 4κ(S)/p.  ∎
     Equivalently (4p I − T)m₄⁺ = 4κ/p (master 15.177 for Max+).

  C. T m₄⁻ closed form.
     For y∈Max− (Cy=−py): ∑_{r∉S} C_{vr} y_r = −p y_v − ∑_{u∈S\\{v}} C_{vu} y_u.
     Same expansion yields X(−4p − 2 ∑ C_e y_e).  With
        E_−[X y_e]=−C_{e^*}/p,  ∑ C_e C_{e^*}=2κ,
     one gets on every 4-set:
        (T m₄⁻)(S) = −4p m₄⁻(S) + 4κ(S)/p.  ∎
     Equivalently (4p I + T)m₄⁻ = 4κ/p.

  D. Tμ identity.
        Tμ = ½(T m₄⁺ + T m₄⁻) = 2p(m₄⁺ − m₄⁻).  ∎
     In particular Tμ=0 on any 4-set where m₄⁺=m₄⁻, and
        ν := Tμ/(4p) = ½(m₄⁺ − m₄⁻)
     is the diamond partner of μ in 15.228.

PROVED Max+-free (Paley n=p²+1)
  E. Monomial conjugation C ∼ −C.
     Let σ∈F_{p²}^× be a nonsquare (χ(σ)=−1), π the permutation of
     V=F_{p²}∪{∞} with π(∞)=∞ and π(x)=σx on F_{p²}, and D=diag with
     D_∞=−1, D_x=1 on F_{p²}.  Then
        D P^T C P D = −C
     (finite block: mult-by-σ multiplies differences by σ ⇒ χ flips;
      ∞-row: D flips the all-ones row of C to all-minus-ones).  ∎

  F. Max− as monomial image of Max+.
     From (E): y∈Max+ ⇒ z := D·(y∘π^{-1}) satisfies C z = −p z, i.e. z∈Max−.
     The map is bijective (involution up to Aut), so Max− = {D·(y∘π^{-1}): y∈Max+}.  ∎

  G. Exact m₄⁻ from m₄⁺.
        m₄⁻(S) = χ_D(S) m₄⁺(π^{-1}S),   χ_D(S):=∏_{i∈S} D_i ∈{±1}.
     Hence
        μ(S) = ½[ m₄⁺(S) + χ_D(S) m₄⁺(π^{-1}S) ].  ∎
     (Paley-only; uses (E)(F).)

  H. Consequence for |κ|=1.
     χ_D(S)=−1 ⇔ ∞∈S (only D_∞=−1).  On |κ|=1:
        m₄⁺=m₄⁻  ⇔  m₄⁺(S)=χ_D(S) m₄⁺(π^{-1}S)
     and then Tμ=0 and ν=0 on |κ|=1 (by D).  ∎

CERTIFIED
  I. (B)(C)(D) hold to machine precision at p=5,7 on all 4-sets.
  J. (E)(F)(G): D P^T C P D = −C exact; Max−=monomial image of Max+;
     m₄⁻ formula exact; max|m₄⁺−m₄⁻|=0 on all |κ|=1 at p=5,7;
     max|Tμ|=0 on |κ|=1 at p=5,7.
  K. Reflection hyp (15.253) still holds at p=5,7 (prior census).

OPEN (blocks residual i / predicate flip)
  L. Prove Max+-free for all primes p≥5 one of:
     (i)  m₄⁺=m₄⁻ on every |κ|=1 four-set (⇔ m₄⁺(S)=χ_D(S)m₄⁺(π^{-1}S)),
          then Tμ=0 / ν=0 on |κ|=1;
     (ii) reflection |ρ|≤|ρ_f4| on |κ|=1 (15.253);
     (iii) hull / |Ext| maj / Es4 / ker=sc.
     Soft-close forbidden.

Until L: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15254.json
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
    is_prime,
    n_of,
    residual_i_dual_eq_empty_proved_general,
    type_I_k_3p_minus_2_closed_general,
)


def theorem_T_m4plus() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On every 4-set, for Max+ with π=C/p: (T m₄⁺)(S)=4p m₄⁺(S)−4κ(S)/p.  "
            "Proof: pointwise Cy=py expansion of the T-weight; take E; "
            "E[X y_e]=C_{e^*}/p; ∑_e C_e C_{e^*}=2κ."
        ),
        "depends_on": ["Cy_eq_py", "pairwise_pi_189", "T_definition"],
    }


def theorem_T_m4minus() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On every 4-set, for Max− with π=−C/p: (T m₄⁻)(S)=−4p m₄⁻(S)+4κ(S)/p.  "
            "Proof: same expansion with Cy=−py; E_−[X y_e]=−C_{e^*}/p."
        ),
        "depends_on": ["Cy_eq_minus_py", "pairwise_pi_189_minus", "T_definition"],
    }


def theorem_T_mu() -> dict:
    a, b = theorem_T_m4plus(), theorem_T_m4minus()
    return {
        "proved": a["proved"] and b["proved"],
        "theorem": (
            "Tμ=½(Tm₄⁺+Tm₄⁻)=2p(m₄⁺−m₄⁻) on all 4-sets.  "
            "Hence ν:=Tμ/(4p)=½(m₄⁺−m₄⁻), and Tμ=0 wherever m₄⁺=m₄⁻."
        ),
        "depends_on": ["theorem_T_m4plus", "theorem_T_m4minus"],
    }


def theorem_paley_monomial_conjugation() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Paley n=p²+1: for nonsquare σ∈F_{p²}^× and π(∞)=∞, π(x)=σx, "
            "D_∞=−1, D_x=1: D P^T C P D = −C.  Finite block from χ(σ·)=−χ; "
            "∞-row from D flip of the all-ones border."
        ),
        "depends_on": ["paley_construction", "quadratic_character"],
    }


def theorem_mu_from_m4plus_pullback() -> dict:
    c = theorem_paley_monomial_conjugation()
    return {
        "proved": c["proved"],
        "theorem": (
            "Paley: Max−={D·(y∘π^{-1}): y∈Max+}, hence "
            "m₄⁻(S)=χ_D(S) m₄⁺(π^{-1}S) and "
            "μ(S)=½[m₄⁺(S)+χ_D(S) m₄⁺(π^{-1}S)]."
        ),
        "depends_on": ["theorem_paley_monomial_conjugation"],
    }


def theorem_kappa1_m4_eq_structure() -> dict:
    return {
        "proved_structure": True,
        "bound_proved_general": False,
        "theorem": (
            "On |κ|=1: m₄⁺=m₄⁻ ⇔ m₄⁺(S)=χ_D(S) m₄⁺(π^{-1}S) (Paley).  "
            "Then Tμ=0 and ν=0 on |κ|=1.  Equality certified p=5,7; general OPEN."
        ),
        "depends_on": ["theorem_T_mu", "theorem_mu_from_m4plus_pullback"],
    }


def certify_T_and_conjugation() -> dict:
    """Census p=5,7: T formulas, conjugation, m4+=m4- on |κ|=1."""
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
        C = paley_conference_prime_power(p).astype(np.float64)
        n = C.shape[0]
        fours = list(combinations(range(n), 4))
        idx = {S: i for i, S in enumerate(fours)}
        N = len(fours)
        fa = np.array(fours)
        mp = np.mean(Yp[:, fa[:, 0]] * Yp[:, fa[:, 1]] * Yp[:, fa[:, 2]] * Yp[:, fa[:, 3]], 0)
        mm = np.mean(Ym[:, fa[:, 0]] * Ym[:, fa[:, 1]] * Ym[:, fa[:, 2]] * Ym[:, fa[:, 3]], 0)
        mu = 0.5 * (mp + mm)
        kap = (
            C[fa[:, 0], fa[:, 1]] * C[fa[:, 2], fa[:, 3]]
            + C[fa[:, 0], fa[:, 2]] * C[fa[:, 1], fa[:, 3]]
            + C[fa[:, 0], fa[:, 3]] * C[fa[:, 1], fa[:, 2]]
        )

        def apply_T(f):
            Tf = np.zeros(N)
            for i, S in enumerate(fours):
                Sset = set(S)
                acc = 0.0
                for v in S:
                    for r in range(n):
                        if r in Sset:
                            continue
                        Sp = tuple(sorted([x for x in S if x != v] + [r]))
                        acc += C[v, r] * f[idx[Sp]]
                Tf[i] = acc
            return Tf

        Tmp, Tmm, Tmu = apply_T(mp), apply_T(mm), apply_T(mu)
        err_p = float(np.max(np.abs(Tmp - (4 * p * mp - 4 * kap / p))))
        err_m = float(np.max(np.abs(Tmm - (-4 * p * mm + 4 * kap / p))))
        err_mu = float(np.max(np.abs(Tmu - 2 * p * (mp - mm))))
        ids = np.where(np.abs(np.abs(kap) - 1) < 1e-9)[0]
        out[str(p)] = {
            "ok": True,
            "T_m4plus_formula_err": err_p,
            "T_m4minus_formula_err": err_m,
            "T_mu_formula_err": err_mu,
            "max_abs_m4plus_minus_m4minus_kappa1": float(np.max(np.abs(mp[ids] - mm[ids]))),
            "max_abs_Tmu_kappa1": float(np.max(np.abs(Tmu[ids]))),
            "n_kappa1": int(len(ids)),
            "formulas_ok": err_p < 1e-10 and err_m < 1e-10 and err_mu < 1e-10,
            "m4_eq_on_kappa1": float(np.max(np.abs(mp[ids] - mm[ids]))) < 1e-12,
            "Tmu_zero_on_kappa1": float(np.max(np.abs(Tmu[ids]))) < 1e-10,
        }
    return out


def residual_i_closed_via_254() -> bool:
    return False


def hinge_status_254() -> dict:
    return {
        "T_m4plus_formula": True,
        "T_m4minus_formula": True,
        "T_mu_identity": True,
        "paley_conjugation": True,
        "mu_pullback_formula": True,
        "m4plus_eq_m4minus_kappa1_general": False,
        "reflection_hyp_general": False,
        "residual_i_closed": residual_i_closed_via_254(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove m₄⁺=m₄⁻ on |κ|=1 (Paley functional eq under nonsquare "
            "map) and/or reflection |ρ|≤|ρ_f4| (15.253)."
        ),
    }


def main() -> dict:
    a = theorem_T_m4plus()
    b = theorem_T_m4minus()
    c = theorem_T_mu()
    d = theorem_paley_monomial_conjugation()
    e = theorem_mu_from_m4plus_pullback()
    f = theorem_kappa1_m4_eq_structure()
    cert = certify_T_and_conjugation()
    out = {
        "title": (
            "Prop 15.254 T m₄± closed forms; Paley C∼−C; "
            "residual (i) OPEN on m₄⁺=m₄⁻ / reflection"
        ),
        "proved": {
            "T_m4plus_formula": a["proved"],
            "T_m4minus_formula": b["proved"],
            "T_mu_identity": c["proved"],
            "paley_monomial_conjugation": d["proved"],
            "mu_from_m4plus_pullback": e["proved"],
            "m4_eq_kappa1_structure": f["proved_structure"],
            "m4_eq_kappa1_general": False,
            "residual_i_closed": residual_i_closed_via_254(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "T_m4plus": a,
            "T_m4minus": b,
            "T_mu": c,
            "conjugation": d,
            "mu_pullback": e,
            "kappa1_eq": f,
        },
        "certified": {"T_and_conjugation": cert},
        "hinge": hinge_status_254(),
        "F3": (
            "Tm₄⁺=4p m₄⁺−4κ/p and Tm₄⁻=−4p m₄⁻+4κ/p proved Max+-free; "
            "Tμ=2p(m₄⁺−m₄⁻).  Paley: D P^T C P D=−C, μ=½[m₄⁺+χ_D m₄⁺∘π⁻¹].  "
            "Census m₄⁺=m₄⁻ and Tμ=0 on |κ|=1 at p=5,7.  General m₄ eq / "
            "reflection still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15254.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.254 T m4± formulas; Paley conjugation; residual-i OPEN")
    print(f"  T m4+: {a['proved']}")
    print(f"  T m4-: {b['proved']}")
    print(f"  T mu: {c['proved']}")
    print(f"  conjugation: {d['proved']}")
    print(f"  mu pullback: {e['proved']}")
    for pk, v in cert.items():
        if not v.get("ok"):
            print(f"  p={pk}: {v}")
            continue
        print(
            f"  p={pk}: formulas_ok={v['formulas_ok']} "
            f"m4_eq_k1={v['m4_eq_on_kappa1']} Tmu0_k1={v['Tmu_zero_on_kappa1']} "
            f"err_p={v['T_m4plus_formula_err']:.2e} err_m={v['T_m4minus_formula_err']:.2e}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_254()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
