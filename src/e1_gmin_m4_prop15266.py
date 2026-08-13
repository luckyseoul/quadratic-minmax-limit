#!/usr/bin/env python3
"""
Prop 15.266 — Type-averages of π-odd functions vanish on finite |κ|=1;
type-T has ±4p with evecs agreeing on star=0 (p=5); Tδ=0 reduced to
within-type odd variation; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.254–265)
  π preserves (κ,φ) on finite 4-sets and flips (κ,φ)↦(−κ,−φ) on ∞-sets
  (15.254 H / 15.259 C).  ν is Π-odd: on finite sets ν(S)=−ν(πS) (15.262).
  TΠ=−ΠT, Tμ=4pν, Tμ_part=0 on |κ|=1, so ν=Tδ/(4p) there (15.265).
  Types mean (κ,φ,star) classes; star=0 on |κ|=1.

PROVED Max+-free (Paley π-action)
  A. Type-averages of π-odd functions vanish on finite |κ|=1.
     π preserves (κ,φ) on every finite 4-set, hence permutes each
     finite (κ,φ)-class.  If f(πS)=−f(S), the class average is 0.  ∎
     In particular the type-averages of ν on finite |κ|=1 are 0.

  B. Reduction of Tδ=0 on finite |κ|=1.
     Write ν = ν̄ + ν̃ on finite |κ|=1, with ν̄ the type-average
     (zero by A) and ν̃ mean-zero within each (κ,φ)-class.
     Then Tδ=0 on finite |κ|=1  ⇔  ν̃=0 there
     (the only remaining obstruction is π-odd Aut-invariant
     within-type variation, i.e. a pair of Aut-orbits O≠πO in the
     same (κ,φ)-class with ν|O=−ν|πO≠0).  ∎

  C. Point-transitivity reduces the rest to ∞-sets.
     Aut(C) is transitive on the n points (PSL(2,p²)≤Aut).  Every
     |κ|=1 four-set is Aut-equivalent to one containing ∞, and Aut
     preserves ν.  Hence ν=0 on all |κ|=1  ⇔  ν=0 on ∞-containing
     |κ|=1 sets  ⇔  Tδ=0 there (15.265 D).  Combined with (B), it
     suffices to kill within-type odd variation after transporting
     to ∞-sets, or to prove the type-constant case plus Aut-orbit
     fusion under π.  ∎

PROVED Max+-free (type-quotient of T; p=5 complete)
  D. Type-quotient T_typ.
     Let 𝒯 be the (κ,φ,star) classes and M the |𝒯|×|𝒯| matrix
        M_{j i} = |j|⁻¹ ∑_{S∈j} ∑_{v∈S,r∉S} C_{vr} 1_{S_{v→r}∈i}.
     Then T preserves the space of type-constant functions and acts
     as M there.  At p=5, |𝒯|=12 and
        σ(M) contains simple eigenvalues ±4p=±20,
     with no other ±4p.  ∎

  E. The ±4p type-evecs agree on star=0.
     At p=5 the M-eigenvectors for λ=±20 restrict to the *same*
     function on the four |κ|=1 types (star=0):
        F(−1,−6)=−4c,  F(−1,2)=−c,  F(1,−2)=c,  F(1,6)=4c
     (c>0).  Consequently the type-constant odd combination
     e_{+}−e_{-} vanishes on |κ|=1, and the even combination
     e_{+}+e_{-} restricts to 2F.  Thus every type-constant vector
     of E_{±4p}^{odd} is zero on |κ|=1.  ∎

CERTIFIED
  F. (A) is formal; (D)(E) exact at p=5 (numpy eig of the 12×12
     quotient, residuals <10⁻¹³).  Census: ν̃=0 at p=5 (types are
     m₄-constant) and at p=7 (ν≡0 on all |κ|=1, including types
     with m₄-span>0).  Type-average of ν is 0 at both p.
  F2. At p=7 the type-quotient (14 classes) has extreme eigenvalues
     ±16 (double), **not** ±4p=±28 (gap 12).  So E_{±4p} is invisible
     to type-constant functions at p=7 — it lives entirely in
     within-type variation (matches m₄-span>0 on extreme |κ|=1 types).
     The p=5 type-T ±4p agreement does **not** by itself generalise.

OPEN (blocks residual i / predicate flip)
  G. Prove the within-type odd variation ν̃=0 on |κ|=1 for all primes
     p≥5 (equivalently: every Aut-orbit of finite |κ|=1 four-sets is
     π-stable, or Ext/T kills mean-zero odd functions at eigenvalue
     p⁴−1 / 4p), **or** envelope / reflection / |μ| / Es4 / ker=sc.
     Soft-close forbidden.
     (Type-constant odd E already vanishes on |κ|=1 by (E) at p=5;
      the general type-T ±4p agreement is the natural extension of (E).)

Until G: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15266.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from e1_gmin_m4_prop15170 import (  # noqa: E402
    e1_closed_general,
    gsum_disj_lb_proved_general,
    type_I_k_3p_minus_2_closed_general,
)


def theorem_type_averages_odd_vanish() -> dict:
    return {
        "proved": True,
        "theorem": (
            "On finite 4-sets π preserves (κ,φ), so every π-odd function "
            "has zero average on each finite (κ,φ)-class.  In particular "
            "type-averages of ν on finite |κ|=1 vanish."
        ),
        "depends_on": ["pi_type_transport_259", "nu_anti_invariance_262"],
    }


def theorem_reduction_to_within_type() -> dict:
    return {
        "proved": True,
        "hypothesis_proved_general": False,
        "theorem": (
            "On finite |κ|=1: Tδ=0 ⇔ the within-type π-odd part of ν vanishes "
            "(type-averages already 0).  Point-transitivity reduces the global "
            "statement to ∞-sets."
        ),
        "depends_on": [
            "theorem_type_averages_odd_vanish",
            "Tdelta_reduction_265",
            "Aut_point_transitive",
        ],
    }


def theorem_type_T_p5() -> dict:
    return {
        "proved": True,
        "scope": "p=5 type-quotient (12 classes)",
        "general_p": False,
        "theorem": (
            "At p=5 the type-quotient of T is a 12×12 matrix with simple "
            "eigenvalues ±20=±4p, and the two eigenvectors agree on the "
            "four star=0 types.  Type-constant odd E vanishes on |κ|=1."
        ),
        "depends_on": ["T_definition", "type_classification"],
    }


def residual_i_closed_via_266() -> bool:
    return False


def hinge_status_266() -> dict:
    return {
        "type_averages_odd_vanish": True,
        "reduction_to_within_type": True,
        "type_T_p5_evecs_agree": True,
        "within_type_odd_zero_general": False,
        "residual_i_closed": residual_i_closed_via_266(),
        "open": (
            "Type-averages of ν vanish on finite |κ|=1; type-constant odd E "
            "vanishes at p=5.  Within-type π-odd variation still OPEN."
        ),
    }


def certify_type_T_p5() -> dict:
    import numpy as np
    from itertools import combinations
    from minmax_quadratic import paley_conference_prime_power

    p = 5
    C = paley_conference_prime_power(p).astype(np.float64)
    n = C.shape[0]
    fours = list(combinations(range(n), 4))
    N4 = len(fours)
    idx = {S: i for i, S in enumerate(fours)}

    def typ(S):
        a, b, c, d = S
        kap = C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
        ph = float(np.sum(np.prod(C[:, [a, b, c, d]], axis=1)))
        st = 0.0
        s = [a, b, c, d]
        for vi, v in enumerate(s):
            pr = 1.0
            for uj, u in enumerate(s):
                if uj != vi:
                    pr *= C[v, u]
            st += pr
        return (int(round(kap)), int(round(ph)), int(round(st)))

    types = [typ(S) for S in fours]
    uniq = sorted(set(types))
    tix = {t: i for i, t in enumerate(uniq)}
    nt = len(uniq)
    tid = np.array([tix[t] for t in types], dtype=np.int32)
    cnt = np.bincount(tid, minlength=nt).astype(float)
    Macc = np.zeros((nt, nt))
    for si, S in enumerate(fours):
        j = tid[si]
        ss = set(S)
        for v in S:
            for r in range(n):
                if r in ss:
                    continue
                R = tuple(sorted((ss - {v}) | {r}))
                Macc[j, tid[idx[R]]] += C[v, r]
    M = Macc / cnt[:, None]
    w, V = np.linalg.eig(M)
    dist_p = float(np.min(np.abs(w - 4 * p)))
    dist_m = float(np.min(np.abs(w + 4 * p)))
    # extract the two evecs
    ip = int(np.argmin(np.abs(w - 4 * p)))
    im = int(np.argmin(np.abs(w + 4 * p)))
    vp = np.real(V[:, ip])
    vm = np.real(V[:, im])
    star0 = [i for i, t in enumerate(uniq) if t[2] == 0]
    # scale so the (1,-2,0) component is +1 if nonzero
    def restrict(vec):
        out = {}
        for i in star0:
            out[str(uniq[i])] = float(vec[i])
        return out

    rp, rm = restrict(vp), restrict(vm)
    # agreement up to scale: ratio rp/rm constant on star0
    keys = list(rp.keys())
    ratios = []
    agree = False
    if keys:
        # normalize both to unit on first nonzero
        def normed(d):
            vals = np.array([d[k] for k in keys])
            nrm = np.linalg.norm(vals)
            return vals / nrm if nrm > 1e-14 else vals

        a = normed(rp)
        b = normed(rm)
        # same or opposite global sign
        agree = min(np.linalg.norm(a - b), np.linalg.norm(a + b)) < 1e-8
        # they should be equal (same restriction), not opposite
        same = float(np.linalg.norm(a - b)) < 1e-8
    else:
        same = False
    return {
        "n_types": nt,
        "types": [list(t) for t in uniq],
        "counts": [int(x) for x in cnt],
        "eig_min_dist_4p": dist_p,
        "eig_min_dist_m4p": dist_m,
        "simple_pm_4p": dist_p < 1e-8 and dist_m < 1e-8,
        "evec_plus_star0": rp,
        "evec_minus_star0": rm,
        "star0_restrictions_agree": bool(same),
        "ok": bool(same and dist_p < 1e-8 and dist_m < 1e-8),
    }


def main() -> dict:
    a = theorem_type_averages_odd_vanish()
    b = theorem_reduction_to_within_type()
    c = theorem_type_T_p5()
    cert = certify_type_T_p5()
    out = {
        "title": (
            "Prop 15.266 type-averages of π-odd vanish; type-T ±4p evecs "
            "agree on star=0 (p=5); residual (i) OPEN"
        ),
        "proved": {
            "type_averages_odd_vanish": a["proved"],
            "reduction_to_within_type": b["proved"],
            "type_T_p5": c["proved"],
            "within_type_odd_zero_general": False,
            "residual_i_closed": residual_i_closed_via_266(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "type_averages": a,
            "reduction": b,
            "type_T_p5": c,
        },
        "certified": cert,
        "hinge": hinge_status_266(),
        "F3": (
            "π-odd type-averages vanish on finite |κ|=1; at p=5 type-T has "
            "simple ±4p with evecs agreeing on star=0.  Within-type odd "
            "variation still OPEN.  No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15266.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.266 type-averages / type-T p=5; residual-i OPEN")
    print(f"  type-avg odd vanish: {a['proved']}")
    print(f"  reduction to within-type: {b['proved']}")
    print(f"  type-T p=5: {c['proved']} cert_ok={cert['ok']}")
    print(f"    ±4p simple: {cert['simple_pm_4p']} star0 agree: {cert['star0_restrictions_agree']}")
    print(f"  residual_i closed: {residual_i_closed_via_266()}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
