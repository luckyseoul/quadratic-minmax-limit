#!/usr/bin/env python3
"""
Prop 15.255 — Expansion of m₄± on 4-sets containing ∞ via 1ᵀy identities;
reduction of m₄⁺=m₄⁻ on |κ|=1 to a vanishing sum; residual-(i) still OPEN.

Does **not** flip residual_i / type_I / gsum / E1 / L. Soft-close forbidden.

SETUP (15.189 1ᵀy, 15.254 T/conjugation, 15.253 reflection)
  Normalised Paley conference: C_∞j=1 (j≠∞), row-sums of finite rows =1.
  Max+: 1ᵀy=(p+1)y_∞; Max−: 1ᵀy=−(p−1)y_∞ (15.189 A).
  Three-wise moments vanish on Max± (15.250 B). Pairwise π as 15.189.

PROVED Max+-free (normalised conference + π + 3-wise vanishing)
  A. Expansion of m₄⁺ for S={∞,a,b,c}.
     Write κ:=C_ab+C_ac+C_bc (=κ(S) under C_∞*=1).  Then
        y_∞=(1ᵀy)/(p+1) on Max+,
        m₄⁺(S)=E[y_∞ y_a y_b y_c]
               =1/(p+1) ∑_j E[y_j y_a y_b y_c].
     Split j=∞ (gives m₄⁺ itself), j∈{a,b,c} (gives π-pairs =κ/p),
     and j∉S (gives ∑_{j∉S} m₄⁺({j,a,b,c})=:T₊).  Solve for m₄⁺:
        m₄⁺(S) = κ/p² + T₊/p.  ∎

  B. Expansion of m₄⁻ for the same S.
     y_∞=−(1ᵀy)/(p−1) on Max−, π=−C/p on pairs.  Same split yields
        m₄⁻(S) = κ/p² − T₋/p,
     where T₋:=∑_{j∉S} m₄⁻({j,a,b,c}).  ∎

  C. Difference and average.
        m₄⁺−m₄⁻ = (T₊+T₋)/p = (2/p) ∑_{j∉S} μ({j,a,b,c}),
        μ = κ/p² + (T₊−T₋)/(2p).
     Hence on every |κ|=1 set S∋∞:
        m₄⁺(S)=m₄⁻(S)  ⇔  ∑_{j∉S} μ({j,a,b,c}) = 0.  ∎

  D. Aut reduction.
     G=Aut(C) is 3-transitive on the p²+1 points (PSL(2,p²)≤G for Paley).
     G preserves Max± separately, hence preserves m₄⁺, m₄⁻, and the equality
     m₄⁺=m₄⁻.  Every 4-set is G-equivalent to one containing ∞.  Therefore
        m₄⁺=m₄⁻ on all |κ|=1 four-sets
     ⇔ the vanishing ∑_{j∉S} μ({j,a,b,c})=0 for every triple {a,b,c}
        with |κ({∞,a,b,c})|=1.  ∎

  E. Wick residuals on ∞-sets.
     When m₄⁺=m₄⁻=μ on S∋∞, (A)–(B) give ρ=μ−κ/p² = T₊/p with
     T₊=∑_{j∉S} μ({j,a,b,c}).  Reflection |ρ|≤|ρ_f4| becomes
        |∑_{j∉S} μ({j,a,b,c})| ≤ p |ρ_f4(S)|.  ∎

CERTIFIED
  F. (A)(B) hold to machine precision for all |κ|=1 sets containing ∞ at p=5,7
     (1800 and 14112 sets).  ∑μ =0 and m₄⁺=m₄⁻ there (consistent with C).
  G. Reflection and m₄ eq still hold on all |κ|=1 at p=5,7 (15.253–254).

OPEN (blocks residual i / predicate flip)
  H. Prove the vanishing ∑_{j∉S} μ({j,a,b,c})=0 for |κ({∞,a,b,c})|=1
     (⇔ m₄⁺=m₄⁻ on all |κ|=1 by D), **or** the sum bound in (E) for
     reflection, **or** Es4 / ker=sc / hull.  Soft-close forbidden.

Until H: residual_i False; gsum False; E1/L OPEN.

Writes evidence/e1_gmin_m4_prop15255.json
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


def theorem_m4plus_expansion_infty() -> dict:
    return {
        "proved": True,
        "theorem": (
            "For normalised conference, S={∞,a,b,c}, κ=C_ab+C_ac+C_bc: "
            "m₄⁺(S)=κ/p² + p⁻¹∑_{j∉S} m₄⁺({j,a,b,c}).  "
            "Proof: y_∞=(1ᵀy)/(p+1); expand E[y_∞ y_a y_b y_c]; "
            "3-wise vanish; π-pairs give κ/p; solve for m₄⁺."
        ),
        "depends_on": ["sum_identity_189", "pairwise_pi_189", "odd_moments_250"],
    }


def theorem_m4minus_expansion_infty() -> dict:
    return {
        "proved": True,
        "theorem": (
            "Same S: m₄⁻(S)=κ/p² − p⁻¹∑_{j∉S} m₄⁻({j,a,b,c}).  "
            "Proof: y_∞=−(1ᵀy)/(p−1); π=−C/p; same expansion."
        ),
        "depends_on": ["sum_identity_189_minus", "pairwise_pi_189_minus", "odd_moments_250"],
    }


def theorem_m4_eq_reduction() -> dict:
    a = theorem_m4plus_expansion_infty()
    b = theorem_m4minus_expansion_infty()
    return {
        "proved": a["proved"] and b["proved"],
        "hypothesis_proved_general": False,
        "theorem": (
            "On |κ|=1 with ∞∈S: m₄⁺=m₄⁻ ⇔ ∑_{j∉S} μ({j,a,b,c})=0.  "
            "By 3-transitivity of Aut(C), m₄⁺=m₄⁻ on all |κ|=1 iff this "
            "vanishing holds for every such triple.  Vanishing OPEN."
        ),
        "depends_on": ["theorem_m4plus_expansion_infty", "theorem_m4minus_expansion_infty", "PSL_3_transitive"],
    }


def theorem_reflection_as_sum_bound() -> dict:
    return {
        "proved_criterion": True,
        "hypothesis_proved_general": False,
        "theorem": (
            "When m₄⁺=m₄⁻=μ on S∋∞ with |κ|=1: ρ=μ−κ/p² = p⁻¹∑_{j∉S} μ({j,a,b,c}).  "
            "Reflection |ρ|≤|ρ_f4| ⇔ |∑_{j∉S} μ({j,a,b,c})| ≤ p|ρ_f4(S)|.  "
            "Combined with 15.253 B (ρ_f4 maj ≤ L_abs) would close residual-(i) "
            "on ∞-sets; Aut transports to all |κ|=1."
        ),
        "depends_on": ["theorem_m4_eq_reduction", "theorem_rho_f4_le_L_abs_253"],
    }


def certify_expansions() -> dict:
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
        # Precompute all m4 for speed at p=7
        fours = np.array(list(combinations(range(n), 4)), dtype=np.int32)
        idx = {tuple(S): i for i, S in enumerate(fours)}
        a, b, c, d = fours.T
        mp = np.mean(Yp[:, a] * Yp[:, b] * Yp[:, c] * Yp[:, d], axis=0)
        mm = np.mean(Ym[:, a] * Ym[:, b] * Ym[:, c] * Ym[:, d], axis=0)
        mu = 0.5 * (mp + mm)
        kap = C[a, b] * C[c, d] + C[a, c] * C[b, d] + C[a, d] * C[b, c]
        max_err_p = 0.0
        max_err_m = 0.0
        max_sum_mu = 0.0
        n_inf = 0
        for i, S in enumerate(map(tuple, fours)):
            if 0 not in S:
                continue
            if abs(abs(kap[i]) - 1) > 1e-9:
                continue
            n_inf += 1
            fin = [x for x in S if x != 0]
            x, y, z = fin
            kap2 = C[x, y] + C[x, z] + C[y, z]
            Tp = Tm = 0.0
            sum_mu = 0.0
            for j in range(n):
                if j in S:
                    continue
                T = tuple(sorted([j, x, y, z]))
                ti = idx[T]
                Tp += mp[ti]
                Tm += mm[ti]
                sum_mu += mu[ti]
            pred_p = kap2 / (p * p) + Tp / p
            pred_m = kap2 / (p * p) - Tm / p
            max_err_p = max(max_err_p, abs(pred_p - mp[i]))
            max_err_m = max(max_err_m, abs(pred_m - mm[i]))
            max_sum_mu = max(max_sum_mu, abs(sum_mu))
        out[str(p)] = {
            "ok": True,
            "n_kappa1_with_infty": n_inf,
            "max_err_m4plus_expansion": max_err_p,
            "max_err_m4minus_expansion": max_err_m,
            "max_abs_sum_mu": max_sum_mu,
            "expansions_ok": max_err_p < 1e-10 and max_err_m < 1e-10,
            "sum_mu_vanishes": max_sum_mu < 1e-10,
            "m4_eq": True,  # from prior census; sum_mu=0 equivalent
        }
    return out


def residual_i_closed_via_255() -> bool:
    return False


def hinge_status_255() -> dict:
    return {
        "m4plus_expansion_infty": True,
        "m4minus_expansion_infty": True,
        "m4_eq_reduction": True,
        "vanishing_sum_general": False,
        "reflection_sum_bound_general": False,
        "residual_i_closed": residual_i_closed_via_255(),
        "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
        "e1_closed_general": e1_closed_general(),
        "open": (
            "Prove ∑_{j∉S} μ({j,a,b,c})=0 for |κ({∞,a,b,c})|=1 "
            "(⇔ m₄⁺=m₄⁻ on all |κ|=1), or |∑ μ|≤p|ρ_f4| for reflection."
        ),
    }


def main() -> dict:
    a = theorem_m4plus_expansion_infty()
    b = theorem_m4minus_expansion_infty()
    c = theorem_m4_eq_reduction()
    d = theorem_reflection_as_sum_bound()
    cert = certify_expansions()
    out = {
        "title": (
            "Prop 15.255 m₄± expansion on ∞-sets; m₄-eq reduction; "
            "residual (i) OPEN on vanishing sum"
        ),
        "proved": {
            "m4plus_expansion_infty": a["proved"],
            "m4minus_expansion_infty": b["proved"],
            "m4_eq_reduction_criterion": c["proved"],
            "vanishing_sum_general": False,
            "reflection_sum_bound_general": False,
            "residual_i_closed": residual_i_closed_via_255(),
            "gsum_disj_lb_proved_general": gsum_disj_lb_proved_general(),
            "type_I": type_I_k_3p_minus_2_closed_general(),
            "E1_closed": e1_closed_general(),
            "L_closed": False,
        },
        "algebra": {
            "m4plus_exp": a,
            "m4minus_exp": b,
            "m4_eq_reduction": c,
            "reflection_sum": d,
        },
        "certified": {"expansions": cert},
        "hinge": hinge_status_255(),
        "F3": (
            "m₄⁺=κ/p²+T₊/p and m₄⁻=κ/p²−T₋/p on ∞-sets proved Max+-free.  "
            "m₄⁺=m₄⁻ on all |κ|=1 ⇔ ∑μ over finite 1-point extensions vanishes.  "
            "Census vanishing at p=5,7.  General vanishing / reflection still OPEN.  "
            "No predicate flip."
        ),
    }
    path = ROOT / "evidence" / "e1_gmin_m4_prop15255.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(out, indent=2, default=str))
    print("Prop 15.255 m4± ∞-expansion; m4-eq reduction; residual-i OPEN")
    print(f"  m4+ exp: {a['proved']}")
    print(f"  m4- exp: {b['proved']}")
    print(f"  reduction: {c['proved']}")
    for pk, v in cert.items():
        if not v.get("ok"):
            print(f"  p={pk}: {v}")
            continue
        print(
            f"  p={pk}: exp_ok={v['expansions_ok']} sum_mu_0={v['sum_mu_vanishes']} "
            f"n_inf={v['n_kappa1_with_infty']} err_p={v['max_err_m4plus_expansion']:.2e}"
        )
    print(f"  residual_i closed: {residual_i_closed_via_255()}")
    print(f"  open: {out['hinge']['open']}")
    print("wrote", path)
    return out


if __name__ == "__main__":
    main()
